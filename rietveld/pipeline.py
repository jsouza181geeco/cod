from pathlib import Path

from models import (
    CandidateInput, CandidateResult, RietveldResult,
    PhaseFraction, MultiPhaseResult,
)
from data_loader import parse_diffractogram as parse_xye
from pattern_calc import build_icalc_unit, build_icalc_unit_absolute
from linear_fit import linear_fit
from multi_phase_fit import multi_phase_fit
from fom import calc_fom
from qpa import weight_fractions

FIXED_PARAMS = {
    'U': 0.01,
    'V': -0.002,
    'W': 0.005,
    'eta': 0.5,
    'n_bg': 4,
    'wavelength': 1.54056,
}


def run_pipeline(
    xye_path: str | Path,
    candidates: list[CandidateInput],
    db_client=None,
    params: dict | None = None,
) -> RietveldResult:
    """
    Rietveld phase identification pipeline [R69; MCC99 p.36-50].

    For each candidate:
      1. build_icalc_unit  → Icalc_unit (scale=1, no background)
      2. linear_fit        → optimal scale S + background Ibg  [Y93 p.18-22]
      3. calc_fom          → Rwp, Rp, Rexp, chi2               [Y93 p.22-26]

    Candidates ranked by Rwp ascending (lower = better match).
    db_client is optional — metadata enriches output but does not affect FoM.
    """
    p = {**FIXED_PARAMS, **(params or {})}

    tth, Iobs, sigma = parse_xye(xye_path)

    metadata_map = {}
    if db_client is not None:
        metadata_map = db_client.fetch_metadata([c.cod_id for c in candidates])

    results = []
    for cand in candidates:
        Icalc_unit, n_used = build_icalc_unit(
            tth, cand.reflections,
            U=p['U'], V=p['V'], W=p['W'], eta=p['eta'],
        )
        scale, Icalc = linear_fit(tth, Iobs, sigma, Icalc_unit, n_bg=p['n_bg'])
        fom = calc_fom(Iobs, Icalc, sigma, n_params=1 + p['n_bg'])

        results.append(CandidateResult(
            cod_id=cand.cod_id,
            Rwp=fom['Rwp'],
            Rp=fom['Rp'],
            Rexp=fom['Rexp'],
            chi2=fom['chi2'],
            scale=scale,
            n_peaks_used=n_used,
            metadata=metadata_map.get(cand.cod_id),
        ))

    results.sort(key=lambda r: r.Rwp)
    return RietveldResult(xye_file=str(xye_path), n_points=len(tth), candidates=results)


def run_multiphase(
    xye_path: str | Path,
    candidates: list[CandidateInput],
    db_client=None,
    top_k: int = 4,
    params: dict | None = None,
) -> MultiPhaseResult:
    """
    Multi-phase Rietveld + QPA for mixtures (mining residues) [BH88; HH87].

    1. run_pipeline (single-phase) → rank candidates by Rwp
    2. select top-K by Rwp RANK — NOT by single-phase viability threshold
       (crit. 6.7: a real mixture phase has high single-phase Rwp, a
       threshold would drop it)
    3. build_icalc_unit_absolute per phase → multi_phase_fit (S_k >= 0)
    4. calc_fom on combined model (n_params = K + n_bg)
    5. weight_fractions (Hill-Howard) — needs db_client for Z, M, V

    QPA requires metadata: with db_client=None, weight_pct = 0 for all
    phases (scales/Rwp still computed). Rwp_single_best uses the ABSOLUTE
    basis too, so the crit. 6.6 comparison Rwp_combined <= Rwp_single is fair.
    """
    p = {**FIXED_PARAMS, **(params or {})}

    # 1-2. single-phase ranking → top-K by Rwp (no metadata needed here)
    single = run_pipeline(xye_path, candidates, db_client=None, params=params)
    top = single.candidates[:top_k]
    top_cods = [c.cod_id for c in top]
    cand_map = {c.cod_id: c for c in candidates}
    chosen = [cand_map[cod] for cod in top_cods]

    tth, Iobs, sigma = parse_xye(xye_path)

    # 3. absolute unit patterns (NOT intensity_rel — crit. 7.6/7.9)
    units = [
        build_icalc_unit_absolute(
            tth, c.reflections, U=p['U'], V=p['V'], W=p['W'], eta=p['eta'],
        )[0]
        for c in chosen
    ]

    # crit. 6.6: best single-phase Rwp on the SAME absolute basis (fair)
    _, Icalc_best, _ = multi_phase_fit(tth, Iobs, sigma, [units[0]], n_bg=p['n_bg'])
    fom_best = calc_fom(Iobs, Icalc_best, sigma, n_params=1 + p['n_bg'])

    # 4. combined fit + FoM
    scales, Icalc, _ = multi_phase_fit(tth, Iobs, sigma, units, n_bg=p['n_bg'])
    fom = calc_fom(Iobs, Icalc, sigma, n_params=len(units) + p['n_bg'])

    # 5. metadata (top-K only) + QPA
    meta_map = {}
    if db_client is not None:
        meta_map = db_client.fetch_metadata(top_cods)
    metadatas = [meta_map.get(cod) for cod in top_cods]
    qpa_rows = weight_fractions(scales, metadatas)

    phases = [
        PhaseFraction(
            cod_id=r['cod_id'], scale=r['scale'], weight_pct=r['weight_pct'],
            Z=r['Z'], M=r['M'], V=r['V'], ZMV=r['ZMV'], metadata=metadatas[i],
        )
        for i, r in enumerate(qpa_rows)
    ]

    return MultiPhaseResult(
        xye_file=str(xye_path), n_points=len(tth), phases=phases,
        Rwp=fom['Rwp'], Rp=fom['Rp'], Rexp=fom['Rexp'], chi2=fom['chi2'],
        Rwp_single_best=fom_best['Rwp'],
    )


def prefilter_candidates(
    xye_path: str | Path,
    candidates: list[CandidateInput],
    top_n: int,
    db_client=None,
    wavelength: float = 1.54056,
    tol_deg: float = 0.2,
    min_matches: int = 3,
    prominence_sigma: float = 5.0,
) -> list[CandidateInput]:
    """Hanawalt pre-selection before Rietveld fit, in d-space [H38; crit. 8.1-8.13].

    Detects experimental peaks (prominence-based [SV75]), converts them to
    d-spacing via wavelength (crit. 8.9), then counts d-matches per candidate.
    Returns top_n sorted by n_matches DESC.

    wavelength : sample radiation (A). Used for 2theta -> d conversion. Default
                 CuKalpha. Pass the actual sample lambda (CoKa/CrKa for Fe-rich
                 samples) so the d windows are correct.

    Modes:
      db_client=None → match_candidates_csv (in-memory; suitable for N < 10k)
      db_client set  → match_candidates_db via peak_fingerprints MV (d_hkl),
                       filtered against the provided candidates list. Falls back
                       to CSV mode if MV is unavailable (warns on stderr).

    run_pipeline / run_multiphase receive the smaller list unchanged — the
    pipeline API is unaffected by the source of candidates.
    """
    import sys
    from peak_matcher import detect_peaks, match_candidates_csv, match_candidates_db

    tth, Iobs, sigma = parse_xye(xye_path)
    peaks_tth, _ = detect_peaks(tth, Iobs, sigma, prominence_sigma=prominence_sigma)

    if db_client is not None:
        try:
            cod_ids = match_candidates_db(
                peaks_tth, db_client._conn, wavelength=wavelength,
                tol_deg=tol_deg, min_matches=min_matches, top_n=top_n,
            )
            cand_map = {c.cod_id: c for c in candidates}
            filtered = [cand_map[cid] for cid in cod_ids if cid in cand_map]
            if filtered:
                return filtered[:top_n]
            # MV returned no overlap with provided candidates → fall through
        except Exception as e:
            print(f'[WARN] prefilter DB mode failed ({e}); falling back to CSV mode',
                  file=sys.stderr)

    return match_candidates_csv(
        peaks_tth, candidates, wavelength=wavelength,
        tol_deg=tol_deg, min_matches=min_matches,
    )[:top_n]


def candidates_from_db(
    xye_path: str | Path,
    db_client,
    top_n: int = 50,
    wavelength: float = 1.54056,
    tol_deg: float = 0.2,
    min_matches: int = 3,
    prominence_sigma: float = 5.0,
) -> list[CandidateInput]:
    """DB-only candidate discovery: detect peaks → MV search → fetch reflections [crit. 9.1-9.5].

    Replaces load_candidates_csv when no candidate CSV exists. Intrinsic Hanawalt
    filter (match_candidates_db) returns top_n cod_ids; fetch_reflections loads
    full CuKα reflection data from reference_patterns. The rest of the pipeline
    (run_pipeline / run_multiphase) is unchanged — interface is list[CandidateInput].
    """
    import sys
    from peak_matcher import detect_peaks, match_candidates_db

    tth, Iobs, sigma = parse_xye(xye_path)
    peaks_tth, _ = detect_peaks(tth, Iobs, sigma, prominence_sigma=prominence_sigma)

    if peaks_tth.size == 0:
        raise ValueError(
            'No peaks detected — check prominence_sigma or verify XYE file is valid.'
        )

    cod_ids = match_candidates_db(
        peaks_tth, db_client._conn,
        wavelength=wavelength, tol_deg=tol_deg,
        min_matches=min_matches, top_n=top_n,
    )

    if not cod_ids:
        raise ValueError(
            f'{peaks_tth.size} peaks detected but 0 COD phases matched '
            f'(min_matches={min_matches}). '
            'Check: MV peak_fingerprints built with d_hkl schema? '
            'Run migrations/create_peak_fingerprints.sql first.'
        )

    print(
        f'[DB] {peaks_tth.size} picos detectados → {len(cod_ids)} cod_ids '
        f'→ buscando reflexoes no reference_patterns',
        file=sys.stderr,
    )

    candidates = db_client.fetch_reflections(cod_ids)

    if not candidates:
        raise ValueError(
            'fetch_reflections returned empty list — '
            'nenhum padrao CuKα em reference_patterns para esses cod_ids (crit. 9.3/9.4).'
        )

    return candidates


def _main_multiphase(xye_path, csv_path, top_k):
    """T-025v viz: QPA table + pie + overlay. Requires DB for weight fractions."""
    import matplotlib.pyplot as plt
    from data_loader import load_candidates_csv, parse_xye
    from db_client import DBClient

    candidates = load_candidates_csv(csv_path)
    try:
        db = DBClient()
    except Exception as e:
        print(f"[WARN] DB indisponível ({e}) — wt% será 0")
        db = None
    try:
        res = run_multiphase(xye_path, candidates, db_client=db, top_k=top_k)
    finally:
        if db:
            db.close()

    print(f"\nMulti-fase QPA — {res.xye_file}  |  {res.n_points} pts  |  top-{top_k}\n")
    print(f"Rwp combinado={res.Rwp:.4f}  chi²={res.chi2:.2f}  "
          f"(Rwp melhor fase única={res.Rwp_single_best:.4f})\n")
    print(f"{'cod_id':>10} {'scale':>12} {'Z':>4} {'M':>10} {'V':>10} {'wt%':>8}  formula")
    print('-' * 80)
    for ph in res.phases:
        f = (ph.metadata.formula if ph.metadata else None) or '?'
        print(f"{str(ph.cod_id):>10} {ph.scale:>12.4e} {ph.Z:>4.0f} "
              f"{ph.M:>10.2f} {ph.V:>10.2f} {ph.weight_pct:>8.2f}  {f}")
    print(f"\nSoma wt% = {sum(ph.weight_pct for ph in res.phases):.2f}  (base cristalina)")

    # pie of weight fractions (only phases with wt% > 0)
    present = [ph for ph in res.phases if ph.weight_pct > 0.01]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))
    if present:
        labels = [f"{ph.cod_id}\n{ph.weight_pct:.1f}%" for ph in present]
        ax1.pie([ph.weight_pct for ph in present], labels=labels,
                autopct='%1.1f%%', startangle=90)
        ax1.set_title('Frações em peso (QPA Hill-Howard)')
    else:
        ax1.text(0.5, 0.5, 'Sem wt% (DB inativo?)', ha='center')

    tth, Iobs, sigma = parse_xye(xye_path)
    units = [build_icalc_unit_absolute(tth, next(c for c in candidates if c.cod_id == ph.cod_id).reflections)[0]
             for ph in res.phases]
    scales = [ph.scale for ph in res.phases]
    mu = tth.mean(); std = max(float(tth.std()), 1e-6)
    Icalc = sum(s * u for s, u in zip(scales, units))
    ax2.plot(tth, Iobs, 'k-', lw=0.7, label='Iobs')
    for ph, s, u in zip(res.phases, scales, units):
        if ph.weight_pct > 0.01:
            ax2.plot(tth, s * u, lw=0.6, alpha=0.7, label=f'{ph.cod_id} ({ph.weight_pct:.0f}%)')
    ax2.set_xlabel('2 theta (graus)'); ax2.set_ylabel('Intensidade')
    ax2.legend(fontsize=7); ax2.set_title(f'Decomposição  Rwp={res.Rwp:.4f}')
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    import sys
    import matplotlib.pyplot as plt
    from data_loader import load_candidates_csv

    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    flags = [a for a in sys.argv[1:] if a.startswith('--')]

    xye_path = args[0] if len(args) > 0 else 'synthetic_candidate17.xye'
    csv_path = args[1] if len(args) > 1 else 'data-1782394014136.csv'

    if '--mp' in flags:
        top_k = int(args[2]) if len(args) > 2 else 4
        _main_multiphase(xye_path, csv_path, top_k)
        sys.exit(0)

    candidates = load_candidates_csv(csv_path)
    result = run_pipeline(xye_path, candidates, db_client=None)

    print(f"\nRanking — {result.xye_file}  |  {result.n_points} pts  |  {len(result.candidates)} candidatos\n")
    print(f"{'#':>3}  {'cod_id':>10}  {'Rwp':>7}  {'Rp':>7}  {'chi2':>7}  {'scale':>10}  {'picos':>6}")
    print('-' * 65)
    viable = result.viable()
    for i, r in enumerate(result.candidates):
        flag = '★' if r in viable else ' '
        print(f"{flag}{i+1:>2}  {r.cod_id:>10}  {r.Rwp:>7.4f}  {r.Rp:>7.4f}  "
              f"{r.chi2:>7.2f}  {r.scale:>10.4f}  {r.n_peaks_used:>6}")
    print(f"\n★ {len(viable)} candidato(s) viável(is)  (Rwp<0.15, chi2<3)")

    rwps  = [r.Rwp  for r in result.candidates]
    chi2s = [r.chi2 for r in result.candidates]
    ids   = [str(r.cod_id) for r in result.candidates]

    plt.figure(figsize=(8, 6))
    plt.scatter(rwps, chi2s, c='steelblue', s=60, zorder=3)
    for x, y, label in zip(rwps, chi2s, ids):
        plt.annotate(label, (x, y), fontsize=7, xytext=(4, 4),
                     textcoords='offset points')
    plt.axvline(0.15, color='r', ls='--', lw=0.8, label='Rwp=0.15')
    plt.axhline(3.0,  color='g', ls='--', lw=0.8, label='chi²=3')
    plt.xlabel('Rwp')
    plt.ylabel('chi²')
    plt.title('Ranking de candidatos')
    plt.legend()
    plt.tight_layout()
    plt.show()
