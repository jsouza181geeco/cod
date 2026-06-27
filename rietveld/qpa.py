"""
Quantitative Phase Analysis (QPA) — Hill-Howard weight fractions [HH87; BH88].

    W_k = S_k·(ZMV)_k / Σ_j S_j·(ZMV)_j

Z = formula units per cell, M = molar mass of the formula unit, V = unit
cell volume. The scale factor S_k MUST come from a fit on an ABSOLUTE
intensity basis (build_icalc_unit_absolute), not on per-phase-normalised
`intensity_rel` [crit. 7.6/7.9], else the fractions are biased by each
phase's normalisation constant.

Limitations (screening, not certification):
- Amorphous content NOT quantified (no internal/external standard) → W_k
  are fractions of the CRYSTALLINE portion only. [BH88]
- Preferred orientation (March-Dollase) not corrected. [D86]
- Microabsorption (Brindley) not corrected. [B45]
"""
from crystallo_utils import cell_volume, molar_mass


def weight_fractions(scales, metadatas) -> list[dict]:
    """
    Hill-Howard weight fractions from fitted phase scales + DB metadata.

    scales    : sequence of K phase scale factors (from multi_phase_fit on
                the absolute intensity basis).
    metadatas : sequence of K StructureMetadata (same order). A phase whose
                metadata is missing Z, cell params or a parseable formula
                contributes ZMV=0 (excluded from QPA, weight 0%) [crit. 7.5].

    Returns list[dict] (same order) with keys:
        cod_id, scale, Z, M, V, ZMV, weight_pct

    Σ weight_pct = 100.0 (crystalline basis [crit. 7.4]) unless every ZMV=0.
    """
    rows = []
    for s, m in zip(scales, metadatas):
        s = float(s)
        Z = M = V = zmv = 0.0
        cod_id = None
        if m is not None:
            cod_id = m.cod_id
            ok = (m.Z and m.a and m.b and m.c
                  and m.alpha and m.beta and m.gamma and m.formula)
            if ok:
                try:
                    V = cell_volume(m.a, m.b, m.c, m.alpha, m.beta, m.gamma)
                    M = molar_mass(m.formula)
                    Z = float(m.Z)
                    zmv = s * Z * M * V
                except ValueError:
                    # unknown element / bad formula → exclude this phase
                    zmv = 0.0
        # negative scale (shouldn't happen with bounded fit) → no mass
        if zmv < 0.0:
            zmv = 0.0
        rows.append({'cod_id': cod_id, 'scale': s,
                     'Z': Z, 'M': M, 'V': V, 'ZMV': zmv, 'weight_pct': 0.0})

    total = sum(r['ZMV'] for r in rows)
    if total > 0.0:
        for r in rows:
            r['weight_pct'] = 100.0 * r['ZMV'] / total
    return rows


if __name__ == '__main__':
    import sys
    from data_loader import parse_xye, load_candidates_csv
    from db_client import DBClient
    from pattern_calc import build_icalc_unit_absolute
    from multi_phase_fit import multi_phase_fit
    from fom import calc_fom

    xye_path = sys.argv[1] if len(sys.argv) > 1 else 'synthetic_candidate17.xye'
    csv_path = sys.argv[2] if len(sys.argv) > 2 else 'data-1782394014136.csv'
    idxs = [int(x) for x in sys.argv[3:]] if len(sys.argv) > 3 else [17, 0, 5]

    tth, Iobs, sigma = parse_xye(xye_path)
    candidates = load_candidates_csv(csv_path)
    phases = [candidates[i] for i in idxs]

    units = [build_icalc_unit_absolute(tth, c.reflections)[0] for c in phases]
    scales, Icalc, _ = multi_phase_fit(tth, Iobs, sigma, units, n_bg=4)
    fom = calc_fom(Iobs, Icalc, sigma, n_params=len(units) + 4)

    with DBClient() as db:
        meta_map = db.fetch_metadata([c.cod_id for c in phases])
    metadatas = [meta_map.get(c.cod_id) for c in phases]

    qpa = weight_fractions(scales, metadatas)

    print(f"\nQPA Hill-Howard — Rwp={fom['Rwp']:.4f}  chi2={fom['chi2']:.2f}\n")
    print(f"{'cod_id':>10} {'scale':>12} {'Z':>4} {'M':>10} {'V':>10} {'wt%':>8}  formula")
    print('-' * 78)
    for r, m in zip(qpa, metadatas):
        formula = (m.formula if m else None) or '?'
        print(f"{str(r['cod_id']):>10} {r['scale']:>12.4e} {r['Z']:>4.0f} "
              f"{r['M']:>10.2f} {r['V']:>10.2f} {r['weight_pct']:>8.2f}  {formula}")
    print(f"\nSoma wt% = {sum(r['weight_pct'] for r in qpa):.2f}  (base 100% cristalina)")
