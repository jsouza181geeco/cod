"""
Hanawalt-style XRD pre-selection (Stage 0) — d-spacing, radiation-agnostic.

Reduces O(500k COD candidates) -> O(50) before the Rietveld fit.

Matching is done in d-spacing (NOT 2theta): d_hkl is intrinsic to the crystal
(wavelength-independent), so one fingerprint set covers Cu/Co/Cr/Mo. This is the
canonical Hanawalt/ICDD method (index of (d, I) pairs) [H38; K74].

The sample's observed peaks (in 2theta at lambda_sample) are converted to d via
Bragg before matching. The angular tolerance is propagated to a PER-PEAK d window
via Bragg bounds — a fixed d tolerance would be wrong (Delta_d/d ~ cot(theta),
~25x wider at low angle).

Two matching modes
------------------
match_candidates_csv : in-memory, against a pre-loaded list  (CSV/offline/test)
match_candidates_db  : SQL against peak_fingerprints MV (d_hkl) (production)

Science review criteria 8.1-8.13 [H38; SV75; K74].
"""
from __future__ import annotations

import numpy as np
from scipy.signal import find_peaks

from models import CandidateInput

DEFAULT_WAVELENGTH = 1.54056  # CuKalpha1 (A)


# ---------------------------------------------------------------------------
# Peak detection (unchanged — returns 2theta positions)
# ---------------------------------------------------------------------------

def detect_peaks(
    tth: np.ndarray,
    Iobs: np.ndarray,
    sigma: np.ndarray,
    prominence_sigma: float = 5.0,
    min_distance_deg: float = 0.3,
) -> tuple[np.ndarray, np.ndarray]:
    """Detect experimental XRD peaks by prominence criterion [SV75].

    prominence >= prominence_sigma * median(sigma) adapts to local noise
    level and is robust to sloped background / amorphous humps (crit. 8.1).
    min_distance_deg prevents double-detection of a single broad peak (crit. 8.2).

    Returns
    -------
    tth_peaks : ndarray  2theta positions of detected peaks (degrees)
    I_peaks   : ndarray  intensities at those positions
    """
    if len(tth) < 3:
        return np.empty(0), np.empty(0)

    step = float(tth[1] - tth[0])
    min_dist = max(1, int(min_distance_deg / step))
    prom = prominence_sigma * float(np.median(sigma))

    idx, _ = find_peaks(Iobs, prominence=prom, distance=min_dist)
    return tth[idx], Iobs[idx]


# ---------------------------------------------------------------------------
# 2theta -> d conversion and per-peak d windows (crit. 8.9, 8.10)
# ---------------------------------------------------------------------------

def two_theta_to_d(two_theta_deg, wavelength: float) -> np.ndarray:
    """Bragg: d = lambda / (2 sin(theta)), theta = 2theta/2.  Vectorised."""
    tt = np.asarray(two_theta_deg, dtype=float)
    return wavelength / (2.0 * np.sin(np.radians(tt / 2.0)))


def d_windows(
    peaks_tth: np.ndarray,
    wavelength: float,
    tol_deg: float = 0.2,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-peak d match window from an angular tolerance (crit. 8.10).

    d shrinks as 2theta grows, so:
      d_lo = d(2theta + tol)   (upper angle -> smaller d)
      d_hi = d(2theta - tol)   (lower angle -> larger  d)

    The window width is non-uniform in d (Delta_d/d ~ cot(theta)) — this
    exactly preserves the instrument's angular tolerance, which a fixed
    d tolerance or fixed Delta_d/d would not.

    Returns (d_lo, d_hi), each shape (P,).
    """
    p = np.asarray(peaks_tth, dtype=float)
    d_lo = two_theta_to_d(p + tol_deg, wavelength)
    d_hi = two_theta_to_d(np.maximum(p - tol_deg, 1e-6), wavelength)
    return d_lo, d_hi


# ---------------------------------------------------------------------------
# CSV mode (in-memory)
# ---------------------------------------------------------------------------

def match_candidates_csv(
    peaks_tth: np.ndarray,
    candidates: list[CandidateInput],
    wavelength: float = DEFAULT_WAVELENGTH,
    tol_deg: float = 0.2,
    min_matches: int = 3,
) -> list[CandidateInput]:
    """In-memory Hanawalt matching [H38] in d-space against a candidate list.

    Observed 2theta peaks are converted to per-peak d windows via lambda
    (crit. 8.9) then matched against each candidate reflection's 'd_hkl'
    (crit. 8.8). Discards candidates with n_matches < min_matches (crit. 8.6).
    Updates candidate.peak_matches in place.

    Returns candidates sorted by peak_matches DESC.
    """
    peaks = np.asarray(peaks_tth, dtype=float)
    if peaks.size == 0:
        return []

    d_lo, d_hi = d_windows(peaks, wavelength, tol_deg)   # (P,)

    scored: list[tuple[int, CandidateInput]] = []
    for c in candidates:
        ref_d = np.array(
            [r['d_hkl'] for r in c.reflections
             if r.get('intensity_rel', 0) > 0 and r.get('d_hkl')],
            dtype=float,
        )
        if ref_d.size == 0:
            continue
        # obs peak i matches if any ref_d lies in [d_lo[i], d_hi[i]]
        in_win = (ref_d[None, :] >= d_lo[:, None]) & (ref_d[None, :] <= d_hi[:, None])  # (P,R)
        n_matched = int(np.sum(np.any(in_win, axis=1)))
        if n_matched >= min_matches:
            c.peak_matches = n_matched
            scored.append((n_matched, c))

    scored.sort(key=lambda x: -x[0])
    return [c for _, c in scored]


# ---------------------------------------------------------------------------
# DB mode (peak_fingerprints MV, d_hkl)
# ---------------------------------------------------------------------------

def match_candidates_db(
    peaks_tth: np.ndarray,
    conn,
    wavelength: float = DEFAULT_WAVELENGTH,
    tol_deg: float = 0.2,
    min_matches: int = 3,
    top_n: int = 50,
) -> list[int]:
    """Hanawalt search via PostgreSQL peak_fingerprints MV (d_hkl column).

    Per-peak d windows (crit. 8.10) are joined against pf.d_hkl with BETWEEN
    so the B-tree index on d_hkl is used (crit. 8.5). COUNT(DISTINCT obs.oid)
    counts how many observed peaks matched >=1 reference line (consistent with
    CSV mode).

    Returns list[int] cod_ids sorted by n_matched DESC.

    Requires: xrd_analysis.peak_fingerprints MV with d_hkl
              (see migrations/create_peak_fingerprints.sql).
    """
    peaks = np.asarray(peaks_tth, dtype=float)
    if peaks.size == 0:
        return []

    d_lo, d_hi = d_windows(peaks, wavelength, tol_deg)

    with conn.cursor() as cur:
        cur.execute(
            """
            WITH obs AS (
                SELECT dlo, dhi, oid
                FROM unnest(%s::float8[], %s::float8[])
                     WITH ORDINALITY AS t(dlo, dhi, oid)
            )
            SELECT pf.cod_id, COUNT(DISTINCT obs.oid) AS n_matched
            FROM   xrd_analysis.peak_fingerprints pf
            JOIN   obs ON pf.d_hkl BETWEEN obs.dlo AND obs.dhi
            GROUP  BY pf.cod_id
            HAVING COUNT(DISTINCT obs.oid) >= %s
            ORDER  BY n_matched DESC
            LIMIT  %s
            """,
            (d_lo.tolist(), d_hi.tolist(), min_matches, top_n),
        )
        return [int(row[0]) for row in cur.fetchall()]


if __name__ == '__main__':
    # T-028v — visualisation: detected peaks + per-candidate match table (d-space)
    import sys
    import matplotlib.pyplot as plt
    from data_loader import parse_xye, load_candidates_csv

    xye_path   = sys.argv[1] if len(sys.argv) > 1 else 'synthetic_candidate17.xye'
    csv_path   = sys.argv[2] if len(sys.argv) > 2 else 'data-1782394014136.csv'
    prom_sigma = float(sys.argv[3]) if len(sys.argv) > 3 else 5.0
    wl         = float(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_WAVELENGTH
    tol        = 0.2

    tth, Iobs, sigma = parse_xye(xye_path)
    candidates = load_candidates_csv(csv_path)

    tth_peaks, I_peaks = detect_peaks(tth, Iobs, sigma, prominence_sigma=prom_sigma)
    d_lo, d_hi = d_windows(tth_peaks, wl, tol)

    print(f'Arquivo   : {xye_path}')
    print(f'Pontos    : {len(tth)}')
    print(f'Lambda    : {wl} A')
    print(f'Picos det.: {len(tth_peaks)}  '
          f'(prominence >= {prom_sigma:.1f} x median(sigma))')
    if len(tth_peaks):
        d_peaks = two_theta_to_d(tth_peaks, wl)
        print(f'Range 2th : {tth_peaks[0]:.2f}° -- {tth_peaks[-1]:.2f}°')
        print(f'Range d   : {d_peaks.min():.3f} -- {d_peaks.max():.3f} A')

    # n_matches (d-space) para todos os candidatos, min_matches=0
    all_n: dict[int, int] = {}
    for c in candidates:
        ref_d = np.array(
            [r['d_hkl'] for r in c.reflections
             if r.get('intensity_rel', 0) > 0 and r.get('d_hkl')],
            dtype=float,
        )
        if ref_d.size == 0 or tth_peaks.size == 0:
            all_n[c.cod_id] = 0
            continue
        in_win = (ref_d[None, :] >= d_lo[:, None]) & (ref_d[None, :] <= d_hi[:, None])
        all_n[c.cod_id] = int(np.sum(np.any(in_win, axis=1)))

    top_match = max(all_n.values()) if all_n else 0
    print(f'\nMatching d-space (tol={tol}°, {len(tth_peaks)} picos):')
    print(f"{'cod_id':>10}  {'n_matches':>10}")
    print('-' * 24)
    for cid, nm in sorted(all_n.items(), key=lambda x: -x[1]):
        flag = '  <-- melhor' if nm == top_match and nm > 0 else ''
        print(f'{cid:>10}  {nm:>10}{flag}')

    # Plot
    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(13, 8), sharex=True,
        gridspec_kw={'height_ratios': [4, 1]},
    )
    ax1.plot(tth, Iobs, 'k-', lw=0.7, label='Iobs')
    if tth_peaks.size:
        ax1.scatter(tth_peaks, I_peaks, color='red', s=25, zorder=5,
                    label=f'{len(tth_peaks)} picos detectados')
    ax1.set_ylabel('Intensidade')
    ax1.set_title(f'Deteccao de picos — {xye_path}  (d-space match, lambda={wl}A)')
    ax1.legend()

    top5 = sorted(all_n.items(), key=lambda x: -x[1])[:5]
    colours = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    cand_map = {c.cod_id: c for c in candidates}
    for (cid, nm), col in zip(top5, colours):
        if cid not in cand_map:
            continue
        ref_tths = [
            r['two_theta'] for r in cand_map[cid].reflections
            if r.get('intensity_rel', 0) > 0 and tth[0] <= r['two_theta'] <= tth[-1]
        ]
        for rt in ref_tths:
            ax2.axvline(rt, color=col, alpha=0.45, lw=0.7)
        ax2.plot([], [], color=col, lw=1.5, label=f'{cid} ({nm} matches)')
    ax2.set_xlabel('2 theta (graus)')
    ax2.set_yticks([])
    ax2.legend(fontsize=7, ncol=min(5, len(top5)))
    ax2.set_title('Posicoes Bragg dos top-5 candidatos (linhas coloridas)')

    plt.tight_layout()
    plt.show()
