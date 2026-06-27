"""
Multi-phase weighted least squares (Épico 11) [Y93 cap.5; BH88].

Fits K phase scales + a shared background polynomial to the observed
diffractogram. Phase scales are constrained S_k >= 0 (a negative phase
fraction in a mixture is unphysical [crit. 6.2]); background coefficients
are unconstrained [crit. 6.3].

    I_calc(2θ) = Σ_k S_k · I_calc,unit,k(2θ) + I_bg(2θ)
"""
import numpy as np
from scipy.optimize import lsq_linear


def multi_phase_fit(
    tth: np.ndarray,
    Iobs: np.ndarray,
    sigma: np.ndarray,
    Icalc_units: list[np.ndarray],
    n_bg: int = 4,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    WLS fit of K phases + background, with non-negative phase scales.

    Icalc_units : list[np.ndarray (N,)] — one unit pattern per phase
                  (from build_icalc_unit_absolute for QPA; see crit. 7.9).
    Weights wᵢ = 1/σᵢ² applied via Aw = A·√w, bw = y·√w (same as linear_fit).
    Background 2θ normalised by mean/std for numerical stability.

    Bounds [crit. 6.2/6.3]:
        phase columns      → [0, +∞)   (fração de fase não-negativa)
        background columns  → (-∞, +∞)  (background subtraível)

    Returns (scales (K,), Icalc (N,), bg_coeffs (n_bg,)).
    """
    K = len(Icalc_units)
    if K == 0:
        raise ValueError("Icalc_units is empty — need at least one phase")

    w = 1.0 / np.maximum(sigma**2, 1e-10)
    sqrt_w = np.sqrt(w)

    mu = tth.mean()
    std = max(float(tth.std()), 1e-6)
    tth_norm = (tth - mu) / std

    phase_cols = list(Icalc_units)                       # K columns
    bg_cols = [tth_norm**k for k in range(n_bg)]         # n_bg columns
    A = np.column_stack(phase_cols + bg_cols)            # N x (K + n_bg)

    Aw = A * sqrt_w[:, None]
    bw = Iobs * sqrt_w

    lb = np.array([0.0] * K + [-np.inf] * n_bg)
    ub = np.full(K + n_bg, np.inf)
    res = lsq_linear(Aw, bw, bounds=(lb, ub), method='bvls')

    scales = res.x[:K]
    bg_coeffs = res.x[K:]

    Ibg = sum(bg_coeffs[k] * tth_norm**k for k in range(n_bg))
    Icalc = sum(scales[k] * Icalc_units[k] for k in range(K)) + Ibg

    return scales, Icalc, bg_coeffs


if __name__ == '__main__':
    import sys
    import matplotlib.pyplot as plt
    from data_loader import parse_xye, load_candidates_csv
    from pattern_calc import build_icalc_unit_absolute
    from fom import calc_fom

    xye_path = sys.argv[1] if len(sys.argv) > 1 else 'synthetic_candidate17.xye'
    csv_path = sys.argv[2] if len(sys.argv) > 2 else 'data-1782394014136.csv'
    # phase indices to mix (default: candidate 17 = idx 17, plus two others)
    idxs = [int(x) for x in sys.argv[3:]] if len(sys.argv) > 3 else [17, 0, 5]

    tth, Iobs, sigma = parse_xye(xye_path)
    candidates = load_candidates_csv(csv_path)
    phases = [candidates[i] for i in idxs]

    units = [build_icalc_unit_absolute(tth, c.reflections)[0] for c in phases]
    scales, Icalc, bg = multi_phase_fit(tth, Iobs, sigma, units, n_bg=4)
    fom = calc_fom(Iobs, Icalc, sigma, n_params=len(units) + 4)

    print(f"Fases: {[c.cod_id for c in phases]}")
    print(f"Scales: {np.array2string(scales, precision=4)}")
    print(f"Rwp combinado: {fom['Rwp']:.4f}  chi2: {fom['chi2']:.3f}")
    print(f"Scales >= 0? {(scales >= 0).all()}")

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True,
                             gridspec_kw={'height_ratios': [3, 1]})
    axes[0].plot(tth, Iobs, 'k-', lw=0.8, label='Iobs')
    axes[0].plot(tth, Icalc, 'r-', lw=0.8, label='Icalc total', alpha=0.8)
    mu = tth.mean(); std = max(float(tth.std()), 1e-6); tth_norm = (tth - mu) / std
    Ibg = sum(bg[k] * tth_norm**k for k in range(4))
    for c, s, u in zip(phases, scales, units):
        axes[0].plot(tth, s * u + Ibg, lw=0.6, alpha=0.6, label=f'cod={c.cod_id} (S={s:.3e})')
    axes[0].set_ylabel('Intensidade'); axes[0].legend(fontsize=7)
    axes[0].set_title(f'Multi-fase  Rwp={fom["Rwp"]:.4f}  chi2={fom["chi2"]:.2f}')
    axes[1].plot(tth, Iobs - Icalc, 'g-', lw=0.5)
    axes[1].axhline(0, color='k', lw=0.5)
    axes[1].set_xlabel('2 theta (graus)'); axes[1].set_ylabel('residuo')
    plt.tight_layout()
    plt.show()
