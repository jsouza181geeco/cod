import numpy as np
import scipy.linalg


def linear_fit(
    tth: np.ndarray,
    Iobs: np.ndarray,
    sigma: np.ndarray,
    Icalc_unit: np.ndarray,
    n_bg: int = 4,
) -> tuple[float, np.ndarray]:
    """
    Weighted least-squares: minimise Σ wᵢ(Iobs - S·Icalc_unit - Ibg)²  [Y93 p.18-22]

    Design matrix A = [Icalc_unit | 1 | tth_norm | tth_norm² | ...]
    Weights wᵢ = 1/σᵢ²  (not 1/σ — see [Y93] eq. 2.3)
    Background polynomial normalised by tth mean/std for numerical stability.

    Returns (scale S, Icalc = S·Icalc_unit + Ibg)
    """
    w = 1.0 / np.maximum(sigma**2, 1e-10)
    sqrt_w = np.sqrt(w)

    mu = tth.mean()
    std = max(float(tth.std()), 1e-6)
    tth_norm = (tth - mu) / std

    A_cols = [Icalc_unit] + [tth_norm**k for k in range(n_bg)]
    A = np.column_stack(A_cols)

    Aw = A * sqrt_w[:, None]
    bw = Iobs * sqrt_w

    x, _, _, _ = scipy.linalg.lstsq(Aw, bw)
    scale = float(x[0])
    bg_coeffs = x[1:]

    Ibg = sum(bg_coeffs[k] * tth_norm**k for k in range(n_bg))
    Icalc = scale * Icalc_unit + Ibg

    return scale, Icalc


if __name__ == '__main__':
    import sys
    import matplotlib.pyplot as plt
    from data_loader import parse_xye, load_candidates_csv
    from pattern_calc import build_icalc_unit

    xye_path = sys.argv[1] if len(sys.argv) > 1 else 'synthetic_candidate17.xye'
    csv_path = sys.argv[2] if len(sys.argv) > 2 else 'data-1782394014136.csv'
    idx      = int(sys.argv[3]) if len(sys.argv) > 3 else 0

    tth, Iobs, sigma = parse_xye(xye_path)
    candidates = load_candidates_csv(csv_path)
    cand = candidates[idx]

    Icalc_unit, n_used = build_icalc_unit(tth, cand.reflections)
    scale, Icalc = linear_fit(tth, Iobs, sigma, Icalc_unit)
    diff = Iobs - Icalc

    print(f"Candidato : cod_id={cand.cod_id}")
    print(f"Scale     : {scale:.6f}")
    print(f"|diff| max: {np.abs(diff).max():.1f}  mean: {np.abs(diff).mean():.1f}")

    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True,
                             gridspec_kw={'height_ratios': [3, 3, 1]})

    axes[0].plot(tth, Iobs,  'k-', lw=0.8, label='Iobs')
    axes[0].plot(tth, Icalc, 'r-', lw=0.8, label='Icalc', alpha=0.8)
    axes[0].set_ylabel('Intensidade')
    axes[0].legend()
    axes[0].set_title(f'cod_id={cand.cod_id}  |  scale={scale:.4f}  |  picos={n_used}')

    baseline_max = np.percentile(Iobs, 90)
    axes[1].plot(tth, Iobs,  'k-', lw=0.8, label='Iobs')
    axes[1].plot(tth, Icalc, 'r-', lw=0.8, label='Icalc', alpha=0.8)
    axes[1].set_ylim(0, baseline_max * 1.5)
    axes[1].set_ylabel('Intensidade (zoom)')
    axes[1].legend()
    axes[1].set_title('Zoom — picos secundários')

    axes[2].plot(tth, diff, 'g-', lw=0.5)
    axes[2].axhline(0, color='k', lw=0.5)
    axes[2].set_xlabel('2θ (graus)')
    axes[2].set_ylabel('Δ')
    axes[2].set_title('Resíduo (Iobs − Icalc)')

    plt.tight_layout()
    plt.show()
