import numpy as np


def calc_fom(
    Iobs: np.ndarray,
    Icalc: np.ndarray,
    sigma: np.ndarray,
    n_params: int,
) -> dict[str, float]:
    """
    Weighted profile R-factors [Y93 p.22-26 eq.2.5-2.7; T06].

    Rwp  = sqrt( Σ wᵢ(Iobs-Icalc)² / Σ wᵢ·Iobs² )
    Rp   = Σ|Iobs-Icalc| / Σ|Iobs|
    Rexp = sqrt( (N-P) / Σ wᵢ·Iobs² )
    χ²   = (Rwp/Rexp)²   [ideal ≈ 1.0]

    wᵢ = 1/σᵢ²  (NOT 1/σ).  Denominator floors prevent division by zero.
    """
    w = 1.0 / np.maximum(sigma**2, 1e-10)
    diff = Iobs - Icalc

    sum_w_diff2 = np.sum(w * diff**2)
    sum_w_Iobs2 = np.sum(w * Iobs**2)

    Rwp  = float(np.sqrt(sum_w_diff2 / max(sum_w_Iobs2, 1e-20)))
    Rp   = float(np.sum(np.abs(diff)) / max(np.sum(np.abs(Iobs)), 1e-20))
    N    = len(Iobs)
    Rexp = float(np.sqrt(max(N - n_params, 1) / max(sum_w_Iobs2, 1e-20)))
    chi2 = float((Rwp / Rexp)**2) if Rexp > 0 else float('inf')

    return {'Rwp': Rwp, 'Rp': Rp, 'Rexp': Rexp, 'chi2': chi2}


if __name__ == '__main__':
    import sys
    from data_loader import parse_xye, load_candidates_csv
    from pattern_calc import build_icalc_unit
    from linear_fit import linear_fit

    xye_path = sys.argv[1] if len(sys.argv) > 1 else 'synthetic_candidate17.xye'
    csv_path = sys.argv[2] if len(sys.argv) > 2 else 'data-1782394014136.csv'
    idx      = int(sys.argv[3]) if len(sys.argv) > 3 else 0

    tth, Iobs, sigma = parse_xye(xye_path)
    candidates = load_candidates_csv(csv_path)
    cand = candidates[idx]

    Icalc_unit, _ = build_icalc_unit(tth, cand.reflections)
    scale, Icalc  = linear_fit(tth, Iobs, sigma, Icalc_unit, n_bg=4)
    fom = calc_fom(Iobs, Icalc, sigma, n_params=5)

    def tag(rwp, chi2):
        if rwp < 0.10 and chi2 < 2:  return '✓ BOM'
        if rwp < 0.20 and chi2 < 5:  return '⚠ MÉDIO'
        return '✗ RUIM'

    print(f"\nFigures of Merit — cod_id={cand.cod_id}")
    print(f"  Rwp  = {fom['Rwp']:.4f}")
    print(f"  Rp   = {fom['Rp']:.4f}")
    print(f"  Rexp = {fom['Rexp']:.4f}")
    print(f"  chi2 = {fom['chi2']:.3f}")
    print(f"  {tag(fom['Rwp'], fom['chi2'])}")
    print(f"  scale= {scale:.6f}")

    fom_perfeito = calc_fom(Iobs, Iobs,               sigma, n_params=5)
    fom_nulo     = calc_fom(Iobs, np.zeros_like(Iobs), sigma, n_params=5)
    print(f"\n--- Sanidade ---")
    print(f"Icalc=Iobs  → Rwp={fom_perfeito['Rwp']:.2e}  (esperado ≈ 0)")
    print(f"Icalc=zeros → Rwp={fom_nulo['Rwp']:.4f}      (esperado ≈ 1)")
