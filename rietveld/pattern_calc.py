import numpy as np


def pseudo_voigt_profile(
    tth_grid: np.ndarray,
    two_theta_peak: float,
    fwhm: float,
    eta: float,
) -> np.ndarray:
    """
    Pseudo-Voigt peak profile normalised to 1.0 at centre [TCH87 eq.1-3].

    pV(Δ) = η·L(Δ) + (1-η)·G(Δ),  η=1 → pure Lorentz, η=0 → pure Gauss.
    G and L are peak-normalised (value 1.0 at Δ=0), not area-normalised.
    """
    delta = tth_grid - two_theta_peak
    sigma_g = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    G = np.exp(-delta**2 / (2.0 * sigma_g**2))
    L = 1.0 / (1.0 + (delta / (fwhm / 2.0))**2)
    return eta * L + (1.0 - eta) * G


def caglioti_fwhm(
    two_theta_deg: float | np.ndarray,
    U: float,
    V: float,
    W: float,
) -> float | np.ndarray:
    """
    Caglioti et al. (1958) angular dependence of FWHM.

    Γ²(2θ) = U·tan²θ + V·tanθ + W,  θ = 2θ/2  [CPR58 eq. 2]

    Returns FWHM in degrees (same units as two_theta_deg).
    Floor at 1e-8 prevents imaginary FWHM when V very negative.
    """
    theta_rad = np.radians(np.asarray(two_theta_deg) / 2.0)
    tan_t = np.tan(theta_rad)
    fwhm2 = U * tan_t**2 + V * tan_t + W
    return np.sqrt(np.maximum(fwhm2, 1e-8))


def lorentz_polarization(
    two_theta_deg: float | np.ndarray,
    two_theta_mono_deg: float | None = None,
) -> float | np.ndarray:
    """
    Lorentz-polarization factor for a Bragg-Brentano powder diffractometer.

    No monochromator:
        Lp(θ) = (1 + cos²2θ) / (sin²θ · cosθ)

    With a diffracted-beam monochromator (e.g. graphite, 2θ_m):
        Lp(θ) = (1 + cos²2θ_m · cos²2θ) / (sin²θ · cosθ)

    θ = 2θ/2 (Bragg angle). Lp diverges as 2θ → 0 (the 1/sin²θ Lorentz
    term), so it dominates the angular variation of the calculated
    intensity at low angle — this is why CSV `intensity_rel` is NOT a
    cross-phase-comparable basis for QPA [crit. 7.6/7.8].
    """
    tt = np.radians(np.asarray(two_theta_deg, dtype=np.float64))
    th = tt / 2.0
    if two_theta_mono_deg is None:
        num = 1.0 + np.cos(tt) ** 2
    else:
        cm2 = np.cos(np.radians(two_theta_mono_deg)) ** 2
        num = 1.0 + cm2 * np.cos(tt) ** 2
    denom = np.maximum(np.sin(th) ** 2 * np.cos(th), 1e-12)
    return num / denom


def build_icalc_unit(
    tth: np.ndarray,
    reflections: list[dict],
    U: float = 0.01,
    V: float = -0.002,
    W: float = 0.005,
    eta: float = 0.5,
    cutoff_fwhm: float = 10.0,
) -> tuple[np.ndarray, int]:
    """
    Icalc_unit(2θ_i) = Σ_hkl intensity_rel_hkl · pV(2θ_i − 2θ_hkl, Γ_hkl)  [R69 eq.1]

    intensity_rel already includes Lp · DW · M · |F|² (pre-computed in CSV).
    Scale factor S applied externally by linear_fit.

    Returns (Icalc_unit, n_peaks_used).
    """
    Icalc = np.zeros(len(tth), dtype=np.float64)
    n_used = 0
    tth_min, tth_max = tth[0], tth[-1]

    for refl in reflections:
        two_theta_peak = refl['two_theta']
        intensity_rel = refl.get('intensity_rel') or 0.0
        if intensity_rel <= 0.0:
            continue
        fwhm = caglioti_fwhm(two_theta_peak, U, V, W)
        # skip peaks whose entire contribution falls outside the grid
        if two_theta_peak + cutoff_fwhm * fwhm < tth_min:
            continue
        if two_theta_peak - cutoff_fwhm * fwhm > tth_max:
            continue
        profile = pseudo_voigt_profile(tth, two_theta_peak, fwhm, eta)
        Icalc += intensity_rel * profile
        n_used += 1

    return Icalc, n_used


def build_icalc_unit_absolute(
    tth: np.ndarray,
    reflections: list[dict],
    U: float = 0.01,
    V: float = -0.002,
    W: float = 0.005,
    eta: float = 0.5,
    cutoff_fwhm: float = 10.0,
    two_theta_mono_deg: float | None = None,
) -> tuple[np.ndarray, int]:
    """
    Like build_icalc_unit, but weights each peak by the ABSOLUTE calculated
    intensity instead of the per-phase-normalised `intensity_rel`:

        I_abs_hkl = multiplicity · F_sq · Lp(θ_hkl)

    Required for QPA (Épico 11): `intensity_rel` is normalised per phase
    (every phase max = 100, verified), so its scale factor is NOT
    cross-phase comparable for Hill-Howard [crit. 7.6/7.9]. Fitting on this
    absolute basis makes the fitted S_k the true Rietveld scale.

    Debye-Waller (DW) is omitted — it needs U_iso from atomic_sites and is a
    minor, angle-growing, phase-dependent correction (screening only).

    Returns (Icalc_unit_absolute, n_peaks_used).
    """
    Icalc = np.zeros(len(tth), dtype=np.float64)
    n_used = 0
    tth_min, tth_max = tth[0], tth[-1]

    for refl in reflections:
        two_theta_peak = refl['two_theta']
        mult = refl.get('multiplicity') or 0
        f_sq = refl.get('F_sq') or 0.0
        if mult <= 0 or f_sq <= 0.0:
            continue
        lp = float(lorentz_polarization(two_theta_peak, two_theta_mono_deg))
        intensity_abs = mult * f_sq * lp
        if intensity_abs <= 0.0:
            continue
        fwhm = caglioti_fwhm(two_theta_peak, U, V, W)
        if two_theta_peak + cutoff_fwhm * fwhm < tth_min:
            continue
        if two_theta_peak - cutoff_fwhm * fwhm > tth_max:
            continue
        profile = pseudo_voigt_profile(tth, two_theta_peak, fwhm, eta)
        Icalc += intensity_abs * profile
        n_used += 1

    return Icalc, n_used


if __name__ == '__main__':
    import sys
    import matplotlib.pyplot as plt
    from data_loader import parse_xye, load_candidates_csv

    xye_path = sys.argv[1] if len(sys.argv) > 1 else 'Cu_synthetic.xye'
    csv_path = sys.argv[2] if len(sys.argv) > 2 else 'data-1782394014136.csv'
    idx      = int(sys.argv[3]) if len(sys.argv) > 3 else 0

    tth, Iobs, sigma = parse_xye(xye_path)
    candidates = load_candidates_csv(csv_path)
    cand = candidates[idx]

    Icalc_unit, n_used = build_icalc_unit(tth, cand.reflections)
    scale_approx = Iobs.max() / max(Icalc_unit.max(), 1e-10)
    Icalc_scaled = Icalc_unit * scale_approx

    peak_positions = [r['two_theta'] for r in cand.reflections
                      if r.get('intensity_rel', 0) > 0]

    print(f"Candidato : cod_id={cand.cod_id}")
    print(f"Picos     : {n_used}/{len(cand.reflections)} usados")
    print(f"Icalc max : {Icalc_unit.max():.2f}  em 2θ = {tth[Icalc_unit.argmax()]:.3f}°")

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

    axes[0].plot(tth, Iobs, 'k-', lw=0.8, label='Iobs (exp)')
    axes[0].set_ylabel('Intensidade')
    axes[0].legend()
    axes[0].set_title('Experimental')

    axes[1].plot(tth, Icalc_scaled, 'b-', lw=0.8,
                 label=f'Icalc (cod={cand.cod_id}, escala aprox.)')
    for pp in peak_positions:
        if tth[0] <= pp <= tth[-1]:
            axes[1].axvline(pp, color='r', alpha=0.3, lw=0.5)
    axes[1].set_xlabel('2θ (graus)')
    axes[1].set_ylabel('Intensidade')
    axes[1].legend()
    axes[1].set_title('Padrão calculado  (linhas vermelhas = posições Bragg)')

    plt.tight_layout()
    plt.show()
