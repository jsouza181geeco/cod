"""
Tests for pattern_calc — science review criteria 2.1–2.13 [CPR58; TCH87; R69].
"""
import numpy as np
import pytest

from pattern_calc import caglioti_fwhm, pseudo_voigt_profile, build_icalc_unit


TTH = np.linspace(20.0, 120.0, 5001)
U, V, W, ETA = 0.01, -0.002, 0.005, 0.5


# ============================================================
# caglioti_fwhm  — criteria 2.1–2.4
# ============================================================

def test_fwhm_positive_scalar():
    """2.2: FWHM > 0 for a single 2theta value."""
    fwhm = caglioti_fwhm(43.32, U, V, W)
    assert fwhm > 0


def test_fwhm_positive_array():
    """2.2+2.4: FWHM > 0 for all 2theta in [5°, 150°]."""
    tth_test = np.linspace(5.0, 150.0, 2000)
    fwhm = caglioti_fwhm(tth_test, U, V, W)
    assert (fwhm > 0).all(), "FWHM must be positive across full range [2.2]"


def test_fwhm_uses_half_angle():
    """2.1: theta = two_theta/2. Check FWHM(43°) ≈ 0.07–0.10° (lab typical)."""
    fwhm = float(caglioti_fwhm(43.32, U, V, W))
    assert 0.03 < fwhm < 0.30, f"FWHM={fwhm:.4f} out of typical lab range (suggests /2 error)"


def test_fwhm_units_degrees():
    """2.3: FWHM in same units as 2theta (degrees)."""
    fwhm = float(caglioti_fwhm(90.0, U, V, W))
    # if theta not halved, fwhm(90°) would use tan(90°) → huge or NaN
    assert np.isfinite(fwhm) and fwhm < 2.0


def test_fwhm_monotone_above_30():
    """2.4: FWHM monotonically increasing for 2theta > 30° (typical U>0, W>0)."""
    tth_test = np.linspace(31.0, 110.0, 100)
    fwhm = caglioti_fwhm(tth_test, U, V, W)
    diff = np.diff(fwhm)
    # allow a tiny flat region but no decrease
    assert (diff >= -1e-6).all(), "FWHM not monotone above 30° [2.4]"


# ============================================================
# pseudo_voigt_profile  — criteria 2.5–2.8
# ============================================================

def test_pv_peak_is_one():
    """2.5: pV(Δ=0) = 1.0."""
    peak_pos = 43.32
    val = pseudo_voigt_profile(np.array([peak_pos]), peak_pos, fwhm=0.1, eta=ETA)[0]
    assert abs(val - 1.0) < 1e-10


def test_pv_symmetric():
    """2.6: pV(-Δ) = pV(Δ)."""
    peak_pos = 50.0
    fwhm = 0.15
    deltas = np.array([0.05, 0.10, 0.20, 0.50])
    pos_vals = pseudo_voigt_profile(peak_pos + deltas, peak_pos, fwhm, ETA)
    neg_vals = pseudo_voigt_profile(peak_pos - deltas, peak_pos, fwhm, ETA)
    np.testing.assert_allclose(pos_vals, neg_vals, atol=1e-12)


def test_pv_decays_from_peak():
    """Profile decays monotonically from peak."""
    peak_pos = 43.0
    fwhm = 0.1
    tth_test = np.array([peak_pos - 0.2, peak_pos - 0.1, peak_pos,
                          peak_pos + 0.1, peak_pos + 0.2])
    vals = pseudo_voigt_profile(tth_test, peak_pos, fwhm, ETA)
    assert vals[2] > vals[1] > vals[0]
    assert vals[2] > vals[3] > vals[4]


def test_pv_gauss_limit():
    """eta=0 → pure Gaussian: G(fwhm/2) ≈ 0.5."""
    peak_pos = 43.0
    fwhm = 0.10
    val = pseudo_voigt_profile(np.array([peak_pos + fwhm / 2]), peak_pos, fwhm, eta=0.0)[0]
    assert abs(val - 0.5) < 0.01


def test_pv_lorentz_limit():
    """eta=1 → pure Lorentzian: L(fwhm/2) = 0.5."""
    peak_pos = 43.0
    fwhm = 0.10
    val = pseudo_voigt_profile(np.array([peak_pos + fwhm / 2]), peak_pos, fwhm, eta=1.0)[0]
    assert abs(val - 0.5) < 1e-10


def test_pv_lorentz_denominator_factor4():
    """2.8: L uses (fwhm/2)^2 → half-height at fwhm/2 exactly."""
    peak_pos = 40.0
    fwhm = 0.20
    val = pseudo_voigt_profile(np.array([peak_pos + fwhm / 2]), peak_pos, fwhm, eta=1.0)[0]
    assert abs(val - 0.5) < 1e-10, f"L(HWHM) = {val:.6f} (expected 0.5; denominator wrong if >> 0.5)"


# ============================================================
# build_icalc_unit  — criteria 2.10–2.13
# ============================================================

REFL_CU_FCC = [
    {'h': 1, 'k': 1, 'l': 1, 'two_theta': 43.32, 'intensity_rel': 100.0,
     'd_hkl': 2.088, 'multiplicity': 8, 'F_sq': 1000.0},
    {'h': 2, 'k': 0, 'l': 0, 'two_theta': 50.44, 'intensity_rel': 46.0,
     'd_hkl': 1.808, 'multiplicity': 6, 'F_sq': 460.0},
    {'h': 2, 'k': 2, 'l': 0, 'two_theta': 74.13, 'intensity_rel': 27.0,
     'd_hkl': 1.278, 'multiplicity': 12, 'F_sq': 270.0},
]


def test_icalc_unit_single_peak_max():
    """2.13: max of Icalc_unit near position of most intense peak."""
    Icalc, n_used = build_icalc_unit(TTH, REFL_CU_FCC, U=U, V=V, W=W, eta=ETA)
    peak_pos_est = TTH[np.argmax(Icalc)]
    assert abs(peak_pos_est - 43.32) < 0.5, f"Max at {peak_pos_est:.2f}°, expected ~43.32°"
    assert n_used == 3


def test_icalc_unit_nonnegative():
    """2.12: Icalc_unit >= 0 everywhere."""
    Icalc, _ = build_icalc_unit(TTH, REFL_CU_FCC, U=U, V=V, W=W, eta=ETA)
    assert (Icalc >= 0).all(), "Icalc_unit has negative values [2.12]"


def test_intensity_rel_zero_skipped():
    """2.10: intensity_rel = 0 → peak skipped."""
    refls = [{'h': 1, 'k': 0, 'l': 0, 'two_theta': 43.32, 'intensity_rel': 0.0,
               'd_hkl': 2.0, 'multiplicity': 2, 'F_sq': 0.0}]
    Icalc, n_used = build_icalc_unit(TTH, refls)
    assert n_used == 0
    assert Icalc.sum() == 0.0


def test_intensity_rel_negative_skipped():
    """2.10: intensity_rel < 0 → skipped."""
    refls = [{'h': 1, 'k': 0, 'l': 0, 'two_theta': 43.32, 'intensity_rel': -5.0,
               'd_hkl': 2.0, 'multiplicity': 2, 'F_sq': 0.0}]
    Icalc, n_used = build_icalc_unit(TTH, refls)
    assert n_used == 0


def test_peak_outside_range_skipped():
    """2.11: peak far outside grid → n_peaks_used = 0, Icalc = zeros."""
    refls = [{'h': 1, 'k': 0, 'l': 0, 'two_theta': 200.0, 'intensity_rel': 100.0,
               'd_hkl': 1.0, 'multiplicity': 2, 'F_sq': 1000.0}]
    Icalc, n_used = build_icalc_unit(TTH, refls)
    assert n_used == 0
    assert Icalc.sum() == pytest.approx(0.0)


def test_empty_reflections():
    Icalc, n_used = build_icalc_unit(TTH, [])
    assert n_used == 0
    assert Icalc.sum() == 0.0
