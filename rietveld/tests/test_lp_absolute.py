"""
Tests for lorentz_polarization + build_icalc_unit_absolute
— science review criteria 7.8, 7.9.
"""
import numpy as np
import pytest

from pattern_calc import lorentz_polarization, build_icalc_unit_absolute, build_icalc_unit


TTH = np.linspace(20.0, 120.0, 5001)


# --- 7.8: Lp uses θ = 2θ/2 ---

def test_lp_positive():
    tth = np.linspace(5.0, 150.0, 1000)
    lp = lorentz_polarization(tth)
    assert (lp > 0).all()


def test_lp_diverges_low_angle():
    """7.8: Lp grows strongly toward 2θ→0 (1/sin²θ Lorentz term)."""
    assert lorentz_polarization(2.0) > lorentz_polarization(20.0) > lorentz_polarization(60.0)


def test_lp_ratio_matches_csv():
    """7.8: Lp(5.17°)/Lp(12.22°) ≈ 5–6 (matches observed CSV ir/(m·Fsq) ~6×)."""
    ratio = float(lorentz_polarization(5.171)) / float(lorentz_polarization(12.222))
    assert 4.0 < ratio < 7.0, f"ratio={ratio:.2f}"


def test_lp_half_angle_not_full():
    """7.8: θ=2θ/2. Lp(2θ=90°): θ=45°, num=1+cos²(90°)=1,
    denom=sin²(45°)·cos(45°)=0.354 → Lp≈2.828.
    If code wrongly used θ=2θ=90°: denom=sin²(90°)·cos(90°)=0 → divergence."""
    lp90 = float(lorentz_polarization(90.0))
    assert lp90 == pytest.approx(2.828, rel=0.02)


def test_lp_with_monochromator():
    """Monochromator changes numerator → different (smaller) Lp."""
    lp_no = float(lorentz_polarization(40.0))
    lp_mono = float(lorentz_polarization(40.0, two_theta_mono_deg=26.6))
    assert lp_mono < lp_no
    assert lp_mono > 0


# --- 7.9: absolute pattern weights by mult·F_sq·Lp, not intensity_rel ---

REFL = [
    {'two_theta': 30.0, 'intensity_rel': 100.0, 'multiplicity': 4, 'F_sq': 500.0},
    {'two_theta': 60.0, 'intensity_rel': 50.0,  'multiplicity': 2, 'F_sq': 800.0},
]


def test_absolute_differs_from_rel():
    """7.9: absolute basis differs from intensity_rel basis (different weights)."""
    abs_pat, n_abs = build_icalc_unit_absolute(TTH, REFL)
    rel_pat, n_rel = build_icalc_unit(TTH, REFL)
    assert n_abs == 2 and n_rel == 2
    # peak-height ratio differs because absolute uses mult·Fsq·Lp(θ),
    # rel uses the stored intensity_rel
    assert not np.allclose(abs_pat / abs_pat.max(), rel_pat / rel_pat.max(), atol=1e-3)


def test_absolute_weight_formula():
    """7.9: peak weight = mult·F_sq·Lp(θ)."""
    # single peak, check integrated weight proportional to mult·Fsq·Lp
    r = [{'two_theta': 40.0, 'intensity_rel': 1.0, 'multiplicity': 3, 'F_sq': 200.0}]
    pat, n = build_icalc_unit_absolute(TTH, r)
    expected_peak = 3 * 200.0 * float(lorentz_polarization(40.0))
    # pseudo-Voigt peak normalised to 1.0 at centre → max ≈ expected_peak
    assert pat.max() == pytest.approx(expected_peak, rel=0.02)
    assert n == 1


def test_absolute_skips_zero_fsq():
    """F_sq=0 or mult=0 → skipped."""
    r = [{'two_theta': 40.0, 'intensity_rel': 100.0, 'multiplicity': 0, 'F_sq': 200.0},
         {'two_theta': 50.0, 'intensity_rel': 100.0, 'multiplicity': 4, 'F_sq': 0.0}]
    pat, n = build_icalc_unit_absolute(TTH, r)
    assert n == 0
    assert pat.sum() == 0.0


def test_absolute_nonnegative():
    pat, _ = build_icalc_unit_absolute(TTH, REFL)
    assert (pat >= 0).all()
