"""
Tests for multi_phase_fit — science review criteria 6.1–6.4, 6.6 [Y93; BH88].
"""
import numpy as np
import pytest

from multi_phase_fit import multi_phase_fit


def _two_phase(n=600, seed=0):
    tth = np.linspace(20, 120, n)
    uA = np.exp(-0.5*((tth-40)/3)**2) + 0.5*np.exp(-0.5*((tth-70)/3)**2)
    uB = np.exp(-0.5*((tth-55)/3)**2) + 0.3*np.exp(-0.5*((tth-90)/3)**2)
    sigma = np.ones(n)
    return tth, uA, uB, sigma


# --- 6.1 / 6.4: design matrix K phases + background, recovers known scales ---

def test_recover_known_scales():
    """6.1: Iobs = 3·uA + 7·uB + bg → scales [3, 7]."""
    tth, uA, uB, sigma = _two_phase()
    bg = 100 + 2*tth - 0.01*tth**2
    Iobs = 3.0*uA + 7.0*uB + bg
    scales, Icalc, _ = multi_phase_fit(tth, Iobs, sigma, [uA, uB], n_bg=4)
    assert scales[0] == pytest.approx(3.0, abs=0.01)
    assert scales[1] == pytest.approx(7.0, abs=0.01)


def test_combined_fit_matches():
    tth, uA, uB, sigma = _two_phase()
    bg = 50 + tth
    Iobs = 2.0*uA + 5.0*uB + bg
    _, Icalc, _ = multi_phase_fit(tth, Iobs, sigma, [uA, uB], n_bg=4)
    np.testing.assert_allclose(Icalc, Iobs, rtol=1e-3)


# --- 6.2: non-negativity ---

def test_scales_nonnegative_anticorrelated():
    """6.2: phase absent from sample → scale clamped to 0, never negative."""
    tth, uA, uB, sigma = _two_phase()
    bg = 100 + tth
    Iobs = 5.0*uA + bg      # only A present
    scales, _, _ = multi_phase_fit(tth, Iobs, sigma, [uA, uB], n_bg=4)
    assert (scales >= -1e-9).all()
    assert scales[1] == pytest.approx(0.0, abs=1e-6)   # B → 0
    assert scales[0] == pytest.approx(5.0, abs=0.05)


def test_all_phases_nonnegative_random():
    """6.2: random noisy data → no negative scale."""
    rng = np.random.default_rng(1)
    tth, uA, uB, sigma = _two_phase()
    Iobs = rng.uniform(50, 500, len(tth))
    scales, _, _ = multi_phase_fit(tth, Iobs, sigma, [uA, uB], n_bg=4)
    assert (scales >= -1e-9).all()


# --- 6.3: background can be negative ---

def test_background_can_be_negative():
    """6.3: background coeffs unconstrained (lower bound -inf)."""
    tth, uA, uB, sigma = _two_phase()
    bg = -50 + 0.5*tth     # background dips negative at low 2θ
    Iobs = 3.0*uA + bg
    scales, Icalc, bg_coeffs = multi_phase_fit(tth, Iobs, sigma, [uA, uB], n_bg=4)
    # fit must reproduce data including the negative-trending background
    np.testing.assert_allclose(Icalc, Iobs, rtol=1e-2, atol=1.0)


# --- 6.6: Rwp_combined <= Rwp_single (by construction) ---

def test_combined_le_single():
    """6.6: adding phases cannot worsen the WLS fit (single = feasible subset)."""
    from fom import calc_fom
    from linear_fit import linear_fit
    tth, uA, uB, sigma = _two_phase()
    bg = 80 + tth
    Iobs = 4.0*uA + 6.0*uB + bg

    _, Icalc_single = linear_fit(tth, Iobs, sigma, uA, n_bg=4)
    rwp_single = calc_fom(Iobs, Icalc_single, sigma, n_params=5)['Rwp']

    _, Icalc_multi, _ = multi_phase_fit(tth, Iobs, sigma, [uA, uB], n_bg=4)
    rwp_multi = calc_fom(Iobs, Icalc_multi, sigma, n_params=6)['Rwp']

    assert rwp_multi <= rwp_single + 1e-9


# --- edge: empty phase list ---

def test_empty_raises():
    tth, uA, uB, sigma = _two_phase()
    with pytest.raises(ValueError):
        multi_phase_fit(tth, np.ones_like(tth), sigma, [], n_bg=4)


# --- single phase behaves like linear_fit scale ---

def test_single_phase_recovers_scale():
    tth, uA, uB, sigma = _two_phase()
    Iobs = 8.0*uA + 100.0
    scales, _, _ = multi_phase_fit(tth, Iobs, sigma, [uA], n_bg=4)
    assert scales[0] == pytest.approx(8.0, abs=0.05)
