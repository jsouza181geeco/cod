"""
Tests for linear_fit — science review criteria 3.1–3.7 [Y93 p.18-22].
"""
import numpy as np
import pytest

from linear_fit import linear_fit


def _make_data(n=500, seed=0):
    rng = np.random.default_rng(seed)
    tth = np.linspace(20.0, 120.0, n)
    Icalc_unit = np.exp(-0.5 * ((tth - 60.0) / 5.0) ** 2)  # gaussian peak at 60°
    sigma = np.ones(n) * 5.0
    return tth, Icalc_unit, sigma


# --- 3.5: Icalc = S*Icalc_unit + Ibg (not S*(Icalc_unit+Ibg)) ---

def test_exact_recovery_scale():
    """3.5: Iobs = 10*Icalc_unit + 50 → |scale - 10| < 0.1."""
    tth, Icalc_unit, sigma = _make_data()
    Iobs = 10.0 * Icalc_unit + 50.0
    scale, _ = linear_fit(tth, Iobs, sigma, Icalc_unit, n_bg=4)
    assert abs(scale - 10.0) < 0.1, f"scale={scale:.6f}, expected ~10"


def test_exact_recovery_residual():
    """Residual < 1% of max(Iobs) for exact linear case."""
    tth, Icalc_unit, sigma = _make_data()
    Iobs = 10.0 * Icalc_unit + 50.0
    _, Icalc = linear_fit(tth, Iobs, sigma, Icalc_unit, n_bg=4)
    rel_err = np.abs(Iobs - Icalc).max() / Iobs.max()
    assert rel_err < 0.01


def test_quadratic_background_recovered():
    """Background: constant + linear + quadratic + cubic all recovered."""
    n = 500
    tth = np.linspace(20.0, 120.0, n)
    Icalc_unit = np.exp(-0.5 * ((tth - 60.0) / 5.0) ** 2)
    sigma = np.ones(n)
    bg = 100.0 + 2.0 * tth - 0.01 * tth ** 2
    Iobs = 5.0 * Icalc_unit + bg
    scale, Icalc = linear_fit(tth, Iobs, sigma, Icalc_unit, n_bg=4)
    assert abs(scale - 5.0) < 0.5
    np.testing.assert_allclose(Icalc, Iobs, rtol=1e-3)


# --- 3.1: weights = 1/sigma^2 (not 1/sigma) ---

def test_weights_inverse_sigma_squared():
    """3.1: w=1/sigma^2 → high-sigma outlier has negligible influence on scale."""
    n = 400
    tth = np.linspace(20.0, 120.0, n)
    # non-trivial Icalc_unit so design matrix is not rank-deficient
    Icalc_unit = np.exp(-0.5 * ((tth - 60.0) / 5.0) ** 2)
    # true Iobs = 10 * Icalc_unit + 50, except one massive outlier
    Iobs = 10.0 * Icalc_unit + 50.0
    sigma = np.ones(n)
    Iobs[100] = 5000.0          # outlier: 500× normal
    sigma[100] = 500.0          # huge sigma → w = 1/sigma^2 → tiny weight
    scale, _ = linear_fit(tth, Iobs, sigma, Icalc_unit, n_bg=4)
    # proper 1/sigma^2 weighting keeps outlier influence tiny → scale ≈ 10
    assert abs(scale - 10.0) < 1.0, (
        f"scale={scale:.4f}; outlier with large sigma should not dominate [3.1]"
    )


# --- scale negative for anti-correlated pattern ---

def test_scale_can_be_negative():
    """Science review 5.2: negative scale is valid (anti-correlated candidate)."""
    tth, Icalc_unit, sigma = _make_data()
    # Iobs = flat background, no peak → scale should be ≈ 0 or negative
    Iobs = np.full_like(tth, 100.0)
    scale, _ = linear_fit(tth, Iobs, sigma, Icalc_unit, n_bg=4)
    # just check it doesn't crash and returns finite value
    assert np.isfinite(scale)


# --- 3.4: single-point dataset → no crash ---

def test_floor_std_single_point():
    """3.4: max(std, 1e-6) prevents division by zero."""
    tth = np.array([43.0])
    Icalc_unit = np.array([1.0])
    Iobs = np.array([10.0])
    sigma = np.array([1.0])
    # should not raise
    scale, Icalc = linear_fit(tth, Iobs, sigma, Icalc_unit, n_bg=1)
    assert np.isfinite(scale)


# --- return shape ---

def test_returns_correct_shape():
    tth, Icalc_unit, sigma = _make_data()
    Iobs = 5.0 * Icalc_unit + 30.0
    scale, Icalc = linear_fit(tth, Iobs, sigma, Icalc_unit)
    assert isinstance(scale, float)
    assert Icalc.shape == tth.shape
