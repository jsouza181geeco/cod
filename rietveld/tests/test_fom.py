"""
Tests for fom.calc_fom — science review criteria 4.1–4.7 [Y93 p.22-26; T06].
"""
import numpy as np
import pytest

from fom import calc_fom


def _simple(n=100, seed=0):
    rng = np.random.default_rng(seed)
    Iobs = rng.uniform(100.0, 1000.0, n)
    sigma = np.sqrt(np.maximum(Iobs, 1.0))
    return Iobs, sigma


# --- 4.6: sanity: Icalc=Iobs → Rwp ≈ 0 ---

def test_perfect_fit_rwp_zero():
    """4.6: Icalc = Iobs → Rwp < 1e-10."""
    Iobs, sigma = _simple()
    fom = calc_fom(Iobs, Iobs, sigma, n_params=5)
    assert fom['Rwp'] < 1e-10, f"Rwp={fom['Rwp']:.2e} for perfect fit (expected < 1e-10)"


# --- 4.6: sanity: Icalc=0 → Rwp ≈ 1 ---

def test_null_fit_rwp_one():
    """4.6: Icalc = 0 → Rwp ≈ 1.0."""
    Iobs, sigma = _simple()
    fom = calc_fom(Iobs, np.zeros_like(Iobs), sigma, n_params=5)
    assert abs(fom['Rwp'] - 1.0) < 0.01, f"Rwp={fom['Rwp']:.4f} for null fit (expected ≈ 1)"


# --- 4.7: all >= 0 ---

def test_all_nonnegative():
    """4.7: Rwp, Rp, Rexp, chi2 >= 0."""
    Iobs, sigma = _simple()
    Icalc = Iobs * 0.8  # partial fit
    fom = calc_fom(Iobs, Icalc, sigma, n_params=5)
    assert fom['Rwp']  >= 0
    assert fom['Rp']   >= 0
    assert fom['Rexp'] >= 0
    assert fom['chi2'] >= 0


# --- 4.5: chi2 = (Rwp/Rexp)^2 ---

def test_chi2_equals_rwp_over_rexp_squared():
    """4.5: chi2 = (Rwp/Rexp)^2 algebraically."""
    Iobs, sigma = _simple()
    Icalc = Iobs * 0.9
    fom = calc_fom(Iobs, Icalc, sigma, n_params=5)
    expected = (fom['Rwp'] / fom['Rexp']) ** 2
    assert abs(fom['chi2'] - expected) < 1e-8


# --- 4.2: N-P degrees of freedom in Rexp ---

def test_rexp_uses_n_minus_p():
    """4.2: larger n_params → smaller Rexp (more params = fewer DoF)."""
    Iobs, sigma = _simple(n=200)
    Icalc = Iobs * 0.95
    fom_p5  = calc_fom(Iobs, Icalc, sigma, n_params=5)
    fom_p10 = calc_fom(Iobs, Icalc, sigma, n_params=10)
    assert fom_p10['Rexp'] < fom_p5['Rexp'], "Larger n_params must reduce Rexp [4.2]"


# --- 4.3: floor max(N-P, 1) ---

def test_floor_n_minus_p():
    """4.3: n_params > N → no division by zero, Rexp finite."""
    Iobs = np.array([100.0, 200.0])
    sigma = np.array([10.0,  14.0])
    Icalc = Iobs * 0.9
    fom = calc_fom(Iobs, Icalc, sigma, n_params=100)  # P >> N
    assert np.isfinite(fom['Rexp'])
    assert np.isfinite(fom['chi2'])


# --- 4.4: floor sum_w_Iobs2 ---

def test_floor_zero_intensities():
    """4.4: Iobs = 0 → no division by zero."""
    Iobs  = np.zeros(50)
    sigma = np.ones(50)
    Icalc = np.zeros(50)
    fom = calc_fom(Iobs, Icalc, sigma, n_params=5)
    assert np.isfinite(fom['Rwp'])
    assert np.isfinite(fom['Rexp'])


# --- 4.1: Rwp denominator uses Iobs^2, not Icalc^2 ---

def test_rwp_denominator_uses_iobs():
    """4.1: Rwp denominator = sum(w*Iobs^2), not sum(w*Icalc^2)."""
    # If Icalc = 0.5*Iobs, then:
    #   correct:  Rwp = sqrt(sum(w*(0.5*I)^2) / sum(w*I^2)) = 0.5
    #   wrong:    Rwp = sqrt(sum(w*(0.5*I)^2) / sum(w*(0.5*I)^2)) = 1.0
    Iobs  = np.full(100, 100.0)
    sigma = np.full(100, 10.0)
    Icalc = Iobs * 0.5
    fom = calc_fom(Iobs, Icalc, sigma, n_params=5)
    assert abs(fom['Rwp'] - 0.5) < 0.01, (
        f"Rwp={fom['Rwp']:.4f} — if ~1.0, denominator uses Icalc^2 instead of Iobs^2 [4.1]"
    )


# --- chi2 ≈ 1.0 for correct Poisson noise model ---

def test_chi2_near_one_for_correct_model():
    """chi2 ≈ 1 when residual = pure Poisson noise (scale = 1, large N)."""
    rng = np.random.default_rng(42)
    Icalc_true = np.full(5000, 200.0)
    Iobs = rng.poisson(Icalc_true).astype(float)
    sigma = np.sqrt(np.maximum(Iobs, 1.0))
    fom = calc_fom(Iobs, Icalc_true, sigma, n_params=5)
    assert 0.8 < fom['chi2'] < 1.3, f"chi2={fom['chi2']:.3f} far from 1 for Poisson model"
