"""
Tests for pipeline.run_pipeline — science review criteria 5.1–5.5 [R69; MCC99].
"""
import numpy as np
import pytest
from pathlib import Path

from pipeline import run_pipeline, FIXED_PARAMS
from data_loader import load_candidates_csv, parse_xye

DATA_DIR = Path(__file__).parent.parent

XYE_PATH = DATA_DIR / "synthetic_candidate17.xye"
CSV_PATH = DATA_DIR / "data-1782394014136.csv"
BEST_COD_ID = 1569653  # candidate 17, used to generate the synthetic XYE


@pytest.fixture(scope="module")
def result():
    candidates = load_candidates_csv(CSV_PATH)
    return run_pipeline(XYE_PATH, candidates, db_client=None)


# --- 5.3: sorted by Rwp ascending ---

def test_sorted_by_rwp(result):
    """5.3: candidates sorted Rwp ascending (best first)."""
    rwps = [r.Rwp for r in result.candidates]
    assert rwps == sorted(rwps), "Candidates not sorted by Rwp [5.3]"


# --- 5.1: n_params = 1 + n_bg — chi2 consistent ---

def test_chi2_consistent_with_rwp_rexp(result):
    """5.1: chi2 = (Rwp/Rexp)^2 for all candidates."""
    for r in result.candidates:
        expected = (r.Rwp / r.Rexp) ** 2
        assert abs(r.chi2 - expected) < 1e-6, f"chi2 inconsistent for cod_id={r.cod_id} [5.1]"


# --- 5.4: db_client=None works ---

def test_no_db_metadata_is_none(result):
    """5.4: db_client=None → metadata=None for all candidates."""
    assert all(r.metadata is None for r in result.candidates)


# --- 5.2: scale can be negative ---

def test_negative_scale_not_filtered(result):
    """5.2: candidates with negative scale included (not removed)."""
    negative_scale = [r for r in result.candidates if r.scale < 0]
    # we know from manual run that 6 candidates have negative scale
    assert len(negative_scale) >= 1, "Expected some candidates with negative scale [5.2]"


# --- 5.5: n_peaks_used = 0 → no crash, high Rwp ---

def test_zero_peaks_candidate_included():
    """5.5: candidate with all intensity_rel=0 → n_peaks_used=0, no crash, high Rwp (background-only fit)."""
    from models import CandidateInput
    empty_cand = CandidateInput(
        cod_id=9999999,
        reflections=[{'h': 1, 'k': 0, 'l': 0, 'two_theta': 43.0,
                       'intensity_rel': 0.0, 'd_hkl': 2.0,
                       'multiplicity': 2, 'F_sq': 0.0}],
    )
    result = run_pipeline(XYE_PATH, [empty_cand], db_client=None)
    r = result.candidates[0]
    assert r.n_peaks_used == 0
    # background-only fit: Rwp > 0.3 (no phase signal contributes)
    assert r.Rwp > 0.3, f"Rwp={r.Rwp:.4f} for empty candidate (should be > 0.3)"


# --- correct candidate is #1 ---

def test_best_candidate_is_correct(result):
    """Synthetic XYE from 1569653 → should rank first."""
    assert result.best().cod_id == BEST_COD_ID, (
        f"Best = {result.best().cod_id}, expected {BEST_COD_ID}"
    )


def test_best_rwp_below_threshold(result):
    """Correct candidate: Rwp < 0.15 [T06]."""
    assert result.best().Rwp < 0.15, f"Rwp={result.best().Rwp:.4f}"


def test_best_chi2_near_one(result):
    """Correct candidate: chi2 ≈ 1 (Poisson noise, correct model)."""
    chi2 = result.best().chi2
    assert 0.8 < chi2 < 1.5, f"chi2={chi2:.3f} — should be ≈1 for correct phase [T06]"


def test_viable_not_empty(result):
    """result.viable() returns at least the correct candidate."""
    assert len(result.viable()) >= 1


# --- total candidates ---

def test_all_candidates_returned(result):
    candidates = load_candidates_csv(CSV_PATH)
    assert len(result.candidates) == len(candidates)


# --- n_points ---

def test_n_points_correct(result):
    tth, _, _ = parse_xye(XYE_PATH)
    assert result.n_points == len(tth)


# --- gap between correct and second ---

def test_large_rwp_gap(result):
    """Correct candidate Rwp << all others (synthetic data separation)."""
    best_rwp  = result.candidates[0].Rwp
    second_rwp = result.candidates[1].Rwp
    assert second_rwp / best_rwp > 5, (
        f"Gap too small: best={best_rwp:.4f} second={second_rwp:.4f}"
    )


# --- FIXED_PARAMS keys present ---

def test_fixed_params_has_required_keys():
    required = {'U', 'V', 'W', 'eta', 'n_bg', 'wavelength'}
    assert required.issubset(set(FIXED_PARAMS.keys()))
