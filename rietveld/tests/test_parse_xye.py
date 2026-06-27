"""
Tests for data_loader.parse_xye — science review criteria 1.1–1.4 [Y93; MCC99].
"""
import numpy as np
import pytest
from pathlib import Path

from data_loader import parse_xye, parse_asc, parse_diffractogram

DATA_DIR = Path(__file__).parent.parent
N_MIN = 15  # > parse_xye threshold of 10


def _make_xye(tmp_path, lines, fname="test.xye"):
    p = tmp_path / fname
    p.write_text("\n".join(lines))
    return str(p)


def _col2_lines(n=N_MIN, i_val=100.0):
    return [f"{20.0 + 0.02*i:.3f} {i_val}" for i in range(n)]


def _col3_lines(n=N_MIN, sigma_val=10.0):
    return [f"{20.0 + 0.02*i:.3f} 100.0 {sigma_val}" for i in range(n)]


# --- shape and range (real file) ---

def test_shape_real_file():
    """Real XYE: 5001 data points (4 comment lines excluded)."""
    tth, Iobs, sigma = parse_xye(DATA_DIR / "synthetic_candidate17.xye")
    assert tth.shape == (5001,)
    assert Iobs.shape == (5001,)
    assert sigma.shape == (5001,)


def test_range_real_file():
    tth, _, _ = parse_xye(DATA_DIR / "synthetic_candidate17.xye")
    assert abs(tth[0]  - 20.0) < 0.1
    assert abs(tth[-1] - 120.0) < 0.1


# --- 1.2: 2theta sorted ascending ---

def test_tth_sorted_after_unsorted_input(tmp_path):
    """1.2: parse_xye must return 2theta in ascending order."""
    lines = [f"{30.0 - 0.02*i:.3f} 100.0 10.0" for i in range(N_MIN)]  # descending
    tth, _, _ = parse_xye(_make_xye(tmp_path, lines))
    assert (np.diff(tth) > 0).all(), "2theta not sorted ascending [1.2]"


# --- comment/empty lines ignored ---

def test_comments_and_blanks_ignored(tmp_path):
    data_lines = [f"{20.0 + 0.02*i:.3f} 100.0 10.0" for i in range(N_MIN)]
    lines = ["# file header", "", "# another comment"] + data_lines
    tth, _, _ = parse_xye(_make_xye(tmp_path, lines))
    assert len(tth) == N_MIN
    assert tth[0] == pytest.approx(20.0, abs=0.01)


# --- 1.4: sigma fallback for 2-column ---

def test_sigma_fallback_two_col(tmp_path):
    """1.4: 2-col file → sigma = sqrt(max(Iobs, 1))."""
    lines = [f"{20.0 + 0.02*i:.3f} {100.0 + i*10}" for i in range(N_MIN)]
    tth, Iobs, sigma = parse_xye(_make_xye(tmp_path, lines))
    expected = np.sqrt(np.maximum(Iobs, 1.0))
    np.testing.assert_allclose(sigma, expected)


def test_sigma_not_constant_for_varying_intensity(tmp_path):
    """1.4: sigma must vary with Iobs, not be constant."""
    # increasing intensities: 100, 400, 900 ...
    lines = [f"{20.0 + 0.02*i:.3f} {(i+1)**2 * 4.0}" for i in range(N_MIN)]
    _, Iobs, sigma = parse_xye(_make_xye(tmp_path, lines))
    # sigma = sqrt(I) → strictly increasing
    assert sigma[-1] > sigma[0]


# --- 1.1: sigma > 0 always ---

def test_sigma_always_positive_two_col(tmp_path):
    """1.1: sigma > 0 even for I=0 or I<0."""
    base = [f"{20.0 + 0.02*i:.3f} 200.0" for i in range(N_MIN - 3)]
    edge = [
        "20.60   0.0",    # I=0 → sigma floor = 1
        "20.62  -5.0",   # I<0 → sigma floor = 1
        "20.64  25.0",   # I=25 → sigma=5
    ]
    lines = base + edge
    _, Iobs, sigma = parse_xye(_make_xye(tmp_path, lines))
    assert (sigma > 0).all(), "sigma must be > 0 at all points [1.1]"


def test_sigma_zero_replaced(tmp_path):
    """1.1: 3-column with sigma=0 → replaced."""
    base = [f"{20.0 + 0.02*i:.3f} 100.0 10.0" for i in range(N_MIN - 2)]
    edge = ["20.60 100.0 0.0", "20.62 200.0 0.0"]
    tth, Iobs, sigma = parse_xye(_make_xye(tmp_path, base + edge))
    assert (sigma > 0).all()


# --- 3-column: sigma taken directly ---

def test_three_col_sigma_direct(tmp_path):
    """Three-column file: sigma taken as-is from column 3."""
    lines = _col3_lines(n=N_MIN, sigma_val=99.9)
    _, _, sigma = parse_xye(_make_xye(tmp_path, lines))
    np.testing.assert_allclose(sigma, 99.9, rtol=1e-6)


# --- ValueError cases ---

def test_nonexistent_file():
    with pytest.raises((ValueError, FileNotFoundError, OSError)):
        parse_xye("this_file_does_not_exist.xye")


def test_too_few_points_raises(tmp_path):
    lines = ["20.0 100.0 10.0", "20.1 200.0 14.0"]  # 2 < 10
    with pytest.raises(ValueError):
        parse_xye(_make_xye(tmp_path, lines))


# --- parse_asc ---

def _make_asc(tmp_path, lines, fname="test.asc"):
    p = tmp_path / fname
    p.write_text("\n".join(lines))
    return str(p)


def test_parse_asc_basic(tmp_path):
    """parse_asc: 2-col numeric → correct tth/Iobs, sigma=sqrt(Iobs)."""
    lines = [f"{20.0 + 0.01*i:.3f}  {100.0 + i*5}" for i in range(N_MIN)]
    tth, Iobs, sigma = parse_asc(_make_asc(tmp_path, lines))
    assert len(tth) == N_MIN
    np.testing.assert_allclose(sigma, np.sqrt(np.maximum(Iobs, 1.0)))


def test_parse_asc_skips_header(tmp_path):
    """parse_asc: non-numeric header lines silently skipped."""
    header = ["SAMPLE: Jarosita", "DATE: 2026-06-25", "2THETA COUNTS"]
    data = [f"{20.0 + 0.01*i:.3f}  {500.0}" for i in range(N_MIN)]
    tth, _, _ = parse_asc(_make_asc(tmp_path, header + data))
    assert len(tth) == N_MIN


def test_parse_asc_comma_decimal(tmp_path):
    """parse_asc: European comma decimal separator handled."""
    lines = [f"{(20.0 + 0.01*i):.3f}  {(100.0 + i):.3f}".replace('.', ',') for i in range(N_MIN)]
    tth, Iobs, _ = parse_asc(_make_asc(tmp_path, lines))
    assert len(tth) == N_MIN
    assert tth[0] == pytest.approx(20.0, abs=0.01)


def test_parse_asc_sorted(tmp_path):
    """parse_asc: output 2theta sorted ascending even if file is not."""
    lines = [f"{30.0 - 0.01*i:.3f}  100.0" for i in range(N_MIN)]
    tth, _, _ = parse_asc(_make_asc(tmp_path, lines))
    assert (np.diff(tth) > 0).all()


def test_parse_asc_sigma_floor(tmp_path):
    """parse_asc: sigma >= 1 even for I=0."""
    base = [f"{20.0 + 0.01*i:.3f}  200.0" for i in range(N_MIN - 1)]
    edge = ["20.99  0.0"]
    _, _, sigma = parse_asc(_make_asc(tmp_path, base + edge))
    assert (sigma >= 1.0).all()


def test_parse_asc_nonexistent():
    with pytest.raises(ValueError):
        parse_asc("does_not_exist.asc")


def test_parse_asc_too_few_points(tmp_path):
    lines = ["20.0  100.0", "20.1  200.0"]
    with pytest.raises(ValueError):
        parse_asc(_make_asc(tmp_path, lines))


# --- parse_diffractogram dispatcher ---

def test_dispatcher_routes_asc(tmp_path):
    """parse_diffractogram dispatches .asc extension to parse_asc."""
    lines = [f"{20.0 + 0.01*i:.3f}  {100.0}" for i in range(N_MIN)]
    p = tmp_path / "sample.asc"
    p.write_text("\n".join(lines))
    tth, Iobs, sigma = parse_diffractogram(str(p))
    assert len(tth) == N_MIN
    np.testing.assert_allclose(sigma, np.sqrt(np.maximum(Iobs, 1.0)))


def test_dispatcher_routes_xye(tmp_path):
    """parse_diffractogram dispatches .xye extension to parse_xye."""
    lines = _col3_lines(n=N_MIN, sigma_val=7.0)
    p = tmp_path / "sample.xye"
    p.write_text("\n".join(lines))
    _, _, sigma = parse_diffractogram(str(p))
    np.testing.assert_allclose(sigma, 7.0, rtol=1e-6)


def test_dispatcher_case_insensitive(tmp_path):
    """parse_diffractogram: .ASC (uppercase) dispatches to parse_asc."""
    lines = [f"{20.0 + 0.01*i:.3f}  {100.0}" for i in range(N_MIN)]
    p = tmp_path / "SAMPLE.ASC"
    p.write_text("\n".join(lines))
    tth, _, _ = parse_diffractogram(str(p))
    assert len(tth) == N_MIN
