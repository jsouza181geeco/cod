"""
Tests for crystallo_utils — science review criteria 7.1, 7.2, 7.7 [HH87].
"""
import numpy as np
import pytest

from crystallo_utils import cell_volume, molar_mass


# --- 7.1: cell_volume triclinic general ---

def test_cubic_reduces_to_abc():
    """7.1: orthogonal (90,90,90) → V = a·b·c."""
    assert cell_volume(5.0, 5.0, 5.0, 90, 90, 90) == pytest.approx(125.0, abs=1e-6)
    assert cell_volume(3.0, 4.0, 5.0, 90, 90, 90) == pytest.approx(60.0, abs=1e-6)


def test_hexagonal_quartz():
    """7.1: hexagonal γ=120° (quartz) → V = (√3/2)a²c, NOT a·b·c."""
    a, c = 4.913, 5.405
    V = cell_volume(a, a, c, 90, 90, 120)
    expected = np.sqrt(3) / 2 * a**2 * c
    assert V == pytest.approx(expected, rel=1e-4)
    assert V == pytest.approx(112.99, abs=0.1)
    # must NOT equal the orthogonal a·b·c
    assert abs(V - a * a * c) > 10


def test_triclinic_general_positive():
    """7.1: general triclinic angles → finite positive volume."""
    V = cell_volume(6, 7, 8, 80, 85, 75)
    assert np.isfinite(V) and V > 0
    assert V < 6 * 7 * 8  # non-orthogonal cell is smaller than the box


# --- 7.2: molar_mass full parse ---

def test_molar_mass_simple_oxides():
    """7.2: known oxides."""
    assert molar_mass('SiO2') == pytest.approx(60.083, abs=0.01)
    assert molar_mass('Fe2 O3') == pytest.approx(159.687, abs=0.01)
    assert molar_mass('Ca C O3') == pytest.approx(100.086, abs=0.01)


def test_molar_mass_spaced_and_unspaced():
    """7.2: space-separated == concatenated for unambiguous formulas."""
    assert molar_mass('Ca C O3') == pytest.approx(molar_mass('CaCO3'), rel=1e-9)


def test_molar_mass_fractional_counts():
    """7.2: fractional occupancies in formula (e.g. Br0.8)."""
    m = molar_mass('Br0.8 C2')
    expected = 0.8 * 79.904 + 2 * 12.011
    assert m == pytest.approx(expected, abs=0.01)


def test_molar_mass_implicit_one():
    """7.2: element without count → 1 atom."""
    assert molar_mass('H2 O') == pytest.approx(18.015, abs=0.01)


# --- 7.7: unknown element raises ---

def test_unknown_element_raises():
    """7.7: unknown element → ValueError (NOT silent skip)."""
    with pytest.raises(ValueError):
        molar_mass('Xx2 O3')


def test_empty_formula_raises():
    with pytest.raises(ValueError):
        molar_mass('')


def test_no_tokens_raises():
    with pytest.raises(ValueError):
        molar_mass('123')
