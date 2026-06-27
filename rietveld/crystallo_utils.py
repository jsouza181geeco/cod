"""
Crystallographic utilities for QPA (Épico 11).

cell_volume  — triclinic general unit-cell volume [Å³]
molar_mass   — molar mass [g/mol] from a COD formula string

Both feed the Hill-Howard weight fraction W_k ∝ S_k·Z_k·M_k·V_k [HH87].
"""
import re

import numpy as np

# IUPAC standard atomic weights (g/mol), Z = 1..92.
# Conventional values (CIAAW 2021); abridged uncertainties dropped.
_ATOMIC_WEIGHTS = {
    'H': 1.008, 'He': 4.0026, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81,
    'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
    'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974,
    'S': 32.06, 'Cl': 35.45, 'Ar': 39.95, 'K': 39.098, 'Ca': 40.078,
    'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938,
    'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38,
    'Ga': 69.723, 'Ge': 72.630, 'As': 74.922, 'Se': 78.971, 'Br': 79.904,
    'Kr': 83.798, 'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224,
    'Nb': 92.906, 'Mo': 95.95, 'Tc': 98.0, 'Ru': 101.07, 'Rh': 102.91,
    'Pd': 106.42, 'Ag': 107.87, 'Cd': 112.41, 'In': 114.82, 'Sn': 118.71,
    'Sb': 121.76, 'Te': 127.60, 'I': 126.90, 'Xe': 131.29, 'Cs': 132.91,
    'Ba': 137.33, 'La': 138.91, 'Ce': 140.12, 'Pr': 140.91, 'Nd': 144.24,
    'Pm': 145.0, 'Sm': 150.36, 'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.93,
    'Dy': 162.50, 'Ho': 164.93, 'Er': 167.26, 'Tm': 168.93, 'Yb': 173.05,
    'Lu': 174.97, 'Hf': 178.49, 'Ta': 180.95, 'W': 183.84, 'Re': 186.21,
    'Os': 190.23, 'Ir': 192.22, 'Pt': 195.08, 'Au': 196.97, 'Hg': 200.59,
    'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.98, 'Po': 209.0, 'At': 210.0,
    'Rn': 222.0, 'Fr': 223.0, 'Ra': 226.0, 'Ac': 227.0, 'Th': 232.04,
    'Pa': 231.04, 'U': 238.03,
}

# token = element symbol (1 uppercase + optional lowercase) followed by
# an optional (possibly fractional) count. Two-letter elements match first,
# so 'Co' → cobalt; space-separated COD formulas ('Ca C O3') are unambiguous.
_TOKEN_RE = re.compile(r'([A-Z][a-z]?)(\d*\.?\d*)')


def cell_volume(a: float, b: float, c: float,
                alpha: float, beta: float, gamma: float) -> float:
    """
    Unit-cell volume [Å³], general triclinic formula:

        V = a·b·c·√(1 − cos²α − cos²β − cos²γ + 2cosα·cosβ·cosγ)

    Angles in degrees. Reduces to a·b·c for α=β=γ=90°.
    """
    ca, cb, cg = (np.cos(np.radians(x)) for x in (alpha, beta, gamma))
    radicand = 1.0 - ca**2 - cb**2 - cg**2 + 2.0 * ca * cb * cg
    return float(a * b * c * np.sqrt(max(radicand, 1e-12)))


def molar_mass(formula: str) -> float:
    """
    Molar mass [g/mol] from a COD formula string.

    Accepts 'C538 H654 Bi40', 'SiO2', 'Ca C O3', fractional counts ('Br0.8').
    Unknown element → ValueError (NOT a silent skip — silent skip would
    undercount M and corrupt QPA without warning) [crit. 7.7].
    """
    if not formula or not formula.strip():
        raise ValueError("Empty formula")

    M = 0.0
    found = False
    for elem, count in _TOKEN_RE.findall(formula):
        if elem not in _ATOMIC_WEIGHTS:
            raise ValueError(f"Unknown element '{elem}' in formula '{formula}'")
        n = float(count) if count and count != '.' else 1.0
        M += n * _ATOMIC_WEIGHTS[elem]
        found = True

    if not found:
        raise ValueError(f"No element tokens parsed from '{formula}'")
    return M


if __name__ == '__main__':
    # sanity checks (criteria 7.1, 7.2, 7.7)
    print("--- cell_volume ---")
    print(f"cubic 5,5,5,90,90,90 : {cell_volume(5, 5, 5, 90, 90, 90):.4f}  (esperado 125)")
    print(f"quartz a=b=4.913 c=5.405 90,90,120 : "
          f"{cell_volume(4.913, 4.913, 5.405, 90, 90, 120):.3f}  (esperado ~112.9)")
    print(f"triclinic 6,7,8,80,85,75 : {cell_volume(6, 7, 8, 80, 85, 75):.3f}")

    print("\n--- molar_mass ---")
    for f in ['SiO2', 'Ca C O3', 'Fe2 O3', 'C538 H654 Bi40 Mo2 N6 O98', 'Br0.8 C2']:
        print(f"{f:<30} = {molar_mass(f):.3f} g/mol")

    print("\n--- raise on unknown ---")
    try:
        molar_mass('Xx2 O3')
    except ValueError as e:
        print(f"OK raised: {e}")
