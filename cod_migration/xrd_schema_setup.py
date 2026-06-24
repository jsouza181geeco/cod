#!/usr/bin/env python3
"""
XRD Analysis Schema — creates xrd_analysis schema and populates theoretical
XRD patterns (peak positions + intensities) from COD data.

Peak positions : Bragg's law from cell parameters + space group (data table)
Peak intensities: kinematic theory — structure factor F(hkl) from atomic sites
                  (xrd_analysis.atomic_sites, populated by cod_cif_load.py)

Usage:
    python xrd_schema_setup.py --schema-only
    python xrd_schema_setup.py --limit 100
    python xrd_schema_setup.py --cod-ids 1010369 9000088
    python xrd_schema_setup.py --wavelength 1.54056   # CuKa (default)
    python xrd_schema_setup.py --use-cod-wavelength
    python xrd_schema_setup.py --two-theta-max 120

Columns used from COD:
    data.a/b/c/alpha/beta/gamma  → d_hkl via crystal-system formula
    data.sgNumber                → crystal system + centering absences
    data.sg                      → lattice type (P/I/F/C/R) from 1st char
    data.sgHall                  → stored for reference
    data.wavelength              → d → 2θ via Bragg (falls back to CuKα)
    data.radSymbol               → display label
    atomic_sites.fract_x/y/z    → phase factor exp(2πi·hkl·xyz)
    atomic_sites.type_symbol     → atomic scattering factor f(sinθ/λ)
    atomic_sites.occupancy       → partial occupancy weighting
    atomic_sites.u_iso_or_equiv  → Debye-Waller thermal damping

Intensity formula (kinematic powder diffraction):
    I(hkl) = m·|F(hkl)|²·LP(2θ)
    F(hkl) = Σⱼ fⱼ(s)·occⱼ·exp(-8π²Uⱼs²)·exp(2πi(hxⱼ+kyⱼ+lzⱼ))
    LP(2θ)  = (1+cos²2θ)/(sin²θ·cosθ)        [Lorentz-polarization]
    s        = sinθ/λ = 1/(2d)                [Å⁻¹]
"""
import argparse
import json
import asyncio
import math
import os
import re
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from dotenv import load_dotenv

load_dotenv(Path(__file__).parent / '.env')

import asyncpg
import psycopg2          # used only for one-time synchronous schema DDL

# ---------------------------------------------------------------------------
# Known wavelengths (Å)
# ---------------------------------------------------------------------------

WAVELENGTHS = {
    'CuKa': 1.54056, 'CuKα': 1.54056,
    'MoKa': 0.71073, 'MoKα': 0.71073,
    'CoKa': 1.78897, 'CoKα': 1.78897,
    'CrKa': 2.28970, 'CrKα': 2.28970,
    'AgKa': 0.56086, 'AgKα': 0.56086,
    'FeKa': 1.93604, 'FeKα': 1.93604,
}

DEFAULT_WAVELENGTH = 1.54056
DEFAULT_RAD_SYMBOL = 'CuKα'

# ---------------------------------------------------------------------------
# Cromer-Mann atomic scattering factor coefficients
# Source: International Tables for X-ray Crystallography, Vol. IV (1974), Table 2.2B
# f(s) = Σᵢ aᵢ·exp(-bᵢ·s²) + c      s = sinθ/λ [Å⁻¹]
# Format: element → ([a1,a2,a3,a4], [b1,b2,b3,b4], c)
# ---------------------------------------------------------------------------

_CM: dict[str, tuple] = {
    'H':  ([0.489918, 0.262003, 0.196767, 0.049879],
           [20.6593,  7.74039, 49.5519,  2.20159],  0.001305),
    'He': ([0.873400, 0.630900, 0.311200, 0.178000],
           [9.10370,  3.35680, 22.9276,  0.982100],  0.006400),
    'Li': ([1.12820,  0.750800, 0.617500, 0.465300],
           [3.95460,  1.05240, 85.3905, 168.261],    0.037700),
    'Be': ([1.59190,  1.12780,  0.539100, 0.702900],
           [43.6427,  1.86230, 103.483,  0.542000],  0.038500),
    'B':  ([2.05450,  1.33260,  1.09790,  0.706800],
           [23.2185,  1.02100, 60.3498,  0.140300], -0.193200),
    'C':  ([2.31000,  1.02000,  1.58860,  0.865000],
           [20.8439, 10.2075,  0.568700, 51.6512],   0.215600),
    'N':  ([12.2126,  3.13220,  2.01250,  1.16630],
           [0.005700, 9.89330, 28.9975,  0.582600], -11.5290),
    'O':  ([3.04850,  2.28680,  1.54630,  0.867000],
           [13.2771,  5.70110,  0.323900, 32.9089],  0.250800),
    'F':  ([3.53920,  2.64120,  1.51700,  1.02430],
           [10.2825,  4.29440,  0.261500, 26.1476],  0.277600),
    'Ne': ([3.95530,  3.11250,  1.45460,  1.12510],
           [8.40420,  3.42620,  0.230600, 21.7184],  0.351500),
    'Na': ([4.76260,  3.17360,  1.26740,  1.11280],
           [3.28500,  8.84220,  0.313600, 129.424],  0.676000),
    'Mg': ([5.42040,  2.17350,  1.22690,  2.30730],
           [2.82750, 79.2611,  0.380800,  7.19370],  0.858400),
    'Al': ([6.42020,  1.90020,  1.59360,  1.96460],
           [3.03870,  0.742600, 31.5472, 85.0886],   1.11510),
    'Si': ([6.29150,  3.03530,  1.98910,  1.54100],
           [2.43860, 32.3337,  0.678500, 81.6937],   1.14070),
    'P':  ([6.43450,  4.17910,  1.78000,  1.49080],
           [1.90670, 27.1570,  0.526000, 68.1645],   1.11490),
    'S':  ([6.90530,  5.20340,  1.43790,  1.58630],
           [1.46790, 22.2151,  0.253600, 56.1720],   0.866900),
    'Cl': ([11.4604,  7.19640,  6.25560,  1.64550],
           [0.010400, 1.16620, 18.5194, 47.7784],   -9.55740),
    'Ar': ([7.48450,  6.77230,  0.653900, 1.64420],
           [0.907200, 14.8407, 43.8983, 33.3929],    1.44450),
    'K':  ([8.21860,  7.43980,  1.05190,  0.865900],
           [12.7949,  0.774800, 213.187, 41.6841],   1.42280),
    'Ca': ([8.62660,  7.38730,  1.58990,  1.02110],
           [10.4421,  0.659900, 85.7484, 178.437],   1.37510),
    'Sc': ([9.18900,  7.36790,  1.64090,  1.46800],
           [9.02130,  0.572900, 136.108, 51.3531],   1.33290),
    'Ti': ([9.75950,  7.35580,  1.69910,  1.90210],
           [7.85080,  0.500000, 35.6338, 116.105],   1.28070),
    'V':  ([10.2971,  7.35110,  2.07030,  2.05710],
           [6.86570,  0.438500, 26.8938, 102.478],   1.21990),
    'Cr': ([10.6406,  7.35370,  3.32400,  1.49220],
           [6.10380,  0.392000, 20.2626, 98.7399],   1.18320),
    'Mn': ([11.2819,  7.35730,  3.01930,  2.24410],
           [5.34090,  0.343200, 17.8674, 83.7543],   1.08960),
    'Fe': ([11.7695,  7.35730,  3.52220,  2.30450],
           [4.76110,  0.307200, 15.3535, 76.8805],   1.03690),
    'Co': ([12.2841,  7.34090,  4.00340,  2.34880],
           [4.27910,  0.278400, 13.5359, 71.1692],   1.01180),
    'Ni': ([12.8376,  7.29200,  4.44380,  2.38000],
           [3.87850,  0.256500, 12.1763, 66.3421],   1.03410),
    'Cu': ([13.3380,  7.16760,  5.61580,  1.67350],
           [3.58280,  0.247000, 11.3966, 64.8126],   1.19100),
    'Zn': ([14.0743,  7.03180,  5.16520,  2.41000],
           [3.26550,  0.233300, 10.3163, 58.7097],   1.30410),
    'Ga': ([15.2354,  6.70060,  4.35910,  2.96230],
           [3.06690,  0.241200, 10.7805, 61.4135],   1.71890),
    'Ge': ([16.0816,  6.37470,  3.70680,  3.68300],
           [2.85090,  0.251600, 11.4468, 54.7625],   2.13130),
    'As': ([16.6723,  6.07010,  3.43130,  4.27790],
           [2.63450,  0.264700, 12.9479, 47.7972],   2.53100),
    'Se': ([17.0006,  5.81960,  3.97310,  4.35430],
           [2.40980,  0.272600, 15.2372, 43.8163],   2.84090),
    'Br': ([17.1789,  5.23580,  5.63770,  3.98510],
           [2.17230, 16.5796,  0.260900, 41.4328],   2.95570),
    'Kr': ([17.3555,  6.72860,  5.54930,  3.53750],
           [1.93840, 16.5623,  0.226100, 39.3972],   2.82500),
    'Rb': ([17.1784,  9.64350,  5.13990,  1.52920],
           [1.78880, 17.3151,  0.274800, 164.934],   3.48730),
    'Sr': ([17.5663,  9.81840,  5.42200,  2.66940],
           [1.55640, 14.0988,  0.166400, 132.376],   2.50640),
    'Y':  ([17.7760, 10.2946,  5.72629,  3.26588],
           [1.40290, 12.8006,  0.125599, 104.354],   1.91213),
    'Zr': ([17.8765, 10.9480,  5.41732,  3.65721],
           [1.27618, 11.9160,  0.117622, 87.6627],   2.06929),
    'Nb': ([17.6142, 12.0144,  4.04183,  3.53346],
           [1.18865, 11.7660,  0.204785, 69.7957],   3.75591),
    'Mo': ([3.70250, 17.2356, 12.8876,  3.74290],
           [0.277200, 1.09580, 11.0040, 61.6584],    4.38700),
    'Tc': ([19.1301, 11.0948,  4.64901,  2.71263],
           [0.864132, 8.14487, 21.5707, 86.8472],    5.40428),
    'Ru': ([19.2674, 12.9182,  4.86337,  1.56756],
           [0.808520, 8.43467, 24.7997, 94.2928],    5.37874),
    'Rh': ([19.2957, 14.3501,  4.73425,  1.28918],
           [0.751536, 8.21758, 25.8749, 98.6062],    5.32800),
    'Pd': ([19.3319, 15.5017,  5.29537,  0.605844],
           [0.698655, 7.98929, 25.2052, 76.8986],    5.26593),
    'Ag': ([19.2808, 16.6885,  4.80450,  1.04630],
           [0.644600, 7.47260, 24.6605, 99.8156],    5.17900),
    'Cd': ([19.2214, 17.6444,  4.46100,  1.60290],
           [0.594600, 6.90890, 24.7008, 87.4825],    5.06940),
    'In': ([19.1624, 18.5596,  4.29480,  2.03960],
           [0.547600, 6.37760, 25.8499, 92.8029],    4.93910),
    'Sn': ([19.1889, 19.1005,  4.45850,  2.46630],
           [5.83030,  0.503100, 26.8909, 83.9571],   4.78200),
    'Sb': ([19.6418, 19.0455,  5.03710,  2.68270],
           [5.30340,  0.460700, 27.9074, 75.2825],   4.59090),
    'Te': ([19.9644, 19.0138,  6.14487,  2.52390],
           [4.81742,  0.420885, 28.5284, 70.8403],   4.35200),
    'I':  ([20.1472, 18.9949,  7.51380,  2.27350],
           [4.34700,  0.381400, 27.7660, 66.8776],   4.07120),
    'Xe': ([20.2933, 19.0298,  8.97670,  1.99020],
           [3.92820,  0.344000, 26.4659, 64.2658],   3.71180),
    'Cs': ([20.3892, 19.1062, 10.6620,  1.49530],
           [3.56900,  0.310700, 24.3879, 213.904],   3.33520),
    'Ba': ([20.3361, 19.2970, 10.8880,  2.69590],
           [3.21600,  0.275600, 20.2073, 167.202],   2.77310),
    'La': ([20.5780, 19.5990, 11.3727,  3.28719],
           [2.94817,  0.244475, 18.7726, 133.124],   2.14678),
    'Ce': ([21.1671, 19.7695, 11.8513,  3.33049],
           [2.81219,  0.226836, 17.6083, 127.113],   1.86264),
    'Pr': ([22.0440, 19.6697, 12.3856,  2.82428],
           [2.77393,  0.222087, 16.7669, 143.644],   2.05830),
    'Nd': ([22.6845, 19.6847, 12.7740,  2.85137],
           [2.66248,  0.210628, 15.8850, 137.903],   1.98486),
    'Pm': ([23.3405, 19.6095, 13.1235,  2.87516],
           [2.56270,  0.202088, 15.1009, 132.721],   2.02876),
    'Sm': ([24.0042, 19.4258, 13.4396,  2.89604],
           [2.47274,  0.196451, 14.3996, 128.007],   2.20963),
    'Eu': ([24.6274, 19.0886, 13.7603,  2.92270],
           [2.38790,  0.194200, 13.7546, 123.174],   2.57450),
    'Gd': ([25.0709, 19.0798, 13.8518,  3.54545],
           [2.25341,  0.181951, 12.9331, 101.398],   2.41960),
    'Tb': ([25.8976, 18.2185, 14.3167,  2.95354],
           [2.24256,  0.196143, 12.6648, 115.362],   3.58324),
    'Dy': ([26.5070, 17.6383, 14.5596,  2.96577],
           [2.18020,  0.202172, 12.1899, 111.874],   4.29728),
    'Ho': ([26.9049, 17.2940, 14.5583,  3.63837],
           [2.07051,  0.197940, 11.4407, 92.6566],   4.56796),
    'Er': ([27.6563, 16.4285, 14.9779,  2.98233],
           [2.07356,  0.223545, 11.3604, 105.703],   5.92046),
    'Tm': ([28.1819, 15.8851, 15.1542,  2.98706],
           [2.02859,  0.238849, 10.9975, 102.961],   6.75621),
    'Yb': ([28.6641, 15.4345, 15.3087,  2.98963],
           [1.98890,  0.257119, 10.6647, 100.417],   7.56672),
    'Lu': ([28.9476, 15.2208, 15.1000,  3.71601],
           [1.90182,  9.98519,  0.261033, 84.3298],  7.97628),
    'Hf': ([29.1440, 15.1726, 14.7586,  4.30013],
           [1.83262,  9.59990,  0.275116, 72.0290],  8.58154),
    'Ta': ([29.2024, 15.2293, 14.5135,  4.76492],
           [1.77333,  9.37046,  0.295977, 63.3644],  9.24354),
    'W':  ([29.0818, 15.4300, 14.4327,  5.11982],
           [1.72029,  9.22590,  0.321703, 57.0560],  9.88750),
    'Re': ([28.7621, 15.7189, 14.5564,  5.44174],
           [1.67191,  9.09227,  0.350500, 52.0861], 10.4720),
    'Os': ([28.1894, 16.1550, 14.9305,  5.67589],
           [1.62903,  8.97948,  0.382661, 48.1647], 11.0005),
    'Ir': ([27.3049, 16.7296, 15.6115,  5.83377],
           [1.59279,  8.86553,  0.417916, 45.0011], 11.4722),
    'Pt': ([27.0059, 17.7639, 15.7131,  5.78370],
           [1.51293,  8.81174,  0.424593, 38.6103], 11.6883),
    'Au': ([16.8819, 18.5913, 25.5582,  5.86000],
           [0.461100, 8.62160,  1.48260, 36.3956], 12.0658),
    'Hg': ([20.6809, 19.0417, 21.6575,  5.96760],
           [0.545000, 8.44840,  1.57290, 38.3246], 12.6089),
    'Tl': ([27.5446, 19.1584, 15.5380,  5.52593],
           [0.655150, 8.70751,  1.96347, 45.8149], 13.1746),
    'Pb': ([31.0617, 13.0637, 18.4420,  5.96960],
           [0.690200, 2.35760,  8.61800, 47.2579], 13.4118),
    'Bi': ([33.3689, 12.9510, 16.5877,  6.46920],
           [0.704000, 2.92380,  8.79370, 48.0093], 13.5782),
    'Th': ([35.5645, 23.4219, 12.7473,  4.80703],
           [0.563359, 3.46204, 17.8309, 99.1722], 13.4314),
    'U':  ([36.0228, 23.4128, 14.9491,  4.18800],
           [0.102300, 3.31920, 16.0927, 100.613], 13.3966),
}

_ELEM_RE = re.compile(r'^([A-Z][a-z]?)')


def _element_symbol(type_symbol: str) -> str:
    """
    Extract element symbol from CIF _atom_site_type_symbol.
    Handles: 'Fe', 'Fe2+', 'Fe3+', 'O2-', 'C.ar', 'N.am', 'Ca2+'
    """
    if not type_symbol:
        return ''
    m = _ELEM_RE.match(type_symbol.strip())
    return m.group(1) if m else ''


def atomic_scattering_factor(type_symbol: str, s: float) -> float:
    """
    f(s) for element from type_symbol at s = sinθ/λ [Å⁻¹].
    Falls back to atomic number approximation if element not in table.
    """
    elem = _element_symbol(type_symbol)
    if not elem:
        return 1.0

    if elem in _CM:
        a, b, c = _CM[elem]
        s2 = s * s
        return sum(ai * math.exp(-bi * s2) for ai, bi in zip(a, b)) + c

    # Fallback: rough approximation using atomic number
    _Z = {'Ac':89,'Am':95,'At':85,'Bk':97,'Cf':98,'Cm':96,'Es':99,
          'Fm':100,'Fr':87,'Md':101,'No':102,'Np':93,'Pa':91,'Pu':94,
          'Ra':88,'Rn':86,'Tc':43,'Xe':54}
    Z = _Z.get(elem, 1)
    return Z * math.exp(-3.0 * s * s)


# ---------------------------------------------------------------------------
# Crystal system from sgNumber
# ---------------------------------------------------------------------------

def crystal_system(sg_number: int) -> str:
    if 1 <= sg_number <= 2:      return 'triclinic'
    elif 3 <= sg_number <= 15:   return 'monoclinic'
    elif 16 <= sg_number <= 74:  return 'orthorhombic'
    elif 75 <= sg_number <= 142: return 'tetragonal'
    elif 143 <= sg_number <= 167: return 'trigonal'
    elif 168 <= sg_number <= 194: return 'hexagonal'
    elif 195 <= sg_number <= 230: return 'cubic'
    return 'unknown'


def lattice_type(hm_symbol: str) -> str:
    if not hm_symbol:
        return 'P'
    ch = hm_symbol.strip()[0].upper()
    return ch if ch in 'PIFCABR' else 'P'


def centering_allowed(h: int, k: int, l: int, lattice: str) -> bool:
    if lattice == 'P':   return True
    elif lattice == 'I': return (h + k + l) % 2 == 0
    elif lattice == 'F':
        parities = {h % 2, k % 2, l % 2}
        return len(parities) == 1
    elif lattice == 'C': return (h + k) % 2 == 0
    elif lattice == 'A': return (k + l) % 2 == 0
    elif lattice == 'B': return (h + l) % 2 == 0
    elif lattice == 'R': return (-h + k + l) % 3 == 0
    return True


# ---------------------------------------------------------------------------
# d-spacing per crystal system
# ---------------------------------------------------------------------------

def d_spacing(h: int, k: int, l: int,
              a: float, b: float, c: float,
              alpha_deg: float, beta_deg: float, gamma_deg: float,
              system: str) -> float | None:
    try:
        if system == 'cubic':
            inv_d2 = (h*h + k*k + l*l) / (a*a)

        elif system == 'tetragonal':
            inv_d2 = (h*h + k*k) / (a*a) + l*l / (c*c)

        elif system == 'orthorhombic':
            inv_d2 = h*h/(a*a) + k*k/(b*b) + l*l/(c*c)

        elif system in ('hexagonal', 'trigonal'):
            inv_d2 = (4/3) * (h*h + h*k + k*k) / (a*a) + l*l/(c*c)

        elif system == 'monoclinic':
            be = math.radians(beta_deg)
            sin_be = math.sin(be)
            cos_be = math.cos(be)
            if abs(sin_be) < 1e-10:
                return None
            inv_d2 = (1 / (sin_be * sin_be)) * (
                h*h/(a*a) +
                k*k * sin_be*sin_be / (b*b) +
                l*l/(c*c) -
                2*h*l*cos_be / (a*c)
            )

        elif system == 'triclinic':
            al = math.radians(alpha_deg)
            be = math.radians(beta_deg)
            ga = math.radians(gamma_deg)
            cos_al, cos_be, cos_ga = math.cos(al), math.cos(be), math.cos(ga)
            sin_al, sin_be, sin_ga = math.sin(al), math.sin(be), math.sin(ga)
            V2 = (a*b*c)**2 * (
                1 - cos_al**2 - cos_be**2 - cos_ga**2
                + 2*cos_al*cos_be*cos_ga
            )
            if V2 <= 0:
                return None
            V = math.sqrt(V2)
            S11 = (b*c*sin_al)**2
            S22 = (a*c*sin_be)**2
            S33 = (a*b*sin_ga)**2
            S12 = a*b*c*c * (cos_al*cos_be - cos_ga)
            S23 = a*a*b*c * (cos_be*cos_ga - cos_al)
            S13 = a*b*b*c * (cos_al*cos_ga - cos_be)
            inv_d2 = (
                S11*h*h + S22*k*k + S33*l*l
                + 2*S12*h*k + 2*S23*k*l + 2*S13*h*l
            ) / (V*V)
        else:
            return None

        if inv_d2 <= 1e-12:
            return None
        return 1.0 / math.sqrt(inv_d2)

    except (ValueError, ZeroDivisionError, OverflowError):
        return None


def bragg_two_theta(d: float, wavelength: float) -> float | None:
    ratio = wavelength / (2.0 * d)
    if ratio > 1.0:
        return None
    return math.degrees(2.0 * math.asin(ratio))


def lorentz_polarization(two_theta_deg: float) -> float:
    """LP = (1 + cos²2θ) / (sin²θ · cosθ) — standard powder diffraction factor."""
    theta = math.radians(two_theta_deg / 2.0)
    sin_t = math.sin(theta)
    cos_t = math.cos(theta)
    if abs(sin_t) < 1e-10 or abs(cos_t) < 1e-10:
        return 0.0
    cos2t = math.cos(math.radians(two_theta_deg))
    return (1.0 + cos2t * cos2t) / (sin_t * sin_t * cos_t)


# ---------------------------------------------------------------------------
# Structure factor |F(hkl)|²
# ---------------------------------------------------------------------------

_TWO_PI = 2.0 * math.pi
_8PI2   = 8.0 * math.pi * math.pi


def structure_factor_sq(h: int, k: int, l: int,
                        sites: list[dict], d: float) -> float:
    """
    |F(hkl)|² = |Σⱼ fⱼ(s)·occⱼ·DWⱼ·exp(2πi·φⱼ)|²
    s   = 1/(2d)
    DWⱼ = exp(-8π²·Uⱼ·s²)   [Debye-Waller]
    φⱼ  = h·xⱼ + k·yⱼ + l·zⱼ
    """
    s = 1.0 / (2.0 * d)
    s2 = s * s
    F_re = 0.0
    F_im = 0.0

    for site in sites:
        x = site.get('fract_x') or 0.0
        y = site.get('fract_y') or 0.0
        z = site.get('fract_z') or 0.0
        occ = site.get('occupancy') or 1.0
        u   = site.get('u_iso_or_equiv') or 0.0
        sym = site.get('type_symbol') or ''

        f  = atomic_scattering_factor(sym, s)
        dw = math.exp(-_8PI2 * u * s2) if u > 0 else 1.0
        phase = _TWO_PI * (h*x + k*y + l*z)

        amp = f * occ * dw
        F_re += amp * math.cos(phase)
        F_im += amp * math.sin(phase)

    return F_re*F_re + F_im*F_im


# ---------------------------------------------------------------------------
# Approximate multiplicity
# ---------------------------------------------------------------------------

def multiplicity(h: int, k: int, l: int, system: str) -> int:
    zeros   = sum(1 for x in (h, k, l) if x == 0)
    all_eq  = (h == k == l)
    two_eq  = (h == k or k == l or h == l) and not all_eq

    if system == 'cubic':
        if zeros == 0:  return 8 if all_eq else (24 if two_eq else 48)
        elif zeros == 1: return 12 if two_eq else 24
        else:            return 6
    elif system == 'tetragonal':
        if zeros == 0:   return 8 if h == k else 16
        elif l == 0:     return 4 if (h == 0 or k == 0 or h == k) else 8
        else:            return 4 if (h == 0 or k == 0) else 8
    elif system == 'orthorhombic':
        return {0: 8, 1: 4, 2: 2, 3: 1}[zeros]
    elif system == 'hexagonal':
        return 24 if zeros == 0 else 12
    elif system == 'trigonal':
        return 12 if zeros == 0 else 6
    elif system == 'monoclinic':
        return 2 if k == 0 else 4
    elif system == 'triclinic':
        return 2
    return 1


# ---------------------------------------------------------------------------
# Full pattern calculation
# ---------------------------------------------------------------------------

def calculate_pattern(
    a: float, b: float, c: float,
    alpha: float, beta: float, gamma: float,
    sg_number: int,
    sg_symbol: str,
    wavelength: float,
    sites: list[dict] | None = None,
    two_theta_min: float = 5.0,
    two_theta_max: float = 90.0,
    hkl_max: int = 15,
) -> list[dict]:
    """
    Calculate theoretical Bragg reflections.
    If sites provided → includes intensity_rel (0–100, normalised within pattern).
    If sites absent  → intensity_rel = null in output.

    Returns list sorted by two_theta:
        [{h, k, l, d_hkl, two_theta, multiplicity, F_sq, intensity_rel}]
    """
    system  = crystal_system(sg_number)
    lattice = lattice_type(sg_symbol)

    if system == 'unknown':
        return []

    has_sites = bool(sites)
    seen_d: dict[int, tuple] = {}
    reflections = []

    for h in range(0, hkl_max + 1):
        for k in range(-hkl_max, hkl_max + 1):
            for l in range(-hkl_max, hkl_max + 1):
                if h == 0 and k < 0:         continue
                if h == 0 and k == 0 and l <= 0: continue
                if not centering_allowed(h, k, l, lattice): continue

                d = d_spacing(h, k, l, a, b, c, alpha, beta, gamma, system)
                if d is None or d < 0.4:     continue

                two_theta = bragg_two_theta(d, wavelength)
                if two_theta is None:         continue
                if not (two_theta_min <= two_theta <= two_theta_max): continue

                d_key = round(d * 1000)
                if d_key in seen_d:           continue
                seen_d[d_key] = (h, k, l)

                m = multiplicity(h, k, l, system)

                # Intensity (if atomic sites available)
                if has_sites:
                    try:
                        F2  = structure_factor_sq(h, k, l, sites, d)
                        lp  = lorentz_polarization(two_theta)
                        raw = m * F2 * lp
                    except Exception:
                        raw = 0.0
                else:
                    F2  = None
                    raw = None

                reflections.append({
                    'h': h, 'k': k, 'l': l,
                    'd_hkl':     round(d, 5),
                    'two_theta': round(two_theta, 4),
                    'multiplicity': m,
                    'F_sq':      round(F2, 4) if F2 is not None else None,
                    'intensity_raw': round(raw, 4) if raw is not None else None,
                })

    reflections.sort(key=lambda r: r['two_theta'])

    # Normalise intensities 0–100 within pattern
    if has_sites:
        max_I = max((r['intensity_raw'] for r in reflections
                     if r['intensity_raw'] is not None), default=0.0)
        for r in reflections:
            raw = r.pop('intensity_raw')
            r['intensity_rel'] = (
                round(100.0 * raw / max_I, 2) if max_I > 0 and raw is not None else 0.0
            )
    else:
        for r in reflections:
            r.pop('intensity_raw', None)
            r['intensity_rel'] = None

    return reflections


# ---------------------------------------------------------------------------
# Schema DDL
# ---------------------------------------------------------------------------

SCHEMA_DDL = """
CREATE SCHEMA IF NOT EXISTS xrd_analysis;

COMMENT ON SCHEMA xrd_analysis IS
    'Theoretical XRD analysis tables derived from COD crystallographic data.';

CREATE TABLE IF NOT EXISTS xrd_analysis.reference_patterns (
    id              BIGINT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    cod_id          INTEGER NOT NULL REFERENCES data(file) ON DELETE CASCADE,

    a               DOUBLE PRECISION NOT NULL,
    b               DOUBLE PRECISION NOT NULL,
    c               DOUBLE PRECISION NOT NULL,
    alpha           REAL NOT NULL DEFAULT 90.0,
    beta            REAL NOT NULL DEFAULT 90.0,
    gamma           REAL NOT NULL DEFAULT 90.0,

    sg_number       SMALLINT NOT NULL,
    sg_symbol       VARCHAR(32),
    sg_hall         VARCHAR(64),

    wavelength      REAL NOT NULL,
    rad_symbol      VARCHAR(20),
    wavelength_source VARCHAR(20) NOT NULL,

    two_theta_min   REAL NOT NULL DEFAULT 5.0,
    two_theta_max   REAL NOT NULL DEFAULT 90.0,
    hkl_max         SMALLINT NOT NULL DEFAULT 15,
    has_intensities BOOLEAN NOT NULL DEFAULT FALSE,

    reflections     JSONB NOT NULL,
    n_reflections   INTEGER NOT NULL,

    calculated_at   TIMESTAMP NOT NULL DEFAULT NOW(),

    UNIQUE (cod_id, wavelength, two_theta_min, two_theta_max)
);

CREATE INDEX IF NOT EXISTS rp_cod_id
    ON xrd_analysis.reference_patterns (cod_id);
CREATE INDEX IF NOT EXISTS rp_sg_number
    ON xrd_analysis.reference_patterns (sg_number);
CREATE INDEX IF NOT EXISTS rp_wavelength
    ON xrd_analysis.reference_patterns (wavelength);
CREATE INDEX IF NOT EXISTS rp_has_intensities
    ON xrd_analysis.reference_patterns (has_intensities);
-- NOTE: GIN index on reflections created separately (GIN_CREATE_DDL).
-- It serialises concurrent JSONB writes (fastupdate pending-list lock),
-- so the bulk loader drops it before the parallel phase and rebuilds after.

COMMENT ON TABLE xrd_analysis.reference_patterns IS
    'Theoretical XRD patterns per COD structure. '
    'Each reflection: {h,k,l,d_hkl,two_theta,multiplicity,F_sq,intensity_rel}. '
    'intensity_rel is null when atomic sites were unavailable.';

COMMENT ON COLUMN xrd_analysis.reference_patterns.has_intensities IS
    'TRUE when F(hkl) was calculated from atomic_sites. '
    'FALSE when only peak positions are available (no CIF atom data).';

-- Idempotent: add column if table pre-existed from earlier --schema-only run
ALTER TABLE xrd_analysis.reference_patterns
    ADD COLUMN IF NOT EXISTS has_intensities BOOLEAN NOT NULL DEFAULT FALSE;

CREATE TABLE IF NOT EXISTS xrd_analysis.atomic_sites (
    id              BIGINT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    cod_id          INTEGER NOT NULL REFERENCES data(file) ON DELETE CASCADE,
    label           VARCHAR(20) NOT NULL,
    type_symbol     VARCHAR(16),
    fract_x         REAL,
    fract_y         REAL,
    fract_z         REAL,
    occupancy       REAL DEFAULT 1.0,
    u_iso_or_equiv  REAL,
    wyckoff_symbol  VARCHAR(8),
    site_symmetry   VARCHAR(16),
    UNIQUE (cod_id, label)
);

CREATE INDEX IF NOT EXISTS as_cod_id
    ON xrd_analysis.atomic_sites (cod_id);
CREATE INDEX IF NOT EXISTS as_type_symbol
    ON xrd_analysis.atomic_sites (type_symbol);
"""

# GIN index — managed separately so bulk loader can drop/rebuild it.
# fastupdate=off: no pending-list metapage lock → safe(r) under concurrency,
# but we rebuild after load regardless for speed.
GIN_DROP_DDL   = 'DROP INDEX IF EXISTS xrd_analysis.rp_reflections_gin;'
GIN_CREATE_DDL = (
    'CREATE INDEX IF NOT EXISTS rp_reflections_gin '
    'ON xrd_analysis.reference_patterns USING GIN (reflections) '
    'WITH (fastupdate = off);'
)


# ---------------------------------------------------------------------------
# DB config (built from env — works in spawned workers too, .env loaded at import)
# ---------------------------------------------------------------------------

def pg_params() -> dict:
    # kwargs (not DSN string) — avoids URL-parse breakage on special chars in password
    return {
        'host':     os.environ.get('PG_HOST', 'localhost'),
        'port':     int(os.environ.get('PG_PORT', 5432)),
        'database': os.environ.get('PG_DB', 'cod'),
        'user':     os.environ.get('PG_USER', 'cod_admin'),
        'password': os.environ.get('PG_PASSWORD', ''),
    }


def resolve_wavelength(row: dict, force_wl: float | None, use_cod_wl: bool):
    if force_wl is not None:
        return force_wl, f'{force_wl:.5f}Å', 'custom'
    if use_cod_wl and row.get('wavelength') and float(row['wavelength']) > 0.1:
        wl  = float(row['wavelength'])
        sym = row.get('radSymbol') or f'{wl:.5f}Å'
        return wl, sym, 'cod'
    return DEFAULT_WAVELENGTH, DEFAULT_RAD_SYMBOL, 'default_CuKa'


# Columns in INSERT order — shared by worker insert
_INSERT_COLS = (
    'cod_id', 'a', 'b', 'c', 'alpha', 'beta', 'gamma',
    'sg_number', 'sg_symbol', 'sg_hall',
    'wavelength', 'rad_symbol', 'wavelength_source',
    'two_theta_min', 'two_theta_max', 'hkl_max',
    'has_intensities', 'reflections', 'n_reflections',
)

_INSERT_SQL = f"""
    INSERT INTO xrd_analysis.reference_patterns
        ({', '.join(_INSERT_COLS)})
    VALUES ({', '.join(f'${i}' for i in range(1, len(_INSERT_COLS) + 1))})
    ON CONFLICT (cod_id, wavelength, two_theta_min, two_theta_max)
    DO UPDATE SET
        a = EXCLUDED.a, b = EXCLUDED.b, c = EXCLUDED.c,
        alpha = EXCLUDED.alpha, beta = EXCLUDED.beta, gamma = EXCLUDED.gamma,
        sg_number = EXCLUDED.sg_number,
        sg_symbol = EXCLUDED.sg_symbol,
        sg_hall   = EXCLUDED.sg_hall,
        rad_symbol = EXCLUDED.rad_symbol,
        wavelength_source = EXCLUDED.wavelength_source,
        hkl_max = EXCLUDED.hkl_max,
        has_intensities = EXCLUDED.has_intensities,
        reflections   = EXCLUDED.reflections,
        n_reflections = EXCLUDED.n_reflections,
        calculated_at = NOW()
"""

_STRUCT_SQL = """
    SELECT file,
           a, b, c,
           COALESCE(alpha, 90.0) AS alpha,
           COALESCE(beta,  90.0) AS beta,
           COALESCE(gamma, 90.0) AS gamma,
           "sgNumber" AS sg_number, sg, "sgHall" AS sg_hall,
           wavelength, "radSymbol" AS rad_symbol
    FROM data
    WHERE file = ANY($1::int[])
"""

_SITES_SQL = """
    SELECT cod_id, type_symbol, fract_x, fract_y, fract_z,
           occupancy, u_iso_or_equiv
    FROM xrd_analysis.atomic_sites
    WHERE cod_id = ANY($1::int[])
"""


# ---------------------------------------------------------------------------
# CPU side: build insert tuples for a batch of structure rows
# (pure function — no DB, runs inside worker process)
# ---------------------------------------------------------------------------

def build_records(struct_rows, sites_by_id, cfg):
    """
    struct_rows : list of asyncpg.Record (one per structure)
    sites_by_id : dict cod_id -> list[dict] atomic sites
    cfg         : dict of run options
    Returns (records, counts) where records = list of tuples for INSERT.
    """
    records = []
    n_ok = n_I = n_skip = n_err = 0

    for row in struct_rows:
        cod_id = row['file']
        wl, rad_symbol, wl_source = resolve_wavelength(
            {'wavelength': row['wavelength'], 'radSymbol': row['rad_symbol']},
            cfg['wavelength'], cfg['use_cod_wavelength'],
        )

        sites = None
        if not cfg['angles_only']:
            sites = sites_by_id.get(cod_id) or None

        try:
            reflections = calculate_pattern(
                a=float(row['a']), b=float(row['b']), c=float(row['c']),
                alpha=float(row['alpha']), beta=float(row['beta']),
                gamma=float(row['gamma']),
                sg_number=int(row['sg_number']),
                sg_symbol=row['sg'] or '',
                wavelength=wl,
                sites=sites,
                two_theta_min=cfg['two_theta_min'],
                two_theta_max=cfg['two_theta_max'],
                hkl_max=cfg['hkl_max'],
            )
        except Exception:
            n_err += 1
            continue

        if not reflections:
            n_skip += 1
            continue

        has_I = bool(sites) and reflections[0].get('intensity_rel') is not None

        records.append((
            cod_id,
            float(row['a']), float(row['b']), float(row['c']),
            float(row['alpha']), float(row['beta']), float(row['gamma']),
            int(row['sg_number']), row['sg'], row['sg_hall'],
            wl, rad_symbol, wl_source,
            cfg['two_theta_min'], cfg['two_theta_max'], cfg['hkl_max'],
            has_I, reflections, len(reflections),
        ))
        n_ok += 1
        if has_I:
            n_I += 1

    return records, (n_ok, n_I, n_skip, n_err)


# ---------------------------------------------------------------------------
# Async worker: one event loop + asyncpg pool per process, handles its partition
# ---------------------------------------------------------------------------

async def _init_conn(con):
    # JSONB codec: pass/return python objects directly (no manual json.dumps)
    await con.set_type_codec(
        'jsonb',
        encoder=json.dumps, decoder=json.loads,
        schema='pg_catalog',
    )


async def _aprocess_partition(cod_ids: list[int], cfg: dict) -> tuple:
    pool = await asyncpg.create_pool(
        **pg_params(),
        min_size=2, max_size=cfg['pool_size'],
        init=_init_conn,
        command_timeout=cfg['timeout'],
    )

    tot = [0, 0, 0, 0]   # ok, with_I, skipped, errors
    pending: list = []
    batch = cfg['batch']
    ins_chunk = cfg['insert_chunk']

    async def _insert_one(sub):
        # Short transaction → GIN index lock on reflections released fast.
        # Retry on transient TimeoutError under concurrent GIN contention.
        for attempt in range(3):
            try:
                async with pool.acquire() as con:
                    async with con.transaction():
                        await con.executemany(_INSERT_SQL, sub)
                return
            except (asyncio.TimeoutError, TimeoutError):
                if attempt == 2:
                    raise
                await asyncio.sleep(0.5 * (attempt + 1))

    async def flush(records):
        if not records:
            return
        # Split into small sub-chunks so no single transaction holds the
        # GIN index lock long enough to stall the other worker processes.
        for j in range(0, len(records), ins_chunk):
            await _insert_one(records[j:j + ins_chunk])

    try:
        for start in range(0, len(cod_ids), batch):
            chunk = cod_ids[start:start + batch]

            # Bulk fetch structures + their atomic sites (2 queries, kills N+1)
            async with pool.acquire() as con:
                struct_rows = await con.fetch(_STRUCT_SQL, chunk)
                if cfg['angles_only']:
                    sites_rows = []
                else:
                    sites_rows = await con.fetch(_SITES_SQL, chunk)

            sites_by_id: dict[int, list] = {}
            for s in sites_rows:
                sites_by_id.setdefault(s['cod_id'], []).append({
                    'type_symbol':    s['type_symbol'],
                    'fract_x':        s['fract_x'],
                    'fract_y':        s['fract_y'],
                    'fract_z':        s['fract_z'],
                    'occupancy':      s['occupancy'],
                    'u_iso_or_equiv': s['u_iso_or_equiv'],
                })

            # CPU-bound calc
            records, counts = build_records(struct_rows, sites_by_id, cfg)
            for i in range(4):
                tot[i] += counts[i]

            # Insert concurrently — overlaps with next batch fetch+calc
            pending.append(asyncio.create_task(flush(records)))
            if len(pending) >= cfg['pool_size']:
                await asyncio.gather(*pending)
                pending = []

        if pending:
            await asyncio.gather(*pending)
    finally:
        await pool.close()

    return tuple(tot)


def run_partition(cod_ids: list[int], cfg: dict) -> tuple:
    """Process-pool entry point — one asyncio.run per worker process."""
    return asyncio.run(_aprocess_partition(cod_ids, cfg))


# ---------------------------------------------------------------------------
# Main process: schema, fetch ids (resume-safe), partition, dispatch
# ---------------------------------------------------------------------------

def get_connection():
    return psycopg2.connect(
        host=os.environ.get('PG_HOST', 'localhost'),
        port=int(os.environ.get('PG_PORT', 5432)),
        dbname=os.environ.get('PG_DB', 'cod'),
        user=os.environ.get('PG_USER', 'cod_admin'),
        password=os.environ.get('PG_PASSWORD', ''),
    )


def fetch_ids_to_process(conn, cod_ids, limit, reprocess,
                         two_theta_min, two_theta_max) -> list[int]:
    filters = [
        'status IS NULL',
        'a IS NOT NULL', 'b IS NOT NULL', 'c IS NOT NULL',
        '"sgNumber" IS NOT NULL',
    ]
    params: list = []
    if cod_ids:
        filters.append('file = ANY(%s)')
        params.append(cod_ids)

    if not reprocess:
        # Resume-safe: skip entries already in reference_patterns for this 2θ range
        # (wavelength-agnostic — one pattern per structure/range is enough).
        filters.append("""
            NOT EXISTS (
                SELECT 1 FROM xrd_analysis.reference_patterns rp
                WHERE rp.cod_id = data.file
                  AND rp.two_theta_min = %s
                  AND rp.two_theta_max = %s
            )
        """)
        params.extend([two_theta_min, two_theta_max])

    where = ' AND '.join(filters)
    limit_clause = f'LIMIT {limit}' if limit else ''
    query = f'SELECT file FROM data WHERE {where} ORDER BY file {limit_clause}'

    with conn.cursor() as cur:
        cur.execute(query, params or None)
        return [r[0] for r in cur.fetchall()]


def partition(ids: list[int], n: int) -> list[list[int]]:
    """Round-robin split → balanced load (avoids one worker getting all heavy ids)."""
    parts: list[list[int]] = [[] for _ in range(n)]
    for i, cid in enumerate(ids):
        parts[i % n].append(cid)
    return [p for p in parts if p]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Calculate theoretical XRD patterns (positions + intensities) from COD. '
                    'Parallel: multiprocessing calc + async batched inserts.'
    )
    p.add_argument('--schema-only', action='store_true',
                   help='Create schema/tables only.')
    p.add_argument('--cod-ids', type=int, nargs='+', metavar='ID')
    p.add_argument('--limit', type=int, metavar='N')
    p.add_argument('--wavelength', type=float, metavar='Å',
                   help=f'Force wavelength (Å). Default: {DEFAULT_WAVELENGTH} CuKα. '
                        f'Common: CuKα=1.54056  MoKα=0.71073  CoKα=1.78897')
    p.add_argument('--use-cod-wavelength', action='store_true',
                   help='Use wavelength from COD entry when available.')
    p.add_argument('--two-theta-min', type=float, default=5.0)
    p.add_argument('--two-theta-max', type=float, default=90.0)
    p.add_argument('--hkl-max', type=int, default=15)
    p.add_argument('--angles-only', action='store_true',
                   help='Skip intensity calculation even if atomic_sites are available.')
    p.add_argument('--reprocess', action='store_true',
                   help='Recalculate even if already in reference_patterns. '
                        'Default: skip already-calculated entries (resume-safe).')
    # Performance knobs
    p.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) - 1),
                   help='Parallel processes (default: CPU count - 1).')
    p.add_argument('--batch', type=int, default=400,
                   help='Structures fetched per DB round-trip (default 400).')
    p.add_argument('--insert-chunk', type=int, default=100,
                   help='Rows per INSERT transaction (default 100). Smaller = '
                        'shorter GIN-index lock hold, less cross-worker stall.')
    p.add_argument('--pool-size', type=int, default=4,
                   help='asyncpg connections per worker process (default 4).')
    p.add_argument('--timeout', type=float, default=180.0,
                   help='Per-query command timeout seconds (default 180).')
    p.add_argument('--keep-gin', action='store_true',
                   help='Do NOT drop/rebuild the GIN index around the load. '
                        'Slower parallel inserts but index stays live throughout.')
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # Windows console defaults to cp1252 → reconfigure for θ/α/Å in output
    try:
        sys.stdout.reconfigure(encoding='utf-8')
        sys.stderr.reconfigure(encoding='utf-8')
    except (AttributeError, ValueError):
        pass

    args = parse_args()
    conn = get_connection()

    print('Creating xrd_analysis schema...')
    with conn.cursor() as cur:
        cur.execute(SCHEMA_DDL)
    conn.commit()
    print('Schema OK.')

    if args.schema_only:
        print('Building GIN index on reflections...')
        with conn.cursor() as cur:
            cur.execute(GIN_CREATE_DDL)
        conn.commit()
        conn.close()
        print('Done (schema only).')
        return

    print('Fetching COD ids to process...')
    ids = fetch_ids_to_process(
        conn, args.cod_ids, args.limit, args.reprocess,
        args.two_theta_min, args.two_theta_max,
    )
    conn.close()
    total = len(ids)
    print(f'{total:,} structures to process.')

    if not ids:
        print('Nothing to do (all calculated). Use --reprocess to force.')
        return

    # Drop GIN index before parallel load — concurrent JSONB writes serialise
    # on its pending-list lock, killing throughput. Rebuilt once after.
    if not args.keep_gin:
        print('Dropping GIN index for fast parallel load...')
        conn = get_connection()
        with conn.cursor() as cur:
            cur.execute(GIN_DROP_DDL)
        conn.commit()
        conn.close()

    parts = partition(ids, args.workers)
    cfg = {
        'wavelength':        args.wavelength,
        'use_cod_wavelength': args.use_cod_wavelength,
        'angles_only':       args.angles_only,
        'two_theta_min':     args.two_theta_min,
        'two_theta_max':     args.two_theta_max,
        'hkl_max':           args.hkl_max,
        'batch':             args.batch,
        'insert_chunk':      args.insert_chunk,
        'pool_size':         args.pool_size,
        'timeout':           args.timeout,
    }

    print(f'Workers: {len(parts)}  |  batch: {args.batch}  |  '
          f'pool/worker: {args.pool_size}\n')

    ok = ok_I = skipped = errors = 0
    t0 = time.time()

    with ProcessPoolExecutor(max_workers=len(parts)) as ex:
        futures = {ex.submit(run_partition, part, cfg): len(part)
                   for part in parts}
        done = 0
        for fut in as_completed(futures):
            try:
                a, b, c, d = fut.result()
                ok += a; ok_I += b; skipped += c; errors += d
            except Exception as e:
                errors += futures[fut]
                print(f'  worker FAILED: {type(e).__name__}: {e!r}', file=sys.stderr)
            done += 1
            elapsed = time.time() - t0
            rate = ok / elapsed if elapsed > 0 else 0
            print(f'  partition {done}/{len(parts)} done — '
                  f'{ok:,} ok ({rate:.0f}/s), {errors:,} err')

    elapsed = time.time() - t0
    print(f'\nCalc/insert done in {elapsed:.1f}s '
          f'({ok/elapsed if elapsed else 0:.0f} struct/s).')

    # Rebuild GIN index once (single bulk build ≫ faster than incremental)
    if not args.keep_gin:
        print('Rebuilding GIN index on reflections (single pass)...')
        ti = time.time()
        conn = get_connection()
        with conn.cursor() as cur:
            cur.execute(GIN_CREATE_DDL)
        conn.commit()
        conn.close()
        print(f'GIN index built in {time.time() - ti:.1f}s.')

    print(f'\nDone.')
    print(f'  Calculated:          {ok:>8,}')
    print(f'  With intensities:    {ok_I:>8,}')
    print(f'  No peaks in range:   {skipped:>8,}')
    print(f'  Errors:              {errors:>8,}')

    print("""
Example queries:
  -- Pattern with intensities for one structure
  SELECT cod_id, has_intensities, n_reflections,
         jsonb_array_elements(reflections) AS peak
  FROM xrd_analysis.reference_patterns
  WHERE cod_id = 1010369;

  -- Find strongest peaks near 2θ=33.1° (CuKα)
  SELECT p.cod_id, d.formula, d.mineral,
         (peak->>'two_theta')::real   AS two_theta,
         (peak->>'d_hkl')::real       AS d_hkl,
         (peak->>'intensity_rel')::real AS intensity,
         peak->>'h' || peak->>'k' || peak->>'l' AS hkl
  FROM xrd_analysis.reference_patterns p
  JOIN data d ON d.file = p.cod_id,
  LATERAL jsonb_array_elements(p.reflections) AS peak
  WHERE (peak->>'two_theta')::real BETWEEN 32.8 AND 33.4
    AND (peak->>'intensity_rel')::real > 10
    AND p.wavelength = 1.54056
    AND p.has_intensities = TRUE
  ORDER BY (peak->>'intensity_rel')::real DESC;
""")


if __name__ == '__main__':
    main()
