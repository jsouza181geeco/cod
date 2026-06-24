#!/usr/bin/env python3
"""
Worker puro para cálculo XRD via pymatgen.
Sem estado de DB — recebe job dict, retorna result dict.

Uso: importado por xrd_loader.py e xrd_scale_benchmark.py.
Nunca executar diretamente.
"""
import os
import time
import warnings
from pathlib import Path

try:
    import psutil
    _HAS_PSUTIL = True
except ImportError:
    _HAS_PSUTIL = False

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CIF_ROOT = PROJECT_ROOT / "cod_svn" / "cif"

# CuKα d_min≈1.09 Å → ~3.2 pontos recíprocos por Å³ de volume de célula.
# 20 000 Å³ ≈ 64k pontos → pode levar horas. Acima disso → erro imediato.
MAX_CELL_VOLUME_A3 = 20_000

# Globais inicializados 1x por processo worker via init_worker()
_calculator = None


def cif_path(cod_id: int) -> Path:
    s = str(cod_id).zfill(7)
    return CIF_ROOT / s[0] / s[1:3] / s[3:5] / f"{cod_id}.cif"


def init_worker(wavelength: str, two_theta_max: float = 90.0):
    """Chamado 1x por processo worker via ProcessPoolExecutor(initializer=...)."""
    global _calculator
    warnings.filterwarnings("ignore", category=UserWarning, module="pymatgen")

    # Typo no pymatgen instalado: linha 272 de xrd.py referencia
    # AbstractDiSffractionPatternCalculator (capital S) que nunca é definido.
    # Injeta alias antes de importar XRDCalculator para evitar NameError em get_pattern().
    import pymatgen.analysis.diffraction.xrd as _xrd_mod
    if not hasattr(_xrd_mod, "AbstractDiSffractionPatternCalculator"):
        from pymatgen.analysis.diffraction.core import AbstractDiffractionPatternCalculator
        _xrd_mod.AbstractDiSffractionPatternCalculator = AbstractDiffractionPatternCalculator

    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    _calculator = XRDCalculator(wavelength=wavelength)


def process_cif(job: dict) -> dict:
    """
    Processa 1 CIF. Nunca levanta exceção — erros capturados em status/error_msg.

    job keys obrigatórias: cod_id
    job keys opcionais:    sg_number, sg_symbol, sg_hall
    """
    from pymatgen.core import Structure

    cod_id = job["cod_id"]
    # Prefer DB-catalogued path (public.cod_files) over computed path
    # path from DB already includes "cif/" prefix → join from cod_svn/ (CIF_ROOT.parent)
    if "cif_path" in job and job["cif_path"]:
        path = CIF_ROOT.parent / job["cif_path"]
    else:
        path = cif_path(cod_id)
    rss_before = _rss_mb()

    if not path.exists():
        return _error(cod_id, "CIFNotFound", str(path))

    try:
        t0 = time.perf_counter()
        struct = Structure.from_file(str(path))
        t_load = time.perf_counter() - t0

        if struct.volume > MAX_CELL_VOLUME_A3:
            return _error(
                cod_id, "StructureTooLarge",
                f"cell volume {struct.volume:.0f} Å³ > {MAX_CELL_VOLUME_A3} — skipped",
            )

        t1 = time.perf_counter()
        pattern = _calculator.get_pattern(struct, scaled=True, two_theta_range=(0, 90))
        t_calc = time.perf_counter() - t1

        rss_after = _rss_mb()
        rss_delta = round((rss_after - rss_before), 1) if (rss_after and rss_before) else None

        lattice = struct.lattice
        reflections = _build_reflections(pattern)

        return {
            "status": "ok",
            "cod_id": cod_id,
            # cell params from pymatgen (canonical, symmetry-expanded)
            "a": lattice.a,
            "b": lattice.b,
            "c": lattice.c,
            "alpha": lattice.alpha,
            "beta": lattice.beta,
            "gamma": lattice.gamma,
            # sg from job dict (from DB) — avoids SpacegroupAnalyzer overhead
            "sg_number": job.get("sg_number") or 1,
            "sg_symbol": job.get("sg_symbol"),
            "sg_hall": job.get("sg_hall"),
            "formula": struct.formula,
            "n_reflections": len(pattern.x),
            "reflections": reflections,
            "t_load_s": round(t_load, 4),
            "t_calc_s": round(t_calc, 4),
            "rss_delta_mb": rss_delta,
        }

    except Exception as exc:
        return _error(cod_id, type(exc).__name__, str(exc))


def _build_reflections(pattern) -> list[dict]:
    """
    Converte DiffractionPattern → lista de dicts compatível com schema existente.
    Formato: {h, k, l, d_hkl, two_theta, multiplicity, intensity_rel}
    """
    out = []
    for two_theta, intensity, hkl_families, d_hkl in zip(
        pattern.x, pattern.y, pattern.hkls, pattern.d_hkls
    ):
        # hkl_families = [{"hkl": (h,k,l), "multiplicity": m}, ...]
        primary = hkl_families[0]
        hkl = primary["hkl"]
        # soma multiplicidade de todas as famílias equivalentes no mesmo pico
        mult = sum(f["multiplicity"] for f in hkl_families)
        out.append({
            "h": int(hkl[0]),
            "k": int(hkl[1]),
            "l": int(hkl[2]),
            "d_hkl": round(float(d_hkl), 5),
            "two_theta": round(float(two_theta), 4),
            "multiplicity": int(mult),
            "intensity_rel": round(float(intensity), 3),
        })
    return out


def _rss_mb() -> float | None:
    if not _HAS_PSUTIL:
        return None
    try:
        return psutil.Process(os.getpid()).memory_info().rss / 1024**2
    except Exception:
        return None


def _error(cod_id: int, exc_type: str, msg: str) -> dict:
    return {
        "status": "error",
        "cod_id": cod_id,
        "error_msg": f"{exc_type}: {msg}",
    }