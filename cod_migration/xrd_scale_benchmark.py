#!/usr/bin/env python3
"""
Benchmark de escalabilidade e RAM para xrd_loader.py.

Mede:
  1. Throughput (CIF/s) com 1/2/4/8/N workers → detecta plateau de escalabilidade
  2. RSS memory delta nos CIFs mais pesados (1 worker) → define N_LARGE seguro

Resultados orientam:
  - DEFAULT_N_SMALL / N_MEDIUM / N_LARGE em xrd_loader.py
  - Se plateau em 8 workers → setar OMP_NUM_THREADS=1 e retestar

Uso:
  # teste padrão (200 CIFs, workers 1/2/4/8)
  python xrd_scale_benchmark.py

  # personalizado
  python xrd_scale_benchmark.py --sample 500 --workers 1 2 4 8 12 16

  # só RAM (sem teste de escalabilidade)
  python xrd_scale_benchmark.py --no-scale

  # só escalabilidade (sem RAM)
  python xrd_scale_benchmark.py --no-ram
"""
import argparse
import csv
import multiprocessing as mp
import os
import random
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures import TimeoutError as FuturesTimeout
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
CIF_ROOT = PROJECT_ROOT / "cod_svn" / "cif"
DEFAULT_CSV = PROJECT_ROOT / "CIF de exemplo por cada quantidade de picos.csv"

# CIFs pesados conhecidos do benchmark (cod_ids com >3000 reflections)
HEAVY_COD_IDS = [
    7201223,   # 4483 peaks, 2524s
    7122957,   # 3905 peaks, 651s
    1529376,   # 3724 peaks, 142s
    4307814,   # 3748 peaks, 430s
    1518256,   # 3832 peaks, 463s — 62304 pontos recíprocos
]


def _cif_path(cod_id: int) -> Path:
    s = str(cod_id).zfill(7)
    return CIF_ROOT / s[0] / s[1:3] / s[3:5] / f"{cod_id}.cif"


def load_sample_jobs(csv_path: Path, n: int, seed: int = 42) -> list[dict]:
    """Amostra aleatória de N jobs do benchmark CSV (n_reflections;cod_id)."""
    rows = []
    with open(csv_path, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter=";"):
            cod_id = int(row["min"])
            path = _cif_path(cod_id)
            if path.exists():
                rows.append({
                    "cod_id": cod_id,
                    "sg_number": 1,
                    "sg_symbol": None,
                    "sg_hall": None,
                    "n_reflections_est": int(row["n_reflections"]),
                })
    random.seed(seed)
    sample = random.sample(rows, min(n, len(rows)))
    print(f"  Amostra: {len(sample)} CIFs (de {len(rows)} com arquivo presente)")
    return sample


# ---------------------------------------------------------------------------
# Teste de escalabilidade
# ---------------------------------------------------------------------------

def run_scale_test(
    jobs: list[dict],
    n_workers: int,
    wavelength: str,
    timeout: int,
) -> dict:
    from xrd_worker import init_worker, process_cif

    n_ok = n_err = n_timeout = 0
    t_calc_times = []

    t_start = time.perf_counter()
    with ProcessPoolExecutor(
        max_workers=n_workers,
        initializer=init_worker,
        initargs=(wavelength, 90.0),
    ) as pool:
        futures = {pool.submit(process_cif, job): job for job in jobs}
        for future in as_completed(futures):
            try:
                r = future.result(timeout=timeout)
                if r["status"] == "ok":
                    n_ok += 1
                    t_calc_times.append(r["t_calc_s"])
                else:
                    n_err += 1
            except FuturesTimeout:
                n_timeout += 1

    elapsed = time.perf_counter() - t_start
    throughput = n_ok / elapsed if elapsed > 0 else 0
    return {
        "n_workers":      n_workers,
        "n_ok":           n_ok,
        "n_err":          n_err,
        "n_timeout":      n_timeout,
        "elapsed_s":      round(elapsed, 2),
        "throughput":     round(throughput, 3),
        "avg_t_calc_s":   round(sum(t_calc_times) / len(t_calc_times), 3) if t_calc_times else None,
    }


def scale_benchmark(jobs: list[dict], worker_counts: list[int], wavelength: str, timeout: int):
    print(f"\n  === Benchmark de escalabilidade ({len(jobs)} CIFs) ===")
    print(f"  {'workers':>8} | {'throughput':>12} | {'elapsed':>10} | {'speedup':>8} | {'ok':>6} | {'timeout':>7}")
    print("  " + "-" * 67)

    results = []
    baseline = None

    for n_workers in worker_counts:
        r = run_scale_test(jobs, n_workers, wavelength, timeout)
        results.append(r)

        if baseline is None and r["throughput"] > 0:
            baseline = r["throughput"]

        speedup = r["throughput"] / baseline if baseline else 0
        print(
            f"  {r['n_workers']:>8} | {r['throughput']:>9.3f} CIF/s"
            f" | {r['elapsed_s']:>8.1f}s"
            f" | {speedup:>7.2f}×"
            f" | {r['n_ok']:>6}"
            f" | {r['n_timeout']:>7}",
            flush=True,
        )

    # Análise de plateau
    print()
    if len(results) >= 3:
        last = results[-1]["throughput"]
        prev = results[-2]["throughput"]
        gain = (last - prev) / prev if prev > 0 else 0
        if gain < 0.10:
            opt = results[-2]["n_workers"]
            print(f"  PLATEAU detectado em {opt} workers (ganho adicional: {gain*100:.1f}%)")
            print(f"  → Recomendação: N_SMALL={opt} em xrd_loader.py")
            print(f"  → Se gain < 5%: testar OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 e retestar")
        else:
            print(f"  Sem plateau — considerar testar mais workers.")

    # Estimativa para 534k CIFs
    best = max(results, key=lambda r: r["throughput"])
    if best["throughput"] > 0:
        eta_h = 534606 / best["throughput"] / 3600
        print(f"\n  Melhor: {best['n_workers']} workers @ {best['throughput']:.3f} CIF/s")
        print(f"  ETA 534k (s/ LARGE outliers): {eta_h:.1f}h")

    return results


# ---------------------------------------------------------------------------
# Benchmark de RAM
# ---------------------------------------------------------------------------

def ram_benchmark(wavelength: str):
    """Mede RSS delta nos CIFs mais pesados (single worker no processo atual)."""
    try:
        import psutil
    except ImportError:
        print("\n  [RAM] psutil não instalado — instale com: pip install psutil")
        print("  [RAM] Pulando benchmark de RAM.")
        return

    print("\n  === Benchmark de RAM (worker único, CIFs pesados) ===")
    print(f"  {'cod_id':>10} | {'n_peaks':>8} | {'t_calc':>8} | {'rss_delta':>10} | formula")
    print("  " + "-" * 70)

    from xrd_worker import init_worker, process_cif
    init_worker(wavelength, 90.0)

    total_rss = []
    for cod_id in HEAVY_COD_IDS:
        path = _cif_path(cod_id)
        if not path.exists():
            print(f"  {cod_id:>10} | CIF não encontrado: {path}")
            continue

        job = {"cod_id": cod_id, "sg_number": 1, "sg_symbol": None, "sg_hall": None}
        r = process_cif(job)

        if r["status"] == "ok":
            rss = r.get("rss_delta_mb")
            rss_str = f"{rss:.0f} MB" if rss is not None else "N/A (psutil)"
            total_rss.append(rss or 0)
            print(
                f"  {cod_id:>10} | {r['n_reflections']:>8} "
                f"| {r['t_calc_s']:>6.1f}s "
                f"| {rss_str:>10} "
                f"| {r['formula'][:40]}",
                flush=True,
            )
        else:
            print(f"  {cod_id:>10} | ERRO: {r['error_msg']}")

    if total_rss:
        max_rss = max(total_rss)
        print(f"\n  RSS delta máximo por worker: {max_rss:.0f} MB")
        for n_workers in [4, 8, 12, 16]:
            est_gb = max_rss * n_workers / 1024
            flag = " ← CUIDADO" if est_gb > 16 else ""
            print(f"  {n_workers} workers LARGE → ~{est_gb:.1f} GB RAM{flag}")
        print()
        print(f"  → Recomendação N_LARGE: ", end="")
        # Assumindo 32 GB RAM disponível, margem de 60%
        ram_available_gb = 32
        safe_workers = max(1, int(ram_available_gb * 0.6 * 1024 / max_rss))
        print(f"{min(safe_workers, 4)} (baseado em {max_rss:.0f} MB/worker, {ram_available_gb}GB RAM)")
        print(f"  (ajuste --ram-gb se sua máquina tiver diferente)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Benchmark de escalabilidade e RAM para xrd_loader.py",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--sample",    type=int, default=200,
                   help="Nº de CIFs para teste de escalabilidade")
    p.add_argument("--workers",   nargs="+", type=int, default=[1, 2, 4, 8],
                   help="Contagens de workers a testar")
    p.add_argument("--wavelength", default="CuKa")
    p.add_argument("--timeout",   type=int, default=120,
                   help="Timeout por CIF no teste de escala (s)")
    p.add_argument("--seed",      type=int, default=42,
                   help="Seed aleatória para amostragem")
    p.add_argument("--csv",       type=Path, default=DEFAULT_CSV,
                   help="CSV com CIFs de referência")
    p.add_argument("--no-scale",  action="store_true",
                   help="Pular teste de escalabilidade")
    p.add_argument("--no-ram",    action="store_true",
                   help="Pular benchmark de RAM")
    p.add_argument("--ram-gb",    type=float, default=32.0,
                   help="RAM total disponível (GB) para cálculo de N_LARGE seguro")
    return p.parse_args()


def main():
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass

    args = parse_args()

    if not args.no_scale:
        if not args.csv.exists():
            print(f"  ERRO: CSV não encontrado: {args.csv}")
            sys.exit(1)
        jobs = load_sample_jobs(args.csv, args.sample, args.seed)
        scale_benchmark(jobs, args.workers, args.wavelength, args.timeout)

    if not args.no_ram:
        ram_benchmark(args.wavelength)

    print("\n  Concluído. Use os valores acima para ajustar xrd_loader.py.")


if __name__ == "__main__":
    mp.set_start_method("spawn", force=True)
    main()