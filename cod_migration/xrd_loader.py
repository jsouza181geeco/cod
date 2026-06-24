#!/usr/bin/env python3
"""
XRD Loader — wrapper pymatgen + bulk insert PostgreSQL.

Arquitetura:
  DB (data + reference_patterns) → fetch jobs em lotes por tier
  → 3 × ProcessPoolExecutor (SMALL/MEDIUM/LARGE)
  → mp.Queue
  → writer subprocess (asyncpg executemany + JSONL backup)

Tiers (por n_reflections do cálculo anterior, se disponível):
  SMALL  :   0 – 500  reflections → N_SMALL workers
  MEDIUM : 501 – 2000             → N_MEDIUM workers
  LARGE  : > 2000                 → N_LARGE workers, timeout mais rigoroso

Uso:
  # teste (50 CIFs por tier)
  python xrd_loader.py --limit 50

  # tier SMALL completo, 8 workers
  python xrd_loader.py --tiers SMALL --n-small 8

  # run completo (CONFIRMAR benchmark de RAM antes)
  python xrd_loader.py

  # ajuste manual de workers após benchmark de escalabilidade
  python xrd_loader.py --n-small 10 --n-medium 3 --n-large 1
"""
import argparse
import asyncio
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from pathlib import Path

from dotenv import load_dotenv

load_dotenv(Path(__file__).parent / ".env")

import asyncpg

from xrd_worker import init_worker, process_cif
from xrd_writer import writer_main

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Configuração de tiers
# ---------------------------------------------------------------------------
TIER_SMALL_MAX  = 500
TIER_MEDIUM_MAX = 2000

TIMEOUT_S = 600   # 10 min por CIF — outlier de 42 min → TIMEOUT, não trava pool

# Defaults — benchmark 2026-06-23: plateau em 4 workers (i5-12450H, 8GB RAM)
# 4→8 workers = 0.0% ganho (P-cores saturados). RAM: 170 MB/worker, 4×LARGE = 0.7 GB delta.
DEFAULT_N_SMALL  = 4
DEFAULT_N_MEDIUM = 4
DEFAULT_N_LARGE  = 4

WRITE_QUEUE_MAXSIZE = 2000  # backpressure: bloqueia workers se writer atrasado


# ---------------------------------------------------------------------------
# DB helpers
# ---------------------------------------------------------------------------

def _pg_params() -> dict:
    return {
        "host":     os.environ.get("PG_HOST", "localhost"),
        "port":     int(os.environ.get("PG_PORT", 5432)),
        "database": os.environ.get("PG_DB", "cod"),
        "user":     os.environ.get("PG_USER", "cod_admin"),
        "password": os.environ.get("PG_PASSWORD", ""),
    }


async def fetch_tier_jobs(tier: str, limit: int | None = None) -> list[dict]:
    """
    Busca jobs pendentes para o tier especificado.
    Fonte: tabela data (public schema) + reference_patterns para n_reflections_est.
    Exclui: já processados por este loader (wavelength_source='pymatgen') + failed_patterns.
    """
    lo, hi = {
        "SMALL":  (0,              TIER_SMALL_MAX),
        "MEDIUM": (TIER_SMALL_MAX + 1, TIER_MEDIUM_MAX),
        "LARGE":  (TIER_MEDIUM_MAX + 1, 9_999_999),
    }[tier]

    query = f"""
        SELECT
            d.file                             AS cod_id,
            d."sgNumber"                       AS sg_number,
            d.sg                               AS sg_symbol,
            d."sgHall"                         AS sg_hall,
            COALESCE(rp_est.n_reflections, 0)  AS n_reflections_est,
            cf.path                            AS cif_path
        FROM data d

        -- Apenas CIFs com path real catalogado (evita CIFNotFound)
        JOIN public.cod_files cf ON cf.cod_id = d.file

        -- Exclui já processados pela nova tabela pymatgen
        LEFT JOIN xrd_analysis.reference_patterns_pymatgen rp_done
            ON rp_done.cod_id = d.file

        -- Estimativa de n_reflections da tabela original (cálculo antigo)
        LEFT JOIN LATERAL (
            SELECT n_reflections
            FROM xrd_analysis.reference_patterns
            WHERE cod_id = d.file
            ORDER BY calculated_at DESC
            LIMIT 1
        ) rp_est ON TRUE

        -- Exclui falhas conhecidas do loader pymatgen
        LEFT JOIN xrd_analysis.failed_patterns_pymatgen fp ON fp.cod_id = d.file

        WHERE rp_done.cod_id IS NULL
          AND fp.cod_id IS NULL
          AND d.status IS NULL
          AND d.a IS NOT NULL
          AND d."sgNumber" IS NOT NULL
          AND COALESCE(rp_est.n_reflections, 0) >= {lo}
          AND COALESCE(rp_est.n_reflections, 0) <= {hi}

        ORDER BY COALESCE(rp_est.n_reflections, 0) ASC
        {'LIMIT ' + str(limit) if limit else ''}
    """

    conn = await asyncpg.connect(**_pg_params())
    try:
        rows = await conn.fetch(query)
        return [dict(row) for row in rows]
    finally:
        await conn.close()


# ---------------------------------------------------------------------------
# Processamento de tier
# ---------------------------------------------------------------------------

def _make_pool(n_workers: int, wavelength: str) -> ProcessPoolExecutor:
    return ProcessPoolExecutor(
        max_workers=n_workers,
        initializer=init_worker,
        initargs=(wavelength, 90.0),
    )


def _kill_pool(pool: ProcessPoolExecutor) -> None:
    """Termina workers travados e descarta pool — não espera shutdown."""
    try:
        for proc in pool._processes.values():
            proc.terminate()
    except Exception:
        pass
    pool.shutdown(wait=False)


def run_tier(
    tier: str,
    jobs: list[dict],
    n_workers: int,
    wavelength: str,
    result_queue: mp.Queue,
) -> tuple[int, int]:
    """
    Processa todos os jobs do tier em ProcessPoolExecutor via sliding window.

    Submete no máximo n_workers*2 futures por vez — t_submit reflete tempo real
    de início de execução, não tempo de fila. Evita falso-timeout em cascata
    quando todos os jobs são submetidos de uma vez com fila gigante.

    Retorna (n_ok, n_err+n_timeout).
    """
    if not jobs:
        return 0, 0

    print(
        f"\n  [{tier}] {len(jobs):,} jobs | {n_workers} workers | timeout={TIMEOUT_S}s",
        flush=True,
    )

    MAX_QUEUED = n_workers * 2  # janela deslizante: max futures em voo

    n_ok = n_err = 0
    processed = 0
    t_start = time.perf_counter()

    pool = _make_pool(n_workers, wavelength)
    futures: dict = {}
    t_submit: dict = {}
    pending: set = set()
    job_iter = iter(jobs)

    def _enqueue() -> None:
        """Preenche pipeline até MAX_QUEUED com próximos jobs do iterador."""
        while len(pending) < MAX_QUEUED:
            try:
                job = next(job_iter)
                f = pool.submit(process_cif, job)
                futures[f] = job
                t_submit[f] = time.perf_counter()
                pending.add(f)
            except StopIteration:
                break

    _enqueue()

    try:
        while pending:
            done, _ = wait(pending, timeout=5.0, return_when=FIRST_COMPLETED)

            for future in done:
                pending.discard(future)
                job = futures[future]
                try:
                    result = future.result()
                except Exception as exc:
                    result = {
                        "status": "error",
                        "cod_id": job["cod_id"],
                        "error_msg": f"{type(exc).__name__}: {exc}",
                    }
                result_queue.put(result)
                processed += 1
                if result["status"] == "ok":
                    n_ok += 1
                else:
                    n_err += 1

                if processed % 500 == 0 or processed == len(jobs):
                    elapsed = time.perf_counter() - t_start
                    rate = processed / elapsed if elapsed > 0 else 0
                    print(
                        f"  [{tier}] {processed:>7,}/{len(jobs):,} "
                        f"| ok={n_ok:,} err={n_err:,} "
                        f"| {rate:.2f} CIF/s",
                        flush=True,
                    )

            _enqueue()  # repõe vagas liberadas pelos done

            # Detecta futures que excederam TIMEOUT_S (t_submit ≈ início real)
            now = time.perf_counter()
            timed_out = [f for f in pending if now - t_submit[f] > TIMEOUT_S]
            if timed_out:
                for f in timed_out:
                    job = futures[f]
                    print(
                        f"  [{tier}] TIMEOUT cod_id={job['cod_id']} >{TIMEOUT_S}s — reiniciando pool",
                        flush=True,
                    )
                    result_queue.put({
                        "status": "error",
                        "cod_id": job["cod_id"],
                        "error_msg": f"TimeoutError: >{TIMEOUT_S}s",
                    })
                    n_err += 1
                    processed += 1
                    pending.discard(f)

                # Mata workers travados e reinicia pool
                _kill_pool(pool)
                pool = _make_pool(n_workers, wavelength)

                # Re-submete sobreviventes com t_submit fresco no novo pool
                surviving_jobs = [futures[f] for f in list(pending)]
                pending.clear()
                for job in surviving_jobs:
                    f = pool.submit(process_cif, job)
                    futures[f] = job
                    t_submit[f] = time.perf_counter()
                    pending.add(f)

                _enqueue()  # preenche vagas liberadas pelos timeouts

    finally:
        _kill_pool(pool)

    return n_ok, n_err


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="XRD Loader — pymatgen + PostgreSQL bulk insert",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--wavelength",  default="CuKa",
                   help="Radiação XRD (CuKa, MoKa, CoKa...)")
    p.add_argument("--n-small",  type=int, default=DEFAULT_N_SMALL,
                   help="Workers para tier SMALL (≤500 reflections)")
    p.add_argument("--n-medium", type=int, default=DEFAULT_N_MEDIUM,
                   help="Workers para tier MEDIUM (501–2000)")
    p.add_argument("--n-large",  type=int, default=DEFAULT_N_LARGE,
                   help="Workers para tier LARGE (>2000)")
    p.add_argument("--tiers", nargs="+", default=["SMALL", "MEDIUM", "LARGE"],
                   choices=["SMALL", "MEDIUM", "LARGE"],
                   help="Tiers a processar (ordem: SMALL → MEDIUM → LARGE)")
    p.add_argument("--limit", type=int, default=None,
                   help="Limitar N jobs por tier (para testes)")
    p.add_argument("--timeout", type=int, default=TIMEOUT_S,
                   help="Timeout em segundos por CIF")
    p.add_argument("--batch-size", type=int, default=1000,
                   help="Tamanho do batch do writer para executemany")
    p.add_argument("--backup", type=Path,
                   default=PROJECT_ROOT / "xrd_loader_backup.jsonl",
                   help="Arquivo JSONL de backup (durabilidade)")
    return p.parse_args()


def main():
    # UTF-8 no terminal Windows
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass

    args = parse_args()
    global TIMEOUT_S
    TIMEOUT_S = args.timeout

    # Resolve wavelength float
    try:
        from pymatgen.analysis.diffraction.xrd import WAVELENGTHS
        wavelength_val = WAVELENGTHS.get(args.wavelength, float(args.wavelength))
    except Exception:
        wavelength_val = 1.54184

    tier_workers = {
        "SMALL":  args.n_small,
        "MEDIUM": args.n_medium,
        "LARGE":  args.n_large,
    }

    print(f"  [loader] wavelength={args.wavelength} ({wavelength_val:.5f} Å)")
    print(f"  [loader] tiers={args.tiers} | workers={tier_workers}")
    print(f"  [loader] timeout={TIMEOUT_S}s | batch={args.batch_size}")
    print(f"  [loader] backup={args.backup}")

    # Inicia writer como subprocess separado
    result_queue: mp.Queue = mp.Queue(maxsize=WRITE_QUEUE_MAXSIZE)
    writer_proc = mp.Process(
        target=writer_main,
        args=(result_queue, args.backup, args.batch_size, wavelength_val, args.wavelength),
        daemon=False,
        name="xrd-writer",
    )
    writer_proc.start()
    print(f"  [loader] writer PID={writer_proc.pid}", flush=True)

    total_ok = total_err = 0
    t_global = time.perf_counter()

    try:
        for tier in args.tiers:
            print(f"\n  [loader] Buscando jobs tier={tier} no DB...", flush=True)
            jobs = asyncio.run(fetch_tier_jobs(tier, args.limit))
            print(f"  [loader] {tier}: {len(jobs):,} jobs pendentes", flush=True)

            if not jobs:
                print(f"  [loader] {tier}: nenhum job — pulando.", flush=True)
                continue

            n_ok, n_err = run_tier(
                tier=tier,
                jobs=jobs,
                n_workers=tier_workers[tier],
                wavelength=args.wavelength,
                result_queue=result_queue,
            )
            total_ok += n_ok
            total_err += n_err
            print(f"  [loader] {tier} DONE | ok={n_ok:,} err={n_err:,}", flush=True)

    finally:
        # Sinaliza writer para encerrar
        result_queue.put(None)
        writer_proc.join(timeout=300)
        if writer_proc.is_alive():
            print("  [loader] WARN: writer demorou >300s — terminando forçado.", flush=True)
            writer_proc.terminate()

    elapsed = time.perf_counter() - t_global
    print(
        f"\n  [loader] TUDO CONCLUÍDO | ok={total_ok:,} err={total_err:,} "
        f"| {elapsed/3600:.2f}h",
        flush=True,
    )


if __name__ == "__main__":
    # spawn obrigatório no Windows para multiprocessing funcionar corretamente
    mp.set_start_method("spawn", force=True)
    main()