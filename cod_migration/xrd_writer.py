#!/usr/bin/env python3
"""
Processo writer dedicado para xrd_loader.py.
Recebe resultados via mp.Queue, acumula em buffer e faz bulk insert via asyncpg.
Mantém backup JSONL em disco para durabilidade.

Nunca executar diretamente — iniciado como subprocess por xrd_loader.py.
"""
import asyncio
import json
import multiprocessing as mp
import os
import time
from pathlib import Path

from dotenv import load_dotenv

load_dotenv(Path(__file__).parent / ".env")

# SQL de inserção — ON CONFLICT DO NOTHING garante idempotência em reprocessamento
_INSERT_OK_SQL = """
    INSERT INTO xrd_analysis.reference_patterns_pymatgen (
        cod_id, a, b, c, alpha, beta, gamma,
        sg_number, sg_symbol, sg_hall,
        formula, wavelength, rad_symbol,
        reflections, n_reflections,
        t_load_s, t_calc_s, rss_delta_mb
    ) VALUES (
        $1, $2, $3, $4, $5, $6, $7,
        $8, $9, $10,
        $11, $12, $13,
        $14::jsonb, $15,
        $16, $17, $18
    )
    ON CONFLICT (cod_id, wavelength) DO NOTHING
"""

_INSERT_FAIL_SQL = """
    INSERT INTO xrd_analysis.failed_patterns_pymatgen (cod_id, status, error_msg)
    VALUES ($1, $2, $3)
    ON CONFLICT (cod_id) DO UPDATE
        SET status    = EXCLUDED.status,
            error_msg = EXCLUDED.error_msg,
            failed_at = now()
"""


def writer_main(
    queue_in: mp.Queue,
    backup_path: Path,
    batch_size: int = 1000,
    wavelength_val: float = 1.54184,
    rad_symbol: str = "CuKa",
):
    """Entry point para o subprocess writer."""
    try:
        sys_stdout_fix()
        asyncio.run(_writer_loop(queue_in, Path(backup_path), batch_size, wavelength_val, rad_symbol))
    except KeyboardInterrupt:
        pass


async def _writer_loop(
    queue_in: mp.Queue,
    backup_path: Path,
    batch_size: int,
    wavelength_val: float,
    rad_symbol: str,
):
    import asyncpg

    conn = await asyncpg.connect(**_pg_params())
    buf_ok: list[dict] = []
    buf_fail: list[dict] = []
    total_ok = total_fail = 0
    t_start = time.perf_counter()
    last_print = t_start
    last_flush = t_start  # flush por tempo — evita perda de dados em kill

    backup_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        while True:
            # queue.get() é bloqueante — roda em thread para não travar event loop
            item = await asyncio.to_thread(queue_in.get)

            if item is None:  # sentinel: loader encerrou
                break

            if item["status"] == "ok":
                buf_ok.append(item)
                # Backup imediato em disco (linha por linha)
                _append_jsonl(backup_path, item)
            else:
                status = "TIMEOUT" if "TimeoutError" in item.get("error_msg", "") else "FAILED"
                buf_fail.append({**item, "status": status})

            now = time.perf_counter()

            # Flush batch OK — por tamanho ou por tempo (60s)
            if len(buf_ok) >= batch_size or (buf_ok and now - last_flush >= 60.0):
                await _flush_ok(conn, buf_ok, wavelength_val, rad_symbol)
                total_ok += len(buf_ok)
                buf_ok.clear()
                last_flush = now

            # Flush batch FAIL (menor — não acumula muito)
            if len(buf_fail) >= 100:
                await _flush_fail(conn, buf_fail)
                total_fail += len(buf_fail)
                buf_fail.clear()

            # Progresso a cada 30s
            if now - last_print >= 30:
                _print_progress(total_ok, total_fail, t_start)
                last_print = now

        # Flush final (buffer remanescente)
        if buf_ok:
            await _flush_ok(conn, buf_ok, wavelength_val, rad_symbol)
            total_ok += len(buf_ok)
        if buf_fail:
            await _flush_fail(conn, buf_fail)
            total_fail += len(buf_fail)

        elapsed = time.perf_counter() - t_start
        print(
            f"\n  [writer] CONCLUÍDO | ok={total_ok:,} fail={total_fail:,} "
            f"| {elapsed/60:.1f} min",
            flush=True,
        )

    finally:
        await conn.close()


async def _flush_ok(conn, buf: list[dict], wavelength_val: float, rad_symbol: str):
    records = []
    for r in buf:
        records.append((
            r["cod_id"],
            r["a"], r["b"], r["c"],
            r["alpha"], r["beta"], r["gamma"],
            r["sg_number"],
            r.get("sg_symbol"),
            r.get("sg_hall"),
            r.get("formula"),
            wavelength_val,
            rad_symbol,
            json.dumps(r["reflections"]),
            r["n_reflections"],
            r.get("t_load_s"),
            r.get("t_calc_s"),
            r.get("rss_delta_mb"),
        ))
    await conn.executemany(_INSERT_OK_SQL, records)


async def _flush_fail(conn, buf: list[dict]):
    records = [
        (r["cod_id"], r.get("status", "FAILED"), r.get("error_msg", ""))
        for r in buf
    ]
    await conn.executemany(_INSERT_FAIL_SQL, records)


def _append_jsonl(path: Path, item: dict):
    """Backup mínimo em disco: 1 linha por resultado bem-sucedido."""
    entry = {
        "cod_id":        item["cod_id"],
        "n_reflections": item["n_reflections"],
        "formula":       item.get("formula"),
        "t_calc_s":      item.get("t_calc_s"),
    }
    with open(path, "a", encoding="utf-8") as f:
        f.write(json.dumps(entry) + "\n")


def _print_progress(total_ok: int, total_fail: int, t_start: float):
    elapsed = time.perf_counter() - t_start
    rate = total_ok / elapsed if elapsed > 0 else 0
    eta_h = (534606 - total_ok) / rate / 3600 if rate > 0 else float("inf")
    print(
        f"  [writer] ok={total_ok:>7,} fail={total_fail:>5,} "
        f"| {rate:.2f} CIF/s | ETA 534k: {eta_h:.1f}h",
        flush=True,
    )


def _pg_params() -> dict:
    return {
        "host":     os.environ.get("PG_HOST", "localhost"),
        "port":     int(os.environ.get("PG_PORT", 5432)),
        "database": os.environ.get("PG_DB", "cod"),
        "user":     os.environ.get("PG_USER", "cod_admin"),
        "password": os.environ.get("PG_PASSWORD", ""),
    }


def sys_stdout_fix():
    """UTF-8 no terminal Windows (silencia erro se não aplicável)."""
    import sys
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except (AttributeError, ValueError):
        pass