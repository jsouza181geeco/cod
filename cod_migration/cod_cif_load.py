#!/usr/bin/env python3
"""
COD CIF Loader — parses atomic positions from .cif files and stores in
xrd_analysis.atomic_sites table.

CIF files contain the _atom_site_* loop with:
  label, type_symbol, fract_x, fract_y, fract_z, occupancy, U_iso, Wyckoff

These are NOT in the MySQL dumps — only in the .cif files under cod_svn/cif/.

CIF file path: {SVN_LOCAL}/cif/{cod_id[0]}/{cod_id[1:3]}/{cod_id}.cif
Example: COD 1010369 → cif/1/01/1010369.cif

Usage:
    python cod_cif_load.py --schema-only
    python cod_cif_load.py --limit 500
    python cod_cif_load.py --cod-ids 1010369 9000088
    python cod_cif_load.py                      # all entries that have CIF + no sites yet
    python cod_cif_load.py --reprocess          # reprocess even if already loaded

Prerequisites:
    - xrd_schema_setup.py --schema-only already run
    - cod_svn/cif/ directory present (SVN checkout includes cif/ subdir)
"""
import argparse
import os
import re
import sys
from pathlib import Path

from dotenv import load_dotenv
from tqdm import tqdm

load_dotenv(Path(__file__).parent / '.env')

import psycopg2
import psycopg2.extras

SVN_LOCAL = Path(os.environ.get('SVN_LOCAL_PATH', '../cod_svn'))
CIF_DIR   = SVN_LOCAL / 'cif'


# ---------------------------------------------------------------------------
# CIF file path resolution
# ---------------------------------------------------------------------------

def cif_path(cod_id: int) -> Path:
    """
    COD CIF directory structure: cif/{id[0]}/{id[1:3]}/{id[3:5]}/{id}.cif
    Example: 1010369 → cif/1/01/03/1010369.cif
             1000001 → cif/1/00/00/1000001.cif
    """
    s = str(cod_id)
    return CIF_DIR / s[0] / s[1:3] / s[3:5] / f'{s}.cif'


# ---------------------------------------------------------------------------
# Minimal CIF parser for _atom_site_* loop
# ---------------------------------------------------------------------------

_FLOAT_RE = re.compile(r'^[+-]?\d+\.?\d*(?:\(\d+\))?$')
_UNCERTAINTY_RE = re.compile(r'\(\d+\)')


def _parse_float(s: str) -> float | None:
    """
    Parse CIF numeric value, stripping standard uncertainty in parentheses.
    '0.2513(3)' → 0.2513 · '.' or '?' → None
    """
    s = s.strip()
    if s in ('.', '?', ''):
        return None
    s = _UNCERTAINTY_RE.sub('', s)
    try:
        return float(s)
    except ValueError:
        return None


def _parse_str(s: str) -> str | None:
    s = s.strip().strip("'\"")
    return None if s in ('.', '?', '') else s


def _tokenize_cif(text: str) -> list[str]:
    """
    Split CIF text into tokens, handling:
    - quoted strings ('...' or "...")
    - semicolon-delimited text blocks (;...;)
    - regular whitespace-separated tokens
    """
    tokens = []
    i = 0
    n = len(text)
    while i < n:
        ch = text[i]

        # Skip whitespace
        if ch in ' \t\r\n':
            i += 1
            continue

        # CIF comment
        if ch == '#':
            while i < n and text[i] != '\n':
                i += 1
            continue

        # Semicolon text block (must be at start of line)
        if ch == ';' and (i == 0 or text[i-1] == '\n'):
            i += 1
            start = i
            while i < n:
                if text[i] == ';' and (i == 0 or text[i-1] == '\n'):
                    break
                i += 1
            tokens.append(text[start:i].strip())
            i += 1
            continue

        # Single-quoted string
        if ch == "'":
            i += 1
            start = i
            while i < n and not (text[i] == "'" and (i+1 >= n or text[i+1] in ' \t\r\n')):
                i += 1
            tokens.append(text[start:i])
            i += 1
            continue

        # Double-quoted string
        if ch == '"':
            i += 1
            start = i
            while i < n and text[i] != '"':
                i += 1
            tokens.append(text[start:i])
            i += 1
            continue

        # Regular token
        start = i
        while i < n and text[i] not in ' \t\r\n':
            i += 1
        tokens.append(text[start:i])

    return tokens


# CIF _atom_site_ fields we want to extract
_ATOM_FIELDS = {
    '_atom_site_label':                'label',
    '_atom_site_type_symbol':          'type_symbol',
    '_atom_site_fract_x':              'fract_x',
    '_atom_site_fract_y':              'fract_y',
    '_atom_site_fract_z':              'fract_z',
    '_atom_site_occupancy':            'occupancy',
    '_atom_site_u_iso_or_equiv':       'u_iso',
    '_atom_site_b_iso_or_equiv':       'b_iso',    # convert: U = B/(8π²)
    '_atom_site_wyckoff_symbol':       'wyckoff_symbol',
    '_atom_site_site_symmetry_order':  'site_symmetry',
}

_B_TO_U = 1.0 / (8 * 3.141592653589793**2)


def parse_atom_sites(cif_text: str) -> list[dict]:
    """
    Parse _atom_site_* loop from CIF text.
    Returns list of dicts with keys matching xrd_analysis.atomic_sites columns.
    Returns [] if no atom_site loop found.
    """
    tokens = _tokenize_cif(cif_text)
    results = []
    i = 0
    n = len(tokens)

    while i < n:
        tok = tokens[i].lower()

        # Find loop_ containing _atom_site_ fields
        if tok != 'loop_':
            i += 1
            continue

        i += 1
        # Collect field names for this loop
        fields = []
        while i < n and tokens[i].startswith('_'):
            fields.append(tokens[i].lower())
            i += 1

        # Check if any field is an _atom_site_ field
        col_map = {}  # column index → internal key
        for idx, f in enumerate(fields):
            if f in _ATOM_FIELDS:
                col_map[idx] = _ATOM_FIELDS[f]

        if not col_map:
            # Not an atom_site loop — skip data values
            while i < n and not tokens[i].startswith('_') and tokens[i].lower() != 'loop_':
                i += 1
            continue

        # Read data rows
        n_cols = len(fields)
        while i < n:
            # Stop at next keyword
            t = tokens[i]
            if t.startswith('_') or t.lower() in ('loop_', 'data_', 'save_', 'stop_'):
                break

            row_tokens = tokens[i: i + n_cols]
            if len(row_tokens) < n_cols:
                break
            i += n_cols

            site = {}
            for col_idx, key in col_map.items():
                raw = row_tokens[col_idx] if col_idx < len(row_tokens) else '.'

                if key in ('fract_x', 'fract_y', 'fract_z', 'occupancy',
                           'u_iso', 'b_iso'):
                    site[key] = _parse_float(raw)
                else:
                    site[key] = _parse_str(raw)

            # Convert B_iso → U_iso if U_iso absent
            if site.get('u_iso') is None and site.get('b_iso') is not None:
                site['u_iso'] = site['b_iso'] * _B_TO_U
            site.pop('b_iso', None)

            if site.get('label'):
                results.append(site)

        # Done with this loop
        if results:
            return results   # Return first atom_site loop found

    return results


# ---------------------------------------------------------------------------
# DB helpers
# ---------------------------------------------------------------------------

def get_connection():
    return psycopg2.connect(
        host=os.environ.get('PG_HOST', 'localhost'),
        port=int(os.environ.get('PG_PORT', 5432)),
        dbname=os.environ.get('PG_DB', 'cod'),
        user=os.environ.get('PG_USER', 'cod_admin'),
        password=os.environ.get('PG_PASSWORD', ''),
    )


def ensure_schema(conn):
    """Create atomic_sites table if missing (idempotent)."""
    ddl = """
    CREATE SCHEMA IF NOT EXISTS xrd_analysis;
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
    CREATE INDEX IF NOT EXISTS as_cod_id ON xrd_analysis.atomic_sites (cod_id);
    CREATE INDEX IF NOT EXISTS as_type_symbol ON xrd_analysis.atomic_sites (type_symbol);
    """
    with conn.cursor() as cur:
        cur.execute(ddl)
    conn.commit()


def fetch_cod_ids(conn, cod_ids: list[int] | None,
                  limit: int | None, reprocess: bool) -> list[int]:
    """Fetch COD IDs to process from the data table."""
    filters = ["status IS NULL", "a IS NOT NULL"]
    params: list = []

    if cod_ids:
        filters.append('file = ANY(%s)')
        params.append(cod_ids)

    if not reprocess:
        filters.append("""
            file NOT IN (
                SELECT DISTINCT cod_id FROM xrd_analysis.atomic_sites
            )
        """)

    where = ' AND '.join(filters)
    limit_clause = f'LIMIT {limit}' if limit else ''

    query = f'SELECT file FROM data WHERE {where} ORDER BY file {limit_clause}'
    with conn.cursor() as cur:
        cur.execute(query, params or None)
        return [row[0] for row in cur.fetchall()]


def upsert_sites(conn, cod_id: int, sites: list[dict]):
    """Insert/replace all atom sites for a COD entry."""
    # Deduplicate by label — some CIF files have duplicate _atom_site_label
    # (disordered sites split into A/B components but sharing same label).
    # ON CONFLICT DO UPDATE cannot resolve two rows in the same batch with
    # the same constrained key — deduplicate in Python first, keep first hit.
    seen: set[str] = set()
    unique_sites: list[dict] = []
    for s in sites:
        lbl = s.get('label')
        if lbl and lbl not in seen:
            seen.add(lbl)
            unique_sites.append(s)
    sites = unique_sites

    # DELETE first → plain INSERT (no conflict possible from existing rows)
    with conn.cursor() as cur:
        cur.execute(
            'DELETE FROM xrd_analysis.atomic_sites WHERE cod_id = %s',
            (cod_id,)
        )
        psycopg2.extras.execute_values(
            cur,
            """
            INSERT INTO xrd_analysis.atomic_sites
                (cod_id, label, type_symbol,
                 fract_x, fract_y, fract_z,
                 occupancy, u_iso_or_equiv,
                 wyckoff_symbol, site_symmetry)
            VALUES %s
            ON CONFLICT (cod_id, label) DO NOTHING
            """,
            [
                (
                    cod_id,
                    s.get('label'),
                    s.get('type_symbol'),
                    s.get('fract_x'),
                    s.get('fract_y'),
                    s.get('fract_z'),
                    s.get('occupancy', 1.0),
                    s.get('u_iso'),
                    s.get('wyckoff_symbol'),
                    s.get('site_symmetry'),
                )
                for s in sites
            ],
        )
    conn.commit()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Load atomic positions from COD CIF files into xrd_analysis.atomic_sites.'
    )
    p.add_argument('--schema-only', action='store_true',
                   help='Create table only, skip CIF loading.')
    p.add_argument('--cod-ids', type=int, nargs='+', metavar='ID',
                   help='Process only these COD IDs.')
    p.add_argument('--limit', type=int, metavar='N',
                   help='Process at most N entries (testing).')
    p.add_argument('--reprocess', action='store_true',
                   help='Reprocess even if atomic_sites already loaded for entry.')
    return p.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    if not CIF_DIR.exists():
        print(f'ERROR: CIF directory not found: {CIF_DIR}', file=sys.stderr)
        print('SVN checkout must include cif/ subdir. Run:')
        print(f'  svn update --set-depth infinity {SVN_LOCAL}/cif')
        sys.exit(1)

    conn = get_connection()
    ensure_schema(conn)
    print('Schema xrd_analysis.atomic_sites OK.')

    if args.schema_only:
        conn.close()
        print('Done (schema only).')
        return

    print('Fetching COD IDs to process...')
    cod_ids = fetch_cod_ids(conn, args.cod_ids, args.limit, args.reprocess)
    print(f'{len(cod_ids)} entries to process.')

    if not cod_ids:
        print('Nothing to do.')
        conn.close()
        return

    ok = 0
    no_cif = 0
    no_sites = 0
    errors = 0

    for cod_id in tqdm(cod_ids, desc='CIF files', unit='cif'):
        path = cif_path(cod_id)

        if not path.exists():
            no_cif += 1
            continue

        try:
            text = path.read_text(encoding='utf-8', errors='replace')
            sites = parse_atom_sites(text)

            if not sites:
                no_sites += 1
                tqdm.write(f'  COD {cod_id}: no _atom_site_ loop in CIF')
                continue

            upsert_sites(conn, cod_id, sites)
            ok += 1

        except Exception as e:
            errors += 1
            conn.rollback()
            tqdm.write(f'  COD {cod_id}: ERROR — {e}', file=sys.stderr)

    conn.close()

    print(f'\nDone.')
    print(f'  Loaded:       {ok:>8,}')
    print(f'  No CIF file:  {no_cif:>8,}')
    print(f'  No atom loop: {no_sites:>8,}')
    print(f'  Errors:       {errors:>8,}')

    print("""
Example queries:
  -- All atoms in a structure
  SELECT label, type_symbol, fract_x, fract_y, fract_z, occupancy
  FROM xrd_analysis.atomic_sites
  WHERE cod_id = 1010369
  ORDER BY label;

  -- Structures with Fe atoms in specific Wyckoff position
  SELECT DISTINCT cod_id
  FROM xrd_analysis.atomic_sites
  WHERE type_symbol = 'Fe' AND wyckoff_symbol = '4a';

  -- Count atoms per element per structure (sum formula verification)
  SELECT cod_id, type_symbol, COUNT(*) AS n_sites, SUM(occupancy) AS total_occupancy
  FROM xrd_analysis.atomic_sites
  WHERE cod_id = 1010369
  GROUP BY cod_id, type_symbol
  ORDER BY type_symbol;
""")


if __name__ == '__main__':
    main()
