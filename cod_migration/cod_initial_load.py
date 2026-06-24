#!/usr/bin/env python3
"""
COD Initial Load — converts all MySQL dumps in mysql/ to PostgreSQL.

Usage:
    python cod_initial_load.py [--dry-run] [--table TABLE]

Steps:
    1. Reads every *.sql file in SVN_LOCAL_PATH/mysql/
    2. Converts MySQL DDL → PostgreSQL DDL and executes it
    3. Loads corresponding *.txt data files via COPY (streaming, multi-GB safe)
    4. Records the SVN revision in _cod_sync_state

Prerequisites:
    - PostgreSQL database created and reachable (see .env)
    - SVN checkout present at SVN_LOCAL_PATH with at least mysql/ expanded
    - pip install -r requirements.txt
"""
import argparse
import sys
from pathlib import Path

from dotenv import load_dotenv
from tqdm import tqdm

load_dotenv()

from lib import svn_utils
from lib.pg_loader import (
    copy_from_file,
    execute_ddl,
    get_connection,
    init_sync_state,
    save_revision,
)
from lib.schema_converter import convert_file

SVN_LOCAL = Path(__import__('os').environ.get('SVN_LOCAL_PATH', '../cod_svn'))
MYSQL_DIR = SVN_LOCAL / 'mysql'


def parse_args():
    p = argparse.ArgumentParser(description='COD initial PostgreSQL load')
    p.add_argument('--dry-run', action='store_true',
                   help='Convert schemas and print DDL without executing')
    p.add_argument('--table', metavar='NAME',
                   help='Load only this table (e.g. data, journals)')
    p.add_argument('--schema-only', action='store_true',
                   help='Apply schemas but skip data loading')
    return p.parse_args()


def main():
    args = parse_args()

    if not MYSQL_DIR.exists():
        print(f"ERROR: {MYSQL_DIR} not found. Run SVN checkout first.", file=sys.stderr)
        sys.exit(1)

    sql_files = sorted(MYSQL_DIR.glob('*.sql'))
    if args.table:
        sql_files = [f for f in sql_files if f.stem == args.table]
        if not sql_files:
            print(f"ERROR: {args.table}.sql not found in {MYSQL_DIR}", file=sys.stderr)
            sys.exit(1)

    if args.dry_run:
        for sql_file in sql_files:
            print(f"\n{'='*60}\n-- {sql_file.name}\n{'='*60}")
            stmts = convert_file(sql_file)
            print('\n\n'.join(stmts))
        return

    conn = get_connection()
    init_sync_state(conn)

    try:
        revision = svn_utils.get_revision(SVN_LOCAL)
    except Exception as e:
        print(f"WARNING: Could not get SVN revision: {e}", file=sys.stderr)
        revision = 0

    print(f"SVN revision: {revision}")
    print(f"Loading {len(sql_files)} tables from {MYSQL_DIR}\n")

    loaded_tables = []
    errors = []

    for sql_file in tqdm(sql_files, desc='Tables', unit='table'):
        table = sql_file.stem
        txt_file = sql_file.with_suffix('.txt')

        tqdm.write(f'→ {table}')

        # Apply schema
        try:
            stmts = convert_file(sql_file)
            execute_ddl(conn, stmts)
            tqdm.write(f'  schema OK ({len(stmts)} statements)')
        except Exception as e:
            tqdm.write(f'  SCHEMA ERROR: {e}', file=sys.stderr)
            errors.append((table, 'schema', str(e)))
            conn.rollback()
            continue

        if args.schema_only:
            loaded_tables.append(table)
            continue

        # Load data
        if not txt_file.exists():
            tqdm.write(f'  no .txt file — schema only')
            loaded_tables.append(table)
            continue

        # Skip if table already has rows (safe to re-run)
        with conn.cursor() as cur:
            cur.execute(f'SELECT EXISTS (SELECT 1 FROM "{table}" LIMIT 1)')
            if cur.fetchone()[0]:
                tqdm.write(f'  skip (already loaded)')
                loaded_tables.append(table)
                continue

        try:
            rows = copy_from_file(conn, table, txt_file)
            tqdm.write(f'  data OK ({rows:,} rows)')
            loaded_tables.append(table)
        except Exception as e:
            tqdm.write(f'  DATA ERROR: {e}', file=sys.stderr)
            errors.append((table, 'data', str(e)))
            conn.rollback()

    save_revision(conn, revision, loaded_tables,
                  notes=f'initial load — {len(loaded_tables)} tables')
    conn.close()

    print(f'\n{"="*50}')
    print(f'Done. {len(loaded_tables)}/{len(sql_files)} tables loaded at revision {revision}.')
    if errors:
        print(f'\nErrors ({len(errors)}):')
        for table, stage, msg in errors:
            print(f'  {table} [{stage}]: {msg}')
        sys.exit(1)


if __name__ == '__main__':
    main()
