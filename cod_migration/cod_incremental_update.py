#!/usr/bin/env python3
"""
COD Incremental Update — svn update → upsert changed tables into PostgreSQL.

Usage:
    python cod_incremental_update.py [--dry-run]

How it works:
    1. Reads last SVN revision from _cod_sync_state
    2. Runs `svn update` on the working copy
    3. Detects which mysql/*.sql and *.txt files changed
    4. Re-applies schema (ALTER-safe via CREATE TABLE IF NOT EXISTS)
    5. Upserts changed data via temp table + ON CONFLICT DO UPDATE
    6. Records new revision in _cod_sync_state

Add new tables to TABLE_PKS if they need incremental upsert support.
Tables not listed are skipped for data updates (schema still applied).
"""
import argparse
import sys
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()

from lib import svn_utils
from lib.pg_loader import (
    execute_ddl,
    get_connection,
    get_last_revision,
    save_revision,
    upsert_from_file,
)
from lib.schema_converter import convert_file

SVN_LOCAL = Path(__import__('os').environ.get('SVN_LOCAL_PATH', '../cod_svn'))
MYSQL_DIR = SVN_LOCAL / 'mysql'

# Primary key per table — required for upsert.
# Tables omitted here get schema updates but no data upsert.
TABLE_PKS: dict[str, str] = {
    'data':             'file',
    'spacegroups':      'id',
    'journals':         'id',
    'publishers':       'id',
    'jsequences':       'id',
    'jaltnames':        'id',
    'smiles':           'file',
    'fingerprints':     'file',
    'numbers':          'id',
    'news':             'id',
    'databases':        'id',
    'relations':        'id',
    'successors':       'id',
    'rdf_relations':    'id',
}


def parse_args():
    p = argparse.ArgumentParser(description='COD incremental PostgreSQL update')
    p.add_argument('--dry-run', action='store_true',
                   help='Show what would change without modifying the database')
    return p.parse_args()


def main():
    args = parse_args()

    conn = get_connection()

    last_rev = get_last_revision(conn)
    if last_rev is None:
        print("No previous sync found. Run cod_initial_load.py first.", file=sys.stderr)
        sys.exit(1)

    print(f"Last sync: revision {last_rev}")

    if args.dry_run:
        print("[dry-run] Running svn update to detect changes...")

    new_rev, changed_files = svn_utils.update(SVN_LOCAL)

    if new_rev == last_rev:
        print("Already up to date.")
        conn.close()
        return

    print(f"Updating {last_rev} → {new_rev} ({len(changed_files)} file(s) changed)")

    # Determine which tables have changed mysql/ files
    changed_tables: set[str] = set()
    for f in changed_files:
        p = Path(f)
        if p.parent.name == 'mysql' and p.suffix in ('.sql', '.txt'):
            changed_tables.add(p.stem)

    if not changed_tables:
        print("No mysql/ changes detected. Recording new revision.")
        if not args.dry_run:
            save_revision(conn, new_rev, [], notes=f'update {last_rev}→{new_rev} (no mysql changes)')
        conn.close()
        return

    print(f"Changed tables: {', '.join(sorted(changed_tables))}\n")

    if args.dry_run:
        print("[dry-run] Would update:", ', '.join(sorted(changed_tables)))
        conn.close()
        return

    updated = []
    errors = []

    for table in sorted(changed_tables):
        sql_file = MYSQL_DIR / f'{table}.sql'
        txt_file = MYSQL_DIR / f'{table}.txt'
        pk = TABLE_PKS.get(table)

        print(f'→ {table} (pk={pk or "unknown"})')

        # Re-apply schema (CREATE TABLE IF NOT EXISTS is idempotent)
        if sql_file.exists():
            try:
                stmts = convert_file(sql_file)
                execute_ddl(conn, stmts)
                print(f'  schema OK')
            except Exception as e:
                print(f'  SCHEMA ERROR: {e}', file=sys.stderr)
                errors.append((table, 'schema', str(e)))
                conn.rollback()
                continue

        # Upsert data
        if txt_file.exists():
            if not pk:
                print(f'  WARNING: no PK configured — skipping data update. Add to TABLE_PKS.')
                continue
            try:
                rows = upsert_from_file(conn, table, txt_file, pk)
                print(f'  upserted {rows:,} rows')
                updated.append(table)
            except Exception as e:
                print(f'  UPSERT ERROR: {e}', file=sys.stderr)
                errors.append((table, 'data', str(e)))
                conn.rollback()
        else:
            print(f'  no .txt file — schema only')
            updated.append(table)

    save_revision(conn, new_rev, updated,
                  notes=f'incremental {last_rev}→{new_rev}')
    conn.close()

    print(f'\nDone. {len(updated)} tables updated at revision {new_rev}.')
    if errors:
        print(f'\nErrors ({len(errors)}):')
        for table, stage, msg in errors:
            print(f'  {table} [{stage}]: {msg}')
        sys.exit(1)


if __name__ == '__main__':
    main()
