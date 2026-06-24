"""PostgreSQL helpers: connection, COPY bulk load, upsert, sync state."""
import os
from pathlib import Path

import psycopg2
import psycopg2.extras


def get_connection():
    """Connect using env vars: PG_HOST, PG_PORT, PG_DB, PG_USER, PG_PASSWORD."""
    return psycopg2.connect(
        host=os.environ['PG_HOST'],
        port=int(os.environ.get('PG_PORT', 5432)),
        dbname=os.environ['PG_DB'],
        user=os.environ['PG_USER'],
        password=os.environ['PG_PASSWORD'],
    )


def execute_ddl(conn, statements: list[str]) -> None:
    """Execute DDL statements in a single transaction."""
    with conn.cursor() as cur:
        for stmt in statements:
            if stmt.strip():
                cur.execute(stmt)
    conn.commit()


def copy_from_file(conn, table: str, txt_path: Path) -> int:
    """
    Bulk-load a MySQL-format tab-separated .txt file into a PG table using COPY.
    Handles \\N as NULL. Streams file — safe for multi-GB files.
    Returns number of rows loaded.
    """
    copy_sql = (
        f'COPY "{table}" FROM STDIN '
        f"WITH (FORMAT TEXT, NULL '\\N', DELIMITER E'\\t')"
    )
    with conn.cursor() as cur:
        with open(txt_path, 'rb') as f:
            cur.copy_expert(copy_sql, f)
        count = cur.rowcount
    conn.commit()
    return count


def upsert_from_file(conn, table: str, txt_path: Path, pk_col: str) -> int:
    """
    Upsert data from a .txt file into table.
    Loads into a temp table first, then INSERT ... ON CONFLICT DO UPDATE.
    Safe to call repeatedly — idempotent on pk_col.
    Returns number of rows upserted.
    """
    with conn.cursor() as cur:
        # Resolve column list from PG catalog
        cur.execute(
            """
            SELECT column_name FROM information_schema.columns
            WHERE table_schema = 'public' AND table_name = %s
            ORDER BY ordinal_position
            """,
            (table,),
        )
        columns = [row[0] for row in cur.fetchall()]
        if not columns:
            raise ValueError(f"Table '{table}' not found in public schema")

        non_pk = [c for c in columns if c != pk_col]

        # Create temp table mirroring the target
        temp = f'_tmp_{table}'
        cur.execute(
            f'CREATE TEMP TABLE "{temp}" (LIKE "{table}" INCLUDING DEFAULTS) ON COMMIT DROP;'
        )

        # Load into temp
        copy_sql = (
            f'COPY "{temp}" FROM STDIN '
            f"WITH (FORMAT TEXT, NULL '\\N', DELIMITER E'\\t')"
        )
        with open(txt_path, 'rb') as f:
            cur.copy_expert(copy_sql, f)

        # Upsert into main table
        if non_pk:
            update_set = ', '.join(f'"{c}" = EXCLUDED."{c}"' for c in non_pk)
            cur.execute(
                f'INSERT INTO "{table}" SELECT * FROM "{temp}" '
                f'ON CONFLICT ("{pk_col}") DO UPDATE SET {update_set};'
            )
        else:
            # Table has only PK column — just INSERT IGNORE
            cur.execute(
                f'INSERT INTO "{table}" SELECT * FROM "{temp}" '
                f'ON CONFLICT ("{pk_col}") DO NOTHING;'
            )

        count = cur.rowcount

    conn.commit()
    return count


# ---------------------------------------------------------------------------
# Sync state tracking
# ---------------------------------------------------------------------------

def init_sync_state(conn) -> None:
    """Create the _cod_sync_state bookkeeping table if it doesn't exist."""
    with conn.cursor() as cur:
        cur.execute(
            """
            CREATE TABLE IF NOT EXISTS _cod_sync_state (
                id          SERIAL PRIMARY KEY,
                svn_revision BIGINT NOT NULL,
                synced_at   TIMESTAMP DEFAULT NOW(),
                tables_done TEXT[],
                notes       TEXT
            );
            """
        )
    conn.commit()


def get_last_revision(conn) -> int | None:
    """Return the last recorded SVN revision, or None if never synced."""
    with conn.cursor() as cur:
        cur.execute(
            'SELECT svn_revision FROM _cod_sync_state ORDER BY synced_at DESC LIMIT 1;'
        )
        row = cur.fetchone()
    return row[0] if row else None


def save_revision(conn, revision: int, tables: list[str] = None, notes: str = None) -> None:
    """Record a completed sync event."""
    with conn.cursor() as cur:
        cur.execute(
            'INSERT INTO _cod_sync_state (svn_revision, tables_done, notes) VALUES (%s, %s, %s);',
            (revision, tables or [], notes),
        )
    conn.commit()
