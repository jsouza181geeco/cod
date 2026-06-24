"""Convert MySQL/MariaDB DDL dumps to PostgreSQL DDL."""
import re
from pathlib import Path


# Patterns that match the FULL extracted type string (used with re.fullmatch)
_TYPE_MAP = [
    (r'mediumint(?:\s*\(\s*\d+\s*\))?(?:\s+unsigned)?', 'INTEGER'),
    (r'tinyint(?:\s*\(\s*\d+\s*\))?(?:\s+unsigned)?',   'SMALLINT'),
    (r'smallint(?:\s*\(\s*\d+\s*\))?(?:\s+unsigned)?',  'SMALLINT'),
    (r'bigint(?:\s*\(\s*\d+\s*\))?(?:\s+unsigned)?',    'BIGINT'),
    (r'int(?:\s*\(\s*\d+\s*\))?(?:\s+unsigned)?',       'INTEGER'),
    (r'double(?:\s+unsigned)?',                          'DOUBLE PRECISION'),
    (r'float(?:\s+unsigned)?',                           'REAL'),
    (r'datetime',                                        'TIMESTAMP'),
    (r'(?:long|medium|tiny)?text',                       'TEXT'),
    (r'(?:long|medium|tiny)?blob',                       'BYTEA'),
    (r'year(?:\s*\(\s*\d+\s*\))?',                      'SMALLINT'),
    (r'bool(?:ean)?',                                    'BOOLEAN'),
    # Pass-through (valid in PG): varchar(N), char(N), decimal(M,N), date, time
]


def convert_file(sql_path) -> list[str]:
    """Convert a MySQL .sql dump file to a list of PG SQL statements."""
    content = Path(sql_path).read_text(encoding='utf-8', errors='replace')
    return convert_ddl(content)


def convert_ddl(content: str) -> list[str]:
    """Convert MySQL DDL string → ordered list of PG SQL statements."""
    content = _strip_mysql_noise(content)
    statements = []
    for table_name, body in _find_create_tables(content):
        statements.extend(_convert_table(table_name, body))
    return statements


# ---------------------------------------------------------------------------
# Noise stripping
# ---------------------------------------------------------------------------

def _strip_mysql_noise(sql: str) -> str:
    sql = re.sub(r'/\*M!.*?\*/', '', sql, flags=re.DOTALL)        # MariaDB /*M!...*/
    sql = re.sub(r'/\*!\d+.*?\*/', '', sql, flags=re.DOTALL)      # MySQL /*!NNNNN...*/
    sql = re.sub(r'/\*.*?\*/', '', sql, flags=re.DOTALL)           # remaining block comments
    sql = re.sub(r'^\s*SET\s+@?\w.*?;\s*$', '', sql,
                 flags=re.MULTILINE | re.IGNORECASE)
    sql = re.sub(r'^\s*DROP\s+TABLE\s+IF\s+EXISTS.*?;\s*$', '', sql,
                 flags=re.MULTILINE | re.IGNORECASE)
    return sql


# ---------------------------------------------------------------------------
# CREATE TABLE block finder (depth-counter, handles nested parens in ENUMs)
# ---------------------------------------------------------------------------

def _find_create_tables(content: str) -> list[tuple[str, str]]:
    results = []
    pos = 0
    while True:
        m = re.search(r'CREATE\s+TABLE\s+`?(\w+)`?\s*\(', content[pos:], re.IGNORECASE)
        if not m:
            break
        table_name = m.group(1)
        body_start = pos + m.end()

        depth, j = 1, body_start
        while j < len(content) and depth > 0:
            if content[j] == '(':
                depth += 1
            elif content[j] == ')':
                depth -= 1
            j += 1

        body = content[body_start : j - 1]
        results.append((table_name, body))
        pos = j
    return results


# ---------------------------------------------------------------------------
# Table conversion
# ---------------------------------------------------------------------------

def _convert_table(table: str, body: str) -> list[str]:
    """Convert one CREATE TABLE body to a list of PG statements."""
    pre_stmts = []   # CREATE TYPE etc. — must run before CREATE TABLE
    pg_cols = []
    post_stmts = []  # CREATE INDEX etc. — run after CREATE TABLE
    pk_cols = None

    for col_def in _split_columns(body):
        col_def = col_def.strip()
        if not col_def:
            continue

        # PRIMARY KEY
        if m := re.match(r'PRIMARY\s+KEY\s*\((.+)\)', col_def, re.I):
            pk_cols = _clean_idx_cols(m.group(1))
            continue

        # UNIQUE KEY / UNIQUE INDEX
        if m := re.match(r'UNIQUE\s+(?:KEY|INDEX)\s+`?(\w+)`?\s*\((.+)\)', col_def, re.I):
            cols = _clean_idx_cols(m.group(2))
            post_stmts.append(
                f'CREATE UNIQUE INDEX IF NOT EXISTS "{m.group(1)}" ON "{table}" ({cols});'
            )
            continue

        # FULLTEXT KEY → GIN tsvector index
        if m := re.match(r'FULLTEXT\s+(?:KEY|INDEX)\s+`?(\w+)`?\s*\((.+)\)', col_def, re.I):
            raw = [c.strip().strip('`') for c in m.group(2).replace('`', '').split(',')]
            if len(raw) == 1:
                expr = f'to_tsvector(\'english\', coalesce("{raw[0]}", \'\'))'
            else:
                cat = " || ' ' || ".join(f'coalesce("{c}", \'\')' for c in raw)
                expr = f"to_tsvector('english', {cat})"
            post_stmts.append(
                f'CREATE INDEX IF NOT EXISTS "{m.group(1)}_fts" ON "{table}" USING GIN ({expr});'
            )
            continue

        # Regular KEY / INDEX
        if m := re.match(r'(?:KEY|INDEX)\s+`?(\w+)`?\s*\((.+)\)', col_def, re.I):
            cols = _clean_idx_cols(m.group(2))
            post_stmts.append(
                f'CREATE INDEX IF NOT EXISTS "{m.group(1)}" ON "{table}" ({cols});'
            )
            continue

        # Skip MySQL CONSTRAINT clauses (FK, CHECK) — MyISAM doesn't enforce them
        if re.match(r'CONSTRAINT\b', col_def, re.I):
            continue

        # Column definition: `name` type [modifiers]
        if m := re.match(r'`?(\w+)`?\s+(.*)', col_def, re.DOTALL | re.I):
            col_name = m.group(1)
            rest = m.group(2).strip()
            col_sql, col_pre = _convert_column(col_name, rest, table)
            pre_stmts.extend(col_pre)
            pg_cols.append(f'    {col_sql}')

    if pk_cols:
        pg_cols.append(f'    PRIMARY KEY ({pk_cols})')

    cols_block = ',\n'.join(pg_cols)
    result = pre_stmts + [
        f'CREATE TABLE IF NOT EXISTS "{table}" (\n{cols_block}\n);'
    ] + post_stmts
    return result


# ---------------------------------------------------------------------------
# Column definition conversion
# ---------------------------------------------------------------------------

def _convert_column(name: str, rest: str, table: str) -> tuple[str, list[str]]:
    type_str, after = _extract_type(rest)
    auto_inc = bool(re.search(r'\bAUTO_INCREMENT\b', after, re.I))
    modifiers = _clean_modifiers(after)

    pre_stmts = []

    # MySQL ENUM → PostgreSQL ENUM TYPE
    if em := re.match(r'enum\s*\((.+)\)\s*$', type_str, re.I | re.DOTALL):
        type_name = f"{table}_{name}_type"
        # Use DO block: create only if not exists, never DROP (CASCADE would nuke columns)
        pre_stmts = [
            f'DO $$ BEGIN\n'
            f'    CREATE TYPE "{type_name}" AS ENUM ({em.group(1)});\n'
            f'EXCEPTION WHEN duplicate_object THEN NULL;\n'
            f'END $$;',
        ]
        pg_type = f'"{type_name}"'

    # MySQL SET → TEXT (store comma-separated values as-is)
    elif re.match(r'set\s*\(', type_str, re.I):
        pg_type = 'TEXT'

    else:
        pg_type = _map_type(type_str)

    if auto_inc:
        pg_type += ' GENERATED ALWAYS AS IDENTITY'

    col_sql = f'"{name}" {pg_type}'
    if modifiers:
        col_sql += f' {modifiers}'

    return col_sql, pre_stmts


def _extract_type(rest: str) -> tuple[str, str]:
    """Extract MySQL type token (with parens if any) from start of column rest string."""
    rest = rest.lstrip()

    # Type that starts with parens: varchar(N), enum(...), set(...), int(N), etc.
    if re.match(r'\w+\s*\(', rest, re.I):
        depth, j = 0, 0
        while j < len(rest):
            if rest[j] == '(':
                depth += 1
            elif rest[j] == ')':
                depth -= 1
                if depth == 0:
                    end = j + 1
                    after = rest[end:].lstrip()
                    # Absorb trailing "unsigned" if present
                    if um := re.match(r'unsigned\b', after, re.I):
                        return rest[:end] + ' unsigned', after[um.end():]
                    return rest[:end], after
            j += 1

    # Simple type ± optional "unsigned"
    if m := re.match(r'(\w+)(?:\s+(unsigned)\b)?', rest, re.I):
        return m.group(0), rest[m.end():]

    return rest, ''


def _map_type(t: str) -> str:
    """Map MySQL type name to PostgreSQL equivalent."""
    t = t.strip()
    for pattern, pg in _TYPE_MAP:
        if re.fullmatch(pattern, t, re.I):
            return pg
    return t  # pass-through (varchar, char, decimal, date, time are valid in PG)


def _clean_modifiers(s: str) -> str:
    """Strip MySQL-specific column modifiers and fix PG incompatibilities."""
    s = re.sub(r'\bCHARACTER\s+SET\s+\S+', '', s, flags=re.I)
    s = re.sub(r'\bCOLLATE\s+\S+', '', s, flags=re.I)
    s = re.sub(r'\bAUTO_INCREMENT\b', '', s, flags=re.I)
    s = re.sub(r'\bunsigned\b', '', s, flags=re.I)
    # MySQL current_timestamp() → PG CURRENT_TIMESTAMP (no parens)
    s = re.sub(r'\bcurrent_timestamp\s*\(\s*\)', 'CURRENT_TIMESTAMP', s, flags=re.I)
    return re.sub(r'\s+', ' ', s).strip()


def _clean_idx_cols(cols: str) -> str:
    """Remove backticks and prefix lengths; quote each name to preserve case."""
    cols = cols.replace('`', '')
    cols = re.sub(r'\(\s*\d+\s*\)', '', cols)  # strip prefix lengths e.g. chemname(333)
    parts = [p.strip() for p in cols.split(',') if p.strip()]
    return ', '.join(f'"{p}"' for p in parts)


def _split_columns(body: str) -> list[str]:
    """Split CREATE TABLE body at top-level commas (respects nested parens)."""
    parts, depth, cur = [], 0, []
    for ch in body:
        if ch == '(':
            depth += 1
            cur.append(ch)
        elif ch == ')':
            depth -= 1
            cur.append(ch)
        elif ch == ',' and depth == 0:
            if s := ''.join(cur).strip():
                parts.append(s)
            cur = []
        else:
            cur.append(ch)
    if s := ''.join(cur).strip():
        parts.append(s)
    return parts
