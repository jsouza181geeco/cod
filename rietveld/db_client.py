from pathlib import Path

import psycopg2
import psycopg2.extras
from dotenv import dotenv_values

from models import StructureMetadata, CandidateInput

_HERE = Path(__file__).parent
_DEFAULT_ENV = _HERE.parent / 'cod_migration' / '.env'


def _load_env(env_path: str | Path | None) -> dict:
    p = Path(env_path) if env_path else _DEFAULT_ENV
    if not p.exists():
        raise FileNotFoundError(f".env not found: {p}")
    return dotenv_values(p)


class DBClient:
    def __init__(self, env_path: str | Path | None = None):
        cfg = _load_env(env_path)
        self._conn = psycopg2.connect(
            host=cfg['PG_HOST'],
            port=int(cfg.get('PG_PORT', 5432)),
            dbname=cfg['PG_DB'],
            user=cfg['PG_USER'],
            password=cfg['PG_PASSWORD'],
        )
        self._conn.autocommit = True

    def close(self):
        if self._conn and not self._conn.closed:
            self._conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    def fetch_reflections(self, cod_ids: list[int]) -> list[CandidateInput]:
        """Load CandidateInput (cod_id + reflections JSONB) from reference_patterns.

        Returns only entries with has_intensities=TRUE and CuKa wavelength.
        Used in DB-only mode after match_candidates_db narrows the search space.
        """
        if not cod_ids:
            return []
        ids = list(cod_ids)
        results: list[CandidateInput] = []
        with self._conn.cursor() as cur:
            cur.execute("""
                SELECT DISTINCT ON (cod_id)
                    cod_id, reflections
                FROM xrd_analysis.reference_patterns
                WHERE cod_id = ANY(%s)
                  AND has_intensities = TRUE
                  AND wavelength BETWEEN 1.535 AND 1.546
                ORDER BY cod_id, calculated_at DESC
            """, (ids,))
            for row in cur.fetchall():
                cod_id, reflections = row
                if reflections is None:
                    continue
                if isinstance(reflections, str):
                    import json
                    reflections = json.loads(reflections)
                results.append(CandidateInput(cod_id=cod_id, reflections=reflections))
        return results

    def fetch_metadata(self, cod_ids: list[int]) -> dict[int, StructureMetadata]:
        if not cod_ids:
            return {}

        ids = list(cod_ids)
        meta: dict[int, StructureMetadata] = {}

        with self._conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("""
                SELECT
                    file AS cod_id,
                    formula, calcformula, mineral, chemname,
                    sg AS sg_symbol, "sgNumber" AS sg_number, "sgHall" AS sg_hall,
                    a, b, c,
                    COALESCE(alpha, 90.0) AS alpha,
                    COALESCE(beta,  90.0) AS beta,
                    COALESCE(gamma, 90.0) AS gamma,
                    "Z", "Zprime" AS zprime,
                    wavelength, "radSymbol" AS rad_symbol,
                    method, status, flags,
                    authors, title, journal, year, doi
                FROM data
                WHERE file = ANY(%s)
            """, (ids,))
            for row in cur.fetchall():
                cid = row['cod_id']
                meta[cid] = StructureMetadata(
                    cod_id=cid,
                    formula=row.get('formula'),
                    mineral=row.get('mineral'),
                    chemname=row.get('chemname'),
                    sg_number=row.get('sg_number'),
                    sg_symbol=row.get('sg_symbol'),
                    a=row.get('a'),
                    b=row.get('b'),
                    c=row.get('c'),
                    alpha=row.get('alpha'),
                    beta=row.get('beta'),
                    gamma=row.get('gamma'),
                    Z=row.get('Z') if row.get('Z') is not None else None,
                    wavelength=row.get('wavelength'),
                    rad_symbol=row.get('rad_symbol'),
                    method=str(row['method']) if row.get('method') else None,
                    status=str(row['status']) if row.get('status') else None,
                    authors=row.get('authors'),
                    title=row.get('title'),
                    journal=row.get('journal'),
                    year=row.get('year'),
                    doi=row.get('doi'),
                )

            # supplement from xrd_analysis.reference_patterns (CuKa only)
            cur.execute("""
                SELECT DISTINCT ON (cod_id)
                    cod_id, has_intensities, wavelength AS calc_wavelength
                FROM xrd_analysis.reference_patterns
                WHERE cod_id = ANY(%s)
                  AND wavelength BETWEEN 1.535 AND 1.546
                ORDER BY cod_id, calculated_at DESC
            """, (ids,))
            for row in cur.fetchall():
                cid = row['cod_id']
                if cid in meta:
                    meta[cid].has_intensities = row.get('has_intensities')

        return meta
