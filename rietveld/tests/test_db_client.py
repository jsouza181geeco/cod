"""
Integration tests for db_client.DBClient — requires live PostgreSQL.

Run:  pytest -m integration tests/test_db_client.py -v
Skip: pytest -m "not integration" tests/   (default unit-test run)

T-030c tests require:
  - xrd_analysis.peak_fingerprints MV built with d_hkl schema
    (run migrations/create_peak_fingerprints.sql first)
"""
import pytest
from pathlib import Path

from db_client import DBClient
from data_loader import load_candidates_csv
from models import CandidateInput, StructureMetadata

DATA_DIR = Path(__file__).parent.parent
CSV_PATH = DATA_DIR / "data-1782394014136.csv"

BEST_COD_ID = 1569653  # C538 H654 Bi40 Mo2 N6 O98 — confirmed best candidate


@pytest.fixture(scope="module")
def db():
    """Single DB connection for the whole module (schema not modified)."""
    client = DBClient()
    yield client
    client.close()


@pytest.fixture(scope="module")
def cod_ids():
    return [c.cod_id for c in load_candidates_csv(CSV_PATH)]


@pytest.fixture(scope="module")
def meta(db, cod_ids):
    return db.fetch_metadata(cod_ids)


# --- T-020: basic fetch ---

@pytest.mark.integration
def test_fetch_returns_all_ids(meta, cod_ids):
    """All 20 cod_ids must be returned (all exist in DB)."""
    assert len(meta) == len(cod_ids), (
        f"Expected {len(cod_ids)} entries, got {len(meta)}"
    )


@pytest.mark.integration
def test_fetch_returns_structuremetadata(meta, cod_ids):
    """All values are StructureMetadata instances."""
    assert all(isinstance(v, StructureMetadata) for v in meta.values())


@pytest.mark.integration
def test_any_formula_present(meta):
    """At least one candidate has a formula string."""
    assert any(m.formula for m in meta.values())


@pytest.mark.integration
def test_any_has_intensities_not_none(meta):
    """At least one candidate has has_intensities set (from reference_patterns)."""
    assert any(m.has_intensities is not None for m in meta.values())


@pytest.mark.integration
def test_best_candidate_has_formula(meta):
    """Best candidate (1569653) has a formula in DB."""
    m = meta.get(BEST_COD_ID)
    assert m is not None, f"cod_id={BEST_COD_ID} not found in metadata"
    assert m.formula, f"Expected formula for cod_id={BEST_COD_ID}, got None"


@pytest.mark.integration
def test_best_candidate_has_intensities(meta):
    """Best candidate has_intensities=True (|F|² was calculated)."""
    m = meta[BEST_COD_ID]
    assert m.has_intensities is True, (
        f"Expected has_intensities=True for cod_id={BEST_COD_ID}, got {m.has_intensities}"
    )


@pytest.mark.integration
def test_sg_symbol_present(meta):
    """At least one candidate has sg_symbol (space group)."""
    assert any(m.sg_symbol for m in meta.values())


@pytest.mark.integration
def test_lattice_params_present(meta):
    """At least one candidate has lattice parameter a > 0."""
    assert any(m.a and m.a > 0 for m in meta.values())


@pytest.mark.integration
def test_empty_cod_ids(db):
    """Empty list → empty dict, no crash."""
    result = db.fetch_metadata([])
    assert result == {}


@pytest.mark.integration
def test_nonexistent_cod_id(db):
    """Nonexistent cod_id → not in result, no crash."""
    result = db.fetch_metadata([9999999])
    assert 9999999 not in result


@pytest.mark.integration
def test_context_manager():
    """DBClient as context manager closes connection on exit."""
    with DBClient() as client:
        result = client.fetch_metadata([BEST_COD_ID])
    assert BEST_COD_ID in result
    # connection should be closed after __exit__
    assert client._conn.closed


# --- T-021: end-to-end pipeline with DB metadata ---

@pytest.mark.integration
def test_pipeline_with_db_metadata():
    """T-021: run_pipeline with real DB → best candidate has formula."""
    from pipeline import run_pipeline
    xye_path = DATA_DIR / "synthetic_candidate17.xye"
    candidates = load_candidates_csv(CSV_PATH)

    with DBClient() as db_client:
        result = run_pipeline(xye_path, candidates, db_client=db_client)

    best = result.best()
    assert best.cod_id == BEST_COD_ID
    assert best.Rwp < 0.10
    assert best.metadata is not None
    assert best.metadata.formula is not None
    assert best.metadata.has_intensities is True
    assert best.metadata.sg_symbol is not None


# --- T-030c: DB-only mode (fetch_reflections + candidates_from_db) ---

@pytest.mark.integration
def test_fetch_reflections_single(db):
    """T-030c: fetch_reflections returns CandidateInput with reflections for known cod_id."""
    cands = db.fetch_reflections([BEST_COD_ID])
    assert len(cands) == 1, f"Expected 1, got {len(cands)}"
    assert isinstance(cands[0], CandidateInput)
    assert cands[0].cod_id == BEST_COD_ID
    assert len(cands[0].reflections) > 0, "reflections must be non-empty"
    assert all(isinstance(r, dict) for r in cands[0].reflections), \
        "reflections must be list[dict] (not raw JSON string)"


@pytest.mark.integration
def test_fetch_reflections_empty(db):
    """fetch_reflections([]) → []  (no crash)."""
    assert db.fetch_reflections([]) == []


@pytest.mark.integration
def test_fetch_reflections_nonexistent(db):
    """Nonexistent cod_id → not in result, no crash."""
    cands = db.fetch_reflections([9999999])
    assert all(c.cod_id != 9999999 for c in cands)


@pytest.mark.integration
def test_fetch_reflections_reflection_keys(db):
    """Reflection dicts must have at least d_hkl and intensity_rel (used by pipeline)."""
    cands = db.fetch_reflections([BEST_COD_ID])
    assert cands, "Need at least 1 candidate"
    for r in cands[0].reflections[:5]:
        assert 'd_hkl' in r, f"Missing d_hkl in reflection: {r.keys()}"
        assert 'intensity_rel' in r, f"Missing intensity_rel in reflection: {r.keys()}"


@pytest.mark.integration
def test_candidates_from_db_smoke():
    """T-030c: candidates_from_db + run_pipeline executes end-to-end without error.

    NOTE on fixture choice: synthetic_candidate17.xye is a complex organometallic
    (cod=1569653, 2155 peaks). The MV stores only top-30 peaks per phase, so
    cod=1569653 is NOT expected in the top-50 from a 530k COD search — other
    phases accidentally overlap more of the 114 detected peaks.
    Hanawalt is designed for minerals (30-50 peaks, distinctive d-values).
    The smoke test therefore verifies correctness of the pipeline path only:
    candidates returned, pipeline runs, result is non-degenerate (crit. 9.7).
    """
    from pipeline import candidates_from_db, run_pipeline

    xye_path = DATA_DIR / 'synthetic_candidate17.xye'

    with DBClient() as db_client:
        candidates = candidates_from_db(
            xye_path, db_client, top_n=50, wavelength=1.54056,
        )

    assert candidates, "candidates_from_db returned empty list — MV empty or 0 matches"
    assert all(isinstance(c, CandidateInput) for c in candidates)
    assert len(candidates) <= 50
    assert all(len(c.reflections) > 0 for c in candidates), \
        "fetch_reflections returned candidate with empty reflections"

    with DBClient() as db_client:
        result = run_pipeline(xye_path, candidates, db_client=db_client)

    assert len(result.candidates) > 0
    assert 0.0 < result.best().Rwp < 2.0, \
        f"Best Rwp={result.best().Rwp:.4f} — degenerate result"
    assert result.best().n_peaks_used > 0


@pytest.mark.integration
def test_candidates_from_db_csv_equivalence():
    """crit. 9.7: fetch_reflections + run_pipeline agrees with CSV path on the winner.

    Uses fetch_reflections(csv_cod_ids) so both paths cover the same 20 phases.
    Verifies that the DB fetch path introduces no reordering bug vs the CSV path.
    (Does NOT compare to the MV Hanawalt ranking — that starts from different
    candidates.)
    """
    from pipeline import run_pipeline

    xye_path  = DATA_DIR / 'synthetic_candidate17.xye'
    csv_cands = load_candidates_csv(CSV_PATH)
    csv_cod_ids = [c.cod_id for c in csv_cands]

    with DBClient() as db_client:
        # same 20 cod_ids, loaded via DB path (fetch_reflections) instead of CSV
        db_cands = db_client.fetch_reflections(csv_cod_ids)
        if not db_cands:
            pytest.skip("No CuKa patterns found in reference_patterns for CSV cod_ids")

        result_db  = run_pipeline(xye_path, db_cands,  db_client=db_client)
        result_csv = run_pipeline(xye_path, csv_cands, db_client=db_client)

    assert result_db.best().cod_id == result_csv.best().cod_id, (
        f"DB fetch path best={result_db.best().cod_id} (Rwp={result_db.best().Rwp:.4f}), "
        f"CSV path best={result_csv.best().cod_id} (Rwp={result_csv.best().Rwp:.4f}). "
        "Divergence indicates bug in fetch_reflections reflection content."
    )
