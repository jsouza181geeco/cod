-- peak_fingerprints MV — Epico 12 (T-028) / Epico 13.5 (d-spacing)
--
-- Fingerprint por d-spacing (radiacao-agnostico), NAO por two_theta.
-- d_hkl e intrinseco do cristal (independe de lambda); two_theta = derivado
-- por Bragg. Matching em d-space (Hanawalt classico) → 1 MV cobre Cu/Co/Cr/Mo.
--
-- Pega 1 pattern por fase: o de MAIOR lambda (CuKa > MoKa), tie-break
-- n_reflections. Lambda maior captura d MAIOR (basal de argila ate ~17A em
-- CuKa vs 8A em MoKa) — critico p/ minerio com argilominerais (crit. 8.11).
-- d_hkl independe da lambda; so o range 2theta-cutoff muda por lambda.
--
-- Science review: criterios 8.4 (top-K por intensity_rel), 8.5 (B-tree index),
--                 8.8 (d-spacing radiacao-agnostico).
--
-- ATENCAO: MV ja existe na base de dev (versao two_theta). CREATE ... IF NOT
-- EXISTS seria no-op. Postgres NAO tem CREATE OR REPLACE MATERIALIZED VIEW.
-- DROP + CREATE obrigatorio. MV e dado derivado de reference_patterns →
-- DROP nao perde nada original; reconstroi pela query abaixo.

DROP MATERIALIZED VIEW IF EXISTS xrd_analysis.peak_fingerprints;

CREATE MATERIALIZED VIEW xrd_analysis.peak_fingerprints AS
WITH best_pattern AS (
    -- 1 pattern por fase: maior lambda (capta d alto/basal), tie n_reflections
    SELECT DISTINCT ON (cod_id)
        cod_id, reflections, wavelength, rad_symbol, n_reflections
    FROM xrd_analysis.reference_patterns
    WHERE has_intensities = TRUE
    ORDER BY cod_id, wavelength DESC, n_reflections DESC
),
peaks AS (
    SELECT
        bp.cod_id,
        (r->>'d_hkl')::float8          AS d_hkl,          -- fundamental (A)
        (r->>'intensity_rel')::float8  AS intensity_rel,  -- 0..100 por fase
        (r->>'two_theta')::float8      AS two_theta_src,  -- na lambda de origem
        bp.wavelength,                                     -- proveniencia
        bp.rad_symbol,                                     -- proveniencia
        ROW_NUMBER() OVER (
            PARTITION BY bp.cod_id
            ORDER BY (r->>'intensity_rel')::float8 DESC
        ) AS rank
    FROM best_pattern bp,
         jsonb_array_elements(bp.reflections) AS r
    WHERE (r->>'intensity_rel')::float8 > 0
      AND (r->>'d_hkl')::float8 IS NOT NULL
      AND (r->>'d_hkl')::float8 > 0
)
SELECT cod_id, d_hkl, intensity_rel, two_theta_src, wavelength, rad_symbol, rank
FROM peaks
WHERE rank <= 30;

-- B-tree em d_hkl: matching usa BETWEEN d_obs-tol AND d_obs+tol (crit. 8.5)
-- pipeline converte picos 2theta da amostra → d via lambda DA AMOSTRA:
--   d_obs = lambda_sample / (2 * sin(radians(two_theta_obs / 2)))
CREATE INDEX IF NOT EXISTS pf_d_idx   ON xrd_analysis.peak_fingerprints (d_hkl);
CREATE INDEX IF NOT EXISTS pf_cod_idx ON xrd_analysis.peak_fingerprints (cod_id);

-- Refresh apos atualizar reference_patterns:
-- REFRESH MATERIALIZED VIEW CONCURRENTLY xrd_analysis.peak_fingerprints;
-- (CONCURRENTLY exige indice UNICO; criar se precisar refresh sem lock:)
-- CREATE UNIQUE INDEX IF NOT EXISTS pf_uniq ON xrd_analysis.peak_fingerprints (cod_id, d_hkl);
