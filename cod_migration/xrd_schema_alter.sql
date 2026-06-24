-- xrd_schema_alter.sql
-- Executar UMA VEZ antes de iniciar xrd_loader.py
-- Cria nova tabela reference_patterns_pymatgen (dados fisicamente corretos via pymatgen)
-- A tabela original reference_patterns permanece intacta para comparação.
--
-- Uso:
--   psql -h localhost -U postgres -d "crystallography-open-database" -f cod_migration/xrd_schema_alter.sql
--
-- Após validar os novos dados, para substituir a tabela original:
--   DROP TABLE xrd_analysis.reference_patterns;
--   ALTER TABLE xrd_analysis.reference_patterns_pymatgen RENAME TO reference_patterns;

-- ---------------------------------------------------------------------------
-- 1. Nova tabela de padrões XRD calculados via pymatgen
--    Diferenças da original:
--      + formula        (string da fórmula expandida pela simetria)
--      + t_load_s / t_calc_s / rss_delta_mb (instrumentação)
--      - wavelength_source / two_theta_min / two_theta_max / hkl_max (N/A — fixos)
--    UNIQUE (cod_id, wavelength) → permite recalcular com outra radiação no futuro
-- ---------------------------------------------------------------------------

CREATE TABLE IF NOT EXISTS xrd_analysis.reference_patterns_pymatgen (
    id              BIGINT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
    cod_id          INTEGER NOT NULL REFERENCES data(file) ON DELETE CASCADE,

    -- Parâmetros de cela (extraídos do Structure.from_file + expansão de simetria)
    a               DOUBLE PRECISION NOT NULL,
    b               DOUBLE PRECISION NOT NULL,
    c               DOUBLE PRECISION NOT NULL,
    alpha           REAL NOT NULL DEFAULT 90.0,
    beta            REAL NOT NULL DEFAULT 90.0,
    gamma           REAL NOT NULL DEFAULT 90.0,

    -- Grupo espacial (da tabela data do COD)
    sg_number       SMALLINT,
    sg_symbol       VARCHAR(32),
    sg_hall         VARCHAR(64),

    -- Fórmula química da célula completa (expandida pela simetria)
    formula         TEXT,

    -- Radiação
    wavelength      REAL NOT NULL DEFAULT 1.54184,  -- CuKa médio (pymatgen padrão)
    rad_symbol      VARCHAR(20) DEFAULT 'CuKa',

    -- Padrão XRD
    reflections     JSONB NOT NULL,
    n_reflections   INTEGER NOT NULL,

    -- Instrumentação (timing + memória por CIF)
    t_load_s        REAL,
    t_calc_s        REAL,
    rss_delta_mb    REAL,

    calculated_at   TIMESTAMPTZ NOT NULL DEFAULT now(),

    UNIQUE (cod_id, wavelength)
);

-- Índices
CREATE INDEX IF NOT EXISTS rp_pmg_cod_id
    ON xrd_analysis.reference_patterns_pymatgen (cod_id);

CREATE INDEX IF NOT EXISTS rp_pmg_n_reflections
    ON xrd_analysis.reference_patterns_pymatgen (n_reflections);

CREATE INDEX IF NOT EXISTS rp_pmg_n_refl_t_calc
    ON xrd_analysis.reference_patterns_pymatgen (n_reflections, t_calc_s)
    WHERE t_calc_s IS NOT NULL;

-- GIN index para busca em reflections (criar APÓS carga completa — lento em bulk)
-- CREATE INDEX rp_pmg_reflections_gin
--     ON xrd_analysis.reference_patterns_pymatgen USING GIN (reflections)
--     WITH (fastupdate = off);

COMMENT ON TABLE xrd_analysis.reference_patterns_pymatgen IS
    'Padrões XRD calculados via pymatgen XRDCalculator (fisicamente corretos). '
    'Structure.from_file() expande unidade assimétrica → célula completa. '
    'has_intensities=TRUE sempre. two_theta_range=(0,90). '
    'Reflections: [{h,k,l,d_hkl,two_theta,multiplicity,intensity_rel}]';

-- ---------------------------------------------------------------------------
-- 2. Tabela de falhas do loader pymatgen
--    ON CONFLICT DO UPDATE → reprocessamento atualiza erro anterior
-- ---------------------------------------------------------------------------

CREATE TABLE IF NOT EXISTS xrd_analysis.failed_patterns_pymatgen (
    cod_id    INTEGER PRIMARY KEY REFERENCES data(file) ON DELETE CASCADE,
    status    VARCHAR(20) NOT NULL DEFAULT 'FAILED',  -- FAILED | TIMEOUT
    error_msg TEXT,
    failed_at TIMESTAMPTZ DEFAULT now()
);

CREATE INDEX IF NOT EXISTS fp_pmg_status
    ON xrd_analysis.failed_patterns_pymatgen (status);

-- ---------------------------------------------------------------------------
-- 3. Otimizações PostgreSQL para bulk load
--    Executar ANTES de iniciar xrd_loader.py:
-- ---------------------------------------------------------------------------
-- ALTER SYSTEM SET synchronous_commit = off;
-- ALTER SYSTEM SET work_mem = '256MB';
-- ALTER SYSTEM SET maintenance_work_mem = '1GB';
-- ALTER SYSTEM SET checkpoint_completion_target = 0.9;
-- SELECT pg_reload_conf();
--
-- Executar APÓS carga completa:
-- ALTER SYSTEM RESET synchronous_commit;
-- ALTER SYSTEM RESET work_mem;
-- SELECT pg_reload_conf();
-- VACUUM ANALYZE xrd_analysis.reference_patterns_pymatgen;
-- CREATE INDEX rp_pmg_reflections_gin ON xrd_analysis.reference_patterns_pymatgen
--     USING GIN (reflections) WITH (fastupdate = off);

-- ---------------------------------------------------------------------------
-- 4. Queries de diagnóstico pós-run
-- ---------------------------------------------------------------------------
-- -- Distribuição por tier + tempo médio
-- SELECT
--     CASE
--         WHEN n_reflections <= 500  THEN 'SMALL'
--         WHEN n_reflections <= 2000 THEN 'MEDIUM'
--         ELSE 'LARGE'
--     END AS tier,
--     count(*)                                AS n_patterns,
--     round(avg(t_calc_s)::numeric, 2)        AS avg_t_calc_s,
--     round(max(t_calc_s)::numeric, 2)        AS max_t_calc_s,
--     round(avg(rss_delta_mb)::numeric, 1)    AS avg_rss_mb
-- FROM xrd_analysis.reference_patterns_pymatgen
-- GROUP BY 1 ORDER BY 1;
--
-- -- Correlação n_reflections → t_calc_s (valida heurística de balanceamento)
-- SELECT
--     (n_reflections / 100) * 100             AS n_refl_bucket,
--     count(*)                                AS n,
--     round(avg(t_calc_s)::numeric, 2)        AS avg_t_calc_s
-- FROM xrd_analysis.reference_patterns_pymatgen
-- WHERE t_calc_s IS NOT NULL
-- GROUP BY 1 ORDER BY 1;
--
-- -- Falhas por tipo
-- SELECT status, count(*), left(error_msg, 80) AS sample_error
-- FROM xrd_analysis.failed_patterns_pymatgen
-- GROUP BY 1, 3 ORDER BY 2 DESC LIMIT 20;
--
-- -- Substituição após validação:
-- DROP TABLE xrd_analysis.reference_patterns;
-- ALTER TABLE xrd_analysis.reference_patterns_pymatgen RENAME TO reference_patterns;
-- ALTER TABLE xrd_analysis.failed_patterns_pymatgen RENAME TO failed_patterns;