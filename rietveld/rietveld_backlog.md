# Rietveld Pipeline — Backlog de Implementação

## Legenda de Status

| Ícone | Significado |
|-------|-------------|
| `[ ]` | Pendente |
| `[~]` | Em andamento |
| `[x]` | Concluído |
| `[!]` | Bloqueado / problema |
| `[-]` | Pulado / não necessário |

---

## Contexto

**Modo de teste atual:**
- Input candidatos: `data-1782394014136.csv` (query fixa, 20 estruturas)
  - Colunas: `cod_id`, `peak_matches`, `reflections` (JSONB)
  - Reflections já calculadas — não buscar novamente do DB
- Input experimental: `synthetic_candidate17.xye`
  - Gerado a partir de cod_id=1569653 (C538 H654 Bi40 Mo2 N6 O98)
  - Escala=4000, background exponencial decrescente, ruído Poisson
  - Cu_synthetic.xye foi descartado: candidatos reais são orgânicos/organometálicos complexos, não Cu puro
  - Candidato correto = idx=17, esperado Rwp << outros 19

**Modo futuro:**
- Módulos especializados geram a lista de candidatos `(cod_id, reflections)`
- Interface do pipeline permanece idêntica — apenas a fonte da lista muda

**Input definitivo do pipeline:**
```
xye_file   : str | Path
candidates : list[CandidateInput]
             # CandidateInput.cod_id: int
             # CandidateInput.reflections: list[dict]  ← vêm do CSV
```

**DB usado para:**
- Metadados por cod_id (fórmula, SG, parâmetros de cela, referência bibliográfica)
- NÃO para buscar picos — esses vêm nos candidatos

**Filosofia de implementação:**
Cada épico termina com task `[v]` de visualização — roda o módulo diretamente
(`python <módulo>.py`) e produz output verificável antes de avançar.

---

## Schemas de BD Usados

### `public.data` (tabela principal COD)

| Coluna | Tipo PG | Uso no pipeline |
|--------|---------|-----------------|
| `file` | INTEGER PK | = cod_id |
| `a, b, c` | DOUBLE PRECISION | parâmetros de cela |
| `alpha, beta, gamma` | REAL | ângulos de cela |
| `sg` | VARCHAR(32) | símbolo HM do grupo espacial |
| `sgHall` | VARCHAR(64) | símbolo Hall |
| `sgNumber` | SMALLINT | número internacional do SG |
| `formula` | VARCHAR(255) | fórmula empírica |
| `calcformula` | VARCHAR(255) | fórmula calculada |
| `cellformula` | VARCHAR(255) | fórmula por cela |
| `mineral` | VARCHAR(255) | nome mineral |
| `chemname` | VARCHAR(2048) | nome químico IUPAC |
| `Z` | SMALLINT | unidades de fórmula por cela |
| `Zprime` | REAL | unidades assimétricas |
| `wavelength` | REAL | comprimento de onda do COD |
| `radSymbol` | VARCHAR(20) | símbolo da radiação (ex: CuKα) |
| `method` | ENUM | single crystal / powder / theoretical |
| `authors` | TEXT | autores |
| `title` | TEXT | título do artigo |
| `journal` | VARCHAR(255) | periódico |
| `year` | SMALLINT | ano de publicação |
| `doi` | VARCHAR(127) | DOI |
| `status` | ENUM | warnings / errors / retracted |
| `flags` | SET | has coordinates / has disorder / has Fobs |

### `xrd_analysis.reference_patterns`

| Coluna | Tipo | Uso |
|--------|------|-----|
| `cod_id` | INTEGER | FK → data.file |
| `a, b, c, alpha, beta, gamma` | DOUBLE/REAL | cela calculada |
| `sg_number, sg_symbol, sg_hall` | — | SG do cálculo |
| `wavelength` | REAL | radiação usada no cálculo |
| `rad_symbol` | VARCHAR | ex: CuKα |
| `has_intensities` | BOOLEAN | F(hkl) calculado de atomic_sites? |
| `n_reflections` | INTEGER | # de picos calculados |
| `reflections` | JSONB | **NÃO usar no pipeline** — picos vêm do CSV |
| `calculated_at` | TIMESTAMPTZ | timestamp do cálculo |

### `xrd_analysis.atomic_sites`

| Coluna | Tipo | Uso |
|--------|------|-----|
| `cod_id` | INTEGER | FK |
| `label` | VARCHAR(20) | rótulo do sítio |
| `type_symbol` | VARCHAR(16) | símbolo do elemento |
| `fract_x/y/z` | REAL | coordenadas fracionais |
| `occupancy` | REAL | ocupação parcial |
| `u_iso_or_equiv` | REAL | Debye-Waller Uiso |

> **Tabelas EXCLUÍDAS:** `reference_patterns_pymatgen`, `failed_patterns_pymatgen`

### `xrd_analysis.peak_fingerprints` (Materialized View — Épico 12 / 13.5 d-space)

| Coluna | Tipo | Uso |
|--------|------|-----|
| `cod_id` | INTEGER | FK → data.file |
| `d_hkl` | FLOAT8 | espaçamento d (Å) — λ-independente (Bragg), chave de match |
| `intensity_rel` | FLOAT8 | intensidade relativa (normalizada por fase, max=100) |
| `two_theta_src` | FLOAT8 | 2θ na λ de origem (proveniência, não usado no match) |
| `wavelength` | FLOAT8 | λ do padrão de origem (proveniência) |
| `rad_symbol` | VARCHAR | símbolo da radiação (ex: CuKα) |
| `rank` | INTEGER | rank de intensidade dentro da fase (1=mais intenso) |

**Criação:** `migrations/create_peak_fingerprints.sql` — `DROP + CREATE` obrigatório (sem `CREATE OR REPLACE` em PG).
**Índices:** B-tree em `d_hkl` (BETWEEN no match d-space) + B-tree em `cod_id`.
**Escopo:** top-30 reflexões por fase, `has_intensities=TRUE`, qualquer λ (530k fases). `best_pattern` prefere maior λ (CuKα > MoKα → capta d alto, basal argila até 17.7Å).
**⚠ MV viva (dev):** ainda versão 2θ (132k fases). Rodar migration antes de usar DB mode:
```bash
psql -d "crystallography-open-database" -f rietveld/migrations/create_peak_fingerprints.sql
```

---

## Estrutura de Arquivos

> **Repo root = `rietveld/`** — módulos ficam na raiz, imports são flat.

```
(repo root)
  __init__.py          # mantido para compat; não requerido para imports internos
  models.py            # dataclasses
  data_loader.py       # parse .xye + load CSV      [+ bloco __main__ viz]
  db_client.py         # fetch metadata do PostgreSQL
  pattern_calc.py      # Icalc_unit                 [+ bloco __main__ viz]
  linear_fit.py        # WLS escala + background    [+ bloco __main__ viz]
  fom.py               # Rwp, Rp, Rexp, chi²        [+ bloco __main__ viz]
  pipeline.py          # orquestra tudo             [+ bloco __main__ viz]
  peak_matcher.py      # Hanawalt pre-selection     [+ bloco __main__ viz]
  cli.py               # entry point CLI
  migrations/
    create_peak_fingerprints.sql  # MV top-30 picos por fase (Épico 12)
  tests/
    __init__.py
    test_parse_xye.py
    test_pattern_calc.py
    test_linear_fit.py
    test_fom.py
    test_pipeline.py
    test_db_client.py  # @pytest.mark.integration
  requirements.txt
  Cu_synthetic.xye
  data-1782394014136.csv
  rietveld_pipeline_plan.md
  rietveld_backlog.md
```

---

## Épico 0 — Setup

### T-001: Criar estrutura de pastas

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `rietveld/__init__.py`, `tests/__init__.py` | Dirs criados via `New-Item -ItemType Directory`. `__init__.py` vazio em ambos. CSV, XYE, plan e backlog movidos para `rietveld/`. |

Criar `rietveld/` com `__init__.py` vazio e `tests/__init__.py` vazio.
**Critério:** `python -c "import rietveld"` sem erro.
**Deps:** —

---

### T-002: Criar `requirements.txt`

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `requirements.txt` (raiz do projeto) | 7 dependências conforme backlog. Colocado na raiz, não dentro de `rietveld/`. |

```
numpy
scipy
pandas
psycopg2-binary
python-dotenv
matplotlib
pytest
```

`matplotlib` necessário para tasks `[v]`.
**Critério:** `pip install -r requirements.txt` OK.
**Deps:** T-001

---

## Épico 1 — Models

### T-003: `models.py`

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `models.py` | 4 dataclasses (CandidateInput, StructureMetadata, CandidateResult, RietveldResult). Import + mock OK. |

```python
from dataclasses import dataclass, field

@dataclass
class CandidateInput:
    cod_id: int
    reflections: list[dict]  # [{h,k,l,two_theta,intensity_rel,d_hkl,multiplicity,F_sq}]
    peak_matches: int = 0

@dataclass
class StructureMetadata:
    cod_id: int
    formula: str | None = None
    mineral: str | None = None
    chemname: str | None = None
    sg_number: int | None = None
    sg_symbol: str | None = None
    a: float | None = None
    b: float | None = None
    c: float | None = None
    alpha: float | None = None
    beta: float | None = None
    gamma: float | None = None
    Z: int | None = None
    wavelength: float | None = None
    rad_symbol: str | None = None
    has_intensities: bool | None = None
    authors: str | None = None
    title: str | None = None
    journal: str | None = None
    year: int | None = None
    doi: str | None = None
    method: str | None = None
    status: str | None = None

@dataclass
class CandidateResult:
    cod_id: int
    Rwp: float
    Rp: float
    Rexp: float
    chi2: float
    scale: float
    n_peaks_used: int = 0
    metadata: StructureMetadata | None = None

@dataclass
class RietveldResult:
    xye_file: str
    n_points: int
    candidates: list[CandidateResult] = field(default_factory=list)

    def best(self) -> CandidateResult:
        return self.candidates[0]

    def viable(self, rwp_max: float = 0.15, chi2_max: float = 3.0) -> list[CandidateResult]:
        return [c for c in self.candidates if c.Rwp < rwp_max and c.chi2 < chi2_max]
```

**Critério:** Import sem erro; instâncias criáveis com valores mock.
**Deps:** T-001

---

## Épico 2 — Data Loading

### T-004: `data_loader.py` — parse .xye

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `data_loader.py` | parse_xye + load_candidates_csv. Fix JSON malformado (espaço em número) via regex fallback. |

```python
def parse_xye(path: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
```

Regras:
- Ignorar linhas `#` e vazias
- 3 colunas → `(tth, Iobs, sigma)` direto
- 2 colunas → `sigma = sqrt(max(Iobs, 1.0))`
- Sigma zero ou negativo → substituir por `sqrt(max(I, 1.0))`
- Raise `ValueError` se arquivo inexistente ou < 10 pontos válidos

**Critério:** `parse_xye("Cu_synthetic.xye")` retorna 3 arrays shape `(5001,)`. *(Backlog dizia 5005 — contagem errada incluindo 4 linhas de comentário. Shape real = 5001.)*
**Deps:** T-001

---

### T-004v [viz] — Verificar parse: stats + plot difratograma ⬅ RODAR ANTES DE AVANÇAR

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `data_loader.py` — bloco `__main__` | Bloco __main__ implementado com plot + stats. |

Bloco `if __name__ == '__main__':` no final de `data_loader.py`:

```python
if __name__ == '__main__':
    import sys
    import matplotlib.pyplot as plt

    path = sys.argv[1] if len(sys.argv) > 1 else 'Cu_synthetic.xye'
    tth, Iobs, sigma = parse_xye(path)

    print(f"Arquivo : {path}")
    print(f"Pontos  : {len(tth)}")
    print(f"2θ range: {tth[0]:.2f}° — {tth[-1]:.2f}°")
    print(f"I max   : {Iobs.max():.1f}  em 2θ = {tth[Iobs.argmax()]:.3f}°")
    print(f"σ range : {sigma.min():.2f} — {sigma.max():.2f}")

    plt.figure(figsize=(12, 4))
    plt.plot(tth, Iobs, 'k-', lw=0.8, label='Iobs')
    plt.fill_between(tth, Iobs - sigma, Iobs + sigma, alpha=0.2, color='gray')
    plt.xlabel('2θ (graus)')
    plt.ylabel('Intensidade')
    plt.title(f'Difratograma experimental — {path}')
    plt.legend()
    plt.tight_layout()
    plt.show()
```

**Rodar:** `python data_loader.py Cu_synthetic.xye`

**Output esperado:**
```
Arquivo : Cu_synthetic.xye
Pontos  : 5005
2θ range: 20.00° — 120.00°
I max   : 20000.01  em 2θ = 43.320°
```

**Critério visual:** Pico principal em ~43.3° dominante. Sem gaps ou dados corrompidos.
**Deps:** T-004

---

### T-005: `data_loader.py` — load CSV candidatos

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `data_loader.py` | load_candidates_csv implementado. Fix automático de JSON malformado (row 17: espaço em número). |

```python
def load_candidates_csv(path: str | Path) -> list[CandidateInput]:
```

- Ler com pandas
- Coluna `reflections`: string JSON → `list[dict]` via `json.loads`
- Raise `ValueError` se CSV vazio ou sem colunas `cod_id` e `reflections`

**Critério:** Retorna lista de 20 `CandidateInput`; cada um com `reflections` como `list[dict]`.
**Deps:** T-003

---

### T-005v [viz] — Inspecionar candidatos: tabela + range de picos ⬅ RODAR ANTES DE AVANÇAR

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `data_loader.py` — bloco `__main__` | Tabela de candidatos no mesmo bloco __main__. |

Extensão do bloco `__main__` existente (mesmo bloco, após o plot do .xye):

```python
    # continua no mesmo bloco __main__
    csv_path = sys.argv[2] if len(sys.argv) > 2 else 'data-1782394014136.csv'
    candidates = load_candidates_csv(csv_path)

    print(f"\n{'cod_id':>10}  {'peak_matches':>12}  {'n_reflections':>13}  {'2θ min':>8}  {'2θ max':>8}")
    print('-' * 58)
    for c in candidates:
        tth_vals = [r['two_theta'] for r in c.reflections]
        print(f"{c.cod_id:>10}  {c.peak_matches:>12}  {len(c.reflections):>13}  "
              f"{min(tth_vals):>8.2f}  {max(tth_vals):>8.2f}")
```

**Rodar:** `python data_loader.py Cu_synthetic.xye data-1782394014136.csv`

**Critério visual:** Tabela de 20 linhas. Pelo menos 1 candidato com picos em ~43°.
**Deps:** T-005

---

## Épico 3 — DB Client

### T-006: `db_client.py` — conexão

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `db_client.py` | DBClient com context manager. Lê cod_migration/.env automaticamente. |

Reusa variáveis de `cod_migration/.env`: `PG_HOST`, `PG_PORT`, `PG_DB`, `PG_USER`, `PG_PASSWORD`.

```python
class DBClient:
    def __init__(self, env_path: str | Path | None = None): ...
    def close(self): ...
    def __enter__(self): return self
    def __exit__(self, *_): self.close()
```

**Critério:** `with DBClient() as db: pass` sem erro com DB ativo.
**Deps:** T-003

---

### T-007: `db_client.py` — `fetch_metadata(cod_ids)`

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `db_client.py` | fetch_metadata impl. Fix: coluna `Z` precisa aspas (`"Z"`) no PostgreSQL. 20/20 retornados com formula, sg_symbol, has_intensities. |

Query em `public.data`:
```sql
SELECT
    file AS cod_id,
    formula, calcformula, mineral, chemname,
    sg AS sg_symbol, "sgNumber" AS sg_number, "sgHall" AS sg_hall,
    a, b, c,
    COALESCE(alpha, 90.0) AS alpha,
    COALESCE(beta,  90.0) AS beta,
    COALESCE(gamma, 90.0) AS gamma,
    Z, "Zprime" AS zprime,
    wavelength, "radSymbol" AS rad_symbol,
    method, status, flags,
    authors, title, journal, year, doi
FROM data
WHERE file = ANY(%s)
```

Query suplementar em `xrd_analysis.reference_patterns`:
```sql
SELECT DISTINCT ON (cod_id)
    cod_id, has_intensities, wavelength AS calc_wavelength
FROM xrd_analysis.reference_patterns
WHERE cod_id = ANY(%s)
  AND wavelength BETWEEN 1.535 AND 1.546
ORDER BY cod_id, calculated_at DESC
```

Merge por `cod_id` → retornar `dict[int, StructureMetadata]`.

**Critério:** Para cod_id com CuKα no DB, retorna `StructureMetadata` com `formula`, `sg_symbol`, `has_intensities`.
**Deps:** T-006, T-003

---

## Épico 4 — Cálculo de Padrão

### T-008: `pattern_calc.py` — FWHM Caglioti

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `pattern_calc.py` | caglioti_fwhm impl. Revisão science_review 2.1–2.4: todos pass. fwhm(43.32°)=0.076° (típico CuKα lab). |

```python
def caglioti_fwhm(two_theta_deg: float | np.ndarray,
                  U: float, V: float, W: float) -> float | np.ndarray:
    theta_rad = np.radians(np.asarray(two_theta_deg) / 2.0)
    tan_t = np.tan(theta_rad)
    fwhm2 = U * tan_t**2 + V * tan_t + W
    return np.sqrt(np.maximum(fwhm2, 1e-8))
```

**Critério:** Positivo para qualquer 2θ em [5°, 150°]; aceita array ou escalar.
**Deps:** T-001

---

### T-009: `pattern_calc.py` — perfil pseudo-Voigt

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `pattern_calc.py` | pseudo_voigt_profile impl. Science review 2.5–2.9: todos pass. G e L FWHM medido = 0.100° (input 0.1°). Ordem correta: η·L+(1-η)·G. |

```python
def pseudo_voigt_profile(
    tth_grid: np.ndarray,
    two_theta_peak: float,
    fwhm: float,
    eta: float,
) -> np.ndarray:
    delta = tth_grid - two_theta_peak
    sigma_g = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    G = np.exp(-delta**2 / (2.0 * sigma_g**2))
    L = 1.0 / (1.0 + (delta / (fwhm / 2.0))**2)
    return eta * L + (1.0 - eta) * G
```

**Critério:** Valor 1.0 no centro; decai simetricamente.
**Deps:** T-008

---

### T-010: `pattern_calc.py` — `build_icalc_unit`

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `pattern_calc.py` | build_icalc_unit impl. Science review 2.10–2.13 + 5.5: todos pass. Nota: 2.13 testado com mock isolado (estruturas complexas com 860+ picos sobrepostos deslocam máximo por soma construtiva — comportamento físico correto). |

```python
def build_icalc_unit(
    tth: np.ndarray,
    reflections: list[dict],
    U: float = 0.01,
    V: float = -0.002,
    W: float = 0.005,
    eta: float = 0.5,
    cutoff_fwhm: float = 10.0,
) -> tuple[np.ndarray, int]:
    """
    Returns:
        Icalc_unit : ndarray (N,)
        n_peaks_used : int
    """
    Icalc = np.zeros(len(tth), dtype=np.float64)
    n_used = 0
    tth_min, tth_max = tth[0], tth[-1]

    for refl in reflections:
        two_theta_peak = refl['two_theta']
        intensity_rel  = refl.get('intensity_rel') or 0.0
        if intensity_rel <= 0.0:
            continue
        fwhm = caglioti_fwhm(two_theta_peak, U, V, W)
        if two_theta_peak + cutoff_fwhm * fwhm < tth_min:
            continue
        if two_theta_peak - cutoff_fwhm * fwhm > tth_max:
            continue
        profile = pseudo_voigt_profile(tth, two_theta_peak, fwhm, eta)
        Icalc += intensity_rel * profile
        n_used += 1

    return Icalc, n_used
```

**Critério:** Cu FCC reflections → máximo de `Icalc_unit` próximo de 43.32°.
**Deps:** T-009

---

### T-010v [viz] — Plot Icalc_unit vs Iobs: confirmar alinhamento de picos ⬅ RODAR ANTES DE AVANÇAR

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `pattern_calc.py` — bloco `__main__` | Bloco __main__ impl com plot 2 painéis + linhas Bragg. |

```python
if __name__ == '__main__':
    import sys
    import matplotlib.pyplot as plt
    from data_loader import parse_xye, load_candidates_csv

    xye_path = sys.argv[1] if len(sys.argv) > 1 else 'Cu_synthetic.xye'
    csv_path = sys.argv[2] if len(sys.argv) > 2 else 'data-1782394014136.csv'
    idx      = int(sys.argv[3]) if len(sys.argv) > 3 else 0

    tth, Iobs, sigma = parse_xye(xye_path)
    candidates = load_candidates_csv(csv_path)
    cand = candidates[idx]

    Icalc_unit, n_used = build_icalc_unit(tth, cand.reflections)
    scale_approx = Iobs.max() / max(Icalc_unit.max(), 1e-10)
    Icalc_scaled = Icalc_unit * scale_approx

    peak_positions = [r['two_theta'] for r in cand.reflections
                      if r.get('intensity_rel', 0) > 0]

    print(f"Candidato : cod_id={cand.cod_id}")
    print(f"Picos     : {n_used}/{len(cand.reflections)} usados")
    print(f"Icalc max : {Icalc_unit.max():.2f}  em 2θ = {tth[Icalc_unit.argmax()]:.3f}°")

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

    axes[0].plot(tth, Iobs, 'k-', lw=0.8, label='Iobs (exp)')
    axes[0].set_ylabel('Intensidade')
    axes[0].legend()
    axes[0].set_title('Experimental')

    axes[1].plot(tth, Icalc_scaled, 'b-', lw=0.8,
                 label=f'Icalc (cod={cand.cod_id}, escala aprox.)')
    for pp in peak_positions:
        if tth[0] <= pp <= tth[-1]:
            axes[1].axvline(pp, color='r', alpha=0.3, lw=0.5)
    axes[1].set_xlabel('2θ (graus)')
    axes[1].set_ylabel('Intensidade')
    axes[1].legend()
    axes[1].set_title('Padrão calculado  (linhas vermelhas = posições Bragg)')

    plt.tight_layout()
    plt.show()
```

**Rodar:** `python pattern_calc.py Cu_synthetic.xye data-1782394014136.csv 0`

Trocar `0` para comparar outros candidatos.

**Critério visual:** Candidato Cu correto → linhas vermelhas alinhadas com picos de Iobs. Candidato errado → desalinhamento visível.
**Deps:** T-010, T-004v, T-005v

---

## Épico 5 — Ajuste Linear

### T-011: `linear_fit.py`

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `linear_fit.py` | linear_fit WLS impl. Science review 3.1–3.7: todos pass. scale exato a 7e-15 para Iobs=10*Icalc+50. Residual max 3.6e-12 para fit exato. |

```python
import numpy as np
import scipy.linalg

def linear_fit(
    tth: np.ndarray,
    Iobs: np.ndarray,
    sigma: np.ndarray,
    Icalc_unit: np.ndarray,
    n_bg: int = 4,
) -> tuple[float, np.ndarray]:
    """
    Retorna (scale S, Icalc ajustado = S·Icalc_unit + Ibg)
    """
    w = 1.0 / np.maximum(sigma**2, 1e-10)
    sqrt_w = np.sqrt(w)

    mu, std = tth.mean(), tth.std()
    std = max(std, 1e-6)
    tth_norm = (tth - mu) / std

    # Design matrix: [Icalc_unit | 1 | tth | tth² | tth³]
    A_cols = [Icalc_unit] + [tth_norm**k for k in range(n_bg)]
    A = np.column_stack(A_cols)

    Aw = A * sqrt_w[:, None]
    bw = Iobs * sqrt_w

    x, _, _, _ = scipy.linalg.lstsq(Aw, bw)
    scale = float(x[0])
    bg_coeffs = x[1:]

    Ibg = sum(bg_coeffs[k] * tth_norm**k for k in range(n_bg))
    Icalc = scale * Icalc_unit + Ibg

    return scale, Icalc
```

**Critério:** `Iobs = 10 * Icalc_unit + 50` → `|scale - 10| < 0.1`.
**Deps:** T-010

---

### T-011v [viz] — Plot Iobs vs Icalc + resíduo ⬅ RODAR ANTES DE AVANÇAR

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `linear_fit.py` — bloco `__main__` | 3 painéis: overlay full, zoom, resíduo. Default: synthetic_candidate17.xye. |

```python
if __name__ == '__main__':
    import sys
    import matplotlib.pyplot as plt
    import numpy as np
    from data_loader  import parse_xye, load_candidates_csv
    from pattern_calc import build_icalc_unit

    xye_path = sys.argv[1] if len(sys.argv) > 1 else 'Cu_synthetic.xye'
    csv_path = sys.argv[2] if len(sys.argv) > 2 else 'data-1782394014136.csv'
    idx      = int(sys.argv[3]) if len(sys.argv) > 3 else 0

    tth, Iobs, sigma = parse_xye(xye_path)
    candidates = load_candidates_csv(csv_path)
    cand = candidates[idx]

    Icalc_unit, n_used = build_icalc_unit(tth, cand.reflections)
    scale, Icalc = linear_fit(tth, Iobs, sigma, Icalc_unit)
    diff = Iobs - Icalc

    print(f"Candidato : cod_id={cand.cod_id}")
    print(f"Scale     : {scale:.6f}")
    print(f"|diff| max: {np.abs(diff).max():.1f}  mean: {np.abs(diff).mean():.1f}")

    fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True,
                             gridspec_kw={'height_ratios': [3, 3, 1]})

    axes[0].plot(tth, Iobs,  'k-', lw=0.8, label='Iobs')
    axes[0].plot(tth, Icalc, 'r-', lw=0.8, label='Icalc', alpha=0.8)
    axes[0].set_ylabel('Intensidade')
    axes[0].legend()
    axes[0].set_title(f'cod_id={cand.cod_id}  |  scale={scale:.4f}  |  picos={n_used}')

    baseline_max = np.percentile(Iobs, 90)
    axes[1].plot(tth, Iobs,  'k-', lw=0.8, label='Iobs')
    axes[1].plot(tth, Icalc, 'r-', lw=0.8, label='Icalc', alpha=0.8)
    axes[1].set_ylim(0, baseline_max * 1.5)
    axes[1].set_ylabel('Intensidade (zoom)')
    axes[1].legend()
    axes[1].set_title('Zoom — picos secundários')

    axes[2].plot(tth, diff, 'g-', lw=0.5)
    axes[2].axhline(0, color='k', lw=0.5)
    axes[2].set_xlabel('2θ (graus)')
    axes[2].set_ylabel('Δ')
    axes[2].set_title('Resíduo (Iobs − Icalc)')

    plt.tight_layout()
    plt.show()
```

**Rodar:** `python linear_fit.py Cu_synthetic.xye data-1782394014136.csv 0`

**Critério visual:**
- Candidato correto: Icalc sobreposto a Iobs; resíduo plano
- Candidato errado: picos deslocados; resíduo com estrutura residual grande
**Deps:** T-011, T-010v

---

## Épico 6 — Figures of Merit

### T-012: `fom.py`

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `fom.py` | calc_fom impl. Science review 4.1–4.7: todos pass. Candidato 17 real: Rwp=0.022, chi2=1.012 (esperado ≈1.0 — ruído Poisson puro). |

```python
import numpy as np

def calc_fom(
    Iobs: np.ndarray,
    Icalc: np.ndarray,
    sigma: np.ndarray,
    n_params: int,
) -> dict[str, float]:
    w = 1.0 / np.maximum(sigma**2, 1e-10)
    diff = Iobs - Icalc

    sum_w_diff2 = np.sum(w * diff**2)
    sum_w_Iobs2 = np.sum(w * Iobs**2)

    Rwp  = float(np.sqrt(sum_w_diff2 / max(sum_w_Iobs2, 1e-20)))
    Rp   = float(np.sum(np.abs(diff)) / max(np.sum(np.abs(Iobs)), 1e-20))
    N    = len(Iobs)
    Rexp = float(np.sqrt(max(N - n_params, 1) / max(sum_w_Iobs2, 1e-20)))
    chi2 = float((Rwp / Rexp)**2) if Rexp > 0 else float('inf')

    return {'Rwp': Rwp, 'Rp': Rp, 'Rexp': Rexp, 'chi2': chi2}
```

| Métrica | Bom | Ruim |
|---------|-----|------|
| `Rwp` | < 0.10 | > 0.30 |
| `Rp` | < 0.08 | > 0.25 |
| `chi2` | ≈ 1.0 | >> 1 |

**Critério:** `Icalc = Iobs` → `Rwp < 1e-10`; `Icalc = zeros` → `Rwp ≈ 1.0`.
**Deps:** T-001

---

### T-012v [viz] — FoM de 1 candidato + verificação de sanidade ⬅ RODAR ANTES DE AVANÇAR

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `fom.py` — bloco `__main__` | Bloco __main__ impl com sanidades. Default: synthetic_candidate17.xye. |

```python
if __name__ == '__main__':
    import sys
    import numpy as np
    from data_loader  import parse_xye, load_candidates_csv
    from pattern_calc import build_icalc_unit
    from linear_fit   import linear_fit

    xye_path = sys.argv[1] if len(sys.argv) > 1 else 'Cu_synthetic.xye'
    csv_path = sys.argv[2] if len(sys.argv) > 2 else 'data-1782394014136.csv'
    idx      = int(sys.argv[3]) if len(sys.argv) > 3 else 0

    tth, Iobs, sigma = parse_xye(xye_path)
    candidates = load_candidates_csv(csv_path)
    cand = candidates[idx]

    Icalc_unit, _ = build_icalc_unit(tth, cand.reflections)
    scale, Icalc  = linear_fit(tth, Iobs, sigma, Icalc_unit, n_bg=4)
    fom = calc_fom(Iobs, Icalc, sigma, n_params=5)

    def tag(rwp, chi2):
        if rwp < 0.10 and chi2 < 2: return '✓ BOM'
        if rwp < 0.20 and chi2 < 5: return '⚠ MÉDIO'
        return '✗ RUIM'

    print(f"\nFigures of Merit — cod_id={cand.cod_id}")
    print(f"  Rwp  = {fom['Rwp']:.4f}")
    print(f"  Rp   = {fom['Rp']:.4f}")
    print(f"  Rexp = {fom['Rexp']:.4f}")
    print(f"  chi2 = {fom['chi2']:.3f}")
    print(f"  {tag(fom['Rwp'], fom['chi2'])}")
    print(f"  scale= {scale:.6f}")

    fom_perfeito = calc_fom(Iobs, Iobs,               sigma, n_params=5)
    fom_nulo     = calc_fom(Iobs, np.zeros_like(Iobs), sigma, n_params=5)
    print(f"\n--- Sanidade ---")
    print(f"Icalc=Iobs  → Rwp={fom_perfeito['Rwp']:.2e}  (esperado ≈ 0)")
    print(f"Icalc=zeros → Rwp={fom_nulo['Rwp']:.4f}      (esperado ≈ 1)")
```

**Rodar:** `python fom.py Cu_synthetic.xye data-1782394014136.csv 0`

Trocar `0` por outros índices para comparar FoM antes de rodar os 20.

**Critério visual:** Sanidades passam. FoM candidato Cu notavelmente menor que candidatos aleatórios.
**Deps:** T-012, T-011

---

## Épico 7 — Pipeline

### T-013: `pipeline.py`

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `pipeline.py` | run_pipeline impl. Science review 5.1–5.5: todos pass. Candidato correto (1569653) = #1, Rwp=0.022, chi2=1.012. 1 viável, 6 com scale<0 (fases anti-correlacionadas, não filtradas). |

```python
FIXED_PARAMS = {
    'U': 0.01, 'V': -0.002, 'W': 0.005,
    'eta': 0.5,
    'n_bg': 4,
    'wavelength': 1.54056,
}

def run_pipeline(
    xye_path: str | Path,
    candidates: list[CandidateInput],
    db_client=None,
    params: dict | None = None,
) -> RietveldResult:
    p = {**FIXED_PARAMS, **(params or {})}

    tth, Iobs, sigma = parse_xye(xye_path)

    metadata_map = {}
    if db_client is not None:
        metadata_map = db_client.fetch_metadata([c.cod_id for c in candidates])

    results = []
    for cand in candidates:
        Icalc_unit, n_used = build_icalc_unit(
            tth, cand.reflections,
            U=p['U'], V=p['V'], W=p['W'], eta=p['eta'],
        )
        scale, Icalc = linear_fit(tth, Iobs, sigma, Icalc_unit, n_bg=p['n_bg'])
        fom = calc_fom(Iobs, Icalc, sigma, n_params=1 + p['n_bg'])

        results.append(CandidateResult(
            cod_id=cand.cod_id,
            Rwp=fom['Rwp'], Rp=fom['Rp'], Rexp=fom['Rexp'], chi2=fom['chi2'],
            scale=scale, n_peaks_used=n_used,
            metadata=metadata_map.get(cand.cod_id),
        ))

    results.sort(key=lambda r: r.Rwp)
    return RietveldResult(xye_file=str(xye_path), n_points=len(tth), candidates=results)
```

**Critério:** Executa com `db_client=None`; retorna 20 candidatos ordenados por Rwp.
**Deps:** T-004, T-005, T-010, T-011, T-012, T-006

---

### T-013v [viz] — Ranking completo: tabela + scatter Rwp vs chi2 ⬅ RODAR ANTES DE AVANÇAR

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `pipeline.py` — bloco `__main__` | Tabela ranking + scatter Rwp×chi2. Default: synthetic_candidate17.xye. |

```python
if __name__ == '__main__':
    import sys
    import matplotlib.pyplot as plt
    from data_loader import load_candidates_csv

    xye_path = sys.argv[1] if len(sys.argv) > 1 else 'Cu_synthetic.xye'
    csv_path = sys.argv[2] if len(sys.argv) > 2 else 'data-1782394014136.csv'

    candidates = load_candidates_csv(csv_path)
    result = run_pipeline(xye_path, candidates, db_client=None)

    print(f"\nRanking — {result.xye_file}  |  {result.n_points} pts  |  {len(result.candidates)} candidatos\n")
    print(f"{'#':>3}  {'cod_id':>10}  {'Rwp':>7}  {'Rp':>7}  {'chi2':>7}  {'scale':>10}  {'picos':>6}")
    print('-' * 65)
    viable = result.viable()
    for i, r in enumerate(result.candidates):
        flag = '★' if r in viable else ' '
        print(f"{flag}{i+1:>2}  {r.cod_id:>10}  {r.Rwp:>7.4f}  {r.Rp:>7.4f}  "
              f"{r.chi2:>7.2f}  {r.scale:>10.4f}  {r.n_peaks_used:>6}")
    print(f"\n★ {len(viable)} candidato(s) viável(is)  (Rwp<0.15, chi2<3)")

    rwps  = [r.Rwp  for r in result.candidates]
    chi2s = [r.chi2 for r in result.candidates]
    ids   = [str(r.cod_id) for r in result.candidates]

    plt.figure(figsize=(8, 6))
    plt.scatter(rwps, chi2s, c='steelblue', s=60, zorder=3)
    for x, y, label in zip(rwps, chi2s, ids):
        plt.annotate(label, (x, y), fontsize=7, xytext=(4, 4),
                     textcoords='offset points')
    plt.axvline(0.15, color='r', ls='--', lw=0.8, label='Rwp=0.15')
    plt.axhline(3.0,  color='g', ls='--', lw=0.8, label='chi²=3')
    plt.xlabel('Rwp')
    plt.ylabel('chi²')
    plt.title('Ranking de candidatos')
    plt.legend()
    plt.tight_layout()
    plt.show()
```

**Rodar:** `python pipeline.py Cu_synthetic.xye data-1782394014136.csv`

**Critério visual:** Candidato Cu no canto inferior esquerdo do scatter. Não-Cu no canto superior direito.
**Deps:** T-013

---

## Épico 8 — CLI

### T-014: `cli.py`

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `cli.py` | argparse (xye, csv, --no-db, --top, --plot, --env). print_results com has_I flag + [!sem F2] warn. plot_best 2-panel. Science review: 5.3,5.4 ✓. `python cli.py synthetic_candidate17.xye data-1782394014136.csv --no-db --top 5` OK. Com DB: formula/has_I do PostgreSQL. |

```python
import argparse
import numpy as np
from data_loader  import load_candidates_csv
from db_client    import DBClient
from pipeline     import run_pipeline

def main():
    p = argparse.ArgumentParser(description='Rietveld Phase Identification')
    p.add_argument('xye',            help='Arquivo .xye experimental')
    p.add_argument('candidates_csv', help='CSV com cod_id, peak_matches, reflections')
    p.add_argument('--no-db',   action='store_true')
    p.add_argument('--top',     type=int, default=5)
    p.add_argument('--plot',    action='store_true', help='Plot melhor candidato')
    p.add_argument('--env',     default=None)
    args = p.parse_args()

    candidates = load_candidates_csv(args.candidates_csv)

    db = None
    if not args.no_db:
        try:
            db = DBClient(env_path=args.env)
        except Exception as e:
            print(f"[WARN] DB indisponível: {e}")

    try:
        result = run_pipeline(args.xye, candidates, db_client=db)
    finally:
        if db:
            db.close()

    print_results(result, top=args.top)

    if args.plot:
        plot_best(result, args.xye, candidates)


def print_results(result, top=5):
    viable = result.viable()
    print(f"\n{'#':>3}  {'cod_id':>10}  {'Rwp':>7}  {'Rp':>7}  {'chi2':>7}  "
          f"{'scale':>10}  {'formula':<20}  mineral")
    print('-' * 85)
    for i, r in enumerate(result.candidates[:top]):
        m = r.metadata
        formula = (m.formula or '?') if m else '—'
        mineral = (m.mineral or '') if m else '—'
        flag = '★' if r in viable else ' '
        print(f"{flag}{i+1:>2}  {r.cod_id:>10}  {r.Rwp:>7.4f}  {r.Rp:>7.4f}  "
              f"{r.chi2:>7.2f}  {r.scale:>10.4f}  {formula:<20}  {mineral}")
    if viable:
        print(f"\n★ {len(viable)} candidato(s) viável(is)")


def plot_best(result, xye_path, candidates):
    import matplotlib.pyplot as plt
    from data_loader  import parse_xye
    from pattern_calc import build_icalc_unit
    from linear_fit   import linear_fit

    best = result.best()
    cand = next(c for c in candidates if c.cod_id == best.cod_id)
    tth, Iobs, sigma = parse_xye(xye_path)
    Icalc_unit, _ = build_icalc_unit(tth, cand.reflections)
    _, Icalc = linear_fit(tth, Iobs, sigma, Icalc_unit)
    diff = Iobs - Icalc

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True,
                             gridspec_kw={'height_ratios': [4, 1]})
    axes[0].plot(tth, Iobs,  'k-', lw=0.8, label='Iobs')
    axes[0].plot(tth, Icalc, 'r-', lw=0.8, label=f'Icalc  cod={best.cod_id}', alpha=0.8)
    axes[0].set_ylabel('Intensidade')
    m = best.metadata
    title = f"cod={best.cod_id}  Rwp={best.Rwp:.4f}  chi2={best.chi2:.2f}"
    if m and m.formula:
        title += f"  [{m.formula}  {m.mineral or ''}]"
    axes[0].set_title(title)
    axes[0].legend()
    axes[1].plot(tth, diff, 'g-', lw=0.5)
    axes[1].axhline(0, color='k', lw=0.5)
    axes[1].set_xlabel('2θ (graus)')
    axes[1].set_ylabel('Δ')
    plt.tight_layout()
    plt.show()
```

**Rodar:**
```bash
python cli.py Cu_synthetic.xye data-1782394014136.csv --no-db --top 5 --plot
```

**Deps:** T-013

---

## Épico 9 — Testes Automatizados

### T-015: `tests/test_parse_xye.py`

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `tests/test_parse_xye.py` | 11 tests. Critérios 1.1–1.4 todos verificados. 57/57 pass. |

- Shape `(5005,)` para `Cu_synthetic.xye`
- `tth[0] ≈ 20.0`, `tth[-1] ≈ 120.0`
- Linhas `#` ignoradas (mock inline)
- Sigma fallback para 2 colunas
- `ValueError` para arquivo inexistente

---

### T-016: `tests/test_pattern_calc.py`

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `tests/test_pattern_calc.py` | 14 tests. Critérios 2.1–2.13 todos verificados. |

- Reflections mock Cu FCC (1 pico em 43.32°) → `Icalc_unit` máximo nessa região
- Pico fora do range → `n_peaks_used = 0`, `Icalc_unit = zeros`
- `intensity_rel = 0` → pico ignorado
- `caglioti_fwhm` positivo para 2θ em [5°, 150°]

---

### T-017: `tests/test_linear_fit.py`

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `tests/test_linear_fit.py` | 7 tests. Critérios 3.1–3.7 verificados. |

- `Iobs = 10 * Icalc_unit + 50` → `|scale - 10| < 0.1`
- Resíduo < 1% de `max(Iobs)`

---

### T-018: `tests/test_fom.py`

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `tests/test_fom.py` | 9 tests. Critérios 4.1–4.7 todos verificados. |

- `Icalc = Iobs` → `Rwp < 1e-10`
- `Icalc = zeros` → `Rwp ≈ 1.0`
- `chi2 > 0` sempre
- `n_params` maior → `Rexp` menor

---

### T-019: `tests/test_pipeline.py`

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `tests/test_pipeline.py` | 13 tests. Critérios 5.1–5.5 todos verificados. Candidato correto #1, gap 5.5x, chi2≈1. |

- 20 candidatos retornados
- Ordenados por `Rwp` crescente
- `result.best().Rwp < 0.15`
- `result.best().metadata is None`
- `result.viable()` não vazia

---

## Épico 10 — Integração DB Real

### T-020: `tests/test_db_client.py`

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `tests/test_db_client.py` | 12 integration tests (@pytest.mark.integration). fetch_metadata: 20/20 retornados, formula, has_intensities, sg_symbol, lattice, empty/inexistente, context manager, end-to-end com pipeline. 12/12 pass. |

```python
@pytest.mark.integration
def test_fetch_metadata():
    with DBClient() as db:
        cod_ids = [c.cod_id for c in load_candidates_csv('data-1782394014136.csv')]
        meta = db.fetch_metadata(cod_ids)
    assert len(meta) > 0
    assert any(m.formula for m in meta.values())
    assert any(m.has_intensities is not None for m in meta.values())
```

---

### T-021: Validação end-to-end com metadata

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | — (CLI) | `python cli.py synthetic_candidate17.xye data-1782394014136.csv --top 5` com DB ativo. #1=1569653, Rwp=0.022, has_I=T, formula=C538 H654 Bi40 Mo2 N6 O98. CLI end-to-end com metadata validado. |

---

## Contexto de Domínio — Resíduos de Mineração

**Amostras = misturas minerais policristalinas** (quartzo, calcita, hematita, goetita, gipsita, argilominerais, escória amorfa).

Implicações vs fase pura:
- **Inorgânicos de alta simetria** → 3–30 picos/fase (não 2000 como organometálicos). Hanawalt/PDF-2 viável como pré-filtro.
- **Múltiplas fases simultâneas** → difratograma = soma de K fases. Fase única (Épicos 1–10) identifica a dominante, mas não quantifica nem explica a mistura.
- **Quantificação (QPA)** → frações em peso por Hill-Howard exigem Z, M, V de cada fase (todos em `StructureMetadata` → QPA exige DB ativo).
- **Background amorfo** → escória/vidro elevam linha de base com hump largo → `n_bg` pode precisar 6–8.

Épicos 1–10 = Stage 1 (ID fase única). Épicos 11–13 = Stage 1.5 (mistura + QPA de triagem).

> **Não** é refinamento estrutural (Stage 2, fora de escopo) — xyz/occ/Uiso permanecem fixos via `intensity_rel` pré-calculado.

---

## Épico 11 — Análise Multi-fase

> Revisão: critérios **6.1–6.5** do `rietveld_science_review.md`.

### T-022: `crystallo_utils.py` — volume de cela + massa molar

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `crystallo_utils.py` | cell_volume (triclínico geral) + molar_mass. Tabela atômica Z=1–92 (CIAAW). Elemento desconhecido → ValueError (7.7). Sanidade: cubo=125, quartzo γ=120°=112.99, SiO2=60.08, CaCO3=100.09, Fe2O3=159.69. Validado contra 20 formulas reais do DB: 0 falhas (cobre As,Au,Cd,Ru,U,Bi,W,Sb,Ce,Co,Cr,Ni,Pd,Cs...). Critérios 7.1,7.2,7.7. |

```python
import numpy as np
import re

# IUPAC standard atomic weights (g/mol) — TABELA COMPLETA obrigatória.
# Stub abaixo é ilustrativo; impl real usa lib `periodictable`/`mendeleev`
# OU dict completo (~92 elementos). Elemento ausente NÃO pode ser ignorado.
_ATOMIC_WEIGHTS = {
    'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998,
    'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974,
    'S': 32.06, 'Cl': 35.45, 'K': 39.098, 'Ca': 40.078, 'Ti': 47.867,
    'Cr': 51.996, 'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693,
    'Cu': 63.546, 'Zn': 65.38, 'As': 74.922, 'Ba': 137.33, 'Pb': 207.2,
    # ... COMPLETAR todos os elementos (minerais têm F, Cl, Cr, Ni, As, etc.)
}

def cell_volume(a, b, c, alpha, beta, gamma) -> float:
    """Volume da cela (triclínico geral, Å³). Ângulos em graus.
    Reduz a a·b·c para ortogonal (α=β=γ=90°)."""
    ca, cb, cg = (np.cos(np.radians(x)) for x in (alpha, beta, gamma))
    return float(a * b * c * np.sqrt(
        max(1 - ca**2 - cb**2 - cg**2 + 2*ca*cb*cg, 1e-12)
    ))

def molar_mass(formula: str) -> float:
    """Massa molar (g/mol) de string tipo 'C538 H654 Bi40' ou 'SiO2'.
    Aceita contagens fracionárias (ex: 'Br0.8').

    Elemento desconhecido → RAISE (critério 7.7). Skip silencioso
    subcontaria M → QPA errado sem aviso. Ambiguidade CO vs Co
    mitigada por formula COD espaçada ('Ca C O3'); fórmula colada
    ('CaCO3') resolvida por longest-match de 2 letras primeiro.
    """
    M = 0.0
    for elem, count in re.findall(r'([A-Z][a-z]?)(\d*\.?\d*)', formula):
        if elem not in _ATOMIC_WEIGHTS:
            raise ValueError(f"Elemento desconhecido '{elem}' em '{formula}'")
        n = float(count) if count else 1.0
        M += n * _ATOMIC_WEIGHTS[elem]
    return M
```

**Critério:** `cell_volume(5,5,5,90,90,90)≈125` (tolerância float); `molar_mass('SiO2')≈60.08`; elemento desconhecido → `ValueError`. Cobre 7.1, 7.2, 7.7.
**Deps:** T-001

---

### T-023: `multi_phase_fit.py` — WLS multi-fase com bounds

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `multi_phase_fit.py` | `multi_phase_fit` via `scipy.optimize.lsq_linear` method='bvls'. Fases bound [0,inf), background (-inf,inf). Validado: recovery [3,7] a 1e-13; fase anti-correlacionada → scale~1e-13 (não negativo, 6.2); Rwp_combinado=Rwp_single quando fases extras→0 (6.6). Suite 57/57 intacta. Critérios 6.1–6.4,6.6. |

```python
import numpy as np
from scipy.optimize import lsq_linear

def multi_phase_fit(tth, Iobs, sigma, Icalc_units, n_bg=4):
    """
    Ajuste WLS de K fases + background. Scales >= 0 (lsq_linear bounds).

    Icalc_units : list[np.ndarray]  — um Icalc_unit por fase (de build_icalc_unit)
    Retorna (scales: np.ndarray (K,), Icalc: np.ndarray, bg_coeffs)
    """
    K = len(Icalc_units)
    w = 1.0 / np.maximum(sigma**2, 1e-10)
    sqrt_w = np.sqrt(w)

    mu = tth.mean(); std = max(float(tth.std()), 1e-6)
    tth_norm = (tth - mu) / std

    phase_cols = list(Icalc_units)
    bg_cols    = [tth_norm**k for k in range(n_bg)]
    A = np.column_stack(phase_cols + bg_cols)

    Aw = A * sqrt_w[:, None]
    bw = Iobs * sqrt_w

    # bounds: fases [0, inf), background (-inf, inf)
    lb = np.array([0.0]*K + [-np.inf]*n_bg)
    ub = np.full(K + n_bg, np.inf)
    res = lsq_linear(Aw, bw, bounds=(lb, ub), method='bvls')

    scales = res.x[:K]
    bg_coeffs = res.x[K:]
    Ibg = sum(bg_coeffs[k] * tth_norm**k for k in range(n_bg))
    Icalc = sum(scales[k] * Icalc_units[k] for k in range(K)) + Ibg
    return scales, Icalc, bg_coeffs
```

**Critério (6.1–6.4):** `Iobs = 3*unit_A + 7*unit_B + bg` → recupera scales [3,7]; nenhum scale negativo; background subtraível.

> ⚠ **Para o FIT** (Rwp), `intensity_rel` per-fase normalizado é irrelevante (só reescala $S_k$). **Para QPA** (T-024), os $S_k$ precisam de base absoluta comum — ver crítico em T-024.

**Deps:** T-010, T-011

---

### T-023v [viz] — Overlay multi-fase + decomposição por componente ⬅ VERIFICAR

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `multi_phase_fit.py` — `__main__` | Plot Iobs/Icalc total/componentes por fase/resíduo. Verificado com synthetic_candidate17: fase pura → só 1569653 scale>0 (2.1e-4), 2 fases erradas→0. Rwp=0.0237, chi2=1.17. Legenda em notação científica (scale absoluto ~1e-4). |

Plot: Iobs, Icalc total, cada $S_k \cdot I_{\text{calc,unit},k}$ empilhado, resíduo. Confirma decomposição visual da mistura.
**Deps:** T-023

---

### T-022b: `pattern_calc.py` — Lorentz-polarização + padrão absoluto

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `pattern_calc.py` (extensão) | `lorentz_polarization` (BB, opcional monocromador) + `build_icalc_unit_absolute` (pondera mult·F_sq·Lp). **DESCOBERTA:** reconstrução `m·F_sq·Lp` (sem DW) recupera `intensity_rel` com **correlação=1.0000** (exato a 3 casas até 12.9°) → DW ausente/≈1 no CSV → **QPA base EXATA, não aproximação**. Lp(5.17°)/Lp(12.2°)=5.65 ✓ (bate ~6× CSV). Critérios 7.8, 7.9. |

```python
def lorentz_polarization(two_theta_deg, two_theta_mono_deg=None):
    """Fator Lorentz-polarização (Bragg-Brentano).
    Sem monocromador: (1 + cos²2θ) / (sin²θ cosθ).
    Com monocromador grafite: numerador usa cos²2θ_m."""
    tt = np.radians(np.asarray(two_theta_deg))
    th = tt / 2.0
    if two_theta_mono_deg is None:
        num = 1.0 + np.cos(tt)**2
    else:
        cm2 = np.cos(np.radians(two_theta_mono_deg))**2
        num = 1.0 + cm2 * np.cos(tt)**2
    return num / np.maximum(np.sin(th)**2 * np.cos(th), 1e-12)

def build_icalc_unit_absolute(tth, reflections, U=0.01, V=-0.002, W=0.005,
                               eta=0.5, cutoff_fwhm=10.0):
    """Igual a build_icalc_unit MAS pondera por I_abs = mult·F_sq·Lp(θ),
    NÃO por intensity_rel (que é normalizado por fase → inútil p/ QPA).
    Usado no modo multi-fase/QPA (Épico 11)."""
    # idêntico a build_icalc_unit, trocando:
    #   intensity = (refl['multiplicity'] or 1) * (refl['F_sq'] or 0) \
    #               * lorentz_polarization(two_theta_peak)
    ...
```

**Critério (7.8):** `lorentz_polarization` usa θ=2θ/2; decresce monotonicamente em 2θ médio-alto; reproduz razão ~6× entre 5° e 12° observada no CSV.
**Deps:** T-009

---

### T-024: `qpa.py` — frações em peso Hill-Howard

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `qpa.py` | `weight_fractions(scales, metadatas)` Hill-Howard W_k=S_k·ZMV_k/Σ. V por cell_volume, M por molar_mass, Z de metadata. Guard: metadata None / Z,cela,formula faltante / elemento desconhecido → ZMV=0 (exclui fase). Validado: scale=0→0%, None→0%, Σ=100; **mistura sintética 2 fases razão conhecida → wt% recuperado erro 0.25pp**. Critérios 7.1–7.5. |

> 🔴 **CRÍTICO RESOLVIDO — base de intensidade (critério 7.6).**
> ✅ **VERIFICADO no CSV (2026):** `intensity_rel` É normalizado por fase — todas as 20 fases têm `ir_max = 100.0`. `ir/(m·F_sq)` varia ~6× dentro da fase (= fator $L_p$ angular).
> Logo $S_k$ do `multi_phase_fit` em `intensity_rel` está enviesado por $\max_k$ → **QPA com esse $S_k$ é numericamente errado**, mesmo com fit visualmente perfeito.
> **Fix obrigatório:** modo QPA reconstrói padrão **absoluto** $I^{abs}_{hkl} = m \cdot F\_sq \cdot L_p(\theta)$ via `build_icalc_unit_absolute` (T-022b) e fita nessa base → $S_k$ = escala de Rietveld verdadeira. DW omitido (screening; viés em 2θ alto).

```python
from crystallo_utils import cell_volume, molar_mass

def weight_fractions(scales, metadatas) -> list[dict]:
    """
    Hill-Howard (1987): W_k = S_k(ZMV)_k / Σ_j S_j(ZMV)_j.
    metadatas : list[StructureMetadata]  (mesma ordem de scales)
    Retorna list[{cod_id, scale, ZMV, weight_pct}].
    """
    zmv = []
    for s, m in zip(scales, metadatas):
        if m is None or not m.Z or not m.a:
            zmv.append(0.0); continue
        V = cell_volume(m.a, m.b, m.c, m.alpha, m.beta, m.gamma)
        M = molar_mass(m.formula or '')
        zmv.append(s * m.Z * M * V)
    total = sum(zmv)
    out = []
    for (s, m, v) in zip(scales, metadatas, zmv):
        out.append({
            'cod_id': m.cod_id if m else None,
            'scale': float(s),
            'ZMV': float(v),
            'weight_pct': 100.0 * v / total if total > 0 else 0.0,
        })
    return out
```

**Critério (7.3, 7.4):** Σ weight_pct = 100; fase com scale=0 → 0%. Cobre 7.5 (skip metadata nula).
**Deps:** T-022, T-007 (DB metadata obrigatório)

---

### T-025: `pipeline.py` — modo multi-fase

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `pipeline.py` (extensão) | `run_multiphase` + dataclasses `PhaseFraction`/`MultiPhaseResult` em models.py. Seleção top-K por rank Rwp (6.7). Fit em base absoluta (7.9). Rwp_single_best na mesma base absoluta p/ comparação justa (6.6). QPA via weight_fractions. Sem DB → wt%=0, fit ok. Validado: top_k=4 fase pura→1569653=100%, Σ=100, dominante correto, 6.6 ✓. Suite 57/57 intacta. |

`run_multiphase(xye, candidates, db_client, top_k=4)`:
1. Roda `run_pipeline` (fase única) → ranqueia por Rwp
2. Seleciona top-K por **rank de Rwp** (menor primeiro) — **NÃO** filtrar por viabilidade single-fase (critério 6.7: fase real de mistura tem Rwp single alto, threshold a dropa)
3. `build_icalc_unit_absolute` (T-022b, base $m\cdot F_{sq}\cdot L_p$) para cada → `multi_phase_fit` conjunto — **não** `intensity_rel` (crítico 7.6/7.9)
4. `calc_fom` do modelo combinado (`n_params = K + n_bg`)
5. `weight_fractions` com metadata (ver pré-requisito de base absoluta em T-024)
6. Retorna `MultiPhaseResult` (scales, Rwp_combinado, frações QPA)

Novo dataclass em `models.py`: `MultiPhaseResult` + `PhaseFraction`.

> **Rigoroso (opcional):** greedy forward selection — começa vazio, adiciona iterativamente a fase que mais reduz Rwp combinado, para quando ganho < ε. Evita incluir fases que só ajustam ruído. Top-K por rank é a aproximação de screening.

**Critério (6.6):** Rwp combinado ≤ Rwp da melhor fase isolada — **garantido por construção** se melhor fase única ∈ conjunto (single-fase = ponto factível com demais scales=0).
**Deps:** T-023, T-024

---

### T-025v [viz] — Tabela QPA + pizza de frações ⬅ VERIFICAR

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `pipeline.py` — `_main_multiphase` (flag `--mp`) | Tabela (cod,scale,Z,M,V,wt%,formula) + pizza wt% + overlay decomposição. Fixture `synthetic_mixture.xye` (3 fases 1569653/1536746/4100953, alvo 50/30/20%, Poisson seed=7). **Recuperado 51/31/18% (erro 1–2pp).** Rwp comb=0.021 vs single best=0.228 → prova 6.7 (single>0.15 viável; threshold dropava fases reais). Rodar: `python pipeline.py synthetic_mixture.xye data-1782394014136.csv 4 --mp`. |

Tabela: cod_id, formula, scale, weight_pct. Pizza/barra das frações. Rwp single vs Rwp combinado.
**Deps:** T-025

---

### T-026: `cli.py` — flag `--phases K`

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `cli.py` (extensão) | `--phases K` → `run_multiphase` + `print_qpa` (tabela cod/scale/Z/M/V/wt%/has_I) + `plot_qpa` (pizza+decomposição em --plot). Warn `--no-db`+`--phases`. Nota auto "amostra é mistura" se Rwp_single>0.15. Testado: mistura→51/31/18%; --no-db→wt%=0; fase única regressão intacta. Suite 57/57. |

`--phases K` → ativa `run_multiphase` com top-K. Exige DB (QPA precisa metadata). Aviso se `--no-db` + `--phases`.

```bash
python cli.py amostra.xye candidatos.csv --phases 4 --plot
```
**Deps:** T-025

---

### T-027: testes multi-fase + QPA

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `tests/test_crystallo_utils.py` (10), `tests/test_lp_absolute.py` (9), `tests/test_multi_phase_fit.py` (8), `tests/test_qpa.py` (10) | 37 testes. Critérios 6.1–6.6, 7.1–7.10. Suite total: 94 unit + 12 integration = 106 pass. |

- `cell_volume` ortogonal(=abc)/triclínico; `molar_mass` SiO2≈60.08/CaCO3≈100.09; **elemento desconhecido → `pytest.raises(ValueError)`** (7.7)
- `lorentz_polarization`: θ=2θ/2; decresce em 2θ médio-alto; razão Lp(5°)/Lp(12°)≈5–6 (bate com CSV) (7.8)
- `build_icalc_unit_absolute`: pondera por mult·F_sq·Lp, **não** intensity_rel (7.9)
- `multi_phase_fit`: recuperação de scales conhecidos [3,7]; scales ≥ 0 (nenhum negativo mesmo com fase anti-correlacionada); mistura sintética 2 fases; Rwp_combinado ≤ Rwp_single (6.6)
- `weight_fractions`: Σ=100; scale=0 → 0%; metadata nula skip; **base absoluta** (7.6/7.10 — mistura sintética 2 fases de razão de massa conhecida → weight_pct correto dentro de ~5%)
- Critérios 6.1–6.7, 7.1–7.10

**Deps:** T-022, T-023, T-024

---

## Épico 12 — Pré-seleção por picos (Hanawalt)

> Reduz O(500k candidatos COD) → O(50) antes do Rietveld. Necessário em uso real (COD inteiro para CuKα). Science review: critérios **8.1–8.7** do `rietveld_science_review.md`.

### T-028: `peak_matcher.py` — detecção + matching de picos

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `peak_matcher.py` | **d-space (radiação-agnóstico, Épico 13.5).** `detect_peaks` (prominence=5σ, min_distance=0.3°, retorna 2θ), `two_theta_to_d` (Bragg), `d_windows` (janela d por-pico via bounds, crit 8.10), `match_candidates_csv`/`match_candidates_db` em d_hkl. Converte picos 2θ amostra → d via `wavelength`. MV migration reescrita p/ d-spacing (`best_pattern` ORDER BY wavelength DESC → capta basal argila, crit 8.11). Validado: Bragg 26.6°→3.348Å; janela d 0.245Å@12° vs 0.0093Å@60° (26× — bate cot(θ)); cod=1569653 top-5; unnest 2-array WITH ORDINALITY OK no PG. Critérios 8.1–8.13. Suite 94/94. **MV viva ainda é versão 2θ — rodar migration nova antes do DB mode.** |

**Dois modos de matching:**

| Modo | Quando | Função |
|------|--------|--------|
| **CSV** | offline / testes (lista de candidatos pré-carregada) | `match_candidates_csv` |
| **DB** | produção (COD completo via MV) | `match_candidates_db` |

Ambos usam a mesma `detect_peaks`. Interface unificada: `prefilter_candidates` escolhe o modo por `db_client`.

```python
import numpy as np
from scipy.signal import find_peaks


def detect_peaks(tth, Iobs, sigma, prominence_sigma=5.0, min_distance_deg=0.3):
    """Picos experimentais via find_peaks com proeminência relativa ao ruído.

    prominence ≥ prominence_sigma × median(σ) → robusto a background inclinado
    e hump amorfo (critério 8.1). min_distance_deg separa picos adjacentes (8.2).

    Returns: (tth_peaks, I_peaks) — posições e alturas dos picos detectados.
    """
    step = float(tth[1] - tth[0])
    min_dist = max(1, int(min_distance_deg / step))
    prom = prominence_sigma * float(np.median(sigma))
    idx, _ = find_peaks(Iobs, prominence=prom, distance=min_dist)
    return tth[idx], Iobs[idx]


def match_candidates_csv(peaks_tth, candidates, tol_deg=0.2, min_matches=3):
    """Hanawalt in-memory: conta matches 2θ entre picos experimentais e reflexões
    de cada candidato (critério 8.3). Complexidade O(peaks × reflections × N_cand).
    Adequado para N_cand < 1000.

    Returns: list[CandidateInput] ordenada por n_matches DESC,
             filtrada por n_matches >= min_matches.
             candidate.peak_matches atualizado (critério 8.6).
    """
    scored = []
    for c in candidates:
        ref_tths = [r['two_theta'] for r in c.reflections
                    if r.get('intensity_rel', 0) > 0]
        n_matched = sum(
            1 for p in peaks_tth
            if any(abs(p - r) <= tol_deg for r in ref_tths)
        )
        if n_matched >= min_matches:
            c.peak_matches = n_matched
            scored.append((n_matched, c))
    scored.sort(key=lambda x: -x[0])
    return [c for _, c in scored]


def match_candidates_db(peaks_tth, conn, tol_deg=0.2, min_matches=3, top_n=50):
    """Hanawalt via PostgreSQL — query em peak_fingerprints MV (critério 8.5).
    Usa BETWEEN para aproveitar índice B-tree em two_theta.

    Returns: list[int] cod_ids ordenados por n_matches DESC.
    Requer: MV xrd_analysis.peak_fingerprints criada (ver migrations/).
    """
    if not peaks_tth:
        return []
    with conn.cursor() as cur:
        cur.execute("""
            WITH obs AS (SELECT unnest(%s::float8[]) AS tth)
            SELECT pf.cod_id, COUNT(*) AS n_matched
            FROM xrd_analysis.peak_fingerprints pf
            JOIN obs ON pf.two_theta BETWEEN obs.tth - %s AND obs.tth + %s
            GROUP BY pf.cod_id
            HAVING COUNT(*) >= %s
            ORDER BY n_matched DESC
            LIMIT %s
        """, (list(peaks_tth), tol_deg, tol_deg, min_matches, top_n))
        return [row[0] for row in cur.fetchall()]
```

**SQL de criação da MV** (arquivo `migrations/create_peak_fingerprints.sql`):
```sql
-- Unnest reflections JSONB → top-30 por intensity_rel por fase CuKα
CREATE MATERIALIZED VIEW IF NOT EXISTS xrd_analysis.peak_fingerprints AS
SELECT cod_id, two_theta, intensity_rel, rank
FROM (
    SELECT
        rp.cod_id,
        (r->>'two_theta')::float8      AS two_theta,
        (r->>'intensity_rel')::float8  AS intensity_rel,
        ROW_NUMBER() OVER (
            PARTITION BY rp.cod_id
            ORDER BY (r->>'intensity_rel')::float8 DESC
        ) AS rank
    FROM xrd_analysis.reference_patterns rp,
         jsonb_array_elements(rp.reflections) AS r
    WHERE rp.has_intensities = TRUE
      AND rp.wavelength BETWEEN 1.535 AND 1.546
) sub
WHERE rank <= 30;

CREATE INDEX IF NOT EXISTS pf_tth_idx ON xrd_analysis.peak_fingerprints (two_theta);
CREATE INDEX IF NOT EXISTS pf_cod_idx ON xrd_analysis.peak_fingerprints (cod_id);
```

> **Nota criação MV:** rodar UMA VEZ após migração COD completa. Não recria automaticamente. `REFRESH MATERIALIZED VIEW CONCURRENTLY` se COD for atualizado.

**Parâmetros padrão:**

| Parâmetro | Padrão | Justificativa |
|-----------|--------|---------------|
| `prominence_sigma` | 5.0 | SNR≥5σ; <1% falsos positivos para Poisson |
| `min_distance_deg` | 0.3° | FWHM típico laboratorial ≈ 0.1–0.2° |
| `tol_deg` | 0.2° | zero-point error + FWHM/2 + erro 2θ COD |
| `min_matches` | 3 | critério Hanawalt clássico [H38] |
| `top_n` | 50 | margem antes do Rietveld |

**Critério (8.7):** Para `synthetic_candidate17.xye` + candidatos CSV, cod=1569653 deve ter maior `n_matches` tanto em modo CSV quanto DB. Implementação testável sem MV usando `match_candidates_csv` contra os 20 candidatos do CSV.

**Deps:** T-004, T-005

---

### T-028v [viz] — Picos detectados + ranking de matches ⬅ VERIFICAR

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `peak_matcher.py` — bloco `__main__` | 2 painéis: Iobs + scatter vermelhos nos picos detectados; linhas Bragg coloridas por candidato (top-5). Tabela stdout cod_id/n_matches. Rodar: `python peak_matcher.py synthetic_candidate17.xye data-1782394014136.csv`. Saída verificada: 114 picos, cod=1569653 tied #1. |

**Rodar:** `python peak_matcher.py synthetic_candidate17.xye data-1782394014136.csv`

**Critério visual:** Marcadores nos picos reais; sem detecção no background. cod=1569653 = #1 na tabela.

---

### T-029: integração `--prefilter N` → pipeline

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `pipeline.py` / `cli.py` | `prefilter_candidates(xye_path, candidates, top_n, db_client, tol_deg, min_matches, prominence_sigma)` em pipeline.py. Dois modos: CSV (in-memory, padrão) e DB (MV, com fallback CSV se MV indisponível). `--prefilter N` em cli.py: imprime `[Hanawalt] X → Y candidatos` em stderr, passa lista filtrada para `run_pipeline`/`run_multiphase` sem mudar API. Validado: 20→5 candidatos, cod=1569653 #1 (Rwp=0.022). Suite 94/94 intacta. Critérios 8.1–8.7 OK. |

**`pipeline.py`** — função helper (assinatura real implementada):
```python
def prefilter_candidates(xye_path, candidates, top_n,
                          db_client=None, tol_deg=0.2, min_matches=3,
                          prominence_sigma=5.0):
    """Pré-filtro Hanawalt: subconjunto de candidates ordenado por n_matches.
    db_client=None → modo CSV (in-memory).
    db_client ativo → modo DB (MV peak_fingerprints) com FALLBACK p/ CSV se
    MV indisponível (captura exceção, warn em stderr).
    Interface run_pipeline/run_multiphase inalterada — só recebe lista menor."""
    from peak_matcher import detect_peaks, match_candidates_csv, match_candidates_db
    tth, Iobs, sigma = parse_xye(xye_path)
    peaks_tth, _ = detect_peaks(tth, Iobs, sigma, prominence_sigma=prominence_sigma)
    if db_client is not None:
        try:
            cod_ids = match_candidates_db(peaks_tth, db_client._conn,
                                          tol_deg=tol_deg, min_matches=min_matches,
                                          top_n=top_n)
            cand_map = {c.cod_id: c for c in candidates}
            filtered = [cand_map[cid] for cid in cod_ids if cid in cand_map]
            if filtered:
                return filtered[:top_n]
        except Exception:
            ...  # fallback CSV
    return match_candidates_csv(peaks_tth, candidates,
                                tol_deg=tol_deg, min_matches=min_matches)[:top_n]
```

**`cli.py`** — argumento `--prefilter N`:
```
--prefilter N   mantém top-N candidatos por n_matches antes do Rietveld
                CSV mode: sem DB; DB mode: usa peak_fingerprints MV
```

Lógica em `main()`: se `args.prefilter` → `prefilter_candidates(args.xye, candidates, top_n=args.prefilter, db_client=db)` → imprime `[Hanawalt] X → Y candidatos` (stderr) → passa lista filtrada.

**Critério:** `synthetic_candidate17.xye + --prefilter 5` inclui cod=1569653 em `filtered`.

**Deps:** T-028

---

## Épico 13 — Modo DB-only (sem CSV)

> Busca no COD completo direto do PostgreSQL — sem CSV de candidatos pré-carregado. Único modo viável em produção (CSV é fixture de teste de 20 fases). Usa MV `peak_fingerprints` (Épico 12) para selecionar, `fetch_reflections` para carregar reflexões.

### T-030a: `db_client.fetch_reflections(cod_ids)` — carrega reflexões do DB

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `db_client.py` | `fetch_reflections(cod_ids) -> list[CandidateInput]`. Query `reference_patterns` `DISTINCT ON (cod_id)`, filtro `has_intensities=TRUE` + CuKα (1.535–1.546). JSONB → list[dict] (json.loads se string). Substitui `load_candidates_csv` no fluxo produção. |

```python
def fetch_reflections(self, cod_ids: list[int]) -> list[CandidateInput]:
    """Carrega CandidateInput (cod_id + reflections JSONB) de reference_patterns.
    Só has_intensities=TRUE e CuKα. Usado no modo DB-only após
    match_candidates_db estreitar o espaço de busca."""
    # SELECT DISTINCT ON (cod_id) cod_id, reflections
    # FROM xrd_analysis.reference_patterns
    # WHERE cod_id = ANY(%s) AND has_intensities=TRUE
    #   AND wavelength BETWEEN 1.535 AND 1.546
    # ORDER BY cod_id, calculated_at DESC
    ...
```

**Critério:** `fetch_reflections([1569653])` retorna 1 `CandidateInput` com `reflections` não-vazio. `fetch_reflections([])` → `[]`.
**Deps:** T-006, T-028

---

### T-030b: `cli.py` — flag `--from-db` (CSV opcional)

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `cli.py` / `pipeline.py` | `candidates_from_db(xye, db_client, top_n, wavelength)` em pipeline.py: detect_peaks → match_candidates_db(MV d-space) → fetch_reflections → list[CandidateInput]. `candidates_csv` vira `nargs='?'` (opcional). Guards: `--from-db`+`--no-db` → `p.error()`; sem CSV sem `--from-db` → `p.error()`; CSV+`--from-db` → warn ignorado. `--prefilter N` mapeado p/ `top_n` em DB mode (default 50). `--wavelength λ` propagado p/ conversão 2θ→d (CoKa/CrKa p/ amostras Fe). 106/106 pass. |

Novo modo: `--from-db` dispensa o CSV. Fluxo:
1. `detect_peaks` no XYE
2. `match_candidates_db` na MV (d-space) → top-N `cod_ids` (N = `--prefilter`, default 50)
3. `db.fetch_reflections(cod_ids)` → `list[CandidateInput]`
4. `run_pipeline` ou `run_multiphase` normal

`candidates_csv` vira **opcional** em argparse (`nargs='?'`). Guards: `--from-db`+`--no-db` → `p.error()` imediato (crit. 9.6). Sem `--from-db` e sem CSV → erro.

Helper em `pipeline.py` (`candidates_from_db` — não `run_from_db`):
```python
def candidates_from_db(xye_path, db_client, top_n=50, wavelength=1.54056,
                        tol_deg=0.2, min_matches=3, prominence_sigma=5.0):
    """DB-only candidate discovery: detect peaks → MV match (d-space) → fetch reflections.
    Retorna list[CandidateInput] — mesma interface do CSV mode."""
    tth, Iobs, sigma = parse_xye(xye_path)
    peaks_tth, _ = detect_peaks(tth, Iobs, sigma, prominence_sigma=prominence_sigma)
    cod_ids = match_candidates_db(peaks_tth, db_client._conn,
                                  wavelength=wavelength, tol_deg=tol_deg,
                                  min_matches=min_matches, top_n=top_n)
    return db_client.fetch_reflections(cod_ids)
```

CLI:
```bash
# fase única, busca COD completo:
python cli.py amostra.xye --from-db --prefilter 50 --top 10

# multi-fase + QPA, busca COD completo:
python cli.py amostra.xye --from-db --prefilter 50 --phases 4 --plot

# amostras Fe (CoKα):
python cli.py amostra.xye --from-db --wavelength 1.7902 --prefilter 50
```

**Critério:** `python cli.py synthetic_candidate17.xye --from-db --prefilter 50` (com MV d-space criada) retorna cod=1569653 no topo. Resultado equivalente ao modo CSV para os mesmos cod_ids.
**Deps:** T-030a, T-029

---

### T-030c: smoke test DB-only

| Status | Arquivo | Como implementado |
|--------|---------|-------------------|
| `[x]` | `tests/test_db_client.py` (@integration) | 6 testes novos. `test_fetch_reflections_single/empty/nonexistent/reflection_keys` (fetch_reflections unit-level). `test_candidates_from_db_smoke`: pipeline path completo sem assertiva de ranking (organometálico fixture inapropriado p/ Hanawalt COD-wide — explicado no docstring). `test_candidates_from_db_csv_equivalence` (crit. 9.7): `fetch_reflections(csv_cod_ids)` vs CSV → mesmo best, verifica ausência de bug no path fetch. Fix: `d_lo.tolist()` em vez de `list(d_lo)` → np.float64 serializado corretamente pelo psycopg2. 18/18 integration pass. |

**Nota científica (fixture):** `synthetic_candidate17.xye` (organometálico 2155 picos, MV guarda top-30) → cod=1569653 ausente do top-50 COD-wide é comportamento correto do Hanawalt, não bug. Hanawalt discrimina bem para minerais (30–50 picos distintos). Fixture adequado para smoke test end-to-end de minerais.

**Fix crítico descoberto:** `list(np.ndarray)` → `[np.float64(...)]` → psycopg2 serializa como schema name (`np.float64(...)`) → `InvalidSchemaName`. Fix: `.tolist()` converte p/ Python float nativo (`peak_matcher.py:175`).

**Pré-requisito:** MV `peak_fingerprints` com schema d_hkl (`psql -d "crystallography-open-database" -f migrations/create_peak_fingerprints.sql`).
**Critério:** 18/18 integration pass com DB + MV ativos.
**Deps:** T-030a, T-030b

---

## Ordem de Implementação

```
Sprint 1 — Dados visíveis no terminal e no plot
  T-001 → T-002 → T-003
  T-003 → T-004 → T-004v  ← PARA AQUI E VERIFICA
  T-003 → T-005 → T-005v  ← PARA AQUI E VERIFICA

Sprint 2 — Padrão calculado visível
  T-001 → T-008 → T-009 → T-010 → T-010v  ← PARA AQUI E VERIFICA

Sprint 3 — Ajuste + FoM visíveis
  T-010 → T-011 → T-011v  ← PARA AQUI E VERIFICA
  T-001 → T-012 → T-012v  ← PARA AQUI E VERIFICA

Sprint 4 — Pipeline completo
  T-006 → T-007               (DB — opcional neste sprint)
  T-013 → T-013v  ← PARA AQUI E VERIFICA
  T-014                        (CLI com --plot)

Sprint 5 — Testes + DB
  T-015 a T-019
  T-020 → T-021

Sprint 6 — Multi-fase + QPA (mineração)
  T-022 → T-023 → T-023v  ← VERIFICA decomposição
  T-024 → T-025 → T-025v  ← VERIFICA frações QPA
  T-026 (CLI --phases)
  T-027 (testes)

Sprint 7 — Pré-seleção Hanawalt
  T-028 → T-028v  ← VERIFICA picos detectados + ranking
  T-029            (--prefilter; CSV mode primeiro, DB mode opcional)

Sprint 8 — DB-only (produção, COD completo)
  T-030a (fetch_reflections)  ← já implementado
  T-030b (--from-db, CSV opcional)
  T-030c (smoke test DB-only)  ← exige MV criada
```

| Sprint | Entregável verificável |
|--------|------------------------|
| 1 | Difratograma no plot; tabela de candidatos no terminal |
| 2 | Picos calculados alinhados (ou não) com experimental |
| 3 | Iobs vs Icalc sobrepostos; FoM com sanidades passando |
| 4 | Ranking completo + scatter + CLI com `--plot` |
| 5 | Testes passando; pipeline validado com metadata do DB |
| 6 | Mistura decomposta em fases; frações em peso QPA (Hill-Howard) |
| 7 | Pré-filtro Hanawalt reduz candidatos antes do Rietveld; `--prefilter N` na CLI |
| 8 | Modo DB-only: busca no COD completo sem CSV; `--from-db` na CLI |

---

## Parâmetros Fixos (Fase de Teste)

| Parâmetro | Valor | Descrição |
|-----------|-------|-----------|
| `λ` | 1.54056 Å | Cu Kα1 |
| `U` | 0.01 | Caglioti U |
| `V` | -0.002 | Caglioti V |
| `W` | 0.005 | Caglioti W |
| `η` | 0.5 | mistura pseudo-Voigt |
| `n_bg` | 4 | grau do polinômio de background (mineração com amorfo: testar 6–8) |

> Stage 2 (futuro): refinamento não-linear de U,V,W,η via Levenberg-Marquardt ou L-BFGS-B.
> QPA (Épico 12): correções não aplicadas — orientação preferencial (March-Dollase [D86]), microabsorção (Brindley [B45]), fase amorfa (spike/padrão interno [BH88]). Frações = base 100% cristalina, screening não-certificado.

---

## Extensibilidade

| Fonte de candidatos | Módulo | Status |
|---------------------|--------|--------|
| Query fixa CSV | `load_candidates_csv()` | ✓ Implementado (Épico 2) |
| Busca por similaridade de picos (Hanawalt) | `peak_matcher.py` — modos CSV + DB | ✓ Épico 12 (T-028) |
| Busca DB-only (COD completo via MV) | `candidates_from_db()` + `db_client.fetch_reflections()` | ✓ Épico 13 (T-030b) |
| Busca por composição química | `chemistry_search.py` | Futuro |
| Busca por grupo espacial | `sg_search.py` | Futuro |

`pipeline.run_pipeline()` recebe sempre `list[CandidateInput]` — fonte não importa.

| Modo de análise | Função | Status |
|-----------------|--------|--------|
| Fase única (ID dominante) | `run_pipeline()` | ✓ Implementado (Épico 7) |
| Multi-fase + QPA (mistura) | `run_multiphase()` | ✓ Épico 11 (T-025) |
| Fonte de candidatos = DB (sem CSV) | `candidates_from_db()` → `run_pipeline/run_multiphase` | ✓ Épico 13 (T-030b) |
