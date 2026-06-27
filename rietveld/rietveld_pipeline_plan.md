# Plano de Implementação: Refinamento de Rietveld

## Visão Geral

O refinamento de Rietveld compara um padrão de difração experimental (`.xye`) contra padrões
calculados a partir de estruturas cristalinas (CIF/DB), minimizando a diferença via mínimos
quadrados ponderados. O resultado são os **Figures of Merit (FoM)** que quantificam a qualidade
do ajuste.

Pipeline:
```
Cu_synthetic.xye
      +
candidates.csv (20 estruturas do PostgreSQL)
      |
      v
[1] Parse dados experimentais
      |
      v
[2] Construir Icalc por candidato
      |    ├─ [2a] Peak shape (pseudo-Voigt)
      |    ├─ [2b] FWHM Caglioti
      |    └─ [2c] Somar contribuições de todos os picos
      |
      v
[3] Ajuste linear: escala + background
      |
      v
[4] Calcular FoM: Rwp, Rp, Rexp, χ²
      |
      v
[5] Ranking de candidatos por Rwp
```

---

## Etapa 1 — Parse dos Dados Experimentais

### O que faz
Lê o arquivo `.xye` e extrai os três arrays necessários para o refinamento.

### Input
- **Arquivo:** `Cu_synthetic.xye`
- **Formato:** texto, 3 colunas por linha: `2θ  I  σ`
- **Origem:** padrão XRD experimental (ou sintético para teste)
- **Linhas com `#`:** comentários, ignorados

### Output
| Variável | Tipo | Descrição |
|----------|------|-----------|
| `tth` | `ndarray (N,)` | ângulos 2θ em graus |
| `Iobs` | `ndarray (N,)` | intensidade observada |
| `sigma` | `ndarray (N,)` | desvio padrão de cada ponto (≈ √Iobs para estatística de contagem Poisson) |

### Fórmula
Nenhuma. Leitura direta de arquivo.

---

## Etapa 2 — Construção de Icalc

### O que faz
Para cada candidato, calcula o padrão teórico `Icalc(2θ)` somando a contribuição de cada pico
de Bragg sobre o grid experimental.

### Input
| Variável | Origem |
|----------|--------|
| `tth` (grid experimental) | Etapa 1 |
| `reflections` (lista de picos) | coluna `reflections` do `candidates.csv` (JSONB do PostgreSQL) |
| `U, V, W` (parâmetros Caglioti) | estimativa inicial fixa: `U=0.01, V=-0.002, W=0.005` |
| `eta` (mistura Gaussiana/Lorentziana) | estimativa inicial fixa: `0.5` |

Cada pico em `reflections` contém:
```json
{
  "h": 1, "k": 1, "l": 1,
  "two_theta": 43.317,
  "intensity_rel": 100.0,
  "multiplicity": 8,
  "F_sq": 12500.0,
  "d_hkl": 2.087
}
```
> **Nota:** `intensity_rel` já incorpora `Lp · DW · M · |F|²` normalizado — calculado pelo pipeline `xrd_schema_setup.py`.

### Output
| Variável | Tipo | Descrição |
|----------|------|-----------|
| `Icalc_unit` | `ndarray (N,)` | padrão calculado com escala=1, sem background |

### Sub-etapas e Fórmulas

#### 2a — FWHM via Caglioti
Largura a meia altura do pico em função do ângulo. Modela o alargamento instrumental.

```
FWHM²(θ) = U·tan²θ + V·tanθ + W
```

- `θ` = metade de `two_theta` (em radianos)
- `U, V, W` = parâmetros de Caglioti (refinaveis; fixos na fase de teste)
- Resultado: `FWHM` em graus

#### 2b — Função de Perfil: Pseudo-Voigt
Aproximação da convolução entre Gaussiana (alargamento instrumental) e Lorentziana
(alargamento por tamanho de cristalito/microdeformação).

```
pV(Δ2θ) = η · L(Δ2θ, FWHM) + (1-η) · G(Δ2θ, FWHM)
```

Onde:
```
G(Δ2θ) = exp(-Δ2θ² / (2σ²))         σ = FWHM / (2√(2·ln2))

L(Δ2θ) = 1 / (1 + (Δ2θ / (FWHM/2))²)

Δ2θ = 2θ_grid - 2θ_pico
```

- `η ∈ [0,1]` = fração Lorentziana (`η=0` puro Gaussiano, `η=1` puro Lorentziano)

#### 2c — Soma de picos
```
Icalc_unit(2θ_i) = Σ_hkl  intensity_rel_hkl · pV(2θ_i - 2θ_hkl, FWHM_hkl, η)
```

---

## Etapa 3 — Ajuste Linear: Escala + Background

### O que faz
Determina o **fator de escala** `S` e os **coeficientes de background** `b_k` que minimizam
a diferença ponderada entre `Iobs` e `S·Icalc_unit + Ibg`. Como ambos entram linearmente,
o problema se resolve com **mínimos quadrados ponderados** (uma única operação matricial,
sem iteração).

### Input
| Variável | Origem |
|----------|--------|
| `tth` | Etapa 1 |
| `Iobs` | Etapa 1 |
| `sigma` | Etapa 1 |
| `Icalc_unit` | Etapa 2 |

### Output
| Variável | Tipo | Descrição |
|----------|------|-----------|
| `S` | float | fator de escala ótimo |
| `Icalc` | `ndarray (N,)` | padrão ajustado = `S·Icalc_unit + Ibg` |

### Fórmulas

#### Background polinomial
```
Ibg(2θ) = Σ_{k=0}^{n-1}  b_k · ((2θ - μ) / std)^k
```
- Polinômio de grau `n-1` (padrão: `n=4` → grau 3)
- `2θ` normalizado: `(2θ - média) / desvio_padrão` (estabilidade numérica)

#### Mínimos quadrados ponderados
Monta matriz de design `A` (N × (1 + n)):
```
A = [ Icalc_unit | 1 | tth_norm | tth_norm² | tth_norm³ ]
```

Vetor de pesos:
```
w_i = 1 / σ_i²
```

Resolve:
```
(Aᵀ W A) · x = Aᵀ W · Iobs
```

Onde `W = diag(w)`. Resultado: `x = [S, b_0, b_1, b_2, b_3]`.

Implementado via `scipy.linalg.lstsq` com pré-multiplicação por `√W`.

---

## Etapa 4 — Figures of Merit (FoM)

### O que faz
Quantifica a qualidade do ajuste entre `Iobs` e `Icalc` final.

### Input
| Variável | Origem |
|----------|--------|
| `Iobs` | Etapa 1 |
| `sigma` | Etapa 1 |
| `Icalc` | Etapa 3 |
| `n_params` | número de parâmetros livres (1 escala + n_bg coeficientes) |

### Output
| Métrica | Bom ajuste | Péssimo ajuste |
|---------|-----------|---------------|
| `Rwp` | < 0.10 | > 0.30 |
| `Rp` | < 0.08 | > 0.25 |
| `χ²` | ≈ 1.0 | >> 1 |

### Fórmulas

#### R_wp — Weighted Profile R-factor (mais importante)
```
Rwp = √( Σ w_i(Iobs_i - Icalc_i)² / Σ w_i·Iobs_i² )
```

#### R_p — Profile R-factor
```
Rp = Σ |Iobs_i - Icalc_i| / Σ |Iobs_i|
```

#### R_exp — Expected R-factor (limite estatístico mínimo)
```
Rexp = √( (N - P) / Σ w_i·Iobs_i² )
```
- `N` = número de pontos
- `P` = número de parâmetros livres

#### χ² — Goodness of Fit
```
χ² = (Rwp / Rexp)²
```
- `χ² = 1.0` → ajuste perfeito (limitado apenas pelo ruído estatístico)
- `χ² >> 1` → modelo não descreve os dados

---

## Etapa 5 — Ranking de Candidatos

### O que faz
Roda as etapas 2–4 para cada um dos 20 candidatos e ordena pelo `Rwp`.

### Input
- `candidates.csv`: colunas `cod_id`, `peak_matches`, `reflections`
- `Cu_synthetic.xye`: dados experimentais

### Output
Tabela ordenada:
```
cod_id    Rwp     Rp      chi2    S
-------   -----   -----   -----   ------
1234567   0.043   0.031   1.12    0.0821
2345678   0.187   0.142   4.31    0.0012
...
```

### Critério de seleção
- `Rwp < 0.15` + `χ² < 3` → candidato viável para refinamento completo (Etapa 2 de Rietveld)
- Candidato com menor `Rwp` = melhor correspondência de fase

---

## Parâmetros Fixos neste Teste

| Parâmetro | Valor | Descrição |
|-----------|-------|-----------|
| `λ` | 1.54056 Å | comprimento de onda Cu Kα1 |
| `U` | 0.01 | Caglioti U |
| `V` | -0.002 | Caglioti V |
| `W` | 0.005 | Caglioti W |
| `η` | 0.5 | mistura pseudo-Voigt |
| `n_bg` | 4 | grau do polinômio de background |

> No refinamento completo (Stage 2), esses parâmetros seriam refinados via otimização
> não-linear (Levenberg-Marquardt ou L-BFGS-B).

---

## Dependências Python

```
numpy
scipy
pandas
```

---

## Arquivos Relevantes

| Arquivo | Conteúdo |
|---------|----------|
| `Cu_synthetic.xye` | padrão experimental sintético Cu FCC, Cu Kα |
| `data-*.csv` | 20 candidatos do PostgreSQL (cod_id + reflections JSONB) |
| `xrd_schema_setup.py` | pipeline que gerou `xrd_analysis.reference_patterns` |
