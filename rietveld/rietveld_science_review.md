# Scientific Review Prompt — Rietveld Phase Identification Pipeline

## Instruções de Uso

Cole este documento no início de uma conversa de revisão científica. Diga ao revisor:

> "Você é um revisor científico de uma pipeline de identificação de fases por difração de raios X (XRD). Revise o código/resultados da **[etapa]** contra a metodologia e referências abaixo. Aponte erros de implementação, violações físicas ou desvios das equações canônicas. Seja preciso: cite a equação, a referência e o desvio encontrado."

---

## Contexto da Pipeline

**Objetivo:** Identificar a fase cristalina que melhor descreve um difratograma experimental de pó (Cu Kα, λ = 1.54056 Å), comparando com uma lista de candidatos (cod_id + reflexões calculadas).

**Abordagem:** Rietveld simplificado — ajuste linear de escala + background, sem refinamento de parâmetros estruturais. Fase ID por ranking de Rwp.

**Dados de entrada:**
- `Cu_synthetic.xye`: 2θ (20°–120°), Iobs, σ; 5005 pontos; passo 0.02°
- `data-1782394014136.csv`: 20 candidatos com `cod_id`, `reflections` [{h,k,l, two_theta, intensity_rel, d_hkl, multiplicity, F_sq}]

**Limitações conhecidas desta implementação:**
- Parâmetros de perfil (U, V, W, η) fixos — não refinados
- Factor de Lorentz-polarização já incluso em `intensity_rel` (calculado externamente)
- Background por polinômio de Chebyshev/potência, não Legendre
- Sem correção de assimetria axial de pico (Finger-Cox-Jephcoat)

---

## Referências Canônicas

| ID | Referência completa | Relevância |
|----|---------------------|------------|
| **[R69]** | Rietveld, H.M. (1969). "A profile refinement method for nuclear and magnetic structures." *J. Appl. Crystallogr.* **2**(2), 65–71. doi:[10.1107/S0021889869006558](https://doi.org/10.1107/S0021889869006558) | Fundação do método |
| **[Y93]** | Young, R.A. (ed.) (1993). *The Rietveld Method*. IUCr/Oxford University Press. ISBN 0-19-855577-6 | Referência mestre: equações de R-factors, WLS, background |
| **[CPR58]** | Caglioti, G., Paoletti, A. & Ricci, F.P. (1958). "Choice of collimators for a crystal spectrometer for neutron diffraction." *Nucl. Instrum.* **3**(4), 223–228. doi:[10.1016/0369-643X(58)90029-X](https://doi.org/10.1016/0369-643X(58)90029-X) | FWHM angular dependence |
| **[TCH87]** | Thompson, P., Cox, D.E. & Hastings, J.B. (1987). "Rietveld refinement of Debye–Scherrer synchrotron X-ray data from Al₂O₃." *J. Appl. Crystallogr.* **20**(2), 79–83. doi:[10.1107/S0021889887087090](https://doi.org/10.1107/S0021889887087090) | Pseudo-Voigt com FWHM separados G e L |
| **[T06]** | Toby, B.H. (2006). "R factors in Rietveld analysis: How good is good enough?" *Powder Diffraction* **21**(1), 67–70. doi:[10.1154/1.2179804](https://doi.org/10.1154/1.2179804) | Interpretação de R-factors; thresholds |
| **[MCC99]** | McCusker, L.B., Von Dreele, R.B., Cox, D.E., Louër, D. & Scardi, P. (1999). "Rietveld refinement guidelines." *J. Appl. Crystallogr.* **32**(1), 36–50. doi:[10.1107/S0021889898009856](https://doi.org/10.1107/S0021889898009856) | Guidelines práticos IUCr |

---

## Etapa 1 — Parse do Difratograma (`.xye`)

### Especificação

Arquivo `.xye` = dados de difração de pó no formato passo-a-passo (step-scan):

```
2θ   Iobs   σ
```

**Regras de leitura:**
- Ignorar linhas `#` (comentários) e vazias
- 3 colunas → `(2θ, Iobs, σ)` direto
- 2 colunas → estimar `σ = √max(Iobs, 1.0)` (aproximação de Poisson para contagens)
- Substituir `σ ≤ 0` por `√max(Iobs, 1.0)` (evita divisão por zero em weights)
- Ordenar por 2θ crescente (defensivo)

### Base Física

σ como raiz de Poisson: para detectores de contagem, variância ≈ N_counts (estatística de Poisson). Portanto σ ≈ √I quando I em contagens brutas. Se I já normalizado, a aproximação é válida apenas como floor de incerteza.

**Ref:** [Y93] Cap. 1, p. 5; [MCC99] p. 38 (discussão sobre σ em dados de pó)

### Critérios de Revisão

| # | Verificar | Erro esperado se violado |
|---|-----------|--------------------------|
| 1.1 | `σ > 0` em todos os pontos | weights `wᵢ = 1/σᵢ² → ∞`, Rwp diverge |
| 1.2 | 2θ estritamente crescente | build_icalc_unit assume grid ordenado para cutoff |
| 1.3 | Nenhum `Iobs < 0` após parse | R-factors fisicamente sem sentido |
| 1.4 | `σ = √max(I, 1)` quando 2 colunas | não `σ = 1` constante (enviesa weights) |

---

## Etapa 2 — Cálculo do Padrão (`pattern_calc.py`)

### 2a. FWHM Caglioti

**Equação [CPR58], eq. 2:**

$$\Gamma^2(2\theta) = U \tan^2\!\theta + V \tan\!\theta + W$$

onde θ = 2θ/2 (ângulo de Bragg em **radianos**, mas 2θ em graus).

```python
theta_rad = np.radians(two_theta_deg / 2.0)
tan_t = np.tan(theta_rad)
fwhm2 = U * tan_t**2 + V * tan_t + W
fwhm = np.sqrt(np.maximum(fwhm2, 1e-8))   # floor para evitar FWHM imaginário
```

**Domínio de validade:** U, V, W devem satisfazer `fwhm² > 0` para todo 2θ no range experimental. Com U > 0, W > 0 e V² < 4UW isso é garantido (discriminante negativo).

**Valores típicos laboratoriais (Cu Kα, difratômetro Bragg-Brentano):**
- U ≈ 0.002–0.05 graus²
- V ≈ −0.001 a −0.010 graus²
- W ≈ 0.002–0.015 graus²

**Ref:** [CPR58]; [Y93] p. 14 eq. 2.7; [MCC99] p. 40

### Critérios de Revisão — FWHM

| # | Verificar | |
|---|-----------|--|
| 2.1 | `theta_rad = np.radians(two_theta_deg / 2.0)` — dividir por 2 antes de converter | erro de fator 2 no argumento de tan |
| 2.2 | floor `fwhm² ≥ ε > 0` com `np.maximum(fwhm2, ε)` antes de `sqrt` | `fwhm² < 0` para V muito negativo → NaN |
| 2.3 | FWHM em graus (mesmo sistema de unidades que 2θ) | inconsistência de unidades corrompe perfil |
| 2.4 | FWHM cresce monotonamente para 2θ > 30° (tipicamente) | check básico de sanidade física |

---

### 2b. Perfil Pseudo-Voigt

**Equação TCH simplificada [TCH87], usada nesta implementação:**

$$\text{pV}(\Delta) = \eta \cdot L(\Delta) + (1 - \eta) \cdot G(\Delta), \quad 0 \leq \eta \leq 1$$

**Gaussiana normalizada ao pico:**
$$G(\Delta) = \exp\!\left(-\frac{\Delta^2}{2\sigma_G^2}\right), \quad \sigma_G = \frac{\Gamma}{2\sqrt{2\ln 2}}$$

equivalente a:

$$G(\Delta) = \exp\!\left(-\frac{4\ln 2 \cdot \Delta^2}{\Gamma^2}\right)$$

**Lorentziana normalizada ao pico:**
$$L(\Delta) = \frac{1}{1 + \left(\frac{2\Delta}{\Gamma}\right)^2} = \frac{\Gamma^2}{\Gamma^2 + 4\Delta^2}$$

onde Δ = 2θ − 2θ_hkl (em graus), Γ = FWHM em graus.

**Importante:** As funções G e L são normalizadas a **1.0 no centro** (Δ = 0), não à área unitária. Isso é correto para esta implementação porque `intensity_rel` já é absoluto relativo ao pico mais intenso.

**Implementação:**
```python
delta = tth_grid - two_theta_peak
sigma_g = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
G = np.exp(-delta**2 / (2.0 * sigma_g**2))
L = 1.0 / (1.0 + (delta / (fwhm / 2.0))**2)
return eta * L + (1.0 - eta) * G
```

**Nota sobre TCH original:** Thompson et al. (1987) derivam η a partir dos componentes Γ_G e Γ_L separados via convolução aproximada. A presente implementação usa η fixo global — simplificação aceitável para identificação de fase, mas não para refinamento estrutural de precisão.

**Ref:** [TCH87] eq. 1–3; [Y93] p. 14; [MCC99] p. 40

### Critérios de Revisão — Pseudo-Voigt

| # | Verificar | |
|---|-----------|--|
| 2.5 | `pV(0) = 1.0` (η·1 + (1-η)·1 = 1) | normalização pico |
| 2.6 | `pV` simétrico: `pV(-Δ) = pV(Δ)` | propriedade física fundamental |
| 2.7 | G usa `fwhm/(2√(2ln2))` para σ_G — **não** `fwhm/2.355` (numericamente equivalente mas verificar) | erro de fator se confundir σ com HWHM |
| 2.8 | L usa `(fwhm/2)²` no denominador — **não** `fwhm²` | erro de fator 4 |
| 2.9 | `0 ≤ η ≤ 1` | fora do range → componentes negativos |

---

### 2c. Intensidade Calculada Unitária

**Equação:**

$$I_{\text{calc,unit}}(2\theta_i) = \sum_{hkl} I_{\text{rel},hkl} \cdot \text{pV}(2\theta_i - 2\theta_{hkl},\, \Gamma_{hkl})$$

onde `intensity_rel` já incorpora: multiplicidade m_hkl, |F_hkl|², fator de Lorentz-polarização L_p. O fator de escala S é aplicado externamente pelo ajuste linear.

**Cutoff:** Picos com `|2θ_i - 2θ_hkl| > cutoff_fwhm · Γ_hkl` ignorados para eficiência (contribuição negligenciável).

**Ref:** [R69] eq. 1–2; [Y93] p. 13 eq. 2.1

### Critérios de Revisão — Icalc_unit

| # | Verificar | |
|---|-----------|--|
| 2.10 | Picos com `intensity_rel ≤ 0` pulados | evita contribuição negativa |
| 2.11 | Picos fora de [2θ_min − cutoff·Γ, 2θ_max + cutoff·Γ] ignorados | não apenas dentro do grid |
| 2.12 | `Icalc_unit ≥ 0` em todos os pontos | fisicamente obrigatório |
| 2.13 | Posição do máximo de Icalc_unit corresponde ao pico mais intenso do candidato | sanidade física |

---

## Etapa 3 — Ajuste Linear WLS (`linear_fit.py`)

### Especificação

Minimizar o funcional de mínimos quadrados ponderados:

$$\chi^2_{\text{fit}} = \sum_{i=1}^{N} w_i \left[ y_i^{\text{obs}} - S \cdot I_{\text{calc,unit}}(2\theta_i) - I_{\text{bg}}(2\theta_i) \right]^2$$

**Pesos:** $w_i = 1/\sigma_i^2$

**Background polinomial normalizado:**

$$I_{\text{bg}}(2\theta) = \sum_{k=0}^{n_{bg}-1} b_k \cdot \tilde{\theta}^k, \quad \tilde{\theta} = \frac{2\theta - \mu_{2\theta}}{\sigma_{2\theta}}$$

A normalização de 2θ evita ill-conditioning numérico da matriz de design para polinômios de grau > 2.

**Design matrix:**

$$\mathbf{A} = \left[ I_{\text{calc,unit}} \;\Big|\; \mathbf{1} \;\Big|\; \tilde{\boldsymbol{\theta}} \;\Big|\; \tilde{\boldsymbol{\theta}}^2 \;\Big|\; \cdots \;\Big|\; \tilde{\boldsymbol{\theta}}^{n_{bg}-1} \right] \in \mathbb{R}^{N \times (1 + n_{bg})}$$

**Solução:**

$$\mathbf{x} = \left[ S,\, b_0,\, b_1,\, \ldots,\, b_{n_{bg}-1} \right]^T = (\mathbf{A}^T \mathbf{W} \mathbf{A})^{-1} \mathbf{A}^T \mathbf{W} \mathbf{y}$$

Na prática via `scipy.linalg.lstsq` com `Aw = A * √w` e `bw = y * √w`.

**Importante:** O scale factor S pode ser negativo se o candidato não se correlacionar positivamente com o difratograma. Isso é um sinal de mau candidato, não erro de implementação.

**Ref:** [Y93] Cap. 2 p. 18–22; [MCC99] p. 37–38

### Critérios de Revisão — Linear Fit

| # | Verificar | |
|---|-----------|--|
| 3.1 | `w = 1/σ²` — **não** `w = 1/σ` | weights errados comprometem Rwp |
| 3.2 | `sqrt_w[:, None]` aplicado a A (não apenas a b) | WLS via pre-multiplicação de ambos lados |
| 3.3 | `tth_norm` calculado com `std` de tth (não range) | normalização estatisticamente correta |
| 3.4 | `std = max(std, 1e-6)` — floor para evitar divisão por zero | dataset com 1 ponto |
| 3.5 | `Icalc = S * Icalc_unit + Ibg` — not `S * (Icalc_unit + Ibg)` | background não escala junto com sinal |
| 3.6 | n_params = 1 (scale) + n_bg (background coeffs) | usado em Rexp; sub-contagem → χ² inflado |
| 3.7 | `scipy.linalg.lstsq` — não `np.linalg.solve` (que exige matriz quadrada) | robustez numérica |

---

## Etapa 4 — Figures of Merit (`fom.py`)

### Especificação

**Rwp — Weighted Profile R-factor:**

$$R_{wp} = \sqrt{\frac{\sum_i w_i (y_i^{\text{obs}} - y_i^{\text{calc}})^2}{\sum_i w_i (y_i^{\text{obs}})^2}}$$

**Rp — Profile R-factor:**

$$R_p = \frac{\sum_i |y_i^{\text{obs}} - y_i^{\text{calc}}|}{\sum_i |y_i^{\text{obs}}|}$$

**Rexp — Expected R-factor:**

$$R_{\text{exp}} = \sqrt{\frac{N - P}{\sum_i w_i (y_i^{\text{obs}})^2}}$$

onde N = número de pontos, P = número de parâmetros livres (= 1 + n_bg nesta pipeline).

**χ² — Goodness of Fit:**

$$\chi^2 = \left(\frac{R_{wp}}{R_{\text{exp}}}\right)^2 = \frac{\sum_i w_i (y_i^{\text{obs}} - y_i^{\text{calc}})^2}{N - P}$$

**Valor ideal:** χ² = 1.0 (estatisticamente esperado se σᵢ corretos e modelo perfeito).

**Ref:** [Y93] p. 22–26 eq. 2.5–2.7; [T06]; [MCC99] p. 40–41

### Thresholds de Qualidade [T06]

| Métrica | Excelente | Aceitável (fase ID) | Ruim |
|---------|-----------|---------------------|------|
| Rwp | < 0.08 | < 0.15 | > 0.25 |
| Rp | < 0.06 | < 0.12 | > 0.20 |
| χ² | 1.0–1.5 | < 3.0 | > 5.0 |

**Nota:** Para identificação de fase (não refinamento estrutural), thresholds são menos restritivos. O ranking relativo de Rwp entre candidatos é mais informativo que o valor absoluto. [T06] p. 69.

### Critérios de Revisão — FoM

| # | Verificar | |
|---|-----------|--|
| 4.1 | Rwp usa `w_i(y_obs)²` no denominador — **não** `w_i(y_calc)²` | definição canônica [Y93] eq. 2.6 |
| 4.2 | `N - P` no denominador de Rexp — **não** N | graus de liberdade corretos |
| 4.3 | `max(N - P, 1)` — floor para evitar Rexp = 0 | datasets minúsculos em testes |
| 4.4 | `max(sum_w_Iobs2, ε)` — floor em denominador | evita Rwp = ∞ se Iobs = 0 |
| 4.5 | χ² = (Rwp/Rexp)² — equivalente a Σw(diff²)/(N-P) | verificar algebricamente |
| 4.6 | Sanidade: `Icalc = Iobs → Rwp ≈ 0` e `Icalc = 0 → Rwp ≈ 1` | testes obrigatórios |
| 4.7 | Rwp, Rp, Rexp, χ² são sempre ≥ 0 | propriedade matemática |

---

## Etapa 5 — Pipeline Completa (`pipeline.py`)

### Especificação do Fluxo

```
parse_xye → [Icalc_unit, n_used] → [scale, Icalc] → FoM → rank
```

Para cada candidato `c` com `reflections`:
1. `build_icalc_unit(tth, c.reflections, U, V, W, η)` → `(Icalc_unit, n_used)`
2. `linear_fit(tth, Iobs, σ, Icalc_unit, n_bg)` → `(scale, Icalc)`
3. `calc_fom(Iobs, Icalc, σ, n_params=1+n_bg)` → `{Rwp, Rp, Rexp, χ²}`
4. Armazenar `CandidateResult`
5. Ordenar por `Rwp` crescente

**Critério de viabilidade:**
- Rwp < 0.15 **e** χ² < 3.0 [T06]

**Ref:** [R69]; [MCC99] p. 36–50

### Parâmetros Fixos da Implementação

| Parâmetro | Valor usado | Referência |
|-----------|-------------|------------|
| λ | 1.54056 Å (Cu Kα1) | NIST SRD; [MCC99] p. 37 |
| U | 0.01 grau² | Típico laboratorial [CPR58] |
| V | −0.002 grau² | Típico laboratorial [CPR58] |
| W | 0.005 grau² | Típico laboratorial [CPR58] |
| η | 0.5 | Mixing 50% Lorentz/Gauss [TCH87] |
| n_bg | 4 | Grau do polinômio de background [Y93] |

### Critérios de Revisão — Pipeline

| # | Verificar | |
|---|-----------|--|
| 5.1 | `n_params = 1 + n_bg` passado para `calc_fom` — consistente com os parâmetros do fit | |
| 5.2 | `scale` pode ser negativo — não filtrar | candidato com correlação negativa é Rwp ruim, não erro |
| 5.3 | Candidatos ordenados por Rwp crescente (melhor primeiro) | |
| 5.4 | `db_client=None` funciona (metadados opcionais) | modo teste sem DB |
| 5.5 | Candidatos com `n_peaks_used = 0` incluídos mas com Rwp ≈ 1 (Icalc_unit = 0) | não crashar, reportar como ruim |

---

## Etapa 6 — Análise Multi-fase (resíduos de mineração)

### Contexto

Resíduos de mineração = misturas policristalinas (quartzo + calcita + hematita + argilominerais + ...). O difratograma é a **soma** das contribuições de K fases mais background. Identificação de fase única (Etapas 1–5) é insuficiente: a melhor fase isolada nunca explica todos os picos de uma mistura.

### Modelo

$$I_{\text{calc}}(2\theta_i) = \sum_{k=1}^{K} S_k \cdot I_{\text{calc,unit},k}(2\theta_i) + I_{\text{bg}}(2\theta_i)$$

Cada $I_{\text{calc,unit},k}$ é construído por `build_icalc_unit` independentemente (reuso da Etapa 2c). Background é **único e compartilhado** — não há um background por fase.

### Ajuste — WLS com restrição de não-negatividade

Os fatores de escala $S_k$ **devem ser ≥ 0** — uma fração de fase negativa é não-física numa mistura. O background pode ser negativo. Isso impede `scipy.linalg.lstsq` irrestrito (que permite $S_k < 0$). Usar **`scipy.optimize.lsq_linear`** com bounds por variável:

- colunas de fase ($S_1 ... S_K$): bound $[0, +\infty)$
- colunas de background ($b_0 ... b_{n_{bg}-1}$): bound $(-\infty, +\infty)$

$$\mathbf{A} = \left[ I_{\text{calc,unit},1} \,|\, \cdots \,|\, I_{\text{calc,unit},K} \,|\, \mathbf{1} \,|\, \tilde{\theta} \,|\, \cdots \right]$$

Pesos $w_i = 1/\sigma_i^2$ aplicados via $A_w = A \sqrt{w}$, $b_w = y\sqrt{w}$ (idêntico à Etapa 3).

**Ref:** [Y93] cap. 5; [BH88]

### Critérios de Revisão — Multi-fase

| # | Verificar | |
|---|-----------|--|
| 6.1 | design matrix tem K colunas de fase + n_bg colunas de background | |
| 6.2 | $S_k \geq 0$ imposto via `lsq_linear` bounds — **não** `lstsq` irrestrito | fração negativa não-física |
| 6.3 | background com lower bound $-\infty$ — não forçar $\geq 0$ | background subtraível |
| 6.4 | cada $I_{\text{calc,unit},k}$ via `build_icalc_unit` (reuso, mesmos U,V,W,η) | consistência de perfil |
| 6.5 | Rwp/χ² calculados sobre o modelo combinado, não por fase | FoM da mistura |
| 6.6 | Rwp_combinado ≤ Rwp da melhor fase isolada — **garantido por construção** se a melhor fase única ∈ conjunto (single-fase = ponto factível com demais scales=0) | sanidade, não estritamente < |
| 6.7 | Seleção das K fases por **rank de Rwp** (menor primeiro), **NÃO** por threshold de viabilidade single-fase | fase real de mistura tem Rwp single alto (explica só parte) → threshold a dropa |

---

## Etapa 7 — Quantificação de Fases (QPA)

### Hill-Howard (1987)

Fração em peso da fase $k$:

$$W_k = \frac{S_k (ZMV)_k}{\sum_{j=1}^{K} S_j (ZMV)_j}$$

onde $Z_k$ = unidades de fórmula por cela, $M_k$ = massa molar da unidade de fórmula, $V_k$ = volume da cela.

### Volume de cela (fórmula triclínica geral)

$$V = abc\sqrt{1 - \cos^2\alpha - \cos^2\beta - \cos^2\gamma + 2\cos\alpha\cos\beta\cos\gamma}$$

Reduz automaticamente para casos de maior simetria ($\alpha=\beta=\gamma=90° \Rightarrow V=abc$).

### Massa molar

$$M = \sum_{\text{átomos}} n_i \cdot A_i$$

$A_i$ = peso atômico padrão (IUPAC). Parse da string `formula` (ex: `"SiO2"`, `"Ca C O3"`, contagens fracionárias `"Br0.8"` permitidas).

### Critérios de Revisão — QPA

| # | Verificar | |
|---|-----------|--|
| 7.1 | $V$ pela fórmula triclínica geral — **não** assumir cela ortogonal | minerais não-cúbicos |
| 7.2 | $M$ de parse completo da formula (todos elementos, contagens fracionárias) | |
| 7.3 | $W_k = S_k(ZMV)_k / \sum_j$ — normalizado | Hill-Howard eq. 5 |
| 7.4 | $\sum_k W_k = 1.0$ | base 100% cristalina |
| 7.5 | $Z$, $a,b,c$, ângulos de `StructureMetadata` — checar não-nulos antes de usar | QPA exige DB ativo |
| 7.6 | **Base de intensidade comparável entre fases.** ✅ **VERIFICADO (2026): `intensity_rel` É normalizado por fase** — todas as 20 fases do CSV têm `ir_max = 100.0`. Logo $S_k$ do fit absorve o fator $\max_k$ → **QPA com $S_k$ direto está ERRADO**. QPA **DEVE** reconstruir padrão absoluto de $m \cdot \|F\|^2 \cdot L_p(\theta)$ (`multiplicity`, `F_sq` no CSV; $L_p$ analítico) e fitar nessa base. | confirmado via inspeção do CSV; ver T-024 |
| 7.7 | `molar_mass`: elemento desconhecido → **`raise`**, não skip silencioso | skip → M subcontado → QPA errado sem aviso |

### Reconstrução do padrão absoluto (modo QPA)

Como `intensity_rel` é normalizado por fase (verificado, 7.6), o modo QPA **não** usa `intensity_rel`. Reconstrói intensidade absoluta por reflexão:

$$I^{\text{abs}}_{hkl} = m_{hkl} \cdot |F|^2_{hkl} \cdot L_p(\theta_{hkl})$$

**Lorentz-polarização (Bragg-Brentano, sem monocromador):**

$$L_p(\theta) = \frac{1 + \cos^2 2\theta}{\sin^2\theta \cos\theta}$$

(Com monocromador de grafite: numerador $1 + \cos^2 2\theta_m \cos^2 2\theta$, $2\theta_m$ = ângulo do monocromador.)

`build_icalc_unit_absolute(tth, reflections)` soma perfis pseudo-Voigt ponderados por $I^{\text{abs}}_{hkl}$ em vez de `intensity_rel`. O `multi_phase_fit` nesse padrão → $S_k$ = escala de Rietveld verdadeira → Hill-Howard válido.

**DW (Debye-Waller) omitido** nesta reconstrução (exigiria $U_{iso}$ de `atomic_sites`). ✅ **VERIFICADO (2026):** reconstrução $m\cdot F_{sq}\cdot L_p$ **sem DW** recupera `intensity_rel` do CSV com **correlação = 1.0000** (exato a 3 casas até $2\theta\approx13°$) → DW ausente ou $\approx 1$ neste CSV → **base absoluta EXATA**, não aproximação. Caveat: se um CSV futuro incluir DW em `intensity_rel`, haveria viés dependente da fase em $2\theta$ alto — re-verificar a correlação ao trocar a fonte de dados.

### Critérios de Revisão — padrão absoluto

| # | Verificar | |
|---|-----------|--|
| 7.8 | $L_p(\theta)$ usa $\theta = 2\theta/2$ em radianos | mesmo erro de fator 2 que Caglioti |
| 7.9 | modo QPA fita em $m\cdot F_{sq}\cdot L_p$, **não** em `intensity_rel` | senão $S_k$ enviesado por $\max_k$ |
| 7.10 | mesmo padrão absoluto usado no fit E no ZMV (consistência) | DW omitido uniformemente |

### Limitações conhecidas (QPA)

- **Amorfo NÃO quantificado** sem padrão interno (método do spike) ou externo. $W_k$ são frações da **porção cristalina apenas**. Resíduos de mineração frequentemente têm fase amorfa (escória, vidro). [BH88]
- **Debye-Waller omitido** na reconstrução absoluta → viés em $2\theta$ alto, dependente da fase.
- **Orientação preferencial** (calcita, micas, argilas, gipsita) enviesa `intensity_rel` — correção March-Dollase **não aplicada**. [D86]
- **Microabsorção**: contraste de coeficiente de absorção $\mu$ entre fases (ex: hematita vs quartzo) enviesa QPA — correção Brindley **não aplicada**. [B45]
- `intensity_rel` com U,V,W fixos → QPA de **screening/triagem**, não certificação metrológica.

**Refs adicionais:**

- **[HH87]** Hill, R.J. & Howard, C.J. (1987). Quantitative phase analysis from neutron powder diffraction data using the Rietveld method. *J. Appl. Crystallogr.* **20**(6), 467–474. doi:10.1107/S0021889887086199
- **[BH88]** Bish, D.L. & Howard, S.A. (1988). Quantitative phase analysis using the Rietveld method. *J. Appl. Crystallogr.* **21**(2), 86–91. doi:10.1107/S0021889887009415
- **[D86]** Dollase, W.A. (1986). Correction of intensities for preferred orientation in powder diffractometry: application of the March model. *J. Appl. Crystallogr.* **19**(4), 267–272. doi:10.1107/S0021889886089458
- **[B45]** Brindley, G.W. (1945). The effect of grain or particle size on X-ray reflections from mixed powders. *Phil. Mag.* **36**, 347–369.
- **[C74]** Chung, F.H. (1974). Quantitative interpretation of X-ray diffraction patterns of mixtures (RIR / matrix-flushing). *J. Appl. Crystallogr.* **7**(6), 519–525. doi:10.1107/S0021889874010375

---

## Etapa 8 — Pré-seleção por Picos (Hanawalt Search)

### Contexto e Motivação

O banco COD completo (CuKα) tem ~500 k estruturas × ~80 reflexões cada = ~40 M pares $(2\theta, I_{rel})$. Rodar Rietveld em todos é inviável. A pré-seleção Hanawalt reduz O(500k) → O(50) antes do fit caro.

Método Hanawalt [H38]: dados os 3 picos mais intensos do difratograma experimental, busca no banco as fases cujas 3 reflexões mais intensas coincidem (dentro de tolerância Δ). Generalização moderna: conta todos os picos experimentais que encontram match em qualquer reflexão do candidato; ordena por `n_matches` decrescente; retém top-N.

### Detecção de Picos Experimentais

Entrada: `(tth, Iobs, sigma)`. Saída: `tth_peaks` — posições 2θ dos picos experimentais.

**Critério de proeminência [SV75]:** `scipy.signal.find_peaks` com `prominence ≥ k·median(σ)`. Proeminência mede quanto o pico se destaca do fundo local — robusto a background inclinado (pico de fundo largo não tem proeminência alta). Não usar `height` com threshold global: falha quando background varia ao longo do difratograma.

```python
from scipy.signal import find_peaks
prominence = prominence_sigma * np.median(sigma)  # típico: 5.0
min_dist   = max(1, int(min_distance_deg / step))
idx, _ = find_peaks(Iobs, prominence=prominence, distance=min_dist)
tth_peaks = tth[idx]
```

**Atenção:** fase amorfa gera *hump* largo (~20–35° para silicato amorfo, CuKα). O hump pode ter proeminência alta se o background polinomial não o absorver. Mitigação: aumentar `min_distance_deg` (≥ 2°) para rejeitar picos muito largos com base em `width`.

### Matching In-Memory (modo CSV)

Para cada candidato `c`, conta picos experimentais que encontram match:

```python
n_matched = sum(
    1 for p in peaks_tth
    if any(abs(p - r['two_theta']) <= tol_deg for r in c.reflections
           if r.get('intensity_rel', 0) > 0)
)
```

Complexidade: O(|peaks| × |reflections|) por candidato × N_candidatos. Aceitável para N < 100 k; inviável para COD inteiro.

### Matching via Banco de Dados (modo DB — MV)

Para o COD completo, materialized view `xrd_analysis.peak_fingerprints` pré-computa as top-30 reflexões de cada fase:

```sql
CREATE MATERIALIZED VIEW xrd_analysis.peak_fingerprints AS
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

CREATE INDEX pf_tth_idx ON xrd_analysis.peak_fingerprints (two_theta);
CREATE INDEX pf_cod_idx ON xrd_analysis.peak_fingerprints (cod_id);
```

**Índice B-tree em `two_theta`:** a query de matching usa `BETWEEN obs.tth - tol AND obs.tth + tol` — **não** `ABS(...)` — para habilitar uso de índice. `ABS()` não é indexável por B-tree.

**Query de busca:**
```sql
WITH obs AS (SELECT unnest(%s::float8[]) AS tth)
SELECT pf.cod_id, COUNT(*) AS n_matched
FROM xrd_analysis.peak_fingerprints pf
JOIN obs ON pf.two_theta BETWEEN obs.tth - %s AND obs.tth + %s
GROUP BY pf.cod_id
HAVING COUNT(*) >= %s
ORDER BY n_matched DESC
LIMIT %s;
```

Params: `(peaks_array, tol, tol, min_matches, top_n)`.

### Parâmetros Recomendados

| Parâmetro | Padrão | Justificativa |
|-----------|--------|---------------|
| `prominence_sigma` | 5.0 | SNR ≥ 5 σ → <1% falsos positivos para Poisson |
| `min_distance_deg` | 0.3° | Metade do FWHM típico laboratorial (0.1–0.2°) |
| `tol_deg` | 0.2° | Cobre zero-point error (±0.05°) + FWHM/2 + erro 2θ COD |
| `min_matches` | 3 | Critério Hanawalt clássico (3 linhas mais fortes) [H38] |
| `top_n` | 50 | Margem de segurança antes do Rietveld; reduzir se lento |
| `top_k_mv` | 30 | Top-30 reflexões por fase na MV (picos mais intensos) |

**Tolerância `tol_deg`:** aumentar para 0.3–0.4° em amostras com strain ou desalinhamento de instrumento. Diminuir para 0.1° apenas com calibração rigorosa.

### Referências

- **[H38]** Hanawalt, J.D., Rinn, H.W. & Frevel, L.K. (1938). Chemical analysis by X-ray diffraction. *Ind. Eng. Chem. Anal. Ed.* **10**, 457–512. doi:10.1021/ac50125a001
- **[SV75]** Sonneveld, E.J. & Visser, J.W. (1975). Automatic collection of powder data from photographs. *J. Appl. Crystallogr.* **8**(1), 1–7. doi:10.1107/S0021889875009508

### Critérios de Revisão — Hanawalt

| # | Verificar | Erro esperado se violado |
|---|-----------|--------------------------|
| 8.1 | `detect_peaks` usa `prominence ≥ k·median(σ)`, não `height` com threshold global | falsos positivos no hump amorfo; picos reais perdidos em regiões de alto background |
| 8.2 | `min_distance_deg` imposto para evitar múltipla detecção do mesmo pico | conta o mesmo pico 2× → `n_matched` inflado |
| 8.3 | Matching em 2θ (não d-spacing) — `abs(p - r['two_theta']) ≤ tol` | d-spacing exigiria converter λ, fonte de erros adicionais |
| 8.4 | MV usa top-K por `intensity_rel` (reflexões mais intensas por fase) | picos fracos/ausentes por orientação preferencial → MV com picos aleatórios não discrimina |
| 8.5 | SQL usa `BETWEEN obs.tth-tol AND obs.tth+tol` — **não** `ABS(a-b)<tol` | `ABS()` desativa índice B-tree → varredura sequencial em 15M linhas |
| 8.6 | `min_matches ≥ 3` — não retornar fase com 0–2 matches | matches espúrios; fase aleatória tem ~1–2 matches por acidente |
| 8.7 | modo CSV e modo DB dão mesma ordenação relativa para candidatos em comum | inconsistência sugere bug em um dos caminhos |

---

## Etapa 8b — Matching radiação-agnóstico por d-spacing

### Motivação

`reference_patterns` armazena `two_theta` calculado por Bragg para a **λ específica** de cada entrada COD. Resultado empírico (2026): das 530.793 fases com `has_intensities=TRUE`, só **132.545 são CuKα** — a maioria (344k) é **MoKα** (λ=0.71073), usada em monocristal COD. Uma MV filtrada por CuKα descarta 75% do banco.

`d_hkl` é **intrínseco do cristal** — independe de λ (Bragg: `λ = 2d·sinθ`; `d` vem só de cela+hkl). Fingerprint em d-space → **1 MV cobre Cu/Co/Cr/Mo**, 530k fases. É o método Hanawalt **canônico**: o índice ICDD/PDF é tabela de pares (d, I), nunca 2θ [H38].

**Relevância Saint-Gobain:** resíduo de mineração é rico em Fe (hematita, goethita). Fe **fluoresce sob CuKα** → background alto. Labs usam CoKα/CrKα em amostra ferrosa. MV d-agnóstica é requisito, não luxo.

### Conversão 2θ → d (modo busca)

A amostra mede `two_theta` na sua própria radiação `λ_sample`. Converte cada pico para d antes do matching:

$$d_{obs} = \frac{\lambda_{sample}}{2\sin(\theta)}, \quad \theta = \frac{2\theta_{obs}}{2}$$

Match contra `peak_fingerprints.d_hkl` (B-tree, BETWEEN).

### Propagação de erro — tolerância NÃO é uniforme em d

Diferenciando Bragg ($d = \lambda / (2\sin\theta)$, $\theta = 2\theta/2$):

$$\frac{\Delta d}{d} = -\tfrac{1}{2}\cot(\theta)\,\Delta(2\theta) \quad [\Delta(2\theta)\text{ em rad}]$$

$\cot\theta$ diverge em ângulo baixo → uma tolerância angular **fixa** mapeia para janela d **fortemente variável**:

| 2θ (CuKα) | d | Δd para tol=0.2° | Δd/d |
|-----------|---|------------------|------|
| 12° | 7.37 Å | 0.122 Å | 1.66% |
| 60° | 1.54 Å | 0.0047 Å | 0.30% |

Janela ~25× mais larga em ângulo baixo. Logo:
- `tol_d` **absoluto fixo** (Å) → ERRADO (estreito demais em d alto, largo demais em d baixo)
- `Δd/d` **relativo fixo** → também ERRADO (não constante vs 2θ)

**Correto:** converter o tolerance angular em janela d **por pico**, via Bragg nas bordas:
```python
d_lo = λ_sample / (2*sin(radians((2θ_obs + tol_deg)/2)))
d_hi = λ_sample / (2*sin(radians((2θ_obs - tol_deg)/2)))
# match: pf.d_hkl BETWEEN d_lo AND d_hi
```
Preserva o erro angular físico do difratômetro (≈constante em 2θ) sem aproximar.

### Seleção do pattern por fase (best_pattern) — viés de cobertura

A MV pega 1 pattern por fase. `reflections` foi gerado com corte 2θ∈[5°,90°] **na λ de origem** → o range de d capturado depende de λ:

| λ origem | d_max @ 2θ_min=5° |
|----------|-------------------|
| MoKα 0.71 | 8.15 Å |
| CuKα 1.54 | 17.7 Å |

Reflexão **basal de argilominerais** (frequentes em resíduo mineral) acima de 8 Å some no pattern MoKα:

| Mineral | d(basal) | MoKα capta? |
|---------|----------|-------------|
| Caulinita / gipsita | 7.2 / 7.6 Å | ✓ |
| Ilita, mica | 10.0 Å | ✗ |
| Esmectita | 12–15 Å | ✗ |
| Clorita | 14.1 Å | ✗ |

`ORDER BY n_reflections DESC` escolhe MoKα (mais picos, mas d alto truncado) — **oposto** do que minério com argila precisa. λ maior (CuKα) capta d maior.

**Opções (ordem de rigor):**
1. **Recompute de d direto da cela** (cell params + hkl, λ-independente, range completo) — sem artefato de corte. Rigoroso; custa 1 passe.
2. **Preferir maior λ** no `best_pattern` (`ORDER BY wavelength DESC`, tie-break n_reflections) — capta basal, pragmático.
3. Aceitar truncamento MoKα e documentar perda de basal >8 Å.

### Limitações conhecidas

- **top-30 por intensidade é λ-dependente**: `intensity_rel = m·F²·Lp(θ)`; `Lp` depende de θ→λ → ranking dos 30 mais fortes varia levemente entre λ. F² domina → tolerável para search-match qualitativo (ICDD usa 1 lista por fase em todos instrumentos).
- **Truncamento de d** conforme tabela acima se best_pattern = MoKα.
- d-matching qualitativo (triagem); fit/QPA continua em 2θ na λ da amostra (Etapas 6–7).

### Referências

- **[H38]** Hanawalt et al. (1938) — índice por d-spacing, não 2θ. doi:10.1021/ac50125a001
- **[K74]** Klug, H.P. & Alexander, L.E. (1974). *X-Ray Diffraction Procedures*, 2nd ed., Wiley. Cap. 5 (search-match por d), propagação de erro Bragg.

### Critérios de Revisão — d-spacing

| # | Verificar | Erro esperado se violado |
|---|-----------|--------------------------|
| 8.8 | MV armazena `d_hkl` (λ-independente), não `two_theta`, como chave de match | radiação ≠ CuKα descartada → 75% do COD perdido |
| 8.9 | pipeline converte picos 2θ da amostra → d via `λ_sample` antes do match | usa λ errada → d deslocado → matches falsos |
| 8.10 | tolerância aplicada como janela d **por pico** via Bragg nas bordas (`d_lo`/`d_hi`), NÃO `tol_d` fixo nem `Δd/d` fixo | janela 25× errada em ângulo baixo → falsos pos/neg |
| 8.11 | `best_pattern` prefere maior λ (ou recompute da cela), NÃO `max(n_reflections)` | basal de argila >8 Å truncada em pattern MoKα |
| 8.12 | `rad_symbol`/`wavelength` na MV = proveniência, NÃO filtro de match | filtrar por radiação re-fragmenta cobertura |
| 8.13 | fit/QPA (Etapas 6–7) continuam em 2θ na λ da amostra — d-match é só triagem | misturar bases d/2θ no fit → escala errada |

---

## Etapa 9 — Modo DB-only (busca no COD completo, sem CSV)

### Contexto

O CSV de 20 candidatos é fixture de teste. Em produção a amostra é desconhecida — não há lista pré-definida. O fluxo DB-only descobre os candidatos diretamente do PostgreSQL:

```
detect_peaks(XYE) → match_candidates_db(MV) → top-N cod_ids
    → fetch_reflections(cod_ids) → run_pipeline / run_multiphase
```

Substitui `load_candidates_csv` por `db_client.fetch_reflections`. O resto da pipeline é idêntico — a interface `list[CandidateInput]` não muda (Etapa 5).

### `fetch_reflections` — carga das reflexões

Após `match_candidates_db` retornar os `cod_ids` mais promissores, `fetch_reflections` carrega as reflexões completas de `xrd_analysis.reference_patterns`:

```sql
SELECT DISTINCT ON (cod_id) cod_id, reflections
FROM xrd_analysis.reference_patterns
WHERE cod_id = ANY(%s)
  AND has_intensities = TRUE
  AND wavelength BETWEEN 1.535 AND 1.546
ORDER BY cod_id, calculated_at DESC
```

**`DISTINCT ON (cod_id)` + `ORDER BY calculated_at DESC`:** uma fase pode ter múltiplos cálculos de padrão (re-runs). Sem `DISTINCT ON`, a mesma fase entra N vezes → fit duplicado, QPA enviesado. Pega o cálculo mais recente.

**Filtro `has_intensities = TRUE`:** sem |F|², `intensity_rel` é só multiplicidade × Lp — padrão sem discriminação. Fase entraria com Rwp ruim sem aviso.

**Filtro CuKα (1.535–1.546 Å):** posições 2θ dependem de λ (lei de Bragg). Misturar reflexões calculadas para λ diferente → picos deslocados → matching/fit errados.

### Consistência MV ↔ reference_patterns

A MV `peak_fingerprints` e a query de `fetch_reflections` **devem usar o mesmo filtro** (`has_intensities=TRUE`, CuKα). Senão: a MV indica um `cod_id` que `fetch_reflections` descarta → candidato fantasma (selecionado mas sem reflexões → silenciosamente ausente do fit).

### `--from-db` — modo CLI sem CSV

Flag `--from-db` substitui `candidates_csv` por busca direta no COD. Fluxo:

```
XYE → detect_peaks → match_candidates_db(MV, top_n)
    → fetch_reflections(cod_ids)
    → list[CandidateInput]         (mesma interface do CSV mode)
    → run_pipeline / run_multiphase
```

Implementado em `candidates_from_db()` (`pipeline.py`). CLI recebe `candidates_csv` como argumento opcional (`nargs='?'`); se `--from-db`, CSV ignorado.

**Guards obrigatórios (crit. 9.6):**
- `--from-db` + `--no-db` → `p.error()` imediato (não crash silencioso)
- `--from-db` sem `candidates_csv` → OK
- `candidates_csv` presente + `--from-db` → `[WARN]` CSV ignorado

**`--prefilter N` em modo DB-only:** mapeia para `top_n` em `match_candidates_db` (MV já é o filtro Hanawalt; não há segunda passagem CSV). Default `top_n=50` se `--prefilter` ausente.

**`--wavelength λ`:** propagado para `two_theta_to_d` em `match_candidates_db` e `detect_peaks`. Use `CoKα=1.7902` para amostras com Fe (evita fluorescência); `CrKα=2.2909` para Fe-rico pesado. Reflexões em `reference_patterns` sempre em d-space (λ-agnostico na MV).

**Consistência MV ↔ fetch_reflections:** `match_candidates_db` usa MV (qualquer λ, d-space). `fetch_reflections` filtra CuKα (`wavelength BETWEEN 1.535 AND 1.546`). Se amostra foi medida com CoKα mas `--wavelength` não foi passado, d-space match é correto (intrinseco), mas `fetch_reflections` ainda pega padrão CuKα para o fit — **comportamento correto**: d é λ-independente, o fit em 2θ usa a λ da amostra no `FIXED_PARAMS`.

### Critérios de Revisão — DB-only

| # | Verificar | Erro esperado se violado |
|---|-----------|--------------------------|
| 9.1 | `fetch_reflections` usa `DISTINCT ON (cod_id)` + `ORDER BY calculated_at DESC` | fase duplicada → fit/QPA enviesado |
| 9.2 | `fetch_reflections` filtra `has_intensities = TRUE` | padrão sem |F|² entra sem discriminação |
| 9.3 | `fetch_reflections` filtra CuKα (mesma janela λ da MV e do `fetch_metadata`) | reflexões de λ diferente → 2θ deslocado |
| 9.4 | MV e `fetch_reflections` usam filtro idêntico (consistência) | `cod_id` selecionado mas sem reflexões → candidato fantasma |
| 9.5 | reflexões JSONB convertidas para `list[dict]` (json.loads se string) | `build_icalc_unit` falha em string crua |
| 9.6 | `--from-db` + `--no-db` → `p.error()` imediato (não crash silencioso) | crash obscuro sem DB |
| 9.7 | `candidates_from_db` produz mesmo resultado que modo CSV para os mesmos cod_ids | divergência indica bug no caminho DB |

---

## Checklist de Sanidades Globais

Execute antes de considerar qualquer etapa aprovada:

```python
# 1. Sigma sempre positivo
assert (sigma > 0).all()

# 2. Profile peak = 1.0
assert abs(pseudo_voigt_profile(np.array([peak_pos]), peak_pos, fwhm=0.1, eta=0.5) - 1.0) < 1e-10

# 3. FWHM positivo em todo o range
tth_test = np.linspace(20, 120, 1000)
assert (caglioti_fwhm(tth_test, U=0.01, V=-0.002, W=0.005) > 0).all()

# 4. FoM de Icalc=Iobs
fom_perfeito = calc_fom(Iobs, Iobs, sigma, n_params=5)
assert fom_perfeito['Rwp'] < 1e-10

# 5. FoM de Icalc=0
fom_nulo = calc_fom(Iobs, np.zeros_like(Iobs), sigma, n_params=5)
assert abs(fom_nulo['Rwp'] - 1.0) < 0.01

# 6. Fit exato
scale, Icalc_fit = linear_fit(tth, 10*Icalc_unit + 50, sigma, Icalc_unit, n_bg=4)
assert abs(scale - 10.0) < 0.1

# 7. chi2 = (Rwp/Rexp)^2
fom = calc_fom(Iobs, Icalc, sigma, n_params=5)
assert abs(fom['chi2'] - (fom['Rwp']/fom['Rexp'])**2) < 1e-8
```

---

## Erros Comuns Documentados

| Erro | Sintoma | Diagnóstico |
|------|---------|-------------|
| `theta_rad = radians(two_theta)` sem `/2` | FWHM ~4× estreito | `fwhm(43°) ≈ 0.03°` em vez de `≈ 0.12°` |
| `L(Δ) = 1/(1 + Δ²/Γ²)` — sem fator 4 | Lorentziana 2× mais larga | FWHM medido ≠ Γ |
| `w = 1/σ` em vez de `1/σ²` | Rwp não converge ao mínimo WLS | candidatos não discriminados |
| `Rexp = sqrt(N/ΣwI²)` sem `−P` | χ² sistematicamente < 1 | parece melhor que é |
| `Icalc = S*(Icalc_unit + Ibg)` | Background escalado com fase | Icalc explode nos picos |
| `sigma_G = fwhm` em vez de `fwhm/(2√(2ln2))` | Gaussiana 2.35× mais larga | nenhuma discriminação entre candidatos |
| `pV = eta*G + (1-eta)*L` (invertido) | η=1 → puro Gaussiano (deveria ser Lorentziano) | profiles invertidos |

---

## Referências Bibliográficas Completas

1. **[R69]** Rietveld, H.M. (1969). A profile refinement method for nuclear and magnetic structures. *Journal of Applied Crystallography*, **2**(2), 65–71. https://doi.org/10.1107/S0021889869006558

2. **[CPR58]** Caglioti, G., Paoletti, A. & Ricci, F.P. (1958). Choice of collimators for a crystal spectrometer for neutron diffraction. *Nuclear Instruments*, **3**(4), 223–228. https://doi.org/10.1016/0369-643X(58)90029-X

3. **[TCH87]** Thompson, P., Cox, D.E. & Hastings, J.B. (1987). Rietveld refinement of Debye–Scherrer synchrotron X-ray data from Al₂O₃. *Journal of Applied Crystallography*, **20**(2), 79–83. https://doi.org/10.1107/S0021889887087090

4. **[Y93]** Young, R.A. (ed.) (1993). *The Rietveld Method*. International Union of Crystallography / Oxford University Press. ISBN 0-19-855577-6.

5. **[T06]** Toby, B.H. (2006). R factors in Rietveld analysis: How good is good enough? *Powder Diffraction*, **21**(1), 67–70. https://doi.org/10.1154/1.2179804

6. **[MCC99]** McCusker, L.B., Von Dreele, R.B., Cox, D.E., Louër, D. & Scardi, P. (1999). Rietveld refinement guidelines. *Journal of Applied Crystallography*, **32**(1), 36–50. https://doi.org/10.1107/S0021889898009856

7. **[HH87]** Hill, R.J. & Howard, C.J. (1987). Quantitative phase analysis from neutron powder diffraction data using the Rietveld method. *Journal of Applied Crystallography*, **20**(6), 467–474. https://doi.org/10.1107/S0021889887086199

8. **[BH88]** Bish, D.L. & Howard, S.A. (1988). Quantitative phase analysis using the Rietveld method. *Journal of Applied Crystallography*, **21**(2), 86–91. https://doi.org/10.1107/S0021889887009415

9. **[D86]** Dollase, W.A. (1986). Correction of intensities for preferred orientation in powder diffractometry: application of the March model. *Journal of Applied Crystallography*, **19**(4), 267–272. https://doi.org/10.1107/S0021889886089458

10. **[B45]** Brindley, G.W. (1945). The effect of grain or particle size on X-ray reflections from mixed powders. *Philosophical Magazine*, **36**, 347–369.

11. **[C74]** Chung, F.H. (1974). Quantitative interpretation of X-ray diffraction patterns of mixtures. *Journal of Applied Crystallography*, **7**(6), 519–525. https://doi.org/10.1107/S0021889874010375

12. **[H38]** Hanawalt, J.D., Rinn, H.W. & Frevel, L.K. (1938). Chemical analysis by X-ray diffraction: Classification and use of X-ray diffraction patterns. *Industrial & Engineering Chemistry — Analytical Edition*, **10**(9), 457–512. https://doi.org/10.1021/ac50125a001

13. **[SV75]** Sonneveld, E.J. & Visser, J.W. (1975). Automatic collection of powder data from photographs. *Journal of Applied Crystallography*, **8**(1), 1–7. https://doi.org/10.1107/S0021889875009508

14. **[K74]** Klug, H.P. & Alexander, L.E. (1974). *X-Ray Diffraction Procedures for Polycrystalline and Amorphous Materials*, 2nd ed. Wiley. ISBN 0-471-49369-4. (Cap. 5: search-match por d-spacing; propagação de erro em Bragg.)
