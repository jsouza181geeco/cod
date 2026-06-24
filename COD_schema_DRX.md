# COD — Esquema do Banco de Dados: Significado e Relevância para DRX

> Crystallography Open Database (COD) · 534 593 estruturas · Revisão SVN 306 581  
> Referência CIF: IUCr Core CIF Dictionary v2.4 · Notações: Hermann-Mauguin, Hall, Schoenflies

---

## Visão Geral das 23 Tabelas

| Prioridade | Tabela | Registros | Função principal |
|:---:|---|---:|---|
| ★★★★★ | `data` | 534 593 | Dados cristalográficos completos — núcleo do COD |
| ★★★★☆ | `spacegroups` | 530 | Tabela de grupos espaciais completa |
| ★★★☆☆ | `smiles` | 257 430 | Estrutura química (SMILES) por entrada |
| ★★★☆☆ | `fingerprints` | 257 430 | Fingerprints moleculares (busca por similaridade) |
| ★★☆☆☆ | `amcsd_x_cod` | 20 775 | Cross-ref: AMCSD (minerais) ↔ COD |
| ★★☆☆☆ | `pubchem_x_cod` | 169 593 | Cross-ref: PubChem ↔ COD |
| ★★☆☆☆ | `chemspider_x_cod` | 34 773 | Cross-ref: ChemSpider ↔ COD |
| ★★☆☆☆ | `mpod_x_cod` | 257 | Cross-ref: MPOD (propriedades materiais) ↔ COD |
| ★★☆☆☆ | `relations` | 3 | Tipos de relação entre bancos externos |
| ★★☆☆☆ | `databases` | 9 | Metadados dos bancos externos referenciados |
| ★☆☆☆☆ | `journals` | 1 497 | Periódicos científicos |
| ★☆☆☆☆ | `publishers` | 21 | Editoras |
| ★☆☆☆☆ | `jsequences` | 1 433 | Sequências de periódicos (títulos-mãe) |
| ★☆☆☆☆ | `jaltnames` | 2 001 | Nomes alternativos de periódicos |
| ★☆☆☆☆ | `successors` | 3 | Sucessões entre periódicos |
| ★☆☆☆☆ | `drugbank_x_cod` | 2 | Cross-ref: DrugBank ↔ COD |
| ★☆☆☆☆ | `rod_x_cod` | 1 108 | Cross-ref: ROD (filmes e superfícies) ↔ COD |
| ★☆☆☆☆ | `wikidata_x_cod` | 76 | Cross-ref: Wikidata ↔ COD |
| ★☆☆☆☆ | `wikipedia_x_cod` | 17 | Cross-ref: Wikipedia ↔ COD |
| ★☆☆☆☆ | `xrda_x_cod` | 23 | Cross-ref: XRDA ↔ COD |
| ☆☆☆☆☆ | `numbers` | 25 | Faixas de IDs COD por periódico (administrativa) |
| ☆☆☆☆☆ | `rdf_relations` | 1 | Relações RDF/semânticas |
| ☆☆☆☆☆ | `news` | 31 | Notícias do site COD |

---

## 1. Tabela `data` ★★★★★

**A única tabela necessária para 95% das análises de DRX.** Cada linha = uma estrutura cristalina publicada.  
Mapeada diretamente para campos CIF (Crystallographic Information File) padrão IUCr.

> **Para quem está construindo o software:** pense nessa tabela como o catálogo de "impressões digitais" de materiais. Cada linha é um material conhecido, com todas as suas medidas cristalográficas. Seu software vai comparar os dados do experimento do usuário com essas linhas para encontrar o material.

---

### 1.1 Identificação da Entrada

| Coluna | Tipo | Significado | 🔧 Utilidade no Software |
|---|---|---|---|
| `file` | INTEGER (PK) | **COD ID** — identificador único de 7 dígitos. Usado em citações, URLs e cross-refs. | **Chave primária de tudo.** Use como ID interno nos resultados. Link direto para a página do material: `https://www.crystallography.net/cod/{file}.html`. Mostre para o usuário poder verificar a fonte. |
| `acce_code` | CHAR(6) | Código de depósito no Cambridge Structural Database (CCDC). Permite cruzar COD ↔ CSD. | Secundário. Útil se usuário quiser verificar mesma estrutura em outro banco pago (CSD). Exibe como informação extra nos detalhes do resultado. |
| `svnrevision` | INTEGER | Revisão SVN da última atualização da entrada. | **Controle de atualização do banco local.** Use no script de sync incremental para saber quais entradas foram atualizadas desde o último download. Não exibir ao usuário final. |
| `date`, `time` | DATE, TIME | Data/hora da última modificação no repositório SVN. | Mesma função do `svnrevision` — auditoria interna de quando o dado foi atualizado. Pode exibir como "dados atualizados em X" na interface. |
| `onhold` | DATE | Se preenchido: data de embargo — entrada não pública ainda. NULL = dado público. | **Filtro obrigatório:** `WHERE onhold IS NULL`. Registros com `onhold` preenchido são dados privados ainda não publicados — nunca devem aparecer nos resultados. |
| `duplicateof` | INTEGER | Se esta entrada é duplicata, aponta para o COD ID canônico. | **Deduplicação de resultados.** Se retornar múltiplos resultados para o mesmo material, prefira os que têm `duplicateof IS NULL` (são os canônicos). Ou use `optimal` para pegar o melhor automaticamente. |
| `optimal` | INTEGER | Aponta para a estrutura mais precisa entre duplicatas de um mesmo composto. | **Auto-seleção do melhor resultado.** Se encontrar duplicatas de um material, siga o ponteiro `optimal` → essa é a estrutura mais precisa disponível. Use para "mostrar melhor resultado" na UI. |
| `status` | ENUM | `NULL` = dado limpo; `'warnings'` = avisos; `'errors'` = erros; `'retracted'` = retratado pelos autores. | **Filtro de qualidade mais importante do banco.** Sempre adicione `WHERE status IS NULL` nas queries. Registros com `'errors'` ou `'retracted'` são dados incorretos ou desacreditados — comparar com eles gera resultados errados. |
| `flags` | TEXT | `'has coordinates'` = estrutura completa; `'has disorder'` = desordem estrutural; `'has Fobs'` = dados de intensidade disponíveis. | **Filtro de completude.** Para simulação de difratograma, filtre `flags LIKE '%has coordinates%'` — sem coordenadas atômicas não dá para calcular intensidades dos picos. `'has disorder'` = estrutura com imperfeições → parâmetros menos confiáveis. |

---

### 1.2 Parâmetros de Célula Unitária ★★★★★

> **Analogia para leigo:** imagine que cada cristal é uma caixinha que se repete infinitamente. `a`, `b`, `c` são as dimensões dessa caixinha (largura, profundidade, altura). `alpha`, `beta`, `gamma` são os ângulos entre as paredes (a maioria dos cristais tem ângulos retos = 90°, mas alguns são inclinados). Esses 6 números juntos definem o formato único de cada material.

**Esses são os números que você vai comparar com os dados do experimento para identificar o material.**  
Lei de Bragg: `nλ = 2d·sinθ` — a posição de cada pico no difratograma depende diretamente de `a`, `b`, `c`, `α`, `β`, `γ`.

| Coluna | Tipo | Unidade | Significado | 🔧 Utilidade no Software |
|---|---|---|---|---|
| `a` | DOUBLE PRECISION | Å | Dimensão da célula no eixo **a** (1 Å = 0,0000001 mm). | **Comparação direta com experimento.** Valor medido no DRX do usuário → busca no banco com tolerância: `WHERE a BETWEEN :exp_a * 0.99 AND :exp_a * 1.01`. Tolerância de 1% é boa para início; ajuste conforme precisão do equipamento. |
| `b` | DOUBLE PRECISION | Å | Dimensão da célula no eixo **b**. | Idem `a`. Para sistemas cúbicos `a=b=c` → filtrar só `a` já basta. |
| `c` | DOUBLE PRECISION | Å | Dimensão da célula no eixo **c**. | Idem `a`. |
| `siga` | REAL | Å | Incerteza (precisão) da medição de `a`. Quanto menor, mais preciso. | **Filtro de qualidade opcional.** Se quiser só estruturas bem determinadas: `WHERE siga < 0.005`. Mas ~40% do banco não tem esse valor preenchido (dado antigo) — não use como filtro obrigatório ou vai perder muitos resultados válidos. |
| `sigb` | REAL | Å | Incerteza de `b`. | Idem `siga`. |
| `sigc` | REAL | Å | Incerteza de `c`. | Idem `siga`. |
| `alpha` | REAL | ° | Ângulo entre eixos b e c. 90° em maioria dos cristais. | **Comparação com experimento.** Para sistemas cúbico/tetragonal/ortorrômbico = sempre 90° → não precisa comparar. Para monoclínico/triclínico = diferente de 90° → incluir na busca. |
| `beta` | REAL | ° | Ângulo entre eixos a e c. Para monoclínico: único ângulo ≠ 90°. | Idem `alpha`. Monoclínico: `alpha=90, gamma=90, beta≠90` — só `beta` varia. |
| `gamma` | REAL | ° | Ângulo entre eixos a e b. 120° em hexagonais, 90° nos demais. | Idem `alpha`. |
| `sigalpha` | REAL | ° | Incerteza de `alpha`. | Raramente usado. Mesmo critério de `siga`. |
| `sigbeta` | REAL | ° | Incerteza de `beta`. | Raramente usado. |
| `siggamma` | REAL | ° | Incerteza de `gamma`. | Raramente usado. |
| `vol` | REAL | Å³ | **Volume total da caixinha** — calculado a partir de a,b,c,α,β,γ. Número único que representa o "tamanho" do cristal. | **Busca mais robusta que a+b+c separados.** Use como filtro primário: `WHERE vol BETWEEN :exp_vol * 0.98 AND :exp_vol * 1.02`. Vantagem: configurações alternativas do mesmo material têm a/b/c diferentes mas `vol` igual. Também útil para calcular densidade teórica. |
| `sigvol` | REAL | Å³ | Incerteza do volume. | **Filtro de qualidade relativo (melhor que siga):** `WHERE sigvol/vol < 0.005` (incerteza < 0,5%). Mais justo que comparar valor absoluto pois funciona para cristais grandes e pequenos. NULL em ~40% dos registros antigos — não obrigatório. |

**Sistemas cristalinos e restrições:**

| Sistema | sgNumber | Restrições em a,b,c,α,β,γ | Implicação para busca |
|---|---|---|---|
| Cúbico | 195–230 | a=b=c, α=β=γ=90° | Buscar só por `a` e `sgNumber` basta |
| Tetragonal | 75–142 | a=b≠c, α=β=γ=90° | Buscar `a` e `c` |
| Ortorrômbico | 16–74 | a≠b≠c, α=β=γ=90° | Buscar `a`, `b`, `c` |
| Hexagonal | 168–194 | a=b≠c, α=β=90°, γ=120° | Buscar `a` e `c` |
| Trigonal | 143–167 | varia por setting | Buscar `vol` + `sgNumber` |
| Monoclínico | 3–15 | a≠b≠c, α=γ=90°, β≠90° | Buscar `a`,`b`,`c`,`beta` |
| Triclínico | 1–2 | tudo diferente | Buscar todos os 6 parâmetros |

---

### 1.3 Grupo Espacial ★★★★★

> **Analogia para leigo:** se a célula unitária é a "caixinha", o grupo espacial é o "padrão de simetria" de como os átomos dentro dela se organizam. Existem exatamente 230 padrões possíveis no universo. O `sgNumber` é apenas o número do padrão na lista oficial (1–230). Saber o padrão é fundamental porque ele determina quais picos aparecem e quais somem no difratograma.

| Coluna | Tipo | Significado | 🔧 Utilidade no Software |
|---|---|---|---|
| `sg` | VARCHAR(32) | **Símbolo Hermann-Mauguin** — nome oficial do padrão de simetria. Ex: `P 21/c`, `F m -3 m`. | **Exibição para o usuário.** Mostre esse campo como "Grupo Espacial" na interface — é o que os cientistas reconhecem. **Não use para filtrar** — mesmo grupo pode ter grafias diferentes no banco (`'P 21/c'` vs `'P21/c'`). Para filtrar, use `sgNumber`. |
| `sgHall` | VARCHAR(64) | Notação Hall — versão algébrica não ambígua do grupo espacial. Ex: `-P 2ybc`. | **Uso técnico interno.** Necessário se seu software for calcular posição de picos ou gerar difratograma simulado. Software de difração (SHELX, GSAS, FullProf) usa esse formato. Para UI: não exibir ao usuário comum. |
| `sgNumber` | SMALLINT | **Número do grupo espacial (1–230)** — identificador numérico universal, sem ambiguidade de grafia. | **Filtro principal para grupo espacial.** Se o software do usuário identificou sgNumber = 225 no experimento → `WHERE sgNumber = 225`. Também use para categorizar por sistema cristalino (195–230 = cúbico, etc.). É inteiro → busca instantânea, sem problema de string matching. |

**Grupos espaciais mais comuns no COD:**

| sgNumber | Sistema | Típico em | Filtro rápido |
|---|---|---|---|
| 14 | Monoclínico | Maioria dos fármacos e compostos orgânicos | `sgNumber = 14` |
| 2 | Triclínico | Compostos complexos, baixa simetria | `sgNumber = 2` |
| 62 | Ortorrômbico | Perovskitas, minerais | `sgNumber = 62` |
| 225 | Cúbico FCC | Ouro, prata, cobre, alumínio, NaCl | `sgNumber = 225` |
| 229 | Cúbico BCC | Ferro (α), tungstênio, cromo | `sgNumber = 229` |
| 194 | Hexagonal HCP | Titânio, magnésio, grafite, ZnO | `sgNumber = 194` |

---

### 1.4 Conteúdo Químico ★★★★☆

> **Analogia para leigo:** esses campos descrevem "o que tem dentro" do material — os elementos químicos e suas proporções. São usados para buscar por nome ou composição quando o usuário sabe o que está procurando, mas não tem os parâmetros de célula.

| Coluna | Tipo | Significado | 🔧 Utilidade no Software |
|---|---|---|---|
| `nel` | VARCHAR(4) | **Quantidade de elementos químicos distintos** na fórmula. `'1'`= elemento puro (Fe, Al); `'2'`= binário (Fe₂O₃); `'3'`= ternário (BaTiO₃). | **Filtro de complexidade.** Permite ao usuário refinar: "buscar só metais puros" (`nel='1'`), "só compostos de dois elementos" (`nel='2'`). Rápido pois é string curta. Combine com busca por fórmula. |
| `formula` | VARCHAR(255) | **Fórmula química** como publicada pelo autor. Ex: `- Fe2 O3 -`. Inclui traços como delimitadores. | **Campo de busca por composição.** Use `ILIKE '%Fe%O%'` para buscar todos que contêm ferro e oxigênio. Aviso: formato inconsistente entre entradas antigas e novas — prefira `calcformula` para comparação programática. |
| `calcformula` | VARCHAR(255) | **Fórmula calculada** a partir da estrutura atômica — mais padronizada e confiável que `formula`. | **Melhor campo para comparação programática.** Se seu software precisar verificar se dois registros têm a mesma composição, compare `calcformula`. Quando difere de `formula`, prefira `calcformula`. |
| `cellformula` | VARCHAR(255) | **Fórmula de toda a célula unitária** = `formula × Z`. Ex: Fe₂O₃ com Z=4 → `Fe8 O12`. | **Para cálculo de densidade.** Você precisa da massa total da célula → use `cellformula` para somar os pesos atômicos. Menos útil para identificação de fase. |
| `commonname` | VARCHAR(1024) | Nome comum/trivial do material. Ex: `hematite`, `quartz`, `alumina`, `zirconia`. | **Campo de busca por nome popular.** Permite o usuário digitar "quartzo" → buscar `commonname ILIKE '%quartz%'`. Mostre esse nome no resultado como label principal — mais legível que a fórmula química. |
| `chemname` | VARCHAR(2048) | Nome IUPAC sistemático. Ex: `iron(III) oxide`. | **Campo de busca técnica.** Útil para usuários técnicos que conhecem nomenclatura IUPAC. Use `ILIKE` para busca parcial. Pode ser muito longo — truncar na exibição. |
| `mineral` | VARCHAR(255) | **Nome mineral** (quando aplicável). Ex: `hematite`, `calcite`, `rutile`, `corundum`. | **Campo mais útil para usuários de mineralogia/geologia.** Tem índice FULLTEXT → busca rápida. Mostre como destaque quando preenchido — indica material natural bem conhecido. `WHERE mineral IS NOT NULL` filtra só minerais. |

---

### 1.5 Conteúdo da Célula ★★★★☆

| Coluna | Tipo | Significado | 🔧 Utilidade no Software |
|---|---|---|---|
| `Z` | SMALLINT | **Número de "moléculas" da fórmula dentro de uma caixinha.** Ex: Z=4 para Fe₂O₃ = 4 unidades de Fe₂O₃ dentro de uma célula unitária. | **Necessário para calcular densidade teórica** do material: `densidade = (Z × massa_molar) / (6,022e23 × vol × 1e-30)`. Também valida se a estrutura faz sentido fisicamente. Se `Z IS NULL` → densidade não calculável. |
| `Zprime` | REAL | **Z' = Z dividido pela multiplicidade do grupo espacial.** Z'=1 = estrutura normal. Z'>1 = cristal com estrutura mais complexa (duas moléculas diferentes no mesmo cristal). Z'<1 = estrutura com alta simetria. | Informação avançada. Para o software: Z' muito diferente de 1 pode indicar estrutura mais complexa ou menos confiável. Na UI: pode mostrar como flag "estrutura complexa" se Z'>1. |

---

### 1.6 Condições Experimentais ★★★★☆

> **Para leigo:** esses campos dizem em que temperatura e pressão o material foi medido. Isso é crítico porque o mesmo material se expande com calor — a caixinha fica maior em temperaturas altas. Se o experimento do usuário foi feito a 25°C mas a referência no banco foi medida a -170°C, os parâmetros de célula vão ser diferentes — e o software precisa saber disso.

| Coluna | Tipo | Unidade | Significado | 🔧 Utilidade no Software |
|---|---|---|---|---|
| `celltemp` | REAL | K (Kelvin) | **Temperatura em que a célula foi medida.** 293 K ≈ temperatura ambiente (20°C). 100 K = experimento com resfriamento criogênico. | **Filtro de compatibilidade com experimento do usuário.** Se usuário mediu a 25°C (298 K), prefira referências com `celltemp BETWEEN 280 AND 320`. Estruturas com `celltemp < 200` têm parâmetros de célula menores que à temperatura ambiente — comparação direta vai dar falso negativo. |
| `sigcelltemp` | REAL | K | Incerteza da temperatura de célula. | Raramente usado. Pode ignorar no MVP. |
| `diffrtemp` | REAL | K | Temperatura durante o experimento de difração (pode ser diferente de `celltemp`). | Use `diffrtemp` quando `celltemp` for NULL — às vezes só um dos dois está preenchido. Mesma lógica de filtro. |
| `sigdiffrtemp` | REAL | K | Incerteza da temperatura de difração. | Ignorar no MVP. |
| `cellpressure` | REAL | kPa | Pressão na medição da célula. NULL = pressão ambiente (~101 kPa). | **Filtro para excluir estruturas de alta pressão.** Materiais medidos a alta pressão têm célula comprimida (menor) — não comparáveis com experimento ambiente. Filtro: `WHERE cellpressure IS NULL OR cellpressure < 200` (descarta alta pressão). |
| `sigcellpressure` | REAL | kPa | Incerteza da pressão de célula. | Ignorar no MVP. |
| `diffrpressure` | REAL | kPa | Pressão durante a difração. | Idem `cellpressure` quando esse for NULL. |
| `sigdiffrpressure` | REAL | kPa | Incerteza da pressão de difração. | Ignorar. |
| `thermalhist` | VARCHAR(255) | — | Histórico de tratamento térmico do material. Ex: `annealed at 800°C`, `quenched from melt`. | **Informação contextual para o usuário.** Mostra como o material foi preparado — relevante para entender se é fase estável ou metaestável. Exibir nos detalhes do resultado. Não filtre por esse campo (texto livre inconsistente). |
| `pressurehist` | VARCHAR(255) | — | Histórico de pressão aplicada. Ex: `synthesized at 15 GPa`. | Idem `thermalhist` — informação contextual. Mostra ao usuário para interpretar. |
| `compoundsource` | VARCHAR(255) | — | Origem do material. Ex: `natural`, `commercial`, `synthesized by sol-gel`. | **Contexto do experimento.** `natural` = mineral extraído; `commercial` = produto comprado; `synthesized` = produzido em lab. Útil para o usuário entender a procedência. Exibir nos detalhes. |

---

### 1.7 Método e Radiação ★★★★★

> **Para leigo:** DRX usa raios-X para "fotografar" o cristal. Mas existe mais de um tipo de raio-X (diferentes comprimentos de onda) e mais de um tipo de experimento. Esses campos dizem exatamente como a estrutura foi medida — e é fundamental que o método de referência seja compatível com o método do experimento do usuário.

| Coluna | Tipo | Significado | 🔧 Utilidade no Software |
|---|---|---|---|
| `method` | ENUM | Como a estrutura foi determinada: `'single crystal'` = cristal único grande (mais preciso); `'powder diffraction'` = pó do material (mais comum em laboratório); `'theoretical'` = calculado por computador (sem medição real). | **Filtro de compatibilidade — um dos mais importantes.** Para comparar com DRX de laboratório convencional (pó): `method IN ('single crystal', 'powder diffraction')`. Exclua `'theoretical'` para aplicações que exigem dado experimental. Mostre o método na UI como badge colorido. |
| `radiation` | VARCHAR(32) | Tipo de radiação usada: `'X-ray'`, `'neutron'`, `'electron'`. | **Filtro obrigatório para DRX convencional:** `WHERE radiation = 'X-ray'`. Nêutrons e elétrons geram difratogramas com intensidades completamente diferentes — não comparáveis diretamente. |
| `wavelength` | REAL | **Comprimento de onda dos raios-X usados (em Ångström).** Valores comuns: CuKα = 1,54056 Å (mais comum em lab); MoKα = 0,71073 Å (monocristais); CoKα = 1,78897 Å; CrKα = 2,28970 Å. | **Crítico para simular ou comparar posições de picos.** A posição angular (2θ) de cada pico depende do comprimento de onda: `2θ = 2·arcsin(λ/2d)`. Se usuário usou CuKα e referência usou MoKα, os picos aparecem em posições angulares diferentes — precisa converter. Use `wavelength` para calcular a posição esperada dos picos da referência no equipamento do usuário. |
| `radType` | VARCHAR(80) | Descrição detalhada da fonte e monocromador. Ex: `'Mo K\a'`, `'Cu K\a, graphite monochromator'`. | **Informação técnica complementar ao `wavelength`.** Exibir nos detalhes para usuários avançados. Use para confirmar qual linha espectral foi usada quando `wavelength` estiver NULL. |
| `radSymbol` | VARCHAR(20) | Símbolo abreviado. Ex: `'MoKα'`, `'CuKα'`. | **Para exibição na UI.** Label curto e legível para mostrar a fonte de radiação. Prefira esse para mostrar ao usuário sobre `radType`. |

---

### 1.8 Fatores de Qualidade da Estrutura ★★★☆☆

> **Para leigo:** quando os cientistas determinam uma estrutura cristalina, eles comparam os dados medidos com o modelo calculado. Esses fatores R medem o quanto o modelo coincide com o experimento — é como uma pontuação de erro. Quanto mais baixo, melhor. Pense como a margem de erro do modelo: R=0,05 significa 5% de erro médio.

| Coluna | Tipo | Definição | Intervalo | 🔧 Utilidade no Software |
|---|---|---|---|---|
| `Robs` | REAL | Erro médio do modelo vs. medição (reflexões fortes). | Bom: <0,05 · Aceitável: <0,10 · Ruim: >0,15 | **Principal filtro de qualidade estrutural.** `WHERE Robs < 0.08` retorna estruturas bem determinadas. Mostre como "pontuação de qualidade" na UI (ex: barra de progresso invertida). Use para ordenar resultados: `ORDER BY Robs ASC`. |
| `Rall` | REAL | Mesmo que Robs mas incluindo reflexões fracas. Sempre ≥ Robs. | Bom: <0,08 | **Filtro alternativo quando `Robs` for NULL.** Estruturas antigas às vezes só reportam `Rall`. |
| `Rref` | REAL | R para conjunto de referência interno. | — | Raramente usado sozinho. Fallback quando os outros forem NULL. |
| `wRobs` | REAL | Versão ponderada do Robs — dá mais peso às reflexões mais intensas. | Bom: <0,12 | **Filtro complementar.** Use `wRobs < 0.15` combinado com `Robs`. Mais sensível que Robs para detectar refinamentos problemáticos. |
| `wRall` | REAL | wR para todas as reflexões. | Bom: <0,15 | Alternativo ao `wRobs`. |
| `wRref` | REAL | wR de referência. | — | Fallback. |
| `RFsqd` | REAL | R baseado em intensidades ao quadrado — padrão em software moderno (SHELXL). | Bom: <0,10 | Estruturas modernas (pós-2000) geralmente só reportam esse. Use como alternativa ao Robs em buscas avançadas. |
| `RI` | REAL | Fator R de perfil (Rp) — específico para difração de pó, mede ajuste do perfil Rietveld. | Bom: <0,10 | **Usar como filtro de qualidade quando `method = 'powder diffraction'`** — é o equivalente ao Robs para esse tipo. |
| `gofall` | REAL | Goodness-of-fit: mede se o modelo é consistente internamente. Ideal = 1,0. | Ideal: ~1,0 · Problemático: >2,0 | **Filtro complementar de consistência:** `WHERE gofall < 2.0`. Valores muito altos (>3) indicam refinamento que não convergiu bem. Valores muito baixos (<0,5) indicam dados superestimados. |
| `gofobs` | REAL | GoF só para reflexões fortes. | ~1,0 ideal | Alternativo ao `gofall`. |
| `gofgt` | REAL | GoF para reflexões com I>2σ(I). | ~1,0 ideal | Alternativo. |
| `gofref` | REAL | GoF de referência. | — | Fallback. |

**Query de qualidade recomendada para o software:**
```sql
SELECT file, formula, sg, a, b, c, alpha, beta, gamma, vol, Z, wavelength
FROM data
WHERE status IS NULL                          -- sem erros detectados
  AND onhold IS NULL                          -- dado público
  AND radiation = 'X-ray'                     -- só raios-X
  AND method != 'theoretical'                 -- só experimental
  AND flags LIKE '%has coordinates%'          -- estrutura completa
  AND (Robs < 0.08 OR (Robs IS NULL AND RFsqd < 0.10))  -- qualidade mínima
ORDER BY Robs ASC NULLS LAST;
```

---

### 1.9 Publicação ★★☆☆☆

| Coluna | Tipo | Significado | 🔧 Utilidade no Software |
|---|---|---|---|
| `authors` | TEXT | Lista de autores da publicação original. | **Informação de rastreabilidade.** Exibir nos detalhes do resultado para o usuário verificar a fonte. Não usar para filtrar na busca principal. |
| `title` | TEXT | Título do artigo científico onde a estrutura foi publicada. | **Busca por texto livre** (tem índice FTS): `WHERE to_tsvector('english', title) @@ plainto_tsquery('zirconia phase')`. Exibir nos detalhes. |
| `journal` | VARCHAR(255) | Nome do periódico científico. | **Indicador indireto de qualidade.** Acta Crystallographica = revista especializada de alta credibilidade. Pode usar como fator de ranking de confiança. |
| `year` | SMALLINT | Ano de publicação. | **Filtro de modernidade.** Estruturas recentes (pós-2000) tendem a ter melhor precisão (detectores CCD, software moderno). Filtro útil: `WHERE year > 1990`. Exibir como informação no resultado. |
| `volume` | SMALLINT | Volume do periódico. | Informação bibliográfica. Exibir nos detalhes. Sem utilidade para filtragem. |
| `issue` | VARCHAR(10) | Número do fascículo. | Idem `volume`. |
| `firstpage`, `lastpage` | VARCHAR(20) | Páginas do artigo. | Informação bibliográfica completa. Exibir nos detalhes. |
| `doi` | VARCHAR(127) | **Link direto para o artigo original** (Digital Object Identifier). | **Link clicável nos resultados.** `https://doi.org/{doi}` → abre o artigo. Um dos campos mais valiosos para o usuário verificar a fonte. Sempre exibir quando não NULL. |
| `text` | TEXT | Texto bibliográfico completo concatenado (autores + título + revista + ano). | **Campo para busca textual geral** (tem índice GIN): `WHERE to_tsvector('english', text) @@ plainto_tsquery('hematite')`. Fallback quando busca por `mineral` ou `commonname` não retorna resultado. |

---

## 2. Tabela `spacegroups` ★★★★☆

> **Para leigo:** essa é a tabela auxiliar que contém a definição completa dos 230 padrões de simetria possíveis. Você faz JOIN com ela usando `data.sgNumber = spacegroups.ITCn` para obter informações adicionais sobre o grupo espacial de cada estrutura — principalmente o sistema cristalino (`class`), que é fundamental para categorizar materiais.

| Coluna | Significado | 🔧 Utilidade no Software |
|---|---|---|
| `id` | ID interno (1–530, inclui configurações alternativas do mesmo número ITA) | Chave primária interna. Não usar para JOIN — use `ITCn` em vez disso. |
| `ITCn` | Número ITA (1–230) — **chave de JOIN com `data.sgNumber`** | `JOIN spacegroups s ON d.sgNumber = s.ITCn`. Use para enriquecer resultados com `class`, `HM`, `Schoenflies`. |
| `Hall` | Notação Hall — representação algébrica dos geradores de simetria. | **Para geração de difratograma simulado.** Software de difração precisa dessa notação para calcular fatores de estrutura. Não exibir ao usuário comum. |
| `Schoenflies` | Símbolo Schoenflies — ex: `D4h`, `C2v`, `Oh`. | **Para correlação com espectroscopia** (Raman/IR). Se seu software integrar análise Raman, esse campo mapeia para modos vibracionais permitidos. |
| `HM` | Hermann-Mauguin padrão — mesmo que `data.sg`. | **Exibição.** Use como label "Grupo Espacial" na UI — é o que cientistas leem. |
| `HMu` | HM com eixo único especificado — ex: `P 1 21/c 1` vs `P 21/c`. | Versão mais explícita do `HM`. Útil quando software de difração pede forma não ambígua. |
| `class` | **Sistema cristalino:** cubic, hexagonal, monoclinic, orthorhombic, rhombohedral, tetragonal, triclinic, trigonal. | **Campo mais útil para filtragem por categoria.** Permite ao usuário filtrar: "mostrar só cúbicos" → `WHERE s.class = 'cubic'`. Use para montar dropdown de filtro na UI. Requer JOIN com `data`. |
| `Nau` | Multiplicidade geral = número de posições equivalentes na célula. | **Para calcular Z':** `Z' = data.Z / spacegroups.Nau`. Z'>1 indica estrutura mais complexa. Uso interno. |

**Query útil — categorizar por sistema cristalino:**
```sql
SELECT d.file, d.formula, s.class, s.HM, d.a, d.b, d.c, d.vol, d.Robs
FROM data d
JOIN spacegroups s ON d.sgNumber = s.ITCn
WHERE s.class = 'cubic'
  AND d.status IS NULL
  AND d.radiation = 'X-ray'
ORDER BY d.Robs ASC NULLS LAST;
```

---

## 3. Tabela `smiles` ★★★☆☆

> **Para leigo:** SMILES é uma forma de escrever a estrutura química de um material como texto. Por exemplo, Fe₂O₃ vira `[Fe+3].[Fe+3].[O-2].[O-2].[O-2]`. Serve para buscar por estrutura molecular quando você sabe "parece com isso" mas não sabe o nome exato.

| Coluna | Significado | 🔧 Utilidade no Software |
|---|---|---|
| `id` | PK interno | Chave para JOIN com `fingerprints`. |
| `cod_id` | FK → `data.file` | Chave de JOIN com tabela `data`: `JOIN smiles s ON s.cod_id = d.file`. |
| `value` | String SMILES — representação textual da estrutura química 2D. Ex: `[Fe+3].[Fe+3].[O-2].[O-2].[O-2]` para Fe₂O₃. | **Busca por estrutura química.** Integre com biblioteca RDKit (Python) ou CDK (Java) para busca subescrutural — "encontrar materiais que contêm esse bloco estrutural". Exibir ao usuário como representação visual (usar biblioteca de renderização SMILES → imagem 2D). |

---

## 4. Tabela `fingerprints` ★★★☆☆

> **Para leigo:** fingerprints são uma versão comprimida da estrutura química como um número binário gigante (2048 bits). Servem para comparar materiais rapidamente: dois materiais com fingerprints parecidos têm estrutura química parecida. É como um hash da forma molecular.

| Coluna | Significado | 🔧 Utilidade no Software |
|---|---|---|
| `smiles_id` | FK → `smiles.id` | Chave de JOIN: `JOIN fingerprints f ON f.smiles_id = s.id`. |
| `fp0`–`fp31` | 32 inteiros de 64 bits = 2048 bits totais representando a estrutura molecular em raio de 4 ligações (Morgan fingerprint). | **Busca por similaridade química.** Calcule o fingerprint do composto de interesse → compare com banco via operações de bits (AND, popcount) → Tanimoto similarity. Útil quando usuário quer "materiais estruturalmente parecidos com X". Implementação requer biblioteca de cheminformatics (RDKit). |

---

## 5. Tabelas de Cross-Referência ★★☆☆☆

> **Para leigo:** essas tabelas ligam cada estrutura do COD com a mesma estrutura em outros bancos de dados. Útil quando o usuário quer ver mais informações sobre o material em fontes externas.

Todas têm estrutura: `(id, ext_id, cod_id, relation_id)`.

| Tabela | Banco externo | Utilidade para o software |
|---|---|---|
| `amcsd_x_cod` | AMCSD — banco de minerais | **Alta relevância para mineralogia.** Link para dados de difração de pó de minerais. URL: `https://rruff.geo.arizona.edu/AMS/result.php?ID={ext_id}` |
| `pubchem_x_cod` | PubChem (NCBI) | Dados extras: CAS number, sinônimos, propriedades físicas. URL: `https://pubchem.ncbi.nlm.nih.gov/compound/{ext_id}` |
| `chemspider_x_cod` | ChemSpider (RSC) | Alternativa ao PubChem. URL: `https://www.chemspider.com/Chemical-Structure.{ext_id}.html` |
| `mpod_x_cod` | MPOD — propriedades de materiais | Link para dados de propriedades elásticas/ópticas do material. |
| `rod_x_cod` | ROD — difração de superfícies | Para análise de filmes finos. |
| `drugbank_x_cod` | DrugBank | Para fármacos com estrutura cristalina determinada. |
| `wikidata_x_cod` | Wikidata | Linked data — contexto enciclopédico. |
| `wikipedia_x_cod` | Wikipedia | Link direto para artigo Wikipedia do material. |
| `xrda_x_cod` | XRDA | Cross-ref adicional. |

**Padrão de uso no software:**
```sql
-- Encontrar link PubChem de um resultado COD
SELECT p.ext_id as pubchem_id
FROM pubchem_x_cod p
WHERE p.cod_id = :cod_file_id;
-- → URL: https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem_id}
```

---

## 6. Tabela `relations` ★★☆☆☆

Define os **tipos de relação** entre bancos externos e COD. Usada como lookup pelo campo `relation_id` nas tabelas de cross-referência.

| `name` | Significado | Utilidade no Software |
|---|---|---|
| `is same as` | Mesma estrutura em ambos os bancos | Exibir links externos com confiança |
| `is part of` | Subconjunto estrutural | Exibir com nota "relacionado" |
| `is related to` | Relacionado mas não idêntico | Exibir como "ver também" |

---

## 7. Tabelas Bibliográficas ★☆☆☆☆

> **Para leigo:** essas tabelas descrevem os periódicos científicos onde os artigos foram publicados. Raramente necessárias para análise de DRX — o nome do journal já está direto em `data.journal`. Úteis apenas se quiser exibir informações completas da publicação.

### `journals` (1 497 registros)

| Coluna | Utilidade no Software |
|---|---|
| `name` | JOIN com `data.journal` para enriquecer exibição |
| `abbrev` | Versão curta para exibição compacta |
| `ISSNprint`, `ISSNonline` | Identificador padrão de periódico — pode usar para validar fonte |
| `doiprefix` | Construir URL de artigo: `https://doi.org/{doiprefix}{volume}/{firstpage}` |
| `url` | Link direto para o periódico |
| `firstyear`, `lastyear` | `lastyear IS NULL` = periódico ainda ativo |

### `publishers`, `jsequences`, `jaltnames`, `successors`
Raramente necessárias. Usar só se construir funcionalidade de busca bibliográfica avançada.

---

## 8. Tabelas Administrativas ☆☆☆☆☆

> Não têm utilidade para análise de DRX. Infraestrutura interna do COD.

| Tabela | Ignorar porque |
|---|---|
| `numbers` | Controle interno de numeração de IDs — sem uso na análise |
| `news` | Avisos do site COD — irrelevante para o software |
| `rdf_relations` | Metadados semânticos para linked data — sem uso prático |

---

## Queries Essenciais para Análise DRX

### Busca principal por fase (combina todos os filtros)
```sql
SELECT 
    d.file,
    d.formula,
    d.mineral,
    d.commonname,
    d.sg,
    d.sgNumber,
    s.class AS crystal_system,
    d.a, d.b, d.c,
    d.alpha, d.beta, d.gamma,
    d.vol,
    d.Z,
    d.method,
    d.radiation,
    d.wavelength,
    d.radSymbol,
    d.celltemp,
    d.Robs,
    d.gofall,
    d.year,
    d.doi
FROM data d
JOIN spacegroups s ON d.sgNumber = s.ITCn
WHERE d.status IS NULL
  AND d.onhold IS NULL
  AND d.radiation = 'X-ray'
  AND d.method != 'theoretical'
  AND d.flags LIKE '%has coordinates%'
  AND d.sgNumber = :sgNumber                            -- grupo espacial do experimento
  AND d.vol BETWEEN :vol * 0.98 AND :vol * 1.02        -- ±2% volume
  AND (d.celltemp IS NULL OR d.celltemp BETWEEN 250 AND 320)  -- temperatura ambiente
  AND (d.cellpressure IS NULL OR d.cellpressure < 200)        -- pressão ambiente
ORDER BY d.Robs ASC NULLS LAST
LIMIT 20;
```

### Busca por nome do mineral
```sql
SELECT file, formula, mineral, commonname, sg, sgNumber, a, b, c, vol, Robs, doi
FROM data
WHERE status IS NULL
  AND onhold IS NULL
  AND (mineral ILIKE '%' || :nome || '%'
    OR commonname ILIKE '%' || :nome || '%'
    OR chemname ILIKE '%' || :nome || '%')
  AND radiation = 'X-ray'
ORDER BY Robs ASC NULLS LAST
LIMIT 20;
```

### Busca por grupo espacial + sistema cristalino
```sql
-- Filtro por sistema cristalino (cúbico = sgNumber 195-230)
SELECT d.file, d.formula, d.a, d.b, d.c, d.vol, d.Robs
FROM data d
JOIN spacegroups s ON d.sgNumber = s.ITCn
WHERE s.class = :crystal_system   -- 'cubic', 'hexagonal', 'monoclinic', etc.
  AND d.status IS NULL
  AND d.nel = :num_elements        -- '1', '2', '3' etc.
  AND d.radiation = 'X-ray'
ORDER BY d.Robs ASC NULLS LAST;
```

### Estruturas de pó para referência Rietveld
```sql
SELECT d.file, d.formula, d.sg, d.a, d.b, d.c, d.vol, d.Z,
       d.wavelength, d.Robs, d.RI, d.gofall, d.year, d.doi
FROM data d
WHERE d.method = 'powder diffraction'
  AND d.status IS NULL
  AND d.onhold IS NULL
  AND d.radiation = 'X-ray'
  AND d.flags LIKE '%has coordinates%'
  AND d.Robs < 0.10
ORDER BY d.Robs ASC NULLS LAST;
```

### Calcular densidade teórica (requer massa molar externa)
```sql
-- Retorna dados para calcular: ρ = (Z × M) / (6.022e23 × vol × 1e-30)
SELECT file, formula, calcformula, vol, Z,
       (Z::float / vol) AS Z_per_A3   -- proxy de densidade relativa
FROM data
WHERE status IS NULL
  AND Z IS NOT NULL
  AND vol IS NOT NULL
  AND formula ILIKE '%Fe%O%';
```

---

## Resumo: Colunas Críticas para DRX por Caso de Uso

| Caso de Uso | Colunas Essenciais | Filtros Obrigatórios |
|---|---|---|
| **Identificação de fase** | `formula`, `sg`, `sgNumber`, `a`, `b`, `c`, `alpha`, `beta`, `gamma`, `vol` | `status IS NULL`, `onhold IS NULL`, `radiation = 'X-ray'` |
| **Simulação de difratograma** | + `Z`, `wavelength`, `method`, `flags`, `sgHall` | + `flags LIKE '%has coordinates%'` |
| **Refinamento Rietveld** | + `Robs`, `wRall`, `gofall`, `RI` | + `method = 'powder diffraction'`, `Robs < 0.10` |
| **Polimorfos / fases por T** | + `celltemp`, `diffrtemp`, `thermalhist`, `cellpressure` | + `cellpressure IS NULL` (só ambiente) |
| **Triagem de qualidade** | `status`, `Robs`, `flags`, `sigvol` | `Robs < 0.08`, `sigvol/vol < 0.005` |
| **Busca por nome** | `mineral`, `commonname`, `chemname`, `formula` | `status IS NULL` |
| **Categorização por sistema** | `sgNumber` → JOIN `spacegroups.class` | — |
| **Links externos** | `doi`, `file` → JOIN `pubchem_x_cod`, `amcsd_x_cod` | — |

---

## Dicionário Rápido de Termos para o Desenvolvedor

| Termo | O que é na prática |
|---|---|
| Célula unitária | A menor unidade repetitiva do cristal — como o "tijolo" que forma o material |
| Parâmetros de célula (a,b,c,α,β,γ) | As medidas do "tijolo" — 3 dimensões + 3 ângulos |
| Grupo espacial (sgNumber) | O padrão de simetria dos átomos dentro do tijolo — 230 padrões possíveis |
| Sistema cristalino | Família do padrão: cúbico, hexagonal, monoclínico, etc. |
| Fator R (Robs) | Nota de qualidade do modelo: 0–5% = ótimo, 5–10% = bom, >15% = ruim |
| Z | Quantas "moléculas" cabem dentro de um tijolo |
| Wavelength (λ) | Cor dos raios-X usados — determina a escala angular do difratograma |
| Método | Como foi medido: monocristal (mais preciso) ou pó (mais comum em lab) |
| DRX / Difratograma | Gráfico de intensidade vs. ângulo — a "impressão digital" do material |
| Fase cristalina | Um material com estrutura cristalina específica — mesmo material pode ter fases diferentes |

---

*Documento gerado com base no schema COD SVN rev. 306581 e nas definições CIF IUCr Core Dictionary v2.4.*  
*Referências: Gražulis et al., Nucleic Acids Research 40, D420–D427 (2012) · IUCr CIF Core Dictionary · Spek, Acta Cryst. D65, 148–155 (2009)*
