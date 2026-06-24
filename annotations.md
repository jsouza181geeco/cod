# Tabela "data"
## Colunas para cálculo do DRX: 
- a
- b
- c
- alpha
- beta
- gamma

## Colunas úteis para comparação:
- vol: (Útil para triagem rápida) 
    - Dado difratograma indexado com parâmetros a,b,c,α,β,γ experimentais → calcular vol → buscar no COD com tolerância WHERE vol BETWEEN x*0.99 AND x*1.01 (~1% tolerância). Mais robusto que buscar por a + b + c separados pois variação de configuração (setting) do grupo espacial muda a/b/c individualmente mas preserva vol.
    - Detecção de polimorfos: mesmo composto, mesmo sg → vol diferente = polimorfo ou expansão térmica significativa.
    - Densidade calculada: ρ = (Z × M) / (Nₐ × vol × 1e-30) — confirma se Z e formula fazem sentido fisicamente.

## Possíveis colunas de filtro:
- status: retirando warnings, errors e retracted
- flags: 'has coordinates' = estrutura atômica completa; 'has disorder' = desordem estrutural; 'has Fobs' = fatores de estrutura observados disponíveis.
- siga, sigb, sigc, sigalpha, sigbeta, siggamma, sigvol: σ = NULL → autores não reportaram incerteza. σ anormalmente grande → refinamento não convergiu.
- sg: Grupo espacial determina quais reflexões (hkl) aparecem ou somem no difratograma — as "ausências sistemáticas"
- sgHall
- sgNumber: Definir canditatos a partir do sgNumber do DRX da amostra