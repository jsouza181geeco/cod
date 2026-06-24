# COD → PostgreSQL Migration Plan

## Visão Geral

Baixar Crystallography Open Database (COD) via SVN e migrar para PostgreSQL local, com suporte a atualizações incrementais.

**Fonte:** `svn://www.crystallography.net/cod`  
**Alvo:** PostgreSQL (local)  
**Linguagem:** Python 3

---

## Estrutura do Repositório SVN

```
svn://www.crystallography.net/cod/
├── cif/          # Arquivos .cif (~5.2 GB) — dados cristalográficos
├── hkl/          # Dados de reflexão (difração)
├── mysql/        # Dumps MySQL: data.sql (schema) + data.txt (dados tabulados)
├── retracted/    # Entradas retratadas
├── smi/          # SMILES strings
├── statistics/   # Estatísticas do banco
├── docs/
├── html/
├── trunk/
├── branches/
├── tags/
└── vendor/
```

**Estratégia de download:** Checkout somente dos diretórios necessários (sparse checkout), não o repositório inteiro.

---

## Dependências

### Sistema
- SVN (já instalado)
- PostgreSQL 14+ (servidor local)
- Python 3.9+

### Python
```
psycopg2-binary   # driver PostgreSQL
tqdm              # progress bars
python-dotenv     # config via .env
```

### Instalação
```bash
pip install psycopg2-binary tqdm python-dotenv
```

---

## Fases de Implementação

### Fase 1 — Download Inicial via SVN

**Opção A (recomendada): Sparse checkout apenas dos diretórios necessários**
```bash
svn checkout --depth empty svn://www.crystallography.net/cod cod_svn
cd cod_svn
svn update --set-depth infinity mysql
svn update --set-depth infinity cif
```

**Opção B: rsync apenas CIFs (sem histórico SVN)**
```bash
rsync -av --delete rsync://www.crystallography.net/cif/ ./cif/
```

SVN é preferível: suporta `svn update` incremental por revisão.

**Conteúdo do `mysql/`:**
- `data.sql` — schema MySQL (CREATE TABLE, índices)
- `data.txt` — dados tab-separated, UTF-8

---

### Fase 2 — Conversão de Schema MySQL → PostgreSQL

O `data.sql` usa sintaxe MySQL/MariaDB. Necessita conversão:

| MySQL | PostgreSQL |
|-------|-----------|
| `ENGINE=InnoDB` | remover |
| `AUTO_INCREMENT` | `GENERATED ALWAYS AS IDENTITY` |
| `TINYINT(1)` | `BOOLEAN` |
| `CHARACTER SET utf8` | remover (PostgreSQL usa UTF-8 por padrão) |
| `LOAD DATA LOCAL INFILE` | `COPY FROM STDIN` |
| Comentários `/* mariadb-X.Y */` | remover (quebram parsers) |
| backticks `` `col` `` | aspas duplas `"col"` ou remover |

Script Python fará essa conversão automaticamente via regex antes de executar no PostgreSQL.

---

### Fase 3 — Importação Inicial

**Fluxo do script `cod_initial_load.py`:**

```
1. SVN sparse checkout → mysql/ + cif/
2. Ler data.sql → converter schema MySQL→PG → CREATE TABLE no PG
3. Ler data.txt (tab-separated) → COPY INTO PostgreSQL via psycopg2
4. Criar índices após import (mais rápido que durante)
5. Registrar revisão SVN atual em tabela de controle
```

**Tabela de controle de revisão:**
```sql
CREATE TABLE cod_sync_state (
    id SERIAL PRIMARY KEY,
    svn_revision BIGINT NOT NULL,
    synced_at TIMESTAMP DEFAULT NOW(),
    files_changed INTEGER
);
```

---

### Fase 4 — Atualização Incremental

**Fluxo do script `cod_incremental_update.py`:**

```
1. Ler última revisão salva em cod_sync_state
2. svn update mysql/ cif/
3. svn diff --summarize -r <last_rev>:HEAD → lista de arquivos alterados
4. Para cada arquivo alterado:
   - Se data.sql mudou → reaplicar schema (ALTER TABLE se possível)
   - Se data.txt mudou → re-importar registros alterados
   - Se .cif mudou → atualizar registro correspondente
5. Salvar nova revisão em cod_sync_state
```

**Detecção de mudanças via SVN:**
```bash
svn update --revision HEAD mysql/
svn log --revision <last>:HEAD --verbose
```

---

## Estrutura de Arquivos do Projeto

```
cod_migration/
├── .env                      # credenciais PG (não commitar)
├── requirements.txt
├── cod_initial_load.py       # script carga inicial
├── cod_incremental_update.py # script atualização incremental
├── lib/
│   ├── svn_utils.py          # wrappers SVN (checkout, update, diff)
│   ├── schema_converter.py   # MySQL → PostgreSQL DDL
│   └── pg_loader.py          # COPY, upsert helpers
└── cod_svn/                  # checkout SVN (gitignore)
```

---

## Arquivo `.env`

```env
PG_HOST=localhost
PG_PORT=5432
PG_DB=cod
PG_USER=cod_admin
PG_PASSWORD=senha_aqui
SVN_URL=svn://www.crystallography.net/cod
SVN_LOCAL_PATH=./cod_svn
```

---

## Tabela Principal `data`

Schema MySQL original tem uma tabela `data` com campos cristalográficos. Campos esperados (baseado na documentação):

- `file` — número COD (PK)
- `a`, `b`, `c` — parâmetros de célula unitária (Å)
- `alpha`, `beta`, `gamma` — ângulos (graus)
- `Z` — número de fórmulas por célula
- `vol` — volume da célula
- `sg` — grupo espacial
- `sgHall` — notação Hall
- `formula` — fórmula química
- `calcformula`
- `authors`
- `title`
- `journal`, `year`, `volume`, `issue`, `firstpage`, `lastpage`
- `doi`
- `method` — técnica experimental
- `temperature`
- `status` — ex: `retracted`
- ... (verificar schema completo em `mysql/data.sql`)

---

## Pontos de Atenção

1. **Tamanho:** ~500k+ entradas CIF, data.txt vários GB. Import pode levar horas.
2. **Encoding:** COD usa UTF-8 mas alguns CIFs têm caracteres especiais. Usar `errors='replace'` no Python.
3. **Chave primária:** Campo `file` (integer COD ID) — usar como PK no PG.
4. **Upsert incremental:** `INSERT ... ON CONFLICT (file) DO UPDATE SET ...` para atualizações idempotentes.
5. **Índices:** Criar APÓS bulk import (`DISABLE TRIGGER ALL` não necessário no PG; basta criar índices depois).
6. **SVN no Windows:** Usar `subprocess` com `LC_TIME=en_US.UTF-8` pode não funcionar — usar `encoding='utf-8'` no `subprocess.run`.
7. **Retracted:** `retracted/` contém entradas removidas — decidir se inclui com flag ou exclui.

---

## Ordem de Implementação

- [ ] **1.** Criar banco PG e usuário: `createdb cod`
- [ ] **2.** Fazer sparse SVN checkout de `mysql/`
- [ ] **3.** Implementar `schema_converter.py` — testar com `data.sql`
- [ ] **4.** Implementar `pg_loader.py` — COPY from data.txt
- [ ] **5.** Implementar `cod_initial_load.py` — integrar tudo
- [ ] **6.** Rodar carga inicial e validar contagem de registros
- [ ] **7.** Adicionar sparse checkout de `cif/` se necessário
- [ ] **8.** Implementar `cod_incremental_update.py`
- [ ] **9.** Testar ciclo completo: load → simular mudança → update incremental

---

## Referências

- [How to obtain COD](https://wiki.crystallography.net/howtoobtaincod/)
- [Creating SQL Database](https://wiki.crystallography.net/creatingSQLdatabase/)
- SVN repo: `svn://www.crystallography.net/cod`
- rsync: `rsync://www.crystallography.net/cif/`
