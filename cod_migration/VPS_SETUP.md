# Setup xrd_loader na Oracle VPS

**Target:** Oracle VPS ARM (4 OCPU, 24 GB RAM, Ubuntu)  
**Goal:** rodar `xrd_loader.py` com 8+ workers, sem swap, processando MEDIUM/LARGE tiers

---

## Pré-requisitos locais (Windows)

```powershell
# Dump das tabelas necessárias (sem reference_patterns_pymatgen — recalculada na VPS)
pg_dump -Fc `
  -t public.data `
  -t public.cod_files `
  -t "xrd_analysis.reference_patterns" `
  -h localhost -U cod_admin -d cod `
  > cod_seed.pgc

# Tamanho esperado: ~5-8 GB
```

---

## 1. VPS — Sistema base

```bash
sudo apt update && sudo apt upgrade -y
sudo apt install -y python3.11 python3.11-venv python3-pip \
                    postgresql-15 subversion git \
                    build-essential libpq-dev
```

---

## 2. VPS — PostgreSQL

```bash
# Inicia serviço
sudo systemctl enable --now postgresql

# Cria usuário e banco
sudo -u postgres psql <<'SQL'
CREATE USER cod_admin WITH PASSWORD 'sua_senha_aqui';
CREATE DATABASE cod OWNER cod_admin;
\c cod
CREATE SCHEMA xrd_analysis AUTHORIZATION cod_admin;
SQL
```

---

## 3. VPS — Restaura dump do banco local

```bash
# Transfere dump (do Windows para VPS)
scp cod_seed.pgc user@oracle-vps:/home/user/

# Na VPS: restaura
pg_restore -h localhost -U cod_admin -d cod /home/user/cod_seed.pgc

# Verifica
psql -U cod_admin -d cod -c "SELECT count(*) FROM public.data;"
psql -U cod_admin -d cod -c "SELECT count(*) FROM public.cod_files;"
```

---

## 4. VPS — Clona repositório

```bash
git clone <URL_DO_REPO> ~/projeto
cd ~/projeto
```

---

## 5. VPS — Ambiente Python

```bash
cd ~/projeto/cod_migration
python3.11 -m venv .venv
source .venv/bin/activate
pip install pymatgen asyncpg python-dotenv psutil
```

---

## 6. VPS — .env

```bash
cat > ~/projeto/cod_migration/.env <<'EOF'
PG_HOST=localhost
PG_PORT=5432
PG_DB=cod
PG_USER=cod_admin
PG_PASSWORD=sua_senha_aqui
EOF
chmod 600 ~/projeto/cod_migration/.env
```

---

## 7. VPS — Migration SQL

```bash
psql -h localhost -U cod_admin -d cod -f ~/projeto/cod_migration/xrd_schema_alter.sql
```

---

## 8. VPS — Checkout SVN dos CIFs

```bash
mkdir -p /data/cod_svn
# Roda em background — demora horas (110 GB)
nohup svn checkout svn://www.crystallography.net/cod/trunk/cif /data/cod_svn/cif \
  > /tmp/svn_checkout.log 2>&1 &

echo "PID: $!"
# Monitora: tail -f /tmp/svn_checkout.log
```

> **Aguarda SVN completar antes de continuar.**  
> Verifica: `ls /data/cod_svn/cif/ | head`

---

## 9. VPS — Ajusta CIF_ROOT no worker

O worker usa `CIF_ROOT = PROJECT_ROOT / "cod_svn" / "cif"` por default.  
Se o SVN foi clonado em `/data/cod_svn/cif`, ajusta em `xrd_worker.py`:

```python
# xrd_worker.py linha ~21
CIF_ROOT = Path("/data/cod_svn/cif")
```

Ou cria symlink:
```bash
ln -s /data/cod_svn ~/projeto/cod_svn
```

---

## 10. VPS — Tuning PostgreSQL (bulk load)

```bash
sudo -u postgres psql -d cod <<'SQL'
ALTER SYSTEM SET synchronous_commit = off;
ALTER SYSTEM SET work_mem = '256MB';
ALTER SYSTEM SET maintenance_work_mem = '1GB';
ALTER SYSTEM SET checkpoint_completion_target = 0.9;
SELECT pg_reload_conf();
SQL
```

---

## 11. VPS — Roda o loader

```bash
cd ~/projeto/cod_migration
source .venv/bin/activate

# Roda em background com nohup (sessão SSH pode cair)
nohup python xrd_loader.py \
  --tiers SMALL MEDIUM LARGE \
  --n-small 8 \
  --n-medium 6 \
  --n-large 3 \
  --timeout 300 \
  > /tmp/xrd_loader.log 2>&1 &

echo "PID: $!"
tail -f /tmp/xrd_loader.log
```

---

## 12. Pós-carga — Restaura PostgreSQL e cria índices

```bash
sudo -u postgres psql -d cod <<'SQL'
ALTER SYSTEM RESET synchronous_commit;
ALTER SYSTEM RESET work_mem;
SELECT pg_reload_conf();
VACUUM ANALYZE xrd_analysis.reference_patterns_pymatgen;
-- GIN index para busca em reflections (lento — roda após carga completa)
CREATE INDEX rp_pmg_reflections_gin
    ON xrd_analysis.reference_patterns_pymatgen USING GIN (reflections)
    WITH (fastupdate = off);
SQL
```

---

## 13. Extrai resultados para o banco local

```bash
# Na VPS: dump só da tabela processada
pg_dump -Fc \
  -t "xrd_analysis.reference_patterns_pymatgen" \
  -t "xrd_analysis.failed_patterns_pymatgen" \
  -h localhost -U cod_admin -d cod \
  > xrd_results.pgc

# No Windows: importa
scp user@oracle-vps:/home/user/xrd_results.pgc .
pg_restore -h localhost -U cod_admin -d cod xrd_results.pgc
```

---

## Estimativas de tempo

| Etapa | Tempo estimado |
|-------|----------------|
| SVN checkout (110 GB) | 4-12h (depende da banda) |
| Dump local → VPS (5-8 GB) | 30-90 min |
| SMALL tier (10k CIFs, 8 workers) | ~15 min |
| MEDIUM tier (427k CIFs, 6 workers) | ~15-20h |
| LARGE tier (31k CIFs, 3 workers) | ~5-15h |
| **Total** | **~24-48h de cálculo** |

---

## Monitoramento

```bash
# Progresso em tempo real
tail -f /tmp/xrd_loader.log

# RAM
free -h

# Workers Python ativos
ps aux | grep python | grep -v grep | wc -l

# Contagem no banco
psql -U cod_admin -d cod -c \
  "SELECT count(*) FROM xrd_analysis.reference_patterns_pymatgen;"
```
