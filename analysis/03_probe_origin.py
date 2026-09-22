#!/usr/bin/env python3
"""
Stage 3 — where the cell-type signal lives, and whether it predicts vulnerability.

Stage 2 showed the clocks disagree by up to 35 years across fractions of one
blood draw, and that Hannum and Levine are roughly twice as exposed as the two
Horvath clocks. This asks why: how much of each clock's own probe set is
measuring cell identity rather than time.

The design is paired — every donor contributes all ten fractions — so for each
probe the variance splits three ways: between donors, between cell types, and
residual. The middle term is the one that cannot be ageing.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. the variance fractions must lie in [0, 1] and sum to 1 per probe
  2. the distribution across all 485k probes must be right-skewed — most probes
     are not cell-type markers; if the median probe looked cell-type specific,
     the decomposition would be wrong rather than the biology surprising
  3. CD4+ and CD8+ T cells must be closer to each other, in overall methylation
     distance, than either is to neutrophils. That is settled lineage biology:
     if it fails, the fraction labels are scrambled and nothing else counts.

Usage:  .venv/bin/python analysis/03_probe_origin.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix
from model.clocks import CLOCKS, load_coefficients

MATRIX = ROOT / "reference/data/GSE35069_series_matrix.txt.gz"
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


print("carregando GSE35069 ...", flush=True)
betas, meta = read_series_matrix(MATRIX)
meta = meta.set_index("gsm")
frac = meta["tissue/cell type"].reindex(betas.columns)
donor = meta["title"].reindex(betas.columns).str.extract(r"_(\d+)$")[0]
# COMPLETE CASES ONLY, and the reason matters. The decomposition below relies
# on SS_total = SS_donor + SS_fraction + residual, an identity that holds for a
# balanced design and breaks the moment a cell is missing: group means then come
# from different subsets and the parts stop summing to the whole. The first run
# of this stage skipped that and the sanity check caught it — the components
# exceeded 1. Dropping the incomplete probes restores the balance.
complete = betas.notna().all(axis=1)
dropped = int((~complete).sum())
betas = betas.loc[complete]
X = betas.to_numpy(dtype=np.float32)          # probes x samples
print(f"  {X.shape[0]:,} sondas x {X.shape[1]} amostras "
      f"({dropped:,} sondas com valor ausente, descartadas)", flush=True)


section("CHECAGEM 3 — A BIOLOGIA DE LINHAGEM SE SUSTENTA?")
# mean profile per fraction, then correlation distance between fractions
prof = pd.DataFrame(X, columns=pd.MultiIndex.from_arrays([donor, frac])).T
prof = prof.groupby(level=1).mean()
valid = ~prof.isna().any(axis=0)
corr = np.corrcoef(prof.loc[:, valid].to_numpy())
names = list(prof.index)
def d(a, b):
    return 1 - corr[names.index(a), names.index(b)]
d_t = d("CD4+ T cells", "CD8+ T cells")
d_4n = d("CD4+ T cells", "Neutrophils")
d_8n = d("CD8+ T cells", "Neutrophils")
lineage_ok = d_t < d_4n and d_t < d_8n
print(f"  distancia CD4-CD8        {d_t:.4f}")
print(f"  distancia CD4-neutrofilo {d_4n:.4f}")
print(f"  distancia CD8-neutrofilo {d_8n:.4f}")
print(f"  -> {'ok, os rotulos fazem sentido biologico' if lineage_ok else 'FALHOU: rotulos embaralhados'}")


section("DECOMPOSICAO DE VARIANCIA POR SONDA")
# two-way layout without interaction: donor and fraction, balanced 6 x 10
dn = pd.Categorical(donor).codes
fr = pd.Categorical(frac).codes
nd, nf = dn.max() + 1, fr.max() + 1
grand = X.mean(axis=1, keepdims=True)

def group_means(codes, k):
    M = np.zeros((X.shape[0], k), dtype=np.float32)
    for g in range(k):
        M[:, g] = X[:, codes == g].mean(axis=1)
    return M

Md, Mf = group_means(dn, nd), group_means(fr, nf)
ss_d = nf * ((Md - grand) ** 2).sum(axis=1)
ss_f = nd * ((Mf - grand) ** 2).sum(axis=1)
ss_t = ((X - grand) ** 2).sum(axis=1)
with np.errstate(invalid="ignore", divide="ignore"):
    v_frac = np.where(ss_t > 0, ss_f / ss_t, np.nan)
    v_don = np.where(ss_t > 0, ss_d / ss_t, np.nan)
v = pd.DataFrame({"probe": betas.index, "var_fraction": v_frac, "var_donor": v_don})
v = v.dropna()

in_range = bool(((v.var_fraction >= -1e-6) & (v.var_fraction <= 1 + 1e-6)).all()
                and ((v.var_fraction + v.var_donor) <= 1 + 1e-3).all())
med, p90 = v.var_fraction.median(), v.var_fraction.quantile(0.90)
skewed = med < p90 / 2
print(f"  fracoes em [0,1] e somando <= 1: {'ok' if in_range else 'FALHOU'}")
print(f"  variancia explicada por tipo celular: mediana {med:.3f}, "
      f"p90 {p90:.3f}, p99 {v.var_fraction.quantile(0.99):.3f}")
print(f"  distribuicao assimetrica a direita: {'ok' if skewed else 'FALHOU'}")

if not (lineage_ok and in_range and skewed):
    print("\n  UMA CHECAGEM FALHOU — parando antes de ler o resultado.")
    sys.exit(1)
v.to_csv(OUT / "probe_variance.csv.gz", index=False, compression="gzip")


section("ONDE AS SONDAS DE CADA RELOGIO CAEM NESSA DISTRIBUICAO")
vi = v.set_index("probe")
print(f"  {'relogio':<14}{'sondas':>8}{'mediana':>10}{'percentil':>12}"
      f"{'|peso| x var':>14}")
rows = []
for name in CLOCKS:
    coefs, _, _ = load_coefficients(name)
    present = vi.index.intersection(coefs.index)
    sub = vi.loc[present, "var_fraction"]
    pct = float((v.var_fraction < sub.median()).mean())
    w = coefs.loc[present].abs()
    weighted = float((w * sub).sum() / w.sum())
    rows.append(dict(clock=name, n=len(present), median=float(sub.median()),
                     percentile=pct, weighted=weighted))
    print(f"  {name:<14}{len(present):>8}{sub.median():>10.3f}{pct:>11.0%}"
          f"{weighted:>14.3f}")
res = pd.DataFrame(rows)
res.to_csv(OUT / "clock_probe_variance.csv", index=False)


section("ISSO PREVE A VULNERABILIDADE MEDIDA NA ETAPA 2?")
spread = pd.read_csv(OUT / "purified_variance.csv").set_index("clock")["within"]
res = res.set_index("clock")
joined = res.join(spread.rename("spread_years"))
print(f"  {'relogio':<14}{'|peso| x var':>14}{'desvio etapa 2':>17}")
for c, r in joined.iterrows():
    print(f"  {c:<14}{r.weighted:>14.3f}{r.spread_years:>16.1f}a")
rho = float(joined["weighted"].corr(joined["spread_years"], method="spearman"))
print(f"""
  correlacao de Spearman entre as duas colunas: {rho:+.2f}  (n = 4)

  Com quatro relogios, uma correlacao nao e evidencia de nada sozinha — quatro
  pontos produzem +1.00 por acaso com probabilidade de 1 em 24. Vale como
  consistencia, nao como teste, e esta registrada assim.""")
