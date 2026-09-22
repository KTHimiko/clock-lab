#!/usr/bin/env python3
"""
Stage 2 — what do the clocks say about purified cells from one donor?

Reinius et al. drew blood from six men and split each sample ten ways: whole
blood, PBMC, granulocytes, CD4+ T, CD8+ T, CD14+ monocytes, CD19+ B, CD56+ NK,
neutrophils, eosinophils.

Same person, same day, ten measurements. **Any difference in epigenetic age
between those fractions cannot be ageing.** Its size is a direct measure of how
much a clock reads composition rather than time.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. probe coverage near 100% for every clock, as in validation
  2. exactly 6 donors x 10 fractions parsed from the titles
  3. PBMC must sit closer to whole blood than a random pair does — PBMC IS a
     subset of whole blood, so if the donor parsing or the clock is wrong, that
     relation breaks. This is the check that catches a silent mix-up.

Usage:  .venv/bin/python analysis/02_purified_cells.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix
from model.clocks import CLOCKS, predict

MATRIX = ROOT / "reference/data/GSE35069_series_matrix.txt.gz"
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


print("carregando GSE35069 ...", flush=True)
betas, meta = read_series_matrix(MATRIX)
meta = meta.set_index("gsm")
frac = meta["tissue/cell type"].reindex(betas.columns)
donor = meta["title"].reindex(betas.columns).str.extract(r"_(\d+)$")[0]
print(f"  sondas {betas.shape[0]:,}  amostras {betas.shape[1]}", flush=True)


section("CHECAGENS DE SANIDADE")
n_don, n_frac = donor.nunique(), frac.nunique()
design_ok = (n_don == 6) and (n_frac == 10) and (betas.shape[1] == 60)
print(f"  doadores {n_don}, fracoes {n_frac}, amostras {betas.shape[1]}  "
      f"-> {'ok' if design_ok else 'DESENHO INESPERADO'}")

ages = {}
for name in CLOCKS:
    pred, cov = predict(betas, name)
    ages[name] = pred
    print(f"  {name:<14} cobertura {cov:>6.1%}  "
          f"{'ok' if cov > 0.95 else 'COBERTURA BAIXA'}")
A = pd.DataFrame(ages)
A["donor"], A["fraction"] = donor.values, frac.values
A.to_csv(OUT / "purified_ages.csv")

# check 3: PBMC against whole blood, within donor, versus a random pair
ref = "Horvath2013"
wide = A.pivot_table(index="donor", columns="fraction", values=ref)
d_pbmc_wb = float((wide["PBMC"] - wide["Whole blood"]).abs().mean())
rng = np.random.default_rng(0)
cols = [c for c in wide.columns]
rand = [abs(wide.loc[d, rng.choice(cols)] - wide.loc[d, rng.choice(cols)])
        for d in wide.index for _ in range(200)]
d_rand = float(np.mean(rand))
check3 = d_pbmc_wb < d_rand
print(f"  |PBMC - sangue total| = {d_pbmc_wb:.2f} anos   "
      f"par aleatorio = {d_rand:.2f} anos  -> {'ok' if check3 else 'FALHOU'}")

if not (design_ok and check3):
    print("\n  UMA CHECAGEM FALHOU — parando antes de ler o resultado.")
    sys.exit(1)


section("IDADE EPIGENETICA POR FRACAO CELULAR")
order = (A.groupby("fraction")[ref].mean().sort_values().index)
print(f"  {'fracao':<18}" + "".join(f"{c.replace('2013','').replace('2018',''):>12}" for c in CLOCKS))
for f in order:
    sub = A[A.fraction == f]
    print(f"  {f:<18}" + "".join(f"{sub[c].mean():>12.1f}" for c in CLOCKS))


section("QUANTO DISSO NAO PODE SER ENVELHECIMENTO")
rows = []
for name in CLOCKS:
    w = A.pivot_table(index="donor", columns="fraction", values=name)
    within = float(w.std(axis=1).mean())           # spread across fractions, per donor
    between = float(w.mean(axis=1).std())          # spread across donors
    span = float((w.max(axis=1) - w.min(axis=1)).mean())
    rows.append(dict(clock=name, within=within, between=between, span=span,
                     ratio=within / between if between else np.nan))
    print(f"  {name:<14} entre fracoes do mesmo doador: desvio {within:>5.1f}a, "
          f"amplitude {span:>5.1f}a | entre doadores: {between:>4.1f}a")
res = pd.DataFrame(rows)
res.to_csv(OUT / "purified_variance.csv", index=False)

print(f"""
  O desvio ENTRE FRACOES e medido dentro de uma pessoa, no mesmo dia: e
  impossivel que seja envelhecimento. O desvio ENTRE DOADORES mistura idade real
  com tudo mais. Quando o primeiro e maior que o segundo, a maior parte do que o
  relogio ve numa amostra de sangue nao e a idade de quem doou.

  razao entre-fracoes / entre-doadores:
""" + "".join(f"    {r.clock:<14} {r.ratio:>5.2f}\n" for _, r in res.iterrows()))
