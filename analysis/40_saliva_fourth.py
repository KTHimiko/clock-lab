#!/usr/bin/env python3
"""
Stage 40 — a fourth saliva cohort: which side of the sign reversal is it on?

Stage 39: adjusted for age, Horvath 2018 moves +1.03 (SE 0.29) and +1.90 (0.32)
years per 10 points of immune fraction in the two EPIC saliva cohorts (one
group, GSE232891 / GSE232332) and -0.94 (0.14) in GSE78874 (450k, PEG study).
Array, group and population are confounded in that contrast.

GSE149747 (Methylation Diet and Lifestyle; EPIC; another group) has 44 adults
sampled before an intervention, and again twice after. Too small for a transport
test, big enough to estimate the slope.

PRE-REGISTERED (written before GSE149747 was scored)
  1. Horvath 2018 coverage >= 95% and r with age > 0.5 at baseline. Hard stop.
  2. reading, primary analysis = baseline samples ("Before"), slope of Horvath
     2018 on immune fraction adjusted for age, 95% CI:
       CI above 0 and excluding -0.94  -> sides with the EPIC cohorts: the
                                          reversal follows array (or GSE78874)
       CI below 0 and excluding +1.03  -> sides with GSE78874: the reversal
                                          follows the GSE232891/GSE232332 group
       otherwise                        -> inconclusive
  Secondary, reported: person means over all three timepoints; Levine 2018.

Usage:  .venv/bin/python analysis/40_saliva_fourth.py
"""
import sys, gzip
from pathlib import Path
import numpy as np
import pandas as pd
import rdata

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from model.clocks import predict, load_coefficients
from model.deconvolution import deconvolve

DATA = ROOT / "reference/data"; CACHE = ROOT / "results/cache"; TAG = "GSE149747"
v = rdata.read_rda(str(DATA / "epidish/EpiDISH/data/centEpiFibIC.m.rda"))["centEpiFibIC.m"]
REF1 = pd.DataFrame(v.values, index=[str(x) for x in v.coords[v.dims[0]].values],
                    columns=[str(x) for x in v.coords[v.dims[1]].values])
TWO = ["Levine2018", "Horvath2018"]
needed = set(REF1.index)
for c in TWO:
    needed |= set(load_coefficients(c)[0].index)

meta = {}
with gzip.open(DATA / f"{TAG}_series_matrix.txt.gz", "rt") as fh:
    for line in fh:
        if line.startswith("!series_matrix_table_begin"):
            break
        p = [x.strip('"') for x in line.rstrip("\n").split("\t")]
        if p[0] == "!Sample_title":
            meta["id"] = p[1:]
        elif p[0] == "!Sample_characteristics_ch1":
            k = p[1].partition(":")[0].strip()
            if k in ("age", "timepoint (intervention)"):
                meta[k] = [x.partition(":")[2].strip() for x in p[1:]]
m = pd.DataFrame(meta).set_index("id")
m["age"] = pd.to_numeric(m["age"].str.replace(" years", ""), errors="coerce")
m["person"] = m.index.str.extract(r"MDL(\d+)", expand=False)
m["tp"] = m["timepoint (intervention)"]

parts = []
for ch in pd.read_csv(DATA / f"{TAG}_MDL_Matrix_AverageBetas.csv.gz", index_col=0, chunksize=50_000, low_memory=False,
                      encoding="utf-8-sig"):
    parts.append(ch.loc[ch.index.intersection(needed)])
raw = pd.concat(parts)
bet = raw[[c for c in raw.columns if not c.startswith("Detection")]].apply(pd.to_numeric, errors="coerce")
bet = bet[[c for c in bet.columns if c in m.index]]
print(f"  {bet.shape[1]} amostras, {m.loc[bet.columns, 'person'].nunique()} pessoas; tempos {m.loc[bet.columns, 'tp'].value_counts().to_dict()}")

p1 = [p for p in REF1.index if p in bet.index]
f = deconvolve(REF1.loc[p1].to_numpy(), bet.loc[p1].to_numpy()); f = f / f.sum(1, keepdims=True)
d = pd.DataFrame(f, index=bet.columns, columns=REF1.columns).join(m)
for c in TWO:
    d[c], _ = predict(bet, c, min_coverage=0.0)
print(f"  IC mediano {d.IC.median():.2f}; sondas do painel {len(p1)}/{len(REF1)}")


def slope(df, c):
    X = np.column_stack([np.ones(len(df)), df.age, df.IC]); y = df[c].to_numpy(float)
    b, *_ = np.linalg.lstsq(X, y, rcond=None); e = y - X @ b
    se = np.sqrt(e @ e / (len(y) - 3) * np.linalg.inv(X.T @ X)[2, 2])
    return b[2] / 10, se / 10


base = d[d.tp == "Before"].dropna(subset=["age"])
cov = load_coefficients("Horvath2018")[0].index.isin(bet.index).mean()
r = float(np.corrcoef(base.age, base.Horvath2018)[0, 1])
c1 = cov >= 0.95 and r > 0.5
print(f"CHECAGEM 1 — Horvath2018 cobertura {cov:.1%}, r com idade na linha de base {r:.3f} (n={len(base)}) -> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("parando.")
b, se = slope(base, "Horvath2018"); lo, hi = b - 1.96 * se, b + 1.96 * se
if lo > 0 and lo > -0.94:
    verdict = "lado das coortes EPIC: a inversao acompanha o array (ou o GSE78874)"
elif hi < 0 and hi < 1.03:
    verdict = "lado do GSE78874: a inversao acompanha o grupo do GSE232891/GSE232332"
else:
    verdict = "inconclusivo"
print(f"CHECAGEM 2 — Horvath2018, linha de base: {b:+.2f} anos por +10pp imune, IC95 [{lo:+.2f}, {hi:+.2f}] -> {verdict}")
pm = d.dropna(subset=["age"]).groupby("person")[["age", "IC"] + TWO].mean()
for c in TWO:
    b1, s1 = slope(base, c); b2, s2 = slope(pm, c)
    print(f"  {c:<12} linha de base {b1:+.2f} (EP {s1:.2f}), media por pessoa {b2:+.2f} (EP {s2:.2f}, n={len(pm)})")
d.to_csv(CACHE / f"{TAG}_saliva.csv")
