#!/usr/bin/env python3
"""
Stage 47a — the technical replicates of GSE55763, which stage 35 threw away.

Sehgal et al. (2026) show that adjusting for immune cell fractions *within* a
dataset lowers the biological reliability of nearly every clock. They never
transport coefficients. GSE55763 measured 36 individuals twice for a technical
replication study; stage 35 dropped all 72 arrays so nobody would be counted
twice in the fitting cohort. Those arrays answer the question Sehgal et al. raise
for the transported case: does carrying a correction from another cohort make a
clock less repeatable on the same person's blood, measured twice?

This caches composition, clock ages and DunedinPACE for the replicate arrays,
using the same panels as the rest of the project and the published twelve-type
library alongside. 47b measures reliability.

PRE-REGISTERED
  1. exactly 36 individuals with two arrays each, none of them among the 2,639
     population samples stage 35 fitted on. Hard stop.
  2. every clock tracks age across the 72 arrays at r > 0.5.

Usage:  .venv/bin/python analysis/47a_replicates_cache.py
"""
import sys, gzip
from pathlib import Path
import numpy as np
import pandas as pd
import rdata

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_extended import read_extended, TYPES12
from model.clocks import CLOCKS, predict, load_coefficients
from model.deconvolution import deconvolve, TYPES as TYPES6
from model import dunedinpace

DATA = ROOT / "reference/data"; CACHE = ROOT / "results/cache"; TAG = "GSE55763"

ref12 = pd.read_csv(CACHE / "panel12.csv", index_col=0)[TYPES12]
ref6 = pd.read_csv(CACHE / "panel6.csv", index_col=0)[TYPES6]
v = rdata.read_rda(str(DATA / "epidish/EpiDISH/data/cent12CT450k.m.rda"))["cent12CT450k.m"]
MAP = {"CD4Tnv": "CD4nv", "Baso": "Bas", "CD4Tmem": "CD4mem", "Bmem": "Bmem", "Bnv": "Bnv",
       "Treg": "Treg", "CD8Tmem": "CD8mem", "CD8Tnv": "CD8nv", "Eos": "Eos", "NK": "NK",
       "Neu": "Neu", "Mono": "Mono"}
pub12 = pd.DataFrame(v.values, index=[str(x) for x in v.coords[v.dims[0]].values],
                     columns=[str(x) for x in v.coords[v.dims[1]].values]).rename(columns=MAP)[TYPES12]

meta = {}
with gzip.open(DATA / f"{TAG}_series_matrix.txt.gz", "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("!series_matrix_table_begin"):
            break
        p = [x.strip('"') for x in line.rstrip("\n").split("\t")]
        if p[0] == "!Sample_title":
            meta["barcode"] = [x.split(",")[-1].strip() for x in p[1:]]
        elif p[0] == "!Sample_description":
            meta.setdefault("desc", p[1:])
        elif p[0] == "!Sample_characteristics_ch1" and p[1].partition(":")[0].strip() == "age":
            meta["age"] = [x.partition(":")[2].strip() for x in p[1:]]
m = pd.DataFrame(meta).set_index("barcode")
m["age"] = pd.to_numeric(m["age"], errors="coerce")
rep = m[m.desc.str.startswith("Technical replicate")].copy()
rep["person"] = rep.desc.str.extract(r"sample (\d+)\.")[0]
rep["group"] = rep.desc.str.extract(r"group (\d+)")[0]
pop = set(m.index[m.desc == "Population study sample."])
pairs = rep.groupby("person").size()
c1 = len(pairs) == 36 and (pairs == 2).all() and not (set(rep.index) & pop)
print(f"CHECAGEM 1 — {len(pairs)} pessoas, {rep.shape[0]} arrays, sobreposicao com as 2.639: "
      f"{len(set(rep.index) & pop)} -> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("parando.")

needed = set(ref12.index) | set(ref6.index) | set(pub12.index) | set(dunedinpace.M["gold_probes"])
for c in CLOCKS:
    needed |= set(load_coefficients(c)[0].index)
src = DATA / f"{TAG}_normalized_betas.txt.gz"
hdr = pd.read_csv(src, sep="\t", nrows=0).columns.tolist()
cols = [0] + [i for i, h in enumerate(hdr) if i > 0 and not h.startswith("Detection") and h in rep.index]
print(f"lendo {len(cols) - 1} colunas de replicas ...", flush=True)
parts = []
for ch in pd.read_csv(src, sep="\t", usecols=cols, index_col=0, chunksize=50_000):
    parts.append(ch.loc[ch.index.intersection(needed)].astype(np.float32))
b = pd.concat(parts)
rep = rep.loc[b.columns]
print(f"  {b.shape[0]:,} sondas x {b.shape[1]} arrays")

for name, R, out in (("caseiro", ref12, "comp12"), ("publicado", pub12, "comp12pub")):
    p = [x for x in R.index if x in b.index]
    pd.DataFrame(deconvolve(R.loc[p].to_numpy(), b.loc[p].to_numpy().astype(float)),
                 index=b.columns, columns=TYPES12).to_csv(CACHE / f"{TAG}rep_{out}.csv")
p6 = [x for x in ref6.index if x in b.index]
pd.DataFrame(deconvolve(ref6.loc[p6].to_numpy(), b.loc[p6].to_numpy().astype(float)),
             index=b.columns, columns=TYPES6).to_csv(CACHE / f"{TAG}rep_comp6.csv")
ages, ok2 = {}, True
print("CHECAGEM 2 — relogios:")
for c in CLOCKS:
    s, _ = predict(b.astype(float), c, min_coverage=0.0); ages[c] = s
    r = float(np.corrcoef(rep.age, s)[0, 1]); ok2 &= r > 0.5
    print(f"  {c:<13} r com idade {r:.3f}")
dp, _ = dunedinpace.predict(b.astype(float))
pd.DataFrame(ages).assign(chrono=rep.age, person=rep.person, group=rep.group,
                          DunedinPACE=dp).to_csv(CACHE / f"{TAG}rep_ages.csv")
print(f"  -> {'ok' if ok2 else 'FALHOU'}\ncache das replicas escrito")
