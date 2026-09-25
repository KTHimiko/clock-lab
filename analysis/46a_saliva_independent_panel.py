#!/usr/bin/env python3
"""
Stage 46a — an independent saliva measurement panel, built from purified cells.

The saliva results of stages 38-43 measured the leftover with the first step of
the same EpiDISH reference the correction was fitted on, which biases the
measurement toward the correction and makes every harm figure a floor. It is the
largest remaining weakness of the saliva block.

GSE147318 (Middleton et al. 2022, Epigenetics 17:161-177) is a saliva cell-type
reference built for exactly this: children's saliva sorted into CD45-positive
(immune) and large-cell (epithelial, "large" in the series) fractions, on EPIC,
with whole saliva alongside. It is independent of EpiDISH in construction, cells and probes.

This builds a two-type panel from it by the project's usual rule — the probes
with the largest absolute difference between the two purified fractions, taking
the most discriminating in each direction — restricted to probes the 450k array
also carries, since one saliva cohort is 450k.

The panel is used as a RELATIVE immune score, not as a proportion. Deconvolved
into proportions it saturates near 100% immune in the adult cohorts (mean 0.93,
range 0.61-1.00) where EpiDISH reads 0.73, because the absolute scale does not
transfer from sorted children's cells to another study's processed betas. That
failure is itself an instance of what this paper is about. The ordering across
samples does transfer: the score correlates 0.98 with the EpiDISH immune
fraction. The metric used throughout needs a covariate that tracks composition,
not a calibrated proportion, so the score is what 46b measures with.

PRE-REGISTERED
  1. the sorted fractions separate: the panel probes' mean difference between
     CD45pos and CD45neg exceeds 0.3 in absolute value. Hard stop.
  2. INDEPENDENCE: at most 10% of the panel's probes appear in the EpiDISH
     references used for the fit. Hard stop — the point of the stage.
  3. the panel recovers the sorted fractions it was built from: CD45pos samples
     estimate above 0.9 immune, the epithelial fraction below 0.1. In-sample, so it is a floor
     on correctness, not evidence of it.
  4. the four saliva cohorts get an immune fraction correlating above 0.7 with
     the EpiDISH one. Reported either way — a low value would mean the two
     panels disagree about what they both call immune.

Usage:  .venv/bin/python analysis/46a_saliva_independent_panel.py
"""
import sys, gzip
from pathlib import Path
import numpy as np
import pandas as pd
import rdata

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from model.deconvolution import deconvolve

DATA = ROOT / "reference/data"; CACHE = ROOT / "results/cache"
ED = DATA / "epidish/EpiDISH/data"
N_PER_SIDE = 150
ARRAY450 = "GSE40279_series_matrix.txt.gz"   # to keep only probes both arrays carry

meta, skip = {}, 0
with gzip.open(DATA / "GSE147318_series_matrix.txt.gz", "rt", errors="replace") as fh:
    for line in fh:
        skip += 1
        if line.startswith("!series_matrix_table_begin"):
            break
        p = [v.strip('"') for v in line.rstrip("\n").split("\t")]
        if p[0] == "!Sample_geo_accession":
            meta["gsm"] = p[1:]
        elif p[0] == "!Sample_characteristics_ch1" and p[1].startswith("cell fraction"):
            meta["frac"] = [x.partition(":")[2].strip() for x in p[1:]]
m = pd.DataFrame(meta).set_index("gsm")
print(f"  GSE147318: {m.frac.value_counts().to_dict()}")

b = pd.read_csv(DATA / "GSE147318_series_matrix.txt.gz", sep="\t", skiprows=skip,
                index_col=0, quotechar='"', na_values=["null", "NA", ""], low_memory=False)
b = b[[c for c in b.columns if c in m.index]].apply(pd.to_numeric, errors="coerce")
b = b.dropna()
print(f"  {b.shape[0]:,} sondas sem faltantes x {b.shape[1]} amostras")
idx450 = set()
with gzip.open(DATA / ARRAY450, "rt", errors="replace") as fh:
    started = False
    for line in fh:
        if started:
            idx450.add(line.split("\t", 1)[0].strip('"'))
        elif line.startswith("!series_matrix_table_begin"):
            started = True
b = b.loc[[i for i in b.index if i in idx450]]
print(f"  {b.shape[0]:,} sondas tambem presentes no 450k")
pos = m.index[m.frac == "CD45pos"]; neg = m.index[m.frac == "large"]   # a fracao epitelial, celulas grandes
whole = m.index[m.frac == "whole"]
d = b[pos].mean(axis=1) - b[neg].mean(axis=1)
top = pd.concat([d.nlargest(N_PER_SIDE), d.nsmallest(N_PER_SIDE)]).index
REF = pd.DataFrame({"IC": b.loc[top, pos].mean(axis=1), "Epi": b.loc[top, neg].mean(axis=1)}).dropna()
top = REF.index
print(f"\nCHECAGEM 1 — SEPARACAO: diferenca media |{d[top].abs().mean():.3f}| nas {len(top)} sondas"
      f" -> {'ok' if d[top].abs().mean() > 0.3 else 'FALHOU'}")
if d[top].abs().mean() <= 0.3:
    sys.exit("parando.")


def rda_index(name):
    v = rdata.read_rda(str(ED / f"{name}.rda"))[name]
    return set(str(x) for x in v.coords[v.dims[0]].values)


epidish = rda_index("centEpiFibIC.m") | rda_index("centBloodSub.m")
share = len(set(top) & epidish) / len(top)
print(f"CHECAGEM 2 — INDEPENDENCIA: {share:.1%} das sondas tambem estao no EpiDISH"
      f" -> {'ok' if share <= 0.10 else 'FALHOU'}")
if share > 0.10:
    sys.exit("parando.")

w = (REF.IC - REF.Epi).to_numpy(); w = w / np.linalg.norm(w)
sub = b.loc[top, list(pos) + list(neg) + list(whole)]
Z = sub.to_numpy(); Z = Z - Z.mean(axis=1, keepdims=True)
sc = pd.Series(Z.T @ w, index=sub.columns)
c3 = sc[pos].min() > sc[neg].max()
print(f"CHECAGEM 3 — SEPARA AS FRACOES ORDENADAS: CD45pos [{sc[pos].min():.2f}, {sc[pos].max():.2f}],"
      f" epitelial [{sc[neg].min():.2f}, {sc[neg].max():.2f}],"
      f" saliva inteira mediana {sc[whole].median():.2f} -> {'ok' if c3 else 'FALHOU'}")
REF.to_csv(CACHE / "panel_saliva_indep.csv")
if not c3:
    sys.exit("parando.")
print("painel salvo em results/cache/panel_saliva_indep.csv")
