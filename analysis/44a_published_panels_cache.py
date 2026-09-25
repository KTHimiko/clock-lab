#!/usr/bin/env python3
"""
Stage 44a — recompute blood composition with the published reference libraries.

The panels used since stage 13 were built here, by a simpler probe selection than
the published libraries, and the twelve-type sensitivity analyses scored the
correction with probes the fit had already used. Both are listed as limitations.
The EpiDISH package ships the published alternatives:

  fit          cent12CT450k.m — the Salas et al. 2022 twelve-type reference on
               its IDOL-optimised 450k probe set (600 probes). Same twelve types
               as the panel built here, chosen and weighted by its authors.
  measurement  centDHSbloodDMC.m — the published seven-type blood reference
               (333 probes).

The two share 20 probes out of 600 and 333, so measurement is nearly independent
of the fit, which the home-built pair was not.

This writes {tag}_comp12pub.csv and {tag}_comp7pub.csv for the six blood cohorts;
44b re-runs the central results on them. Clock ages are not recomputed — they do
not depend on the panel.

PRE-REGISTERED
  1. coverage: each cohort carries >= 90% of both probe sets. Hard stop.
  2. the two panels agree with the home-built ones on what they both measure:
     median correlation across cohorts > 0.7 for neutrophils and for naive CD8.
     Reported either way — disagreement would mean the panels measure different
     things, which is itself the answer to the limitation.

Usage:  .venv/bin/python analysis/44a_published_panels_cache.py
"""
import sys, gzip
from pathlib import Path
import numpy as np
import pandas as pd
import rdata

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_extended import TYPES12
from model.deconvolution import deconvolve, TYPES as TYPES6

DATA = ROOT / "reference/data"; CACHE = ROOT / "results/cache"
ED = DATA / "epidish/EpiDISH/data"


def rda(name):
    v = rdata.read_rda(str(ED / f"{name}.rda"))[name]
    return pd.DataFrame(v.values, index=[str(x) for x in v.coords[v.dims[0]].values],
                        columns=[str(x) for x in v.coords[v.dims[1]].values])


REF12 = rda("cent12CT450k.m")
REF7 = rda("centDHSbloodDMC.m")
# the published names, in the order this project uses
MAP12 = {"CD4Tnv": "CD4nv", "Baso": "Bas", "CD4Tmem": "CD4mem", "Bmem": "Bmem",
         "Bnv": "Bnv", "Treg": "Treg", "CD8Tmem": "CD8mem", "CD8Tnv": "CD8nv",
         "Eos": "Eos", "NK": "NK", "Neu": "Neu", "Mono": "Mono"}
REF12 = REF12.rename(columns=MAP12)[TYPES12]
T7 = list(REF7.columns)
NEEDED = set(REF12.index) | set(REF7.index)

SERIES = ["GSE40279", "GSE61151", "GSE50660", "GSE42861"]


def read_rows(path, **kw):
    parts = []
    for ch in pd.read_csv(path, index_col=0, chunksize=50_000, **kw):
        parts.append(ch.loc[ch.index.intersection(NEEDED)].astype(np.float32))
    return pd.concat(parts)


def read_series(tag):
    """beta table of a GEO series matrix, only the rows needed"""
    path = DATA / f"{tag}_series_matrix.txt.gz"
    with gzip.open(path, "rt", errors="replace") as fh:
        skip = 0
        for line in fh:
            skip += 1
            if line.startswith("!series_matrix_table_begin"):
                break
    df = read_rows(path, sep="\t", skiprows=skip, quotechar='"',
                   na_values=["null", "NA", ""], low_memory=False)
    df.index = [str(i).strip('"') for i in df.index]
    df.columns = [str(c).strip('"') for c in df.columns]
    return df.loc[[i for i in df.index if i in NEEDED]]


print("CHECAGEM 1 — COBERTURA E DECONVOLUCAO")
rows = []
for tag in SERIES + ["GSE132203", "GSE55763"]:
    if tag == "GSE132203":
        src = DATA / "GSE132203_Geo_Submission_GTPEpic.csv.gz"
        hdr = pd.read_csv(src, nrows=0).columns.tolist()
        cols = [0] + [i for i, h in enumerate(hdr) if i > 0 and not h.startswith("Detection")]
        b = read_rows(src, usecols=cols)
    elif tag == "GSE55763":
        src = DATA / "GSE55763_normalized_betas.txt.gz"
        hdr = pd.read_csv(src, sep="\t", nrows=0).columns.tolist()
        cols = [0] + [i for i, h in enumerate(hdr) if i > 0 and not h.startswith("Detection")]
        b = read_rows(src, sep="\t", usecols=cols)
    else:
        b = read_series(tag)
    a = pd.read_csv(CACHE / f"{tag}_ages.csv", index_col=0)
    a = a[a.chrono.notna()]
    keep = [c for c in a.index if c in b.columns]
    b = b[keep]
    p12 = [p for p in REF12.index if p in b.index]
    p7 = [p for p in REF7.index if p in b.index]
    cov12, cov7 = len(p12) / len(REF12), len(p7) / len(REF7)
    ok = cov12 >= 0.90 and cov7 >= 0.90 and len(keep) == len(a)
    print(f"  {tag:<11}{len(keep):>6} amostras  12 tipos {cov12:>6.1%}  7 tipos {cov7:>6.1%}"
          f" -> {'ok' if ok else 'FALHOU'}", flush=True)
    if not ok:
        sys.exit("  parando.")
    c12 = deconvolve(REF12.loc[p12].to_numpy(), b.loc[p12].to_numpy().astype(float))
    c7 = deconvolve(REF7.loc[p7].to_numpy(), b.loc[p7].to_numpy().astype(float))
    pd.DataFrame(c12, index=keep, columns=TYPES12).to_csv(CACHE / f"{tag}_comp12pub.csv")
    pd.DataFrame(c7, index=keep, columns=T7).to_csv(CACHE / f"{tag}_comp7pub.csv")
    home = pd.read_csv(CACHE / f"{tag}_comp12.csv", index_col=0).loc[keep, TYPES12]
    rows.append(dict(tag=tag,
                     r_neu=float(np.corrcoef(c12[:, TYPES12.index("Neu")], home.Neu)[0, 1]),
                     r_cd8nv=float(np.corrcoef(c12[:, TYPES12.index("CD8nv")], home.CD8nv)[0, 1])))
    del b

A = pd.DataFrame(rows)
print("\nCHECAGEM 2 — CONCORDANCIA COM O PAINEL CASEIRO")
print(A.round(3).to_string(index=False))
c2 = A.r_neu.median() > 0.7 and A.r_cd8nv.median() > 0.7
print(f"  medianas: Neu {A.r_neu.median():.2f}, CD8nv {A.r_cd8nv.median():.2f} (barra 0,7)"
      f" -> {'ok' if c2 else 'FALHOU'}")
print("cache publicado escrito")
