#!/usr/bin/env python3
"""
Stage 35a — a large fitting cohort: GSE55763, about 2,700 samples on 450k.

Stage 24 found an n-independent floor — the part of the transported correction's
error that more fitting data does not remove — but the largest fitting cohort was
656. With several thousand samples the floor can be measured where estimation
noise is truly negligible, instead of inferred.

GSE55763 (Lehne et al. 2015, peripheral blood, London) carries normalised betas as
a 10.4 GB tab-separated supplement, a detection p-value column after each sample,
columns named by array barcode; the series-matrix title is "Peripheral blood,
<barcode>". Some individuals were measured twice for a technical-replication study;
each individual must appear once.

None of the clocks in use was trained on it (released 2015; Horvath 2013 predates
it; Levine 2018 used InCHIANTI; Horvath 2018's training sets do not include it;
DunedinPACE used the Dunedin Study).

The deconvolution panels are loaded from results/cache/ (saved in stage 33a with
stage 18's fingerprint), not rebuilt.

PRE-REGISTERED
  1. the saved panels load and match their fingerprint on the GSE167998 mixtures
     (r = 0.789, MAE = 0.027). Hard stop.
  2. replicates: 36 individuals were measured twice ("Technical replicate group
     1|2, sample k"; 25 of those 72 arrays are also flagged "population study").
     Only arrays whose description is exactly "Population study sample." are kept,
     which drops every replicate array and so cannot count anyone twice; the kept
     count must be between 2,500 and 2,800. Hard stop.
  3. every age clock tracks age: r > 0.5; coverage >= 95%.
  4. metadata fields matched by EXACT name (stage 33's lesson).

Usage:  .venv/bin/python analysis/35a_large_cohort_cache.py
"""
import sys, gzip
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_extended import read_extended, TYPES12
from model.clocks import CLOCKS, predict, load_coefficients
from model.deconvolution import deconvolve, TYPES as TYPES6
from model import dunedinpace

DATA = ROOT / "reference/data"; CACHE = ROOT / "results/cache"; TAG = "GSE55763"

ref12 = pd.read_csv(CACHE / "panel12.csv", index_col=0)[TYPES12]
ref6 = pd.read_csv(CACHE / "panel6.csv", index_col=0)[TYPES6]
b167, p167 = read_extended(DATA / "GSE167998_matrix_processed.txt.gz", DATA / "BloodExtended_Pheno.csv")
mix = (p167.CellType == "MIX").to_numpy()
est = deconvolve(ref12.to_numpy(), b167.loc[ref12.index].to_numpy()[:, mix])
known = p167.loc[mix, TYPES12].to_numpy()
r = float(np.corrcoef(known.ravel(), est.ravel())[0, 1]); mae = float(np.abs(known - est).mean())
c1 = round(r, 3) == 0.789 and round(mae, 3) == 0.027
print(f"CHECAGEM 1 — paineis salvos: r = {r:.3f}, MAE = {mae:.3f} -> {'ok' if c1 else 'FALHOU'}", flush=True)
if not c1:
    sys.exit("parando.")
del b167

meta = {}
with gzip.open(DATA / f"{TAG}_series_matrix.txt.gz", "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("!series_matrix_table_begin"):
            break
        p = [v.strip('"') for v in line.rstrip("\n").split("\t")]
        if p[0] == "!Sample_title":
            meta["barcode"] = [x.split(",")[-1].strip() for x in p[1:]]
        elif p[0] == "!Sample_description":
            meta["desc"] = p[1:]
        elif p[0] == "!Sample_characteristics_ch1":
            key = p[1].partition(":")[0].strip()
            if key == "age":
                meta[key] = [x.partition(":")[2].strip() for x in p[1:]]
m = pd.DataFrame(meta).set_index("barcode")
m["age"] = pd.to_numeric(m["age"], errors="coerce")
keep = m[m.desc == "Population study sample."].index
c2 = 2500 <= len(keep) <= 2800
print(f"CHECAGEM 2 — amostras do estudo populacional: {len(keep)} de {len(m)} -> {'ok' if c2 else 'FALHOU'}", flush=True)
if not c2:
    sys.exit("parando.")

needed = set(ref12.index) | set(ref6.index) | set(dunedinpace.M["gold_probes"])
for c in CLOCKS:
    needed |= set(load_coefficients(c)[0].index)
src = DATA / f"{TAG}_normalized_betas.txt.gz"
header = pd.read_csv(src, sep="\t", nrows=0).columns.tolist()
cols = [0] + [i for i, h in enumerate(header) if i > 0 and not h.startswith("Detection")]
parts = []
print("lendo (so as linhas necessarias) ...", flush=True)
for ch in pd.read_csv(src, sep="\t", usecols=cols, index_col=0, chunksize=50_000):
    parts.append(ch.loc[ch.index.intersection(needed)].astype(np.float32))
betas = pd.concat(parts)
betas = betas[[c for c in betas.columns if c in keep]]
chrono = m.age.reindex(betas.columns)
print(f"  {betas.shape[0]:,} sondas x {betas.shape[1]} amostras, com idade {int(chrono.notna().sum())}", flush=True)

p12 = [p for p in ref12.index if p in betas.index]; p6 = [p for p in ref6.index if p in betas.index]
pd.DataFrame(deconvolve(ref12.loc[p12].to_numpy(), betas.loc[p12].to_numpy()),
             index=betas.columns, columns=TYPES12).to_csv(CACHE / f"{TAG}_comp12.csv")
pd.DataFrame(deconvolve(ref6.loc[p6].to_numpy(), betas.loc[p6].to_numpy()),
             index=betas.columns, columns=TYPES6).to_csv(CACHE / f"{TAG}_comp6.csv")
ages, usable = {}, []
k = chrono.notna().to_numpy()
print("CHECAGEM 3 — relogios:")
for c in CLOCKS:
    cov = load_coefficients(c)[0].index.isin(betas.index).mean()
    s, _ = predict(betas, c, min_coverage=0.0); ages[c] = s
    rr = float(np.corrcoef(chrono[k], s[k])[0, 1])
    ok = cov >= 0.95 and rr > 0.5
    if ok: usable.append(c)
    print(f"  {c:<13} cobertura {cov:.1%}  r {rr:.3f} -> {'usa' if ok else 'EXCLUIDO'}")
dp, dcov = dunedinpace.predict(betas)
print(f"  DunedinPACE   fundo {dcov:.1%}  media {dp.mean():.3f} dp {dp.std():.3f}")
pd.DataFrame(ages).assign(chrono=chrono).to_csv(CACHE / f"{TAG}_ages.csv")
dp.to_frame().to_csv(CACHE / f"{TAG}_dpace.csv")
pd.Series(usable).to_csv(CACHE / f"{TAG}_usable_clocks.csv", index=False)
print("cache escrito; utilizaveis:", usable)
