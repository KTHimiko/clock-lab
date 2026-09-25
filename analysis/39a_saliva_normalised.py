#!/usr/bin/env python3
"""
Stage 39a — GSE78874 again, quantile-normalised to the EPIC saliva cohorts.

Stage 38 found transports between saliva cohorts harmful in 7 of 8 cells, but
every pair crossed array and preprocessing: GSE78874 is 450k with betas computed
here from raw signal (no normalisation), the other two are processed EPIC betas.
This recomputes GSE78874's composition, clocks and DunedinPACE after
quantile-normalising each of its samples to the mean sorted beta profile of
GSE232891, over the probes both carry among those the panels and clocks need.

It is a crude harmonisation. It equalises the distribution of each sample, not
probe-type differences, but it removes gross scale compression, and that is what
could turn a correct coefficient into a wrong one after transport.

PRE-REGISTERED
  1. after normalisation every usable clock still tracks age r > 0.5. Hard stop.
  2. the normalised immune fraction correlates with the raw one r > 0.8 (the
     normalisation reorders nothing important). Reported either way.

Writes GSE78874qn_* caches alongside the raw ones.
Usage:  .venv/bin/python analysis/39a_saliva_normalised.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import rdata

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from model.clocks import CLOCKS, predict, load_coefficients
from model.deconvolution import deconvolve
from model import dunedinpace

DATA = ROOT / "reference/data"; CACHE = ROOT / "results/cache"


def rda(name):
    v = rdata.read_rda(str(DATA / f"epidish/EpiDISH/data/{name}.rda"))[name]
    return pd.DataFrame(v.values, index=[str(x) for x in v.coords[v.dims[0]].values],
                        columns=[str(x) for x in v.coords[v.dims[1]].values])


REF1 = rda("centEpiFibIC.m"); REF2 = rda("centBloodSub.m")
needed = set(REF1.index) | set(REF2.index) | set(dunedinpace.M["gold_probes"])
for c in CLOCKS:
    needed |= set(load_coefficients(c)[0].index)


def rows(path, cols=None):
    parts = []
    for ch in pd.read_csv(path, index_col=0, chunksize=40_000, usecols=cols):
        parts.append(ch.loc[ch.index.intersection(needed)])
    return pd.concat(parts)


print("lendo GSE232891 (referencia) ...", flush=True)
hdr = pd.read_csv(DATA / "GSE232891_Processed_Beta_Values.csv.gz", nrows=0).columns
ref = rows(DATA / "GSE232891_Processed_Beta_Values.csv.gz",
           [hdr[0]] + [c for c in hdr[1:] if not c.endswith("_detP")]).astype(float)
print("lendo GSE78874 (sinais) ...", flush=True)
raw = rows(DATA / "GSE78874_datSignal.csv.gz")
ids = sorted({c[:-len("Methylated")] for c in raw.columns if c.endswith("Methylated") and not c.endswith("Unmethylated")})
bet = pd.DataFrame(raw[[i + "Methylated" for i in ids]].to_numpy(float)
                   / (raw[[i + "Methylated" for i in ids]].to_numpy(float)
                      + raw[[i + "Unmethylated" for i in ids]].to_numpy(float) + 100),
                   index=raw.index, columns=ids)
common = bet.index.intersection(ref.index)
target = np.sort(ref.loc[common].to_numpy(), axis=0)
target = np.nanmean(target, axis=1)
qn = bet.copy()
X = bet.loc[common].to_numpy()
order = np.argsort(X, axis=0)
Y = np.empty_like(X)
for j in range(X.shape[1]):
    Y[order[:, j], j] = target
qn.loc[common] = Y
print(f"  {len(common):,} sondas em comum normalizadas; {len(bet.index) - len(common):,} ficam brutas")

ages_raw = pd.read_csv(CACHE / "GSE78874_ages.csv", index_col=0)
meas_raw = pd.read_csv(CACHE / "GSE78874_compmeas.csv", index_col=0)
qn = qn[ages_raw.index]
chrono = ages_raw.chrono
p1 = [p for p in REF1.index if p in qn.index]; p2 = [p for p in REF2.index if p in qn.index]
f1 = deconvolve(REF1.loc[p1].to_numpy(), qn.loc[p1].to_numpy()); f1 = f1 / f1.sum(1, keepdims=True)
f2 = deconvolve(REF2.loc[p2].to_numpy(), qn.loc[p2].to_numpy()); f2 = f2 / f2.sum(1, keepdims=True)
fit = np.column_stack([f1[:, 0], f1[:, 1], f2 * f1[:, [2]]])
pd.DataFrame(fit, index=qn.columns, columns=["Epi", "Fib"] + list(REF2.columns)).to_csv(CACHE / "GSE78874qn_compfit.csv")
pd.DataFrame(f1, index=qn.columns, columns=["Epi", "Fib", "IC"]).to_csv(CACHE / "GSE78874qn_compmeas.csv")
r_ic = float(np.corrcoef(f1[:, 2], meas_raw.loc[qn.columns, "IC"])[0, 1])
print(f"CHECAGEM 2 — IC normalizado x bruto: r = {r_ic:.3f}; IC mediano {np.median(f1[:, 2]):.2f} (bruto {meas_raw.IC.median():.2f})"
      f" -> {'ok' if r_ic > 0.8 else 'FALHOU'}")
ages, ok1 = {}, True
print("CHECAGEM 1 — relogios apos normalizacao:")
for c in CLOCKS:
    s, _ = predict(qn, c, min_coverage=0.0); ages[c] = s
    rr = float(np.corrcoef(chrono, s)[0, 1])
    if c in ("Levine2018", "Horvath2018"):
        ok1 &= rr > 0.5
    print(f"  {c:<13} r {rr:.3f} (bruto {np.corrcoef(chrono, ages_raw[c])[0, 1]:.3f})")
if not ok1:
    sys.exit("parando.")
dp, _ = dunedinpace.predict(qn)
pd.DataFrame(ages).assign(chrono=chrono, group=ages_raw.group).to_csv(CACHE / "GSE78874qn_ages.csv")
dp.to_frame().to_csv(CACHE / "GSE78874qn_dpace.csv")
print("cache GSE78874qn escrito")
