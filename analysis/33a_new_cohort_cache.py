#!/usr/bin/env python3
"""
Stage 33a — a fifth cohort on another array and another ancestry: GSE132203.

The Grady Trauma Project (GSE132203): 795 whole-blood samples on the EPIC array,
mostly African American, recruited for trauma exposure. It adds two axes the four
450k cohorts do not have — array generation and ancestry — and none of the clocks
in use was trained on it (Horvath 2013 predates EPIC; Levine 2018 used InCHIANTI;
Horvath 2018's training sets do not include it; DunedinPACE used the Dunedin Study).

GEO carries the betas as a 5.3 GB supplementary CSV, columns named by array
barcode with a detection p-value column after each, and the series matrix title
of each sample IS its barcode. Only the rows the analysis needs are kept: clock
probes, DunedinPACE's 20,000-probe background and the two deconvolution panels.

THE PANELS MUST BE THE SAME PANELS. The stage 18 cache stored compositions, not
the panels that produced them. They are rebuilt here by stage 18's deterministic
procedure — same references, same candidate scoring, same restriction to probes
common to the four original cohorts — and written to results/cache/ so no later
cohort has to do this again.

PRE-REGISTERED
  1. panel identity: the rebuilt twelve-type panel reproduces stage 18's
     fingerprint on the GSE167998 mixtures, r = 0.789 and MAE = 0.027 to three
     decimals. Hard stop: a different panel would make the new cohort
     incomparable with every earlier stage.
  2. panel coverage on EPIC: at least 90% of each panel's probes present in
     GSE132203. The deconvolution uses the present rows.
  3. clock coverage on EPIC: a clock is used on GSE132203 only if at least 95% of
     its probes are present; below that its predictions depend on which probes
     happen to be missing, and it is excluded from every pair with this cohort.
  4. every clock tracks age in GSE132203: r > 0.5 (for DunedinPACE, which is not
     an age estimator, reported only).

Usage:  .venv/bin/python analysis/33a_new_cohort_cache.py
"""
import sys, gzip
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix
from load_extended import read_extended, TYPES12
from model.clocks import CLOCKS, predict, load_coefficients
from model.deconvolution import PANEL, deconvolve, TYPES as TYPES6
from model import dunedinpace

DATA = ROOT / "reference/data"; CACHE = ROOT / "results/cache"
OLD = ["GSE40279", "GSE61151", "GSE50660", "GSE42861"]


def probe_ids(path):
    """Row names of a series matrix without keeping the table."""
    ids, on = [], False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                on = True; next(fh); continue
            if line.startswith("!series_matrix_table_end"):
                break
            if on:
                ids.append(line.split("\t", 1)[0].strip('"'))
    return pd.Index(ids)


def scores_and_cand(betas, labels, types, n_candidate=200):
    sc, cand = {}, set()
    arr = np.asarray(labels)
    for t in types:
        a = betas.loc[:, arr == t]
        c = betas.loc[:, np.isin(arr, [o for o in types if o != t])]
        dd = a.mean(axis=1) - c.mean(axis=1)
        sd = np.sqrt((a.var(axis=1, ddof=1) + c.var(axis=1, ddof=1)) / 2) + 0.01
        s = (dd / sd).replace([np.inf, -np.inf], np.nan).dropna()
        sc[t] = s
        cand |= set(s.nlargest(n_candidate).index) | set(s.nsmallest(n_candidate).index)
    return sc, cand


def final_panel(sc, cand, betas, labels, types, common, n_side=50):
    pool = sorted(cand & set(common)); arr = np.asarray(labels); probes = set()
    for t in types:
        s = sc[t].reindex(pool).dropna()
        probes |= set(s.nlargest(n_side).index) | set(s.nsmallest(n_side).index)
    probes = sorted(probes)
    return pd.DataFrame({t: betas.loc[probes, arr == t].mean(axis=1) for t in types}), probes


print("indices de sondas das quatro coortes originais ...", flush=True)
common = None
for t in OLD:
    ids = probe_ids(DATA / f"{t}_series_matrix.txt.gz")
    common = ids if common is None else common.intersection(ids)
print("referencias ...", flush=True)
b167, p167 = read_extended(DATA / "GSE167998_matrix_processed.txt.gz", DATA / "BloodExtended_Pheno.csv")
is_mix = (p167.CellType == "MIX").to_numpy()
b35, m35 = read_series_matrix(DATA / "GSE35069_series_matrix.txt.gz")
m35 = m35.set_index("gsm").reindex(b35.columns)
lab6 = m35["tissue/cell type"].map({v: k for k, v in PANEL.items()})
common = common.intersection(b167.index).intersection(b35.index)

sc12, cand12 = scores_and_cand(b167.loc[:, ~is_mix], p167.CellType[~is_mix], TYPES12)
sc6, cand6 = scores_and_cand(b35.loc[:, lab6.notna().to_numpy()], lab6.dropna(), TYPES6)
ref12, probes12 = final_panel(sc12, cand12, b167.loc[:, ~is_mix], p167.CellType[~is_mix], TYPES12, common)
ref6, probes6 = final_panel(sc6, cand6, b35.loc[:, lab6.notna().to_numpy()], lab6.dropna(), TYPES6, common)
known = p167.loc[is_mix, TYPES12].to_numpy()
est = deconvolve(ref12.to_numpy(), b167.loc[probes12].to_numpy()[:, is_mix])
r = float(np.corrcoef(known.ravel(), est.ravel())[0, 1]); mae = float(np.abs(known - est).mean())
c1 = round(r, 3) == 0.789 and round(mae, 3) == 0.027
print(f"CHECAGEM 1 — impressao digital do painel: r = {r:.3f}, MAE = {mae:.3f} -> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("parando: painel diferente do da etapa 18.")
ref12.to_csv(CACHE / "panel12.csv"); ref6.to_csv(CACHE / "panel6.csv")
del b167, b35

needed = set(probes12) | set(probes6) | set(dunedinpace.M["gold_probes"])
for c in CLOCKS:
    needed |= set(load_coefficients(c)[0].index)

print("lendo GSE132203 (so as linhas necessarias) ...", flush=True)
src = DATA / "GSE132203_Geo_Submission_GTPEpic.csv.gz"
header = pd.read_csv(src, nrows=0).columns.tolist()
beta_idx = [0] + [i for i, h in enumerate(header) if i > 0 and not h.startswith("Detection")]
parts = []
for chunk in pd.read_csv(src, usecols=beta_idx, index_col=0, chunksize=50_000, dtype=None):
    parts.append(chunk.loc[chunk.index.intersection(needed)].astype(np.float32))
betas = pd.concat(parts)
print(f"  {betas.shape[0]:,} sondas x {betas.shape[1]} amostras", flush=True)

sm = {}
with gzip.open(DATA / "GSE132203_series_matrix.txt.gz", "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("!series_matrix_table_begin"):
            break
        p = [v.strip('"') for v in line.rstrip("\n").split("\t")]
        if p[0] == "!Sample_title":
            sm["title"] = p[1:]
        elif p[0] == "!Sample_geo_accession":
            sm["gsm"] = p[1:]
        # the field name must match EXACTLY: this series also carries an
        # "age acceleration" field, and a prefix match let it overwrite age —
        # caught by check 4 on the first run, every clock at r ~ 0.02 with "age"
        elif p[0] == "!Sample_characteristics_ch1" and p[1].partition(":")[0].strip() == "age":
            sm["age"] = [float(x.partition(":")[2]) if x.partition(":")[2].strip() not in ("", "NA") else np.nan for x in p[1:]]
meta = pd.DataFrame(sm).set_index("title")
betas = betas[[c for c in betas.columns if c in meta.index]]
chrono = meta.age.reindex(betas.columns)
print(f"  com idade: {int(chrono.notna().sum())}")

cov12 = np.mean([p in betas.index for p in probes12]); cov6 = np.mean([p in betas.index for p in probes6])
c2 = cov12 >= 0.9 and cov6 >= 0.9
print(f"CHECAGEM 2 — cobertura dos paineis no EPIC: 12 tipos {cov12:.1%}, 6 tipos {cov6:.1%} -> {'ok' if c2 else 'FALHOU'}")
p12 = [p for p in probes12 if p in betas.index]; p6 = [p for p in probes6 if p in betas.index]
comp12 = deconvolve(ref12.loc[p12].to_numpy(), betas.loc[p12].to_numpy())
comp6 = deconvolve(ref6.loc[p6].to_numpy(), betas.loc[p6].to_numpy())
pd.DataFrame(comp12, index=betas.columns, columns=TYPES12).to_csv(CACHE / "GSE132203_comp12.csv")
pd.DataFrame(comp6, index=betas.columns, columns=TYPES6).to_csv(CACHE / "GSE132203_comp6.csv")

ages, usable = {}, []
print("CHECAGEM 3 e 4 — cobertura dos relogios e correlacao com a idade:")
for c in CLOCKS:
    coefs = load_coefficients(c)[0]; cov = coefs.index.isin(betas.index).mean()
    s, _ = predict(betas, c, min_coverage=0.0)
    k = chrono.notna().to_numpy(); rr = float(np.corrcoef(chrono[k], s[k])[0, 1])
    ok = cov >= 0.95 and rr > 0.5
    if ok: usable.append(c)
    ages[c] = s
    print(f"  {c:<13} cobertura {cov:.1%}  r com idade {rr:.3f} -> {'usa' if ok else 'EXCLUIDO'}")
dp, dcov = dunedinpace.predict(betas)
print(f"  DunedinPACE   cobertura de fundo {dcov:.1%}  media {dp.mean():.3f} dp {dp.std():.3f}")
pd.DataFrame(ages).assign(chrono=chrono).to_csv(CACHE / "GSE132203_ages.csv")
dp.to_frame().to_csv(CACHE / "GSE132203_dpace.csv")
pd.Series(usable).to_csv(CACHE / "GSE132203_usable_clocks.csv", index=False)
print("cache escrito; relogios utilizaveis no GSE132203:", usable)
