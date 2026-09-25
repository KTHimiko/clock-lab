#!/usr/bin/env python3
"""
Stage 38a — saliva: three adult cohorts, cached for a transport test in another tissue.

Every stage so far used blood. The manuscript lists other tissues as untested,
and the one published transport of composition coefficients found (Galkin et al.
2021) is in saliva. In saliva the composition axis is much larger than in blood:
saliva is a mix of buccal epithelial cells and leukocytes, and the immune
fraction ranges widely between people.

Cohorts (GEO, adults with age):
  GSE232891  EPIC, 552: Crohn's disease 227, ulcerative colitis 75, controls 250
  GSE232332  EPIC, 275: oesophageal cancer/dysplasia 98 with 6 technical
             replicates, controls 167 — same group as GSE232891
  GSE78874   450k, 259: PEG study (California), Caucasian and Hispanic; raw
             methylated/unmethylated signal, beta = M / (M + U + 100)

Composition (EpiDISH 2.28.0 references, Bioconductor):
  fit panel, 9 types — HEpiDISH scheme: Epi/Fib/IC from centEpiFibIC.m, then the
    immune fraction split into B, NK, CD4T, CD8T, Mono, Neutro, Eosino with
    centBloodSub.m (fractions renormalised within immune, times IC)
  measurement panel, 3 types — Epi/Fib/IC from centEpiFibIC.m
  The two share the first step, so measurement is biased toward the correction
  (as stage 29's twelve-type measurement was). Unlike blood, no independent
  saliva reference is available here.

Clocks: not trained on any of the three (GSE232891/GSE232332 are 2023; GSE78874
is 2016, after Horvath 2013, and absent from Horvath 2018's training table, whose
saliva set GSE111223 — the same PEG study — was a test set only; Levine 2018 and
DunedinPACE are blood clocks trained elsewhere).

PRE-REGISTERED
  1. parsing: every sample with age maps to exactly one data column; GSE232332's
     replicate columns are dropped; GSE78874 samples fail QC if more than 5% of
     probes have detection p > 0.01. Hard stop if fewer than 200 / 200 / 180
     remain.
  2. OVERLAP between GSE232891 and GSE232332 (same group, both with controls):
     on the genotyping probes (rs*), any cross-cohort pair correlating above 0.9
     is the same person and is removed from GSE232332. Reported; if no rs probes
     are present in either file, the check is reported as not possible and the
     two are treated as possibly overlapping (never both used as fit and target).
  3. saliva sanity: median immune fraction between 0.2 and 0.9 in each cohort.
  4. each clock tracks age r > 0.5 with coverage >= 95%, per cohort.

Usage:  .venv/bin/python analysis/38a_saliva_cache.py
"""
import sys, gzip, ast
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
IMM = list(REF2.columns)
FIT_TYPES = ["Epi", "Fib"] + IMM
MEAS_TYPES = ["Epi", "Fib", "IC"]
needed = set(REF1.index) | set(REF2.index) | set(dunedinpace.M["gold_probes"])
for c in CLOCKS:
    needed |= set(load_coefficients(c)[0].index)


def keep_row(i):
    return i in needed or str(i).startswith("rs")


def series_meta(tag, fields, desc_idx=None):
    out = {}
    with gzip.open(DATA / f"{tag}_series_matrix.txt.gz", "rt", errors="replace") as fh:
        nd = 0
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            p = [v.strip('"') for v in line.rstrip("\n").split("\t")]
            if p[0] in ("!Sample_title", "!Sample_source_name_ch1", "!Sample_geo_accession"):
                out[p[0][8:]] = p[1:]
            elif p[0] == "!Sample_description":
                out[f"desc{nd}"] = p[1:]; nd += 1
            elif p[0] == "!Sample_characteristics_ch1":
                key = p[1].partition(":")[0].strip()
                if key in fields:
                    out[key] = [x.partition(":")[2].strip() for x in p[1:]]
    m = pd.DataFrame(out)
    m["age"] = pd.to_numeric(m["age"], errors="coerce")
    return m


def read_rows(path, sep, names=None, skip=0):
    parts = []
    for ch in pd.read_csv(path, sep=sep, header=None if names else 0, names=names,
                          skiprows=skip, index_col=0, chunksize=40_000, low_memory=False):
        parts.append(ch.loc[[i for i in ch.index if keep_row(i)]])
    return pd.concat(parts)


B, META = {}, {}
section = lambda t: print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)
section("CHECAGEM 1 — LEITURA")

# GSE232891: CSV, columns <id>, <id>_detP; id = source name
tag = "GSE232891"
m = series_meta(tag, ("age",)); m["id"] = m["source_name_ch1"]; m["group"] = m["desc0"].str.replace("Saliva Sample from ", "")
raw = read_rows(DATA / f"{tag}_Processed_Beta_Values.csv.gz", ",")
cols = [c for c in raw.columns if not c.endswith("_detP")]
B[tag] = raw[cols].astype(np.float32); META[tag] = m.set_index("id")
print(f"  {tag}: {B[tag].shape[0]:,} linhas x {B[tag].shape[1]} amostras; metadados {len(m)}; grupos {m.group.value_counts().to_dict()}", flush=True)

# GSE232332: header is a python-list string; data rows are probe + (beta, detP) pairs
tag = "GSE232332"
with gzip.open(DATA / f"{tag}_Matrix_Processed.txt.gz", "rt") as fh:
    head = fh.readline().strip().strip('"'); first = fh.readline()
hdr = ast.literal_eval(head)
nf = len(first.rstrip("\n").split(","))
print(f"  {tag}: cabecalho com {len(hdr)} nomes, primeira linha com {nf} campos")
if nf == len(hdr):
    names = ["probe"] + hdr[1:]
elif nf == len(hdr) - 1:
    names = ["probe"] + hdr[2:]
else:
    sys.exit("  parando: cabecalho e linhas nao batem.")
raw = read_rows(DATA / f"{tag}_Matrix_Processed.txt.gz", ",", names=names, skip=1)
cols = [c for c in raw.columns if not str(c).endswith("_detPValue") and "Replicate" not in str(c) and str(c).strip()]
m = series_meta(tag, ("age", "disease subtype")); m["id"] = m["title"]
m["group"] = np.where(m["desc0"].str.contains("Control"), "Control", "Cancer")
m = m[~m.id.str.contains("Replicate")]
B[tag] = raw[cols].apply(pd.to_numeric, errors="coerce").astype(np.float32); META[tag] = m.set_index("id")
print(f"  {tag}: {B[tag].shape[1]} amostras sem replicas; metadados {len(m)}; grupos {m.group.value_counts().to_dict()}", flush=True)

# GSE78874: signals
tag = "GSE78874"
m = series_meta(tag, ("age", "ethnicity")); m["id"] = m["desc1"]; m["group"] = m["ethnicity"]
raw = read_rows(DATA / f"{tag}_datSignal.csv.gz", ",")
ids = sorted({c[:-len("Methylated")] for c in raw.columns if c.endswith("Methylated") and not c.endswith("Unmethylated")})
Mx = raw[[i + "Methylated" for i in ids]].to_numpy(float); Ux = raw[[i + "Unmethylated" for i in ids]].to_numpy(float)
Px = raw[[i + "DectionPvalue" for i in ids]].to_numpy(float)
fail = (Px > 0.01).mean(axis=0) > 0.05
B[tag] = pd.DataFrame((Mx / (Mx + Ux + 100)).astype(np.float32), index=raw.index, columns=ids).loc[:, ~fail]
META[tag] = m.set_index("id")
print(f"  {tag}: {len(ids)} amostras, {int(fail.sum())} reprovadas no QC, ficam {B[tag].shape[1]}", flush=True)

MIN = {"GSE232891": 200, "GSE232332": 200, "GSE78874": 180}
ok1 = True
for t in B:
    common = [c for c in B[t].columns if c in META[t].index and pd.notna(META[t].loc[c, "age"])]
    dup = META[t].index.duplicated().any()
    B[t] = B[t][common]
    ok = len(common) >= MIN[t] and not dup
    ok1 &= ok
    print(f"  {t}: {len(common)} com idade e dado (minimo {MIN[t]}) -> {'ok' if ok else 'FALHOU'}")
if not ok1:
    sys.exit("  parando.")

section("CHECAGEM 2 — SOBREPOSICAO GSE232891 x GSE232332")
rs = {t: B[t].loc[[i for i in B[t].index if str(i).startswith("rs")]] for t in ("GSE232891", "GSE232332")}
common_rs = rs["GSE232891"].index.intersection(rs["GSE232332"].index)
overlap_checked = len(common_rs) >= 30
if overlap_checked:
    a = rs["GSE232891"].loc[common_rs].to_numpy(float); b = rs["GSE232332"].loc[common_rs].to_numpy(float)
    a = (a - a.mean(0)) / a.std(0); b = (b - b.mean(0)) / b.std(0)
    R = a.T @ b / len(common_rs)
    hits = np.argwhere(R > 0.9)
    print(f"  {len(common_rs)} sondas rs; r maximo cruzado mediano {np.median(R.max(axis=0)):.2f}; pares > 0,9: {len(hits)}")
    drop = sorted({B["GSE232332"].columns[j] for _, j in hits})
    for i, j in hits[:10]:
        print(f"    {B['GSE232891'].columns[i]} ~ {B['GSE232332'].columns[j]}  r = {R[i, j]:.3f}")
    B["GSE232332"] = B["GSE232332"].drop(columns=drop)
    print(f"  removidas do GSE232332: {len(drop)}; ficam {B['GSE232332'].shape[1]}")
else:
    print(f"  so {len(common_rs)} sondas rs em comum: checagem impossivel; as duas nunca serao par ajuste/alvo")
pd.Series({"overlap_checked": overlap_checked}).to_csv(CACHE / "saliva_overlap.csv")

section("CHECAGEM 3 E 4 — COMPOSICAO E RELOGIOS")
usable_all = {}
for t, bet in B.items():
    p1 = [p for p in REF1.index if p in bet.index]; p2 = [p for p in REF2.index if p in bet.index]
    f1 = deconvolve(REF1.loc[p1].to_numpy(), bet.loc[p1].to_numpy().astype(float))
    f1 = f1 / f1.sum(axis=1, keepdims=True)
    f2 = deconvolve(REF2.loc[p2].to_numpy(), bet.loc[p2].to_numpy().astype(float))
    f2 = f2 / f2.sum(axis=1, keepdims=True)
    fit = np.column_stack([f1[:, 0], f1[:, 1], f2 * f1[:, [2]]])
    pd.DataFrame(fit, index=bet.columns, columns=FIT_TYPES).to_csv(CACHE / f"{t}_compfit.csv")
    pd.DataFrame(f1, index=bet.columns, columns=MEAS_TYPES).to_csv(CACHE / f"{t}_compmeas.csv")
    ic = float(np.median(f1[:, 2]))
    print(f"  {t}: sondas painel {len(p1)}/{len(REF1)} e {len(p2)}/{len(REF2)}; IC mediano {ic:.2f}"
          f" [{np.percentile(f1[:, 2], 10):.2f}, {np.percentile(f1[:, 2], 90):.2f}] -> {'ok' if 0.2 <= ic <= 0.9 else 'FALHOU'}")
    chrono = META[t].loc[bet.columns, "age"].astype(float)
    ages, usable = {}, []
    for c in CLOCKS:
        cov = load_coefficients(c)[0].index.isin(bet.index).mean()
        s, _ = predict(bet.astype(float), c, min_coverage=0.0); ages[c] = s
        rr = float(np.corrcoef(chrono, s)[0, 1])
        ok = cov >= 0.95 and rr > 0.5
        if ok: usable.append(c)
        print(f"    {c:<13} cobertura {cov:.1%}  r {rr:.3f} -> {'usa' if ok else 'EXCLUIDO'}")
    dp, dcov = dunedinpace.predict(bet.astype(float))
    print(f"    DunedinPACE   fundo {dcov:.1%}  media {dp.mean():.3f} dp {dp.std():.3f}  r com idade {np.corrcoef(chrono, dp)[0, 1]:.2f}")
    pd.DataFrame(ages).assign(chrono=chrono, group=META[t].loc[bet.columns, "group"]).to_csv(CACHE / f"{t}_ages.csv")
    dp.to_frame().to_csv(CACHE / f"{t}_dpace.csv")
    pd.Series(usable).to_csv(CACHE / f"{t}_usable_clocks.csv", index=False)
    usable_all[t] = usable
print("\ncache escrito;", usable_all)
