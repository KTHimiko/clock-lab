#!/usr/bin/env python3
"""
Stage 37 — why does GSE40279 keep a floor? Ancestry against batch.

Stages 35-36: fitted on all 2,639 samples of GSE55763, where estimation noise is
negligible, the correction transported into GSE40279 still leaves +2.4 (Levine
2018) and +3.5 (Horvath 2018) points of composition, while GSE50660 and GSE61151
get essentially none. Age extrapolation was ruled out in stage 36.

GSE40279 has two features the metadata can test:
  - ancestry: 426 "Caucasian - European" and 230 "Hispanic - Mexican". Horvath
    et al. (2016, Genome Biology 17:171) found Hispanics to have lower intrinsic
    and higher extrinsic epigenetic ageing than Caucasians in blood — so
    composition can relate to the clocks differently by ancestry.
  - batch: nine processing plates. A plate effect that shifts the composition
    estimates and the clocks together would look like a composition effect
    specific to this cohort.

Fit: unpenalised, on all of GSE55763, Levine 2018 and Horvath 2018 (the clocks
clean on GSE40279). Scoring as in stage 35 ("error left" = composition still
visible to the six-type panel after correction, net of its permutation null, as a
share of age-acceleration variance), computed within each target subset.

PRE-REGISTERED (written before any subset was scored)
  1. the metadata joins: every cached GSE40279 sample has an ethnicity and a
     plate, and the counts are 426 / 230. Hard stop.
  2. ANCESTRY: the error left in the Hispanic subset (230) is more than twice the
     median error left in 30 random Caucasian subsets of the same size, for both
     clocks.
  3. BATCH: centring the clocks and both composition matrices within plate in the
     target (fit untouched) removes more than half of the floor, for both clocks.
  Neither, one or both may hold; none are mutually exclusive.

Usage:  .venv/bin/python analysis/37_floor_ancestry_batch.py
"""
import sys, gzip
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_extended import TYPES12
from model.clocks import CLOCKS
from model.deconvolution import TYPES as TYPES6

RES = ROOT / "results"; CACHE = RES / "cache"
RNG = np.random.default_rng(20260924)
N_REPS = 20
COHORTS = ["GSE40279", "GSE61151", "GSE50660", "GSE42861"]
CLOCKS3 = ["Horvath2013", "Levine2018", "Horvath2018"]
ALPHA_FIX = 3.0

# the corrected rule: a clock is out of any pair in which a cohort it trained on
# appears, whichever end of the transport that cohort is on
TRAINED_ON = {"Hannum2013": {"GSE40279"}, "Horvath2013": {"GSE40279"}}


def allowed(clock, *cohorts):
    return not (TRAINED_ON.get(clock, set()) & set(cohorts))


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


def r2(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return 1 - (y - X @ beta).var() / y.var()


def inc_and_base(y, a, C):
    Xa = np.column_stack([np.ones(len(a)), a])
    base = r2(Xa, y)
    return r2(np.column_stack([Xa, C[:, :-1]]), y) - base, base


def null_mean(y, a, C, n_perm=800, rng=RNG):
    Xa = np.column_stack([np.ones(len(a)), a])
    base = r2(Xa, y)
    return float(np.mean([r2(np.column_stack([Xa, C[rng.permutation(len(C))][:, :-1]]), y)
                          - base for _ in range(n_perm)]))


def stratified_draw(ages, n, rng):
    q = pd.qcut(pd.Series(ages), 10, labels=False, duplicates="drop").to_numpy()
    idx = []
    for b in np.unique(q):
        pool = np.where(q == b)[0]
        take = min(max(int(round(n * len(pool) / len(ages))), 1), len(pool))
        idx.extend(rng.choice(pool, take, replace=False))
    idx = np.array(sorted(idx))
    if len(idx) > n:
        idx = np.sort(rng.choice(idx, n, replace=False))
    elif len(idx) < n:
        rest = np.setdiff1d(np.arange(len(ages)), idx)
        idx = np.sort(np.concatenate([idx, rng.choice(rest, n - len(idx), replace=False)]))
    return idx


def partial_out(age, Y):
    Xa = np.column_stack([np.ones(len(age)), age])
    return Y - Xa @ np.linalg.pinv(Xa) @ Y


def sigma(age, C):
    Ct = partial_out(age, C[:, :-1]); Ct = Ct - Ct.mean(axis=0)
    return Ct.T @ Ct / len(Ct)


def transport_index(af, Cf, at, Ct):
    return float(np.trace(np.linalg.inv(sigma(af, Cf)) @ sigma(at, Ct)) / len(af))


def ridge_coefs(age, C, y, alphas):
    Ct = partial_out(age, C[:, :-1])
    yt = partial_out(age, y.reshape(-1, 1)).ravel()
    mu, sd = Ct.mean(axis=0), Ct.std(axis=0) + 1e-12
    Z = (Ct - mu) / sd
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    scale = float((s ** 2).mean()); Uty = U.T @ yt
    return {a: (Vt.T @ (s * Uty / (s ** 2 + a * scale))) / sd for a in alphas}


# ------------------------------------------------------------------ data ----
panel = pd.read_csv(CACHE / "panel.csv", index_col=0).iloc[:, 0]
D = {}
for tag in COHORTS:
    a = pd.read_csv(CACHE / f"{tag}_ages.csv", index_col=0)
    k = a.chrono.notna().to_numpy(); ix = a.index[k]
    D[tag] = dict(chrono=a.chrono[k].to_numpy(),
                  y={c: a[c][k].to_numpy() for c in CLOCKS},
                  C12=pd.read_csv(CACHE / f"{tag}_comp12.csv", index_col=0)
                        .loc[ix, TYPES12].to_numpy(),
                  C6=pd.read_csv(CACHE / f"{tag}_comp6.csv", index_col=0)
                       .loc[ix, TYPES6].to_numpy())
N = {t: len(D[t]["chrono"]) for t in COHORTS}

section("CHECAGEM 1 — CACHE")
c1 = abs(panel.panel_r - 0.789) < 5e-3 and abs(panel.panel_mae - 0.027) < 5e-3
print(f"  painel r = {panel.panel_r:.3f} | MAE = {panel.panel_mae:.3f} -> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("  parando: cache diferente.")




BIG, T = "GSE55763", "GSE40279"
a = pd.read_csv(CACHE / f"{BIG}_ages.csv", index_col=0)
k = a.chrono.notna().to_numpy(); ix = a.index[k]
D[BIG] = dict(chrono=a.chrono[k].to_numpy(), y={c: a[c][k].to_numpy() for c in CLOCKS},
              C12=pd.read_csv(CACHE / f"{BIG}_comp12.csv", index_col=0).loc[ix, TYPES12].to_numpy())

section("CHECAGEM 1 — METADADOS DO GSE40279")
meta = {}
with gzip.open(ROOT / "reference/data" / f"{T}_series_matrix.txt.gz", "rt") as fh:
    for line in fh:
        if line.startswith("!series_matrix_table_begin"):
            break
        p = [v.strip('"') for v in line.rstrip("\n").split("\t")]
        if p[0] == "!Sample_geo_accession":
            meta["gsm"] = p[1:]
        elif p[0] == "!Sample_characteristics_ch1":
            key = p[1].partition(":")[0].strip()
            if key in ("ethnicity", "plate"):
                meta[key] = [x.partition(":")[2].strip() for x in p[1:]]
m = pd.DataFrame(meta).set_index("gsm")
at = pd.read_csv(CACHE / f"{T}_ages.csv", index_col=0)
tix = at.index[at.chrono.notna()]
m = m.reindex(tix)
counts = m.ethnicity.value_counts().to_dict()
c1 = m.notna().all().all() and counts.get("Caucasian - European") == 426 and counts.get("Hispanic - Mexican") == 230
print(f"  {len(m)} amostras; {counts}; placas {m.plate.nunique()} -> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("parando.")
print("  etnia x placa:")
print(pd.crosstab(m.plate, m.ethnicity).T.to_string())

TWO = ["Levine2018", "Horvath2018"]
s = D[BIG]
FIT = {}
for c in TWO:
    FIT[c] = (ridge_coefs(s["chrono"], s["C12"], s["y"][c], [0.0])[0.0], s["C12"].mean(axis=0))


def left(y, age, C12, C6, c):
    b, cbar = FIT[c]
    nul = null_mean(y, age, C6, n_perm=500)
    ib, base = inc_and_base(y, age, C6)
    yc = y - (C12[:, :-1] - cbar[:-1]) @ b
    ia, _ = inc_and_base(yc, age, C6)
    return (ia - nul) / (1 - base), (ib - nul) / (1 - base)


dT = D[T]
eth = m.ethnicity.to_numpy(); plate = m.plate.to_numpy()

section("REFERENCIA — GSE40279 INTEIRO")
base_left = {}
for c in TWO:
    l, b0 = left(dT["y"][c], dT["chrono"], dT["C12"], dT["C6"], c)
    base_left[c] = l
    print(f"  {c:<13} antes {b0:+.1%}  erro restante {l:+.1%}")

section("CHECAGEM 2 — ANCESTRALIDADE")
hisp = np.where(eth == "Hispanic - Mexican")[0]; cauc = np.where(eth == "Caucasian - European")[0]
c2 = True
rows = []
for c in TWO:
    lh, bh = left(dT["y"][c][hisp], dT["chrono"][hisp], dT["C12"][hisp], dT["C6"][hisp], c)
    lc_all, bc = left(dT["y"][c][cauc], dT["chrono"][cauc], dT["C12"][cauc], dT["C6"][cauc], c)
    sub = []
    for r in range(30):
        i = np.sort(np.random.default_rng(3700 + r).choice(cauc, len(hisp), replace=False))
        sub.append(left(dT["y"][c][i], dT["chrono"][i], dT["C12"][i], dT["C6"][i], c)[0])
    med = float(np.median(sub))
    ok = lh > 2 * med and lh > 0
    c2 &= ok
    rows.append(dict(clock=c, hispanic=lh, caucasian_all=lc_all, caucasian_230=med))
    print(f"  {c:<13} hispanicos {lh:+.1%} (antes {bh:+.1%}) | caucasianos todos {lc_all:+.1%} (antes {bc:+.1%}),"
          f" 230 sorteados mediana {med:+.1%} [{np.percentile(sub, 10):+.1%}, {np.percentile(sub, 90):+.1%}]"
          f" -> {'ok' if ok else 'nao'}")
print(f"  2. barra: hispanicos > 2x caucasianos pareados nos dois -> {'ok' if c2 else 'FALHOU'}")

section("CHECAGEM 3 — LOTE (CENTRAR POR PLACA NO ALVO)")


def centre(M, g):
    M = np.asarray(M, float).copy()
    for v in np.unique(g):
        w = g == v
        M[w] = M[w] - M[w].mean(axis=0) + M.mean(axis=0)
    return M


C12c = centre(dT["C12"], plate); C6c = centre(dT["C6"], plate)
c3 = True
for c in TWO:
    yc = centre(dT["y"][c], plate)
    l, b0 = left(yc, dT["chrono"], C12c, C6c, c)
    ok = l < 0.5 * base_left[c]
    c3 &= ok
    rows.append(dict(clock=c, plate_centred=l, raw=base_left[c]))
    print(f"  {c:<13} centrado por placa: antes {b0:+.1%}  erro restante {l:+.1%} (bruto {base_left[c]:+.1%})"
          f" -> {'ok' if ok else 'nao'}")
print(f"  3. barra: centrar por placa tira mais da metade do piso nos dois -> {'ok' if c3 else 'FALHOU'}")
pd.DataFrame(rows).to_csv(RES / "floor_ancestry_batch.csv", index=False)

section("FECHAMENTO")
print(f"  2. ancestralidade: {'ok' if c2 else 'FALHOU'}")
print(f"  3. lote: {'ok' if c3 else 'FALHOU'}")
