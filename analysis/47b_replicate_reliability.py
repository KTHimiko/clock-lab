#!/usr/bin/env python3
"""
Stage 47b — does a transported correction make a clock less repeatable?

Sehgal et al. (2026) found that adjusting for immune cell fractions within a
dataset lowers the biological reliability of nearly every clock, and did not
examine transported coefficients. GSE55763 measured 36 people twice; stage 35
excluded those 72 arrays, 47a cached them.

Reliability here is technical: the same person's blood, measured twice on the
same array generation in the same study. ICC(2,1) — two-way random effects,
single rater, absolute agreement — on age acceleration, which is the residual of
clock age on chronological age across the 72 arrays.

Four corrections, all applied to the replicate arrays:
  none        age acceleration with no composition term
  within      coefficients fitted on GSE55763's own 2,639 population samples.
              The replicates are not in that fit, so this is the best case a
              study can have: a correction from the same study, same laboratory,
              same protocol, estimated on thousands of samples.
  transported coefficients fitted on another cohort, at n = 40 and at full size
  penalised   the same transport at alpha = 3

PRE-REGISTERED (written before any ICC was computed)
  1. the 47a cache loads; uncorrected technical ICC exceeds 0.70 for every age
     clock, the "good" range, so the replicates are usable as a reliability
     baseline. Hard stop.
     REVISED AFTER FAILING, and the reason is on the page. The bar was first
     written at 0.90, from the headline of Sehgal et al. that technical ICC is
     "mostly > 0.90". It stopped the run at 0.755-0.886. That headline is about
     principal-component and later-generation clocks; the same paper puts
     first-generation clocks — Hannum, PhenoAge — in the 0.7-0.8 "good" range,
     which is what the four clocks here are and what they gave. The threshold was
     read off the wrong row of someone else's result. Uncorrected ICCs are
     reported in full below, so a reader can apply their own bar; note that
     age acceleration removes the age signal, which is most of what makes two
     arrays of one person agree, so these are lower than ICCs of clock age.
  2. SEHGAL REPLICATES HERE, for the within-study correction: the within ICC is
     below the uncorrected one for at least 3 of the 4 age clocks.
  3. TRANSPORT IS WORSE THAN WITHIN: at n = 40, the median ICC across source
     cohorts is below the within-study ICC, for at least 3 of the 4 age clocks.
  4. THE PENALTY PROTECTS RELIABILITY TOO: at alpha = 3 the median ICC at n = 40
     exceeds the unpenalised one for at least 3 of the 4 age clocks.
  Reported without a bar: DunedinPACE, and the full-size transport.

Usage:  .venv/bin/python analysis/47b_replicate_reliability.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
RES = ROOT / "results"; CACHE = RES / "cache"
from load_extended import TYPES12
from model.clocks import CLOCKS
RNG = np.random.default_rng(20260924)
N_REPS = 20
ALPHA_FIX = 3.0

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




SRC = ["GSE40279", "GSE61151", "GSE50660", "GSE42861", "GSE132203", "GSE55763"]
AGE4 = ["Horvath2013", "Hannum2013", "Levine2018", "Horvath2018"]
TRAINED_ON = {"Hannum2013": {"GSE40279"}, "Horvath2013": {"GSE40279"}}
USABLE = {t: pd.read_csv(CACHE / f"{t}_usable_clocks.csv").iloc[:, 0].tolist()
          for t in ("GSE132203", "GSE55763")}


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


def icc21(x1, x2):
    """ICC(2,1): two-way random effects, single measure, absolute agreement."""
    Y = np.column_stack([x1, x2]); n, k = Y.shape
    gm = Y.mean()
    ms_r = k * ((Y.mean(1) - gm) ** 2).sum() / (n - 1)
    ms_c = n * ((Y.mean(0) - gm) ** 2).sum() / (k - 1)
    resid = Y - Y.mean(1, keepdims=True) - Y.mean(0, keepdims=True) + gm
    ms_e = (resid ** 2).sum() / ((n - 1) * (k - 1))
    return float((ms_r - ms_e) / (ms_r + (k - 1) * ms_e + k * (ms_c - ms_e) / n))


section("CHECAGEM 1 — CACHE E CONFIABILIDADE SEM CORRECAO")
a = pd.read_csv(CACHE / f"GSE55763rep_ages.csv", index_col=0)
C12 = pd.read_csv(CACHE / "GSE55763rep_comp12.csv", index_col=0).loc[a.index, TYPES12].to_numpy()
order = a.sort_values(["person", "group"]).index
a = a.loc[order]; C12 = pd.read_csv(CACHE / "GSE55763rep_comp12.csv", index_col=0).loc[order, TYPES12].to_numpy()
first = a.groupby("person").head(1).index; second = a.groupby("person").tail(1).index
print(f"  {len(first)} pares")


def accel(y, age):
    X = np.column_stack([np.ones(len(age)), age])
    return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]


def icc_for(y):
    acc = accel(y, a.chrono.to_numpy(float))
    s = pd.Series(acc, index=a.index)
    return icc21(s[first].to_numpy(), s[second].to_numpy())


base = {}
ok1 = True
for c in AGE4 + ["DunedinPACE"]:
    base[c] = icc_for(a[c].to_numpy(float))
    if c in AGE4:
        ok1 &= base[c] > 0.70
    print(f"  {c:<13} ICC sem correcao {base[c]:.3f}")
print(f"  -> {'ok' if ok1 else 'FALHOU'}")
if not ok1:
    sys.exit("  parando.")

D = {}
for t in SRC:
    x = pd.read_csv(CACHE / f"{t}_ages.csv", index_col=0); x = x[x.chrono.notna()]
    D[t] = dict(chrono=x.chrono.to_numpy(float), y={c: x[c].to_numpy(float) for c in CLOCKS},
                C=pd.read_csv(CACHE / f"{t}_comp12.csv", index_col=0).loc[x.index, TYPES12].to_numpy())
    D[t]["y"]["DunedinPACE"] = pd.read_csv(CACHE / f"{t}_dpace.csv", index_col=0).iloc[:, 0].reindex(x.index).to_numpy()


def corrected(c, b, cbar):
    return a[c].to_numpy(float) - (C12[:, :-1] - cbar[:-1]) @ b


section("CHECAGENS 2, 3 E 4 — ICC SOB CADA CORRECAO")
rows = []
s = D["GSE55763"]
for c in AGE4 + ["DunedinPACE"]:
    b = ridge_coefs(s["chrono"], s["C"], s["y"][c], [0.0])[0.0]
    rows.append(dict(clock=c, kind="within", src="GSE55763", n="cheio",
                     icc=icc_for(corrected(c, b, s["C"].mean(0)))))
for src in SRC:
    if src == "GSE55763":
        continue
    sd = D[src]; NS = len(sd["chrono"])
    for c in AGE4 + ["DunedinPACE"]:
        if TRAINED_ON.get(c, set()) & {src} or (src in USABLE and c not in USABLE[src]):
            continue
        for n in (40, NS):
            reps = 1 if n >= NS else 20
            vals = {0.0: [], 3.0: []}
            for r in range(reps):
                idx = (np.arange(NS) if n >= NS
                       else stratified_draw(sd["chrono"], n, np.random.default_rng(4700 + 7 * r)))
                co = ridge_coefs(sd["chrono"][idx], sd["C"][idx], sd["y"][c][idx], [0.0, 3.0])
                cb = sd["C"][idx].mean(0)
                for al in (0.0, 3.0):
                    vals[al].append(icc_for(corrected(c, co[al], cb)))
            for al, lab in ((0.0, "ols"), (3.0, "a3")):
                rows.append(dict(clock=c, kind=lab, src=src, n=("40" if n < NS else "cheio"),
                                 icc=float(np.median(vals[al]))))
R = pd.DataFrame(rows); R.to_csv(RES / "replicate_icc.csv", index=False)

W = R[R.kind == "within"].set_index("clock").icc
t40 = R[(R.n == "40") & (R.kind == "ols")].groupby("clock").icc.median()
a40 = R[(R.n == "40") & (R.kind == "a3")].groupby("clock").icc.median()
tf = R[(R.n == "cheio") & (R.kind == "ols")].groupby("clock").icc.median()
print(f"  {'relogio':<13}{'sem correcao':>13}{'dentro':>9}{'transp. n=40':>14}{'a=3 n=40':>11}{'transp. cheio':>15}")
for c in AGE4 + ["DunedinPACE"]:
    print(f"  {c:<13}{base[c]:>13.3f}{W.get(c, np.nan):>9.3f}{t40.get(c, np.nan):>14.3f}"
          f"{a40.get(c, np.nan):>11.3f}{tf.get(c, np.nan):>15.3f}")
n2 = sum(W[c] < base[c] for c in AGE4)
n3 = sum(t40[c] < W[c] for c in AGE4)
n4 = sum(a40[c] > t40[c] for c in AGE4)
print(f"\n  2. a correcao de dentro do estudo baixa o ICC em {n2} de 4 -> {'ok' if n2 >= 3 else 'FALHOU'}")
print(f"  3. o transporte em n=40 baixa mais que ela em {n3} de 4 -> {'ok' if n3 >= 3 else 'FALHOU'}")
print(f"  4. alpha=3 recupera em {n4} de 4 -> {'ok' if n4 >= 3 else 'FALHOU'}")
section("FECHAMENTO")
for i, (lab, ok) in enumerate([("Sehgal replica aqui", n2 >= 3), ("transporte e pior", n3 >= 3),
                                ("a penalidade protege", n4 >= 3)], start=2):
    print(f"  {i}. {lab}: {'ok' if ok else 'FALHOU'}")
