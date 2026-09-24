#!/usr/bin/env python3
"""
Stage 29 — does the harm survive measurement on the naive/memory axis?

Every stage since 15 fits the correction with the twelve-type panel and measures
what is left with the six-type one, so that a correction is never scored against
its own representation of what it removed. A reviewer pointed out the cost: six
types cannot see the naive/memory split, which is where the introduction says the
largest composition effect lives. The harm might be an artefact of a measurement
blind to the axis that matters — or the benefit might be.

The only panel here that sees the naive/memory axis is the one the correction is
fitted with. Measuring with it is self-scoring: at home it returns zero by
construction, and out of cohort it is biased toward making the correction look
good. So the reading is asymmetric and fixed now:
    the harm PERSISTS measured on twelve types  -> robust to the measurement axis
    the harm VANISHES                           -> uninformative; self-scoring
                                                   could explain it

Two measurements beyond the usual six-type one:
  twelve   the increment in R^2 from all twelve proportions
  nm-axis  the increment from the seven naive/memory columns alone — CD4nv,
           CD4mem, CD8nv, CD8mem, Bnv, Bmem, Treg — the reviewer's concern exactly

PRE-REGISTERED
  1. the cache is the stage 18 cache. Hard stop.
  2. the n-curve at n = 40 (fitting on GSE40279, clean clocks): median net damage
     measured on twelve types still exceeds +5%.
  3. per-clock cells at matched n (stage 27): the share harmful under OLS,
     measured on twelve types, is at least 30% (it is 50% on six).
  4. ridge alpha = 3, measured on twelve types: share of harmful cells reported.
  5. the naive/memory axis alone: n = 40 damage and OLS harmful share, reported.

Usage:  .venv/bin/python analysis/29_measurement_axis.py
"""
import sys
from itertools import permutations
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


NM = [TYPES12.index(t) for t in ("CD4nv", "CD4mem", "CD8nv", "CD8mem", "Bnv", "Bmem", "Treg")]


def inc_cols(y, a, M):
    """R^2 increment from a block of columns used as they are (no sum-to-one
    column to drop — a subset of proportions does not sum to one)."""
    Xa = np.column_stack([np.ones(len(a)), a])
    base = r2(Xa, y)
    return r2(np.column_stack([Xa, M]), y) - base, base


MEAS = {
    "six":    lambda d: d["C6"][:, :-1],
    "twelve": lambda d: d["C12"][:, :-1],
    "nm":     lambda d: d["C12"][:, NM],
}


def perm_null(y, a, M, n_perm=600, rng=RNG):
    Xa = np.column_stack([np.ones(len(a)), a]); base = r2(Xa, y)
    return float(np.mean([r2(np.column_stack([Xa, M[rng.permutation(len(M))]]), y) - base
                          for _ in range(n_perm)]))


CL = ["Levine2018", "Horvath2018"]
print("  nulos por medicao ...", flush=True)
NULL, BEF = {}, {}
for t in COHORTS:
    for c in CL:
        for m, f in MEAS.items():
            M = f(D[t]); y = D[t]["y"][c]; a = D[t]["chrono"]
            NULL[(t, c, m)] = perm_null(y, a, M)
            ib, base = inc_cols(y, a, M)
            BEF[(t, c, m)] = ((ib - NULL[(t, c, m)]) / (1 - base), base)


def run(src, dst, n, reps, alphas=(0.0,)):
    s, d = D[src], D[dst]
    rows = []
    for r in range(reps):
        idx = (np.arange(N[src]) if n >= N[src]
               else stratified_draw(s["chrono"], n, np.random.default_rng(900 + 7 * r + n)))
        af, Cf = s["chrono"][idx], s["C12"][idx]
        cbar = Cf.mean(axis=0)
        for c in CL:
            coefs = ridge_coefs(af, Cf, s["y"][c][idx], list(alphas))
            for alpha, b in coefs.items():
                yc = d["y"][c] - (d["C12"][:, :-1] - cbar[:-1]) @ b
                for m, f in MEAS.items():
                    before, base = BEF[(dst, c, m)]
                    ia, _ = inc_cols(yc, d["chrono"], f(d))
                    rows.append(dict(src=src, dst=dst, n=len(idx), clock=c, alpha=alpha,
                                     meas=m, delta=(ia - NULL[(dst, c, m)]) / (1 - base) - before))
    return rows


section("CHECAGEM 2 — A CURVA EM N=40 E N CHEIO, NAS TRES MEDICOES")
rows = []
for dst in ["GSE61151", "GSE50660", "GSE42861"]:
    rows += run("GSE40279", dst, 40, 30) + run("GSE40279", dst, N["GSE40279"], 1)
C = pd.DataFrame(rows)
tab = C.groupby(["n", "meas"]).delta.median().unstack("meas")[["six", "twelve", "nm"]]
print(tab.to_string(float_format=lambda x: f"{x:+.1%}"))
c2 = tab.loc[40, "twelve"] > 0.05
print(f"\n  n=40 medido em doze: {tab.loc[40, 'twelve']:+.1%} (barra > +5%) -> "
      f"{'o dano PERSISTE' if c2 else 'nao informativo'}")

section("CHECAGEM 3 E 4 — CELULAS POR RELOGIO EM N CASADO, NAS TRES MEDICOES")
rows = []
for src, dst in permutations(COHORTS, 2):
    rows += run(src, dst, min(N[src], N[dst]), 20, alphas=(0.0, 3.0))
    print(f"  {src[3:]} -> {dst[3:]}", flush=True)
P = pd.DataFrame(rows)
P.to_csv(RES / "measurement_axis.csv", index=False)
cell = P.groupby(["src", "dst", "clock", "alpha", "meas"]).delta.median().reset_index()
print(f"\n  {'medicao':<9}{'OLS nocivas':>14}{'a=3 nocivas':>14}")
share = {}
for m in ("six", "twelve", "nm"):
    g = cell[cell.meas == m]
    o = (g[g.alpha == 0.0].delta > 0).mean(); r3 = (g[g.alpha == 3.0].delta > 0).mean()
    share[m] = (o, r3)
    print(f"  {m:<9}{o:>14.0%}{r3:>14.0%}")
c3 = share["twelve"][0] >= 0.30
print(f"\n  OLS nocivas medido em doze: {share['twelve'][0]:.0%} (barra >= 30%) -> "
      f"{'o dano PERSISTE' if c3 else 'nao informativo'}")

section("FECHAMENTO")
print(f"  2. curva n=40 em doze tipos: {'persiste' if c2 else 'nao informativo'}")
print(f"  3. celulas nocivas em doze tipos: {'persiste' if c3 else 'nao informativo'}")
print(f"  4. ridge a=3 em doze tipos: {share['twelve'][1]:.0%} nocivas")
print(f"  5. eixo naive/memoria: n=40 {tab.loc[40, 'nm']:+.1%}, "
      f"OLS nocivas {share['nm'][0]:.0%}, a=3 {share['nm'][1]:.0%}")
