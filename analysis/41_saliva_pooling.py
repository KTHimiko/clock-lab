#!/usr/bin/env python3
"""
Stage 41 — does pooling saliva studies, as Galkin et al. (2021) did, rescue transport?

Stages 38-40: a composition correction carried from one saliva cohort to another
was harmful in 7 of 8 cells, penalised or not, because the clocks' slope on the
immune fraction differs, even in sign, between cohorts. The one published
transport in saliva pooled eight studies before fitting. If model shift varies
around a common value, pooling averages it out; if the target is on the other
side of the split, the pooled coefficient still has the wrong sign.

Leave-one-cohort-out over four saliva cohorts:
  target GSE78874   <- pool GSE232891 + GSE232332 + GSE149747
  target GSE149747  <- pool GSE232891 + GSE232332 + GSE78874
  target GSE232891  <- pool GSE78874 + GSE149747   (GSE232332 left out: may share people)
  target GSE232332  <- pool GSE78874 + GSE149747   (GSE232891 left out, same reason)
GSE149747 enters with its 44 baseline samples. The fit uses the three-type
EpiDISH composition (epithelium, fibroblast, immune; condition number ~1, stage
39) with a fixed intercept per study; the measurement uses the same three types,
so it is biased toward the correction. Levine 2018 and Horvath 2018.

PRE-REGISTERED (written before any pooled fit was computed)
  1. the caches exist and the four cohorts load. Hard stop.
  2. POOLING RESCUES: unpenalised, at most 2 of the 8 (target x clock) cells are
     harmful. The single-source three-type transports of stage 39 left 7 of 8.
  3. POOLING + PENALTY: at alpha = 3, at most 1 of 8.
  Reported: the same without study intercepts, and the pooled immune slope next
  to each target's own.

Usage:  .venv/bin/python analysis/41_saliva_pooling.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"; CACHE = RES / "cache"
RNG = np.random.default_rng(20260924)
N_REPS = 20
ALPHA_FIX = 3.0

# the corrected rule: a clock is out of any pair in which a cohort it trained on
# appears, whichever end of the transport that cohort is on




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




AGE = ["Levine2018", "Horvath2018"]
MEAS = ["Epi", "Fib", "IC"]
section("CHECAGEM 1 — COORTES")
D = {}
try:
    for t in ("GSE232891", "GSE232332", "GSE78874"):
        a = pd.read_csv(CACHE / f"{t}_ages.csv", index_col=0)
        C = pd.read_csv(CACHE / f"{t}_compmeas.csv", index_col=0).loc[a.index, MEAS].to_numpy()
        D[t] = dict(chrono=a.chrono.to_numpy(float), y={c: a[c].to_numpy(float) for c in AGE}, C=C)
    g = pd.read_csv(CACHE / "GSE149747_saliva.csv", index_col=0)
    g = g[g.tp == "Before"].dropna(subset=["age"])
    D["GSE149747"] = dict(chrono=g.age.to_numpy(float), y={c: g[c].to_numpy(float) for c in AGE},
                          C=g[MEAS].to_numpy(float))
except FileNotFoundError as e:
    sys.exit(f"  parando: {e}")
print("  " + ", ".join(f"{t} {len(d['chrono'])}" for t, d in D.items()))

PLAN = {"GSE78874": ["GSE232891", "GSE232332", "GSE149747"],
        "GSE149747": ["GSE232891", "GSE232332", "GSE78874"],
        "GSE232891": ["GSE78874", "GSE149747"],
        "GSE232332": ["GSE78874", "GSE149747"]}


def pooled_fit(srcs, c, alpha, study_fe=True):
    age = np.concatenate([D[s]["chrono"] for s in srcs])
    C = np.vstack([D[s]["C"][:, :-1] for s in srcs])
    y = np.concatenate([D[s]["y"][c] for s in srcs])
    cols = [np.ones(len(age)), age]
    if study_fe:
        k = 0
        for s in srcs[1:]:
            k0 = sum(len(D[x]["chrono"]) for x in srcs[:srcs.index(s)])
            z = np.zeros(len(age)); z[k0:k0 + len(D[s]["chrono"])] = 1; cols.append(z)
    X = np.column_stack(cols)
    H = X @ np.linalg.pinv(X)
    Ct = C - H @ C; yt = y - H @ y
    mu, sd = Ct.mean(0), Ct.std(0) + 1e-12
    Z = (Ct - mu) / sd
    U, sv, Vt = np.linalg.svd(Z, full_matrices=False)
    scale = float((sv ** 2).mean())
    return (Vt.T @ (sv * (U.T @ yt) / (sv ** 2 + alpha * scale))) / sd


def own_slope(t, c):
    d = D[t]; X = np.column_stack([np.ones(len(d["chrono"])), d["chrono"], d["C"][:, 2]])
    return np.linalg.lstsq(X, d["y"][c], rcond=None)[0][2] / 10


def pooled_ic_slope(srcs, c):
    rows = []
    for fe in (True,):
        age = np.concatenate([D[s]["chrono"] for s in srcs]); ic = np.concatenate([D[s]["C"][:, 2] for s in srcs])
        y = np.concatenate([D[s]["y"][c] for s in srcs])
        cols = [np.ones(len(age)), age, ic]
        for s in srcs[1:]:
            k0 = sum(len(D[x]["chrono"]) for x in srcs[:srcs.index(s)])
            z = np.zeros(len(age)); z[k0:k0 + len(D[s]["chrono"])] = 1; cols.append(z)
        rows.append(np.linalg.lstsq(np.column_stack(cols), y, rcond=None)[0][2] / 10)
    return rows[0]


section("RESULTADO")
out = []
for tgt, srcs in PLAN.items():
    d = D[tgt]
    for c in AGE:
        nul = null_mean(d["y"][c], d["chrono"], d["C"], n_perm=500)
        ib, base = inc_and_base(d["y"][c], d["chrono"], d["C"])
        before = (ib - nul) / (1 - base)
        res = {}
        for key, al, fe in (("ols", 0.0, True), ("a3", 3.0, True), ("ols_nofe", 0.0, False)):
            b = pooled_fit(srcs, c, al, fe)
            ia, _ = inc_and_base(d["y"][c] - d["C"][:, :-1] @ b, d["chrono"], d["C"])
            res[key] = (ia - nul) / (1 - base) - before
        out.append(dict(target=tgt, pool="+".join(s[3:] for s in srcs), clock=c, before=before,
                        own=own_slope(tgt, c), pooled=pooled_ic_slope(srcs, c), **res))
R = pd.DataFrame(out); R.to_csv(RES / "saliva_pooling.csv", index=False)
print(f"  {'alvo':<11}{'relogio':<13}{'antes':>8}{'OLS':>9}{'a=3':>8}{'sem EF':>9}{'incl. propria':>15}{'incl. pool':>12}")
for r in R.itertuples():
    print(f"  {r.target:<11}{r.clock:<13}{r.before:>+8.1%}{r.ols:>+9.1%}{r.a3:>+8.1%}{r.ols_nofe:>+9.1%}"
          f"{r.own:>+15.2f}{r.pooled:>+12.2f}")
k0, k3, kn = int((R.ols > 0).sum()), int((R.a3 > 0).sum()), int((R.ols_nofe > 0).sum())
c2, c3 = k0 <= 2, k3 <= 1
print(f"\n  nocivas: OLS {k0} de 8, alpha=3 {k3} de 8, sem efeito de estudo {kn} de 8")
section("FECHAMENTO")
print(f"  2. pooling resgata (OLS <= 2 de 8): {k0} -> {'ok' if c2 else 'FALHOU'}")
print(f"  3. pooling + penalidade (alpha=3 <= 1 de 8): {k3} -> {'ok' if c3 else 'FALHOU'}")
