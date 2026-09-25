#!/usr/bin/env python3
"""
Stage 39b — why stage 38's saliva transports failed: technical shift or collinearity?

Stage 38: between saliva cohorts, 7 of 8 (pair x clock) cells were harmful at
matched n, unpenalised and at alpha = 3 alike; the correction worked at home.
Two readings were written down there, and this stage tests them, on the same
8 cells (4 directed pairs with GSE78874, Levine 2018 and Horvath 2018):

  TECHNICAL — every pair crossed array and preprocessing. 39a re-derived GSE78874
    after quantile-normalising it to GSE232891's beta distribution.
  COLLINEARITY — the nine-type saliva fit has condition numbers near 1,000, and
    ridge leaves the dominant direction nearly untouched. A three-type fit
    (epithelium, fibroblast, immune) is far less collinear.

PRE-REGISTERED (written before any of the variants was scored)
  1. the 38a and 39a caches exist. Hard stop.
  2. TECHNICAL: with GSE78874 normalised, at most 3 of 8 cells are harmful
     unpenalised (38b: 7).
  3. COLLINEARITY: with the three-type fit, at most 1 of 8 cells is harmful at
     alpha = 3 (38b: 7).
  Both, one or neither may hold. The fourth combination (normalised and three-type)
  is reported without a bar.

Usage:  .venv/bin/python analysis/39b_saliva_why.py
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




FIT9 = ["Epi", "Fib", "B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino"]
FIT3 = ["Epi", "Fib", "IC"]
MEAS = ["Epi", "Fib", "IC"]
AGE = ["Levine2018", "Horvath2018"]
need = [CACHE / f"{t}_{w}.csv" for t in ("GSE232891", "GSE232332", "GSE78874", "GSE78874qn")
        for w in ("compfit", "compmeas", "ages")]
section("CHECAGEM 1 — CACHES")
if not all(f.exists() for f in need):
    sys.exit("  parando: cache da 38a/39a ausente.")
print("  ok")


def load(tag, fit_types):
    a = pd.read_csv(CACHE / f"{tag}_ages.csv", index_col=0)
    fit = pd.read_csv(CACHE / f"{tag}_compfit.csv", index_col=0)
    meas = pd.read_csv(CACHE / f"{tag}_compmeas.csv", index_col=0).loc[a.index, MEAS]
    if fit_types == FIT3:
        C12 = meas.to_numpy()
    else:
        C12 = fit.loc[a.index, FIT9].to_numpy()
    return dict(chrono=a.chrono.to_numpy(float), y={c: a[c].to_numpy(float) for c in AGE},
                C12=C12, C6=meas.to_numpy())


def run(variant, p_tag, fit_types):
    D = {t: load(t, fit_types) for t in ("GSE232891", "GSE232332")}
    D["P"] = load(p_tag, fit_types)
    N = {t: len(D[t]["chrono"]) for t in D}
    NUL, BEF = {}, {}
    for t in D:
        for c in AGE:
            NUL[(t, c)] = null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["C6"], n_perm=400)
            ib, base = inc_and_base(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
            BEF[(t, c)] = ((ib - NUL[(t, c)]) / (1 - base), base)
    out = []
    for src, dst in [("GSE232891", "P"), ("P", "GSE232891"), ("GSE232332", "P"), ("P", "GSE232332")]:
        n = min(N[src], N[dst]); s = D[src]
        for c in AGE:
            res = {"ols": [], "a3": []}
            for r in range(20):
                idx = (np.arange(N[src]) if n >= N[src]
                       else stratified_draw(s["chrono"], n, np.random.default_rng(3850 + 7 * r)))
                co = ridge_coefs(s["chrono"][idx], s["C12"][idx], s["y"][c][idx], [0.0, 3.0])
                cb = s["C12"][idx].mean(0); d = D[dst]
                for key, al in (("ols", 0.0), ("a3", 3.0)):
                    yc = d["y"][c] - (d["C12"][:, :-1] - cb[:-1]) @ co[al]
                    ia, _ = inc_and_base(yc, d["chrono"], d["C6"])
                    res[key].append((ia - NUL[(dst, c)]) / (1 - BEF[(dst, c)][1]) - BEF[(dst, c)][0])
                if n >= N[src]:
                    break
            out.append(dict(variant=variant, src=src.replace("P", p_tag) if src == "P" else src,
                            dst=p_tag if dst == "P" else dst, clock=c, before=BEF[(dst, c)][0],
                            ols=float(np.median(res["ols"])), a3=float(np.median(res["a3"]))))
    Z = partial_out(D["P"]["chrono"], D["P"]["C12"][:, :-1]); Z = (Z - Z.mean(0)) / Z.std(0)
    print(f"  {variant}: numero de condicao no GSE78874 {np.linalg.cond(Z):.0f}", flush=True)
    return out


rows = []
rows += run("bruto, 9 tipos", "GSE78874", FIT9)
rows += run("normalizado, 9 tipos", "GSE78874qn", FIT9)
rows += run("bruto, 3 tipos", "GSE78874", FIT3)
rows += run("normalizado, 3 tipos", "GSE78874qn", FIT3)
R = pd.DataFrame(rows); R.to_csv(RES / "saliva_why.csv", index=False)

section("RESULTADO POR VARIANTE")
for v, g in R.groupby("variant", sort=False):
    print(f"\n  {v}: nocivas OLS {(g.ols > 0).sum()} de 8, alpha=3 {(g.a3 > 0).sum()} de 8")
    for r in g.itertuples():
        print(f"    {r.src[3:]+' -> '+r.dst[3:]:<20}{r.clock:<13}antes {r.before:+7.1%}  OLS {r.ols:+8.1%}  a=3 {r.a3:+7.1%}")
g2 = R[R.variant == "normalizado, 9 tipos"]; g3 = R[R.variant == "bruto, 3 tipos"]
c2 = int((g2.ols > 0).sum()) <= 3
c3 = int((g3.a3 > 0).sum()) <= 1

section("FECHAMENTO")
print(f"  2. tecnico (normalizado, OLS <= 3 de 8): {(g2.ols > 0).sum()} -> {'ok' if c2 else 'FALHOU'}")
print(f"  3. colinearidade (3 tipos, alpha=3 <= 1 de 8): {(g3.a3 > 0).sum()} -> {'ok' if c3 else 'FALHOU'}")
