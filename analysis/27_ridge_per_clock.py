#!/usr/bin/env python3
"""
Stage 27 — does the penalty hold clock by clock, or was it averaged too?

Stage 26 withdrew the transport index's safety floor: counted per clock instead
of as a median over clocks, 12 of 40 configurations below the floor were harmful,
because a harmed clock had been outvoted by a helped one. Stage 23's claim that
ridge at alpha = 3 holds all twelve directed pairs at or below zero was counted
the same way — a median over the clocks in each pair. It has to be recounted.

WHAT SHRINKAGE CAN AND CANNOT PROMISE. Under covariate shift alone, penalising
bounds the amplification of estimation noise. Under model shift — beta differs
between cohorts, stage 25 — the fully shrunk correction, b = 0, does exactly
nothing, so its net damage is zero by construction; any penalty moves the
correction toward that point. What a finite penalty cannot promise is to get
there: a cell whose specification error is large may still be harmful at
alpha = 3, and only more shrinkage — giving up more benefit — would fix it.
Patil, Du & Tibshirani (ICML 2024) go further: under shift the optimal penalty
can even be negative. Positive shrinkage is a safe default, not an optimum.

PRE-REGISTERED
  1. the cache is the stage 18 cache. Hard stop.
  2. PRIMARY: at matched n, ridge alpha = 3 leaves the median net damage at or
     below zero in EVERY (directed pair x clock) cell, clocks clean for the pair
     (Levine 2018 and Horvath 2018 everywhere; Horvath 2013 where GSE40279 is
     absent). If any cell fails, the claim becomes "reduces", not "holds".
  3. failing cells, if any, are expected where stage 25's specification term is
     large: reported with it.
  4. the share of harmful cells under OLS and under alpha = 3, per clock.
  5. whether a larger penalty (alpha = 10) holds every cell, and what it costs
     in benefit where OLS was already fine. Reported.

Usage:  .venv/bin/python analysis/27_ridge_per_clock.py
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


print("  nulos ...", flush=True)
NULL, BEFORE = {}, {}
for t in COHORTS:
    for c in CLOCKS3:
        NULL[(t, c)] = null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        ib, base = inc_and_base(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        BEFORE[(t, c)] = ((ib - NULL[(t, c)]) / (1 - base), base)

def transport(src, dst, n, reps, alphas=(0.0,)):
    s, d = D[src], D[dst]
    rows = []
    for r in range(reps):
        idx = (np.arange(N[src]) if n >= N[src]
               else stratified_draw(s["chrono"], n, np.random.default_rng(700 + 7 * r + n)))
        af, Cf = s["chrono"][idx], s["C12"][idx]
        cbar = Cf.mean(axis=0)
        ti = transport_index(af, Cf, d["chrono"], d["C12"])
        for c in CLOCKS3:
            if not allowed(c, src, dst):
                continue
            coefs = ridge_coefs(af, Cf, s["y"][c][idx], list(alphas))
            before, base = BEFORE[(dst, c)]
            for alpha, b in coefs.items():
                yc = d["y"][c] - (d["C12"][:, :-1] - cbar[:-1]) @ b
                ia, _ = inc_and_base(yc, d["chrono"], d["C6"])
                rows.append(dict(src=src, dst=dst, n=len(idx), clock=c, alpha=alpha,
                                 tindex=ti,
                                 delta=(ia - NULL[(dst, c)]) / (1 - base) - before))
    return rows



ALPHAS = (0.0, 3.0, 10.0)
rows = []
for src, dst in permutations(COHORTS, 2):
    n = min(N[src], N[dst])
    rows += transport(src, dst, n, 20, alphas=ALPHAS)
    print(f"  {src[3:]} -> {dst[3:]} (n={n}): pronto", flush=True)
R = pd.DataFrame(rows)
R.to_csv(RES / "ridge_per_clock.csv", index=False)
cell = (R.groupby(["src", "dst", "clock", "alpha"]).delta.median()
          .unstack("alpha").reset_index())
spec = pd.read_csv(RES / "beta_heterogeneity.csv")
SPEC = {(r.src, r.dst, r.clock): r.spec for r in spec.itertuples()}
cell["spec"] = [SPEC.get((r.src, r.dst, r.clock), np.nan) for r in cell.itertuples()]

section("CHECAGEM 2 — ALPHA = 3 SEGURA TODA CELULA PAR x RELOGIO?")
print(f"  {'par':<18}{'relogio':<13}{'OLS':>9}{'a=3':>9}{'a=10':>9}{'espec.':>9}")
for r in cell.sort_values(3.0, ascending=False).itertuples():
    print(f"  {r.src[3:]+' -> '+r.dst[3:]:<18}{r.clock:<13}{r._4:>+9.1%}{r._5:>+9.1%}"
          f"{r._6:>+9.1%}{r.spec:>+9.1%}" + ("   <- nocivo com a=3" if r._5 > 0 else ""))
fails = cell[cell[3.0] > 0]
c2 = len(fails) == 0
print(f"\n  {len(cell) - len(fails)} de {len(cell)} celulas em zero ou abaixo com a=3 -> "
      f"{'SEGURA' if c2 else 'NAO SEGURA TODAS'}")

section("CHECAGEM 3 — AS FALHAS MORAM ONDE O BETA DIFERE?")
if len(fails):
    print(f"  termo de especificacao mediano: falhas {fails.spec.median():+.1%} | "
          f"resto {cell[cell[3.0] <= 0].spec.median():+.1%}")
else:
    print("  nenhuma falha a explicar")

section("CHECAGEM 4 — FRACAO DE CELULAS NOCIVAS, POR RELOGIO")
print(f"  {'relogio':<13}{'celulas':>9}{'OLS':>9}{'a=3':>9}{'a=10':>9}")
for c, g in cell.groupby("clock"):
    print(f"  {c:<13}{len(g):>9}{(g[0.0] > 0).mean():>9.0%}{(g[3.0] > 0).mean():>9.0%}"
          f"{(g[10.0] > 0).mean():>9.0%}")
print(f"  {'todas':<13}{len(cell):>9}{(cell[0.0] > 0).mean():>9.0%}"
      f"{(cell[3.0] > 0).mean():>9.0%}{(cell[10.0] > 0).mean():>9.0%}")

section("CHECAGEM 5 — ALPHA = 10 SEGURA TUDO? E QUANTO CUSTA?")
f10 = cell[cell[10.0] > 0]
print(f"  celulas nocivas com a=10: {len(f10)} de {len(cell)}")
fine = cell[cell[0.0] < 0]
print(f"  onde o OLS ja ajudava ({len(fine)} celulas), beneficio mediano:")
print(f"    OLS {fine[0.0].median():+.1%} | a=3 {fine[3.0].median():+.1%} | "
      f"a=10 {fine[10.0].median():+.1%}")

section("FECHAMENTO")
print(f"  2. alpha=3 segura toda celula: {'sim' if c2 else f'nao ({len(fails)} falham)'}")
