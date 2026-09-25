#!/usr/bin/env python3
"""
Stage 24b — stage 24 rerun with 100 draws per size and fresh seeds.

A reader noted that the n = 40 median varied from +11.8% to +18.9% across three
independent sets of 30 draws, and asked for a pooled figure from more draws. This
is stage 24 unchanged except: 100 draws per fitting size (100 permutations at the
full cohort), new seeds, output to results/null_by_n_r100.csv. Its checks are
stage 24's, re-run.

Stage 24 — the permuted reference, measured at every n.

Stage 17 carried a reference line: the same transport, with composition
coefficients fitted on permuted composition. It came out at +1.8% from four
permutations, was restabilised at +0.5% from 180 for the paper figure, and on
that number three claims were built — "seven times worse than no information",
"twenty times worse", and "the damage comes from real coefficients estimated
badly, not from noise".

A reviewer pointed out that the reference was only ever measured at n = 656, and
was drawn as a horizontal line across a plot whose other curve ran from n = 40.
The comparison that produced every one of those claims set real coefficients at
n = 40 against permuted ones at n = 656. It confounds sample size with
information.

THE ALGEBRA, which is the reviewer's and was checked before this was written.
Let the true composition coefficient be beta. A correction b is applied in the
test cohort and what remains of the composition term is (beta - b)' S (beta - b),
with S the test cohort's composition covariance.
  - real fit:     b = beta + e,  after = e' S e,          delta = e'Se - beta'S beta
  - permuted fit: b = 0 + e_p,   after = beta'S beta + e_p' S e_p,
                                                          delta = e_p' S e_p
The permuted coefficients carry no information, so their error is the whole of
what they carry, with a noise variance at least as large as the real fit's
(composition signal becomes noise once the rows are shuffled). Both error terms
are sigma^2 times the transport index, and the index carries 1/n.

PRE-REGISTERED PREDICTIONS, and what follows from each
  P1  at the same n, the permuted reference does at least as much damage as the
      real fit's error term:  median delta_perm >= median after_real, at no
      fewer than eight of the ten fitting sizes.
      holds  -> "noise is harmless", and every ratio built on it, is retracted
                from the synthesis, the README and both manuscripts. The damage
                IS estimation noise, amplified by covariance mismatch.
  P2  the permuted reference falls with n about as 1/n: over the fitting sizes
      where its median is positive, the slope of log(delta_perm) on log(n) lies
      in [-1.3, -0.7].
  P3  at n = 40, permuted coefficients are MORE harmful than real ones:
      delta_perm > delta_real. The real fit subtracts the signal it correctly
      removes; the permuted one removes nothing. Reported.

CLOCKS. GSE40279 is the fitting cohort, so after stage 23 only Levine 2018 and
Horvath 2018 are clean. Permutation shuffles whole rows of the composition
matrix — every cell type moves together — so collinearity between types is
preserved and only the link to the clock is broken.

SANITY CHECKS
  1. the cache is the stage 18 cache. Hard stop.
  2. at n = 656 the permuted median must land inside stage 17's recomputed
     interquartile range (-0.8% to +2.1%): this stage has to reproduce the one
     point the old reference did measure before its other points mean anything.

Usage:  .venv/bin/python analysis/24_null_by_n.py
"""
import sys
from itertools import permutations, combinations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

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


section("CHECAGEM 1 — CACHE")
c1 = abs(panel.panel_r - 0.789) < 5e-3 and abs(panel.panel_mae - 0.027) < 5e-3
print(f"  -> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("  parando.")

print("  nulos de permutacao ...", flush=True)
CLEAN = [c for c in CLOCKS3 if allowed(c, "GSE40279")]
NULL, BEFORE = {}, {}
for t in COHORTS:
    for c in CLEAN:
        NULL[(t, c)] = null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        ib, base = inc_and_base(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        BEFORE[(t, c)] = ((ib - NULL[(t, c)]) / (1 - base), base)

GRID = [40, 60, 80, 120, 160, 184, 240, 320, 480, 656]
TESTS = ["GSE61151", "GSE50660", "GSE42861"]
REPS = 100
fit = D["GSE40279"]; NF = N["GSE40279"]
rows = []
for n in GRID:
    reps = 1 if n >= NF else REPS
    # the full cohort has one possible draw, so its spread would come from
    # permutations alone; give it the same count of permutations as the others
    perms = REPS if n >= NF else 1
    for r in range(reps):
        idx = (np.arange(NF) if n >= NF
               else stratified_draw(fit["chrono"], n, np.random.default_rng(7000 + 11 * r + n)))
        af, Cf = fit["chrono"][idx], fit["C12"][idx]
        for kind in ("real", "perm"):
            for pr in range(perms if kind == "perm" else 1):
                C_use = Cf if kind == "real" else Cf[np.random.default_rng(
                    150_000 + 101 * r + 7 * pr + n).permutation(len(Cf))]
                cbar = C_use.mean(axis=0)
                for c in CLEAN:
                    b = ridge_coefs(af, C_use, fit["y"][c][idx], [0.0])[0.0]
                    for dst in TESTS:
                        d = D[dst]; before, base = BEFORE[(dst, c)]
                        yc = d["y"][c] - (d["C12"][:, :-1] - cbar[:-1]) @ b
                        ia, _ = inc_and_base(yc, d["chrono"], d["C6"])
                        after = (ia - NULL[(dst, c)]) / (1 - base)
                        rows.append(dict(n=n, rep=r, perm=pr, kind=kind, clock=c,
                                         dst=dst, before=before, after=after,
                                         delta=after - before))
    print(f"  n = {n:>4}: pronto", flush=True)
R = pd.DataFrame(rows)
R.to_csv(RES / "null_by_n_r100.csv", index=False)

real = R[R.kind == "real"]; perm = R[R.kind == "perm"]
tab = pd.DataFrame({
    "delta_real": real.groupby("n").delta.median(),
    "after_real": real.groupby("n").after.median(),
    "delta_perm": perm.groupby("n").delta.median(),
    "perm_q25": perm.groupby("n").delta.quantile(.25),
    "perm_q75": perm.groupby("n").delta.quantile(.75),
})

section("A REFERENCIA EM CADA N")
print(f"  {'n':>5}{'delta real':>12}{'erro do real':>14}{'delta embaralhado':>20}{'IQR':>20}")
for n, r in tab.iterrows():
    print(f"  {n:>5}{r.delta_real:>+12.1%}{r.after_real:>+14.1%}{r.delta_perm:>+20.1%}"
          f"   [{r.perm_q25:+.1%}, {r.perm_q75:+.1%}]")

section("CHECAGEM 2 — REPRODUZ O PONTO QUE A REFERENCIA ANTIGA MEDIU?")
p656 = tab.loc[656, "delta_perm"]
c2 = -0.008 <= p656 <= 0.021
print(f"  n=656: {p656:+.2%}, faixa antiga -0,8% a +2,1% -> {'ok' if c2 else 'FALHOU'}")

section("P1 — NO MESMO N, EMBARALHADO >= ERRO DO REAL?")
ok = (tab.delta_perm >= tab.after_real)
for n, v in ok.items():
    print(f"  n = {n:>4}: {tab.loc[n,'delta_perm']:+.1%} vs {tab.loc[n,'after_real']:+.1%}  "
          f"{'sim' if v else 'NAO'}")
p1 = int(ok.sum()) >= 8
print(f"\n  {int(ok.sum())} de {len(ok)} (barra: >= 8) -> {'CONFIRMADA' if p1 else 'REFUTADA'}")

section("P2 — O EMBARALHADO CAI COMO 1/N?")
pos = tab[tab.delta_perm > 0]
if len(pos) >= 4:
    slope = float(np.polyfit(np.log(pos.index.to_numpy(float)), np.log(pos.delta_perm), 1)[0])
    p2 = -1.3 <= slope <= -0.7
    print(f"  inclinacao log-log sobre {len(pos)} tamanhos com mediana positiva: {slope:.2f}")
    print(f"  faixa [-1.3, -0.7] -> {'CONFIRMADA' if p2 else 'REFUTADA'}")
else:
    p2 = False
    print(f"  so {len(pos)} tamanhos com mediana positiva; inclinacao nao estimavel -> REFUTADA")

section("P3 — EM N=40, EMBARALHADO PIOR QUE REAL? (REPORTADO)")
print(f"  real {tab.loc[40,'delta_real']:+.1%}  |  embaralhado {tab.loc[40,'delta_perm']:+.1%}"
      f"  -> {'sim' if tab.loc[40,'delta_perm'] > tab.loc[40,'delta_real'] else 'nao'}")

section("FECHAMENTO")
print(f"  1. cache: {'ok' if c1 else 'FALHOU'}")
print(f"  2. reproduz n=656: {'ok' if c2 else 'FALHOU'}")
print(f"  P1 embaralhado >= erro do real: {'CONFIRMADA' if p1 else 'REFUTADA'}")
print(f"  P2 queda como 1/n: {'CONFIRMADA' if p2 else 'REFUTADA'}")
