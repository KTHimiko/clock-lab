#!/usr/bin/env python3
"""
Stage 44b — the central blood results, on the published reference libraries.

The panels used since stage 13 were built here. 44a recomputed the composition of
all six blood cohorts with the published alternatives: the Salas et al. 2022
twelve-type reference on its IDOL-optimised 450k probe set for the fit, and the
published seven-type blood reference for the measurement. They share 20 probes of
600 and 333, where the home-built pair shared its whole construction, so this
also removes the "measurement shares the fit's panel" caveat. The two twelve-type
panels agree on what they both estimate (neutrophils r = 0.99, naive CD8 r = 0.90,
median over cohorts).

Everything else is unchanged: the same cohorts, clocks, exclusions, subsampling
and metric. If the results move, the panel was doing the work.

PRE-REGISTERED (written before any published-panel transport was computed)
  1. the 44a cache exists for all six cohorts. Hard stop.
  2. SMALL-N HARM REPLICATES: fitted on 40 samples of GSE40279 and transported to
     the three original test cohorts, the correction is harmful in more than 60%
     of draws. (Home-built panels: 88%.)
  3. THE PENALTY REPLICATES: over all 30 directed pairs at matched n, alpha = 3
     leaves at most 10% of (pair x clock) cells harmful, where unpenalised leaves
     at least 25%. (Home-built: 23 of 72 and 1 of 72.)
  4. THE FLOOR REPLICATES: fitted on all of GSE55763, the median composition left
     across cells is positive.

Usage:  .venv/bin/python analysis/44b_published_panels_results.py
"""
import sys
from itertools import permutations
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
RES = ROOT / "results"; CACHE = RES / "cache"
sys.path.insert(0, str(ROOT))
from load_extended import TYPES12
RNG = np.random.default_rng(20260924)
N_REPS = 20
ALPHA_FIX = 3.0

# the corrected rule: a clock is out of any pair in which a cohort it trained on
# appears, whichever end of the transport that cohort is on






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




BLOOD = ["GSE40279", "GSE61151", "GSE50660", "GSE42861", "GSE132203", "GSE55763"]
AGE = ["Levine2018", "Horvath2018"]
T7 = ["B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino"]
TRAINED_ON = {"Hannum2013": {"GSE40279"}, "Horvath2013": {"GSE40279"}}
USABLE = {t: pd.read_csv(CACHE / f"{t}_usable_clocks.csv").iloc[:, 0].tolist()
          for t in ("GSE132203", "GSE55763")}


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


section("CHECAGEM 1 — CACHE PUBLICADO")
D, N = {}, {}
for t in BLOOD:
    f12, f7 = CACHE / f"{t}_comp12pub.csv", CACHE / f"{t}_comp7pub.csv"
    if not (f12.exists() and f7.exists()):
        sys.exit(f"  parando: falta o cache da 44a para {t}.")
    a = pd.read_csv(CACHE / f"{t}_ages.csv", index_col=0); a = a[a.chrono.notna()]
    D[t] = dict(chrono=a.chrono.to_numpy(float), y={c: a[c].to_numpy(float) for c in AGE},
                C=pd.read_csv(f12, index_col=0).loc[a.index, TYPES12].to_numpy(),
                M=pd.read_csv(f7, index_col=0).loc[a.index, T7].to_numpy())
    N[t] = len(a)
print("  " + ", ".join(f"{t} {N[t]}" for t in BLOOD))

NUL, BEF = {}, {}
for t in BLOOD:
    for c in AGE:
        NUL[(t, c)] = null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["M"], n_perm=400)
        ib, base = inc_and_base(D[t]["y"][c], D[t]["chrono"], D[t]["M"])
        BEF[(t, c)] = ((ib - NUL[(t, c)]) / (1 - base), base)


def after(dst, c, b, cbar):
    d = D[dst]
    yc = d["y"][c] - (d["C"][:, :-1] - cbar[:-1]) @ b
    ia, _ = inc_and_base(yc, d["chrono"], d["M"])
    return (ia - NUL[(dst, c)]) / (1 - BEF[(dst, c)][1])


def allowed(c, *co):
    return not (TRAINED_ON.get(c, set()) & set(co)) and not any(
        t in USABLE and c not in USABLE[t] for t in co)


section("CHECAGEM 2 — DANO EM N=40, AJUSTANDO NO GSE40279")
s = D["GSE40279"]; d40 = []
for r in range(100):
    idx = stratified_draw(s["chrono"], 40, np.random.default_rng(4400 + 11 * r))
    for c in AGE:
        b = ridge_coefs(s["chrono"][idx], s["C"][idx], s["y"][c][idx], [0.0])[0.0]
        cb = s["C"][idx].mean(0)
        for dst in ("GSE61151", "GSE50660", "GSE42861"):
            d40.append(after(dst, c, b, cb) - BEF[(dst, c)][0])
d40 = np.array(d40)
c2 = (d40 > 0).mean() > 0.60
print(f"  {len(d40)} sorteios: nocivo em {(d40 > 0).mean():.0%}, mediana {np.median(d40):+.1%}"
      f" (barra > 60%) -> {'ok' if c2 else 'FALHOU'}")

section("CHECAGEM 3 — A PENALIDADE EM TODOS OS PARES")
cells = []
for src, dst in permutations(BLOOD, 2):
    n = min(N[src], N[dst]); s = D[src]
    for c in AGE:
        if not allowed(c, src, dst):
            continue
        res = {0.0: [], 3.0: []}
        for r in range(20):
            idx = (np.arange(N[src]) if n >= N[src]
                   else stratified_draw(s["chrono"], n, np.random.default_rng(4450 + 7 * r + n)))
            co = ridge_coefs(s["chrono"][idx], s["C"][idx], s["y"][c][idx], [0.0, 3.0])
            cb = s["C"][idx].mean(0)
            for al in (0.0, 3.0):
                res[al].append(after(dst, c, co[al], cb) - BEF[(dst, c)][0])
            if n >= N[src]:
                break
        cells.append(dict(src=src, dst=dst, n=n, clock=c,
                          ols=float(np.median(res[0.0])), a3=float(np.median(res[3.0]))))
P = pd.DataFrame(cells); P.to_csv(RES / "published_panel_cells.csv", index=False)
k0, k3 = int((P.ols > 0).sum()), int((P.a3 > 0).sum())
c3 = k3 <= 0.10 * len(P) and k0 >= 0.25 * len(P)
print(f"  {len(P)} celulas: nocivas sem penalidade {k0}, com alpha=3 {k3}")
print(f"  medianas: sem penalidade {P.ols.median():+.1%}, alpha=3 {P.a3.median():+.1%}"
      f" -> {'ok' if c3 else 'FALHOU'}")

section("CHECAGEM 4 — O PISO AJUSTANDO NO GSE55763 INTEIRO")
s = D["GSE55763"]; left = []
for dst in BLOOD:
    if dst == "GSE55763":
        continue
    for c in AGE:
        if not allowed(c, "GSE55763", dst):
            continue
        b = ridge_coefs(s["chrono"], s["C"], s["y"][c], [0.0])[0.0]
        left.append(after(dst, c, b, s["C"].mean(0)))
left = np.array(left)
c4 = np.median(left) > 0
print(f"  {len(left)} celulas: composicao restante mediana {np.median(left):+.1%}"
      f" [{left.min():+.1%}, {left.max():+.1%}] -> {'ok' if c4 else 'FALHOU'}")
section("FECHAMENTO")
for i, (lab, ok) in enumerate([("dano em n=40", c2), ("penalidade", c3), ("piso", c4)], start=2):
    print(f"  {i}. {lab}: {'ok' if ok else 'FALHOU'}")
