#!/usr/bin/env python3
"""
Stage 38b — does the transport result hold in saliva?

Saliva's composition axis — buccal epithelium against leukocytes — is far larger
than anything in blood, and the one published transport of composition
coefficients (Galkin et al. 2021) is in saliva. 38a cached three adult cohorts:
GSE232891 (EPIC, 552), GSE232332 (EPIC, 265) and GSE78874 (450k, 259). The first
two come from one group and could not be checked for shared people (no rs
probes), so they are never paired: 4 directed pairs, each with GSE78874.

Fit with nine types (HEpiDISH: Epi, Fib and seven immune subtypes); measure with
three (Epi, Fib, IC). The measurement shares the fit's first step, so it is biased
toward the correction: any harm found is a floor. Clocks: Levine 2018 and
Horvath 2018 (usable in all three); DunedinPACE reported, no bar (not validated in
saliva and flat with age here). Counted per (pair x clock) cell.

PRE-REGISTERED (written before any saliva transport was computed)
  1. the 38a cache exists. Hard stop.
  2. SMALL-N HARM: fitted on 40 age-stratified samples and transported, over the
     4 directed pairs, the age clocks' correction is harmful in more than 60% of
     draws.
  3. the shuffled reference at n = 40 is positive.
  4. THE PENALTY: at matched n (259), alpha = 3 leaves at most 1 of the 8 cells
     harmful; the fading penalty (3 x 40 / n) and unpenalised counts reported.
  Reported without a bar: the composition signal before correction (how large
  the saliva axis is for each clock), and error left at matched n.

Usage:  .venv/bin/python analysis/38b_saliva_transport.py
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




SAL = ["GSE232891", "GSE232332", "GSE78874"]
FIT_T = ["Epi", "Fib", "B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino"]
MEAS_T = ["Epi", "Fib", "IC"]
section("CHECAGEM 1 — CACHE DA 38a")
if not all((CACHE / f"{t}_{w}.csv").exists() for t in SAL for w in ("compfit", "compmeas", "ages", "dpace", "usable_clocks")):
    sys.exit("  parando: cache da 38a ausente.")
D, N = {}, {}
for t in SAL:
    a = pd.read_csv(CACHE / f"{t}_ages.csv", index_col=0)
    D[t] = dict(chrono=a.chrono.to_numpy(float),
                y={c: a[c].to_numpy(float) for c in ["Levine2018", "Horvath2018"]},
                C12=pd.read_csv(CACHE / f"{t}_compfit.csv", index_col=0).loc[a.index, FIT_T].to_numpy(),
                C6=pd.read_csv(CACHE / f"{t}_compmeas.csv", index_col=0).loc[a.index, MEAS_T].to_numpy())
    D[t]["y"]["DunedinPACE"] = pd.read_csv(CACHE / f"{t}_dpace.csv", index_col=0).iloc[:, 0].reindex(a.index).to_numpy(float)
    N[t] = len(a)
print("  ok:", N)
AGE = ["Levine2018", "Horvath2018"]
PAIRS = [("GSE232891", "GSE78874"), ("GSE78874", "GSE232891"),
         ("GSE232332", "GSE78874"), ("GSE78874", "GSE232332")]

NUL, BEF = {}, {}
for t in SAL:
    for c in AGE + ["DunedinPACE"]:
        NUL[(t, c)] = null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["C6"], n_perm=500)
        ib, base = inc_and_base(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        BEF[(t, c)] = ((ib - NUL[(t, c)]) / (1 - base), base)
print("\n  sinal de composicao antes da correcao (fracao da variancia da aceleracao):")
for t in SAL:
    print(f"    {t:<10} " + "  ".join(f"{c} {BEF[(t, c)][0]:+.1%}" for c in AGE + ["DunedinPACE"]))


def after(dst, c, b, cbar):
    d = D[dst]
    yc = d["y"][c] - (d["C12"][:, :-1] - cbar[:-1]) @ b
    ia, _ = inc_and_base(yc, d["chrono"], d["C6"])
    return (ia - NUL[(dst, c)]) / (1 - BEF[(dst, c)][1])


section("CHECAGEM 2 E 3 — N = 40")
rows = []
for src, dst in PAIRS:
    s = D[src]
    for r in range(30):
        idx = stratified_draw(s["chrono"], 40, np.random.default_rng(3800 + 11 * r + len(src) + len(dst) * 3))
        Cf = s["C12"][idx]; Cp = Cf[np.random.default_rng(83_000 + r).permutation(40)]
        for c in AGE + ["DunedinPACE"]:
            for kind, Cu in (("real", Cf), ("perm", Cp)):
                b = ridge_coefs(s["chrono"][idx], Cu, s["y"][c][idx], [0.0])[0.0]
                rows.append(dict(src=src, dst=dst, rep=r, clock=c, kind=kind,
                                 delta=after(dst, c, b, Cu.mean(axis=0)) - BEF[(dst, c)][0]))
C = pd.DataFrame(rows); C.to_csv(RES / "saliva_curve40.csv", index=False)
age = C[(C.kind == "real") & C.clock.isin(AGE)].delta
perm = C[(C.kind == "perm") & C.clock.isin(AGE)].delta
c2 = (age > 0).mean() > 0.60; c3 = perm.median() > 0
print(f"  2. relogios de idade: nocivo em {(age > 0).mean():.0%}, mediana {age.median():+.1%} (barra > 60%) -> {'ok' if c2 else 'FALHOU'}")
print(f"  3. embaralhado: mediana {perm.median():+.1%} (barra > 0) -> {'ok' if c3 else 'FALHOU'}")
for (src, dst), g in C[C.kind == "real"].groupby(["src", "dst"]):
    print(f"     {src[3:]} -> {dst[3:]}: " + "  ".join(
        f"{c} {g[g.clock == c].delta.median():+.1%} ({(g[g.clock == c].delta > 0).mean():.0%})" for c in AGE + ["DunedinPACE"]))

section("CHECAGEM 4 — PENALIDADE EM N PAREADO")
cells = []
for src, dst in PAIRS:
    n = min(N[src], N[dst]); s = D[src]; af = 3.0 * 40 / n
    for c in AGE + ["DunedinPACE"]:
        res = {"ols": [], "a3": [], "fade": [], "left": []}
        for r in range(20):
            idx = (np.arange(N[src]) if n >= N[src]
                   else stratified_draw(s["chrono"], n, np.random.default_rng(3850 + 7 * r)))
            co = ridge_coefs(s["chrono"][idx], s["C12"][idx], s["y"][c][idx], [0.0, 3.0, af])
            cbar = s["C12"][idx].mean(axis=0)
            for key, al in (("ols", 0.0), ("a3", 3.0), ("fade", af)):
                res[key].append(after(dst, c, co[al], cbar) - BEF[(dst, c)][0])
            res["left"].append(after(dst, c, co[0.0], cbar))
            if n >= N[src]:
                break
        cells.append(dict(src=src, dst=dst, n=n, clock=c, before=BEF[(dst, c)][0],
                          **{k: float(np.median(v)) for k, v in res.items()}))
P = pd.DataFrame(cells); P.to_csv(RES / "saliva_cells.csv", index=False)
print(f"  {'par':<18}{'relogio':<13}{'antes':>8}{'OLS':>8}{'a=3':>8}{'some':>8}{'restante':>10}")
for r in P.itertuples():
    print(f"  {r.src[3:]+' -> '+r.dst[3:]:<18}{r.clock:<13}{r.before:>+8.1%}{r.ols:>+8.1%}{r.a3:>+8.1%}{r.fade:>+8.1%}{r.left:>+10.1%}")
A = P[P.clock.isin(AGE)]
k0, k3, kf = int((A.ols > 0).sum()), int((A.a3 > 0).sum()), int((A.fade > 0).sum())
c4 = k3 <= 1
print(f"\n  4. celulas nocivas de {len(A)}: OLS {k0}, alpha=3 {k3}, que some {kf} (barra alpha=3 <= 1) -> {'ok' if c4 else 'FALHOU'}")

section("FECHAMENTO")
for i, (lab, ok) in enumerate([("dano em n=40", c2), ("embaralhado positivo", c3), ("penalidade", c4)], start=2):
    print(f"  {i}. {lab}: {'ok' if ok else 'FALHOU'}")
