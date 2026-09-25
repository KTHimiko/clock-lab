#!/usr/bin/env python3
"""
Stage 33b — do the central results hold on a fifth cohort, another array and
another ancestry?

GSE132203 (Grady Trauma Project, EPIC, mostly African American) was cached in
33a with the same deconvolution panels as every earlier stage. Transports to and
from it are also cross-platform — EPIC against 450k — which is the situation any
reuse of published coefficients across array generations is in.

Clocks per pair: a clock is scored on a pair only if it cleared 33a's coverage
and age checks on GSE132203 and is clean for the other cohort (Hannum and Horvath
2013 are out of any pair with GSE40279). DunedinPACE is scored everywhere.
Counted per (directed pair x clock) cell, as always.

PRE-REGISTERED (written before any GSE132203 transport was computed)
  1. the 33a cache exists and 33a's panel fingerprint passed. Hard stop.
  2. THE PENALTY GENERALISES: over the 8 directed pairs between GSE132203 and the
     four original cohorts at matched n, at least 3 cells harmful unpenalised and
     at most 10% of cells harmful at alpha = 3.
  3. SMALL-N HARM, from a new fitting cohort: fitted on 40 age-stratified samples
     of GSE132203 and transported to the four original cohorts, the AGE clocks'
     correction is harmful in more than 60% of draws. DunedinPACE reported
     separately, no bar — stage 32 showed its small-n behaviour depends on the
     fitting cohort.
  4. the permuted reference at n = 40 is positive (the noise component exists
     from this cohort too).

Usage:  .venv/bin/python analysis/33b_new_cohort_transport.py
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


NEW = "GSE132203"
need = [CACHE / f"{NEW}_{w}.csv" for w in ("comp12", "comp6", "ages", "dpace", "usable_clocks")]
if not all(f.exists() for f in need):
    sys.exit("parando: cache da 33a ausente.")
usable = pd.read_csv(CACHE / f"{NEW}_usable_clocks.csv").iloc[:, 0].tolist()
a = pd.read_csv(CACHE / f"{NEW}_ages.csv", index_col=0)
k = a.chrono.notna().to_numpy(); ix = a.index[k]
D[NEW] = dict(chrono=a.chrono[k].to_numpy(),
              y={c: a[c][k].to_numpy() for c in CLOCKS},
              C12=pd.read_csv(CACHE / f"{NEW}_comp12.csv", index_col=0).loc[ix, TYPES12].to_numpy(),
              C6=pd.read_csv(CACHE / f"{NEW}_comp6.csv", index_col=0).loc[ix, TYPES6].to_numpy())
N[NEW] = int(k.sum())
ALL = COHORTS + [NEW]
for t in ALL:
    aa = pd.read_csv(CACHE / f"{t}_ages.csv", index_col=0)
    ixx = aa.index[aa.chrono.notna()]
    D[t]["y"]["DunedinPACE"] = pd.read_csv(CACHE / f"{t}_dpace.csv", index_col=0)["DunedinPACE"].reindex(ixx).to_numpy()
AGE_CLOCKS = ["Horvath2013", "Levine2018", "Horvath2018"]


def clocks_for(src, dst):
    out = [c for c in AGE_CLOCKS if allowed(c, src, dst) and (NEW not in (src, dst) or c in usable)]
    return out + ["DunedinPACE"]


print("  nulos ...", flush=True)
NUL, BEF = {}, {}
for t in ALL:
    for c in AGE_CLOCKS + ["DunedinPACE"]:
        NUL[(t, c)] = null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["C6"], n_perm=500)
        ib, base = inc_and_base(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        BEF[(t, c)] = ((ib - NUL[(t, c)]) / (1 - base), base)


def net(dst, c, b, cbar):
    d = D[dst]
    yc = d["y"][c] - (d["C12"][:, :-1] - cbar[:-1]) @ b
    ia, _ = inc_and_base(yc, d["chrono"], d["C6"])
    before, base = BEF[(dst, c)]
    return (ia - NUL[(dst, c)]) / (1 - base) - before


section("CHECAGEM 2 — A PENALIDADE NOS 8 PARES COM O GSE132203")
cells = []
for src, dst in [(NEW, o) for o in COHORTS] + [(o, NEW) for o in COHORTS]:
    n = min(N[src], N[dst]); s = D[src]
    for c in clocks_for(src, dst):
        res = {0.0: [], 3.0: []}
        for r in range(20):
            idx = (np.arange(N[src]) if n >= N[src]
                   else stratified_draw(s["chrono"], n, np.random.default_rng(1100 + 7 * r + n)))
            coefs = ridge_coefs(s["chrono"][idx], s["C12"][idx], s["y"][c][idx], [0.0, 3.0])
            cbar = s["C12"][idx].mean(axis=0)
            for al, b in coefs.items():
                res[al].append(net(dst, c, b, cbar))
        cells.append(dict(src=src, dst=dst, n=n, clock=c,
                          ols=float(np.median(res[0.0])), r3=float(np.median(res[3.0]))))
    print(f"  {src[3:]} -> {dst[3:]} (n={n}): pronto", flush=True)
P = pd.DataFrame(cells)
P.to_csv(RES / "new_cohort_cells.csv", index=False)
print(f"\n  {'par':<18}{'relogio':<13}{'OLS':>9}{'a=3':>9}")
for r in P.sort_values("ols", ascending=False).itertuples():
    print(f"  {r.src[3:]+' -> '+r.dst[3:]:<18}{r.clock:<13}{r.ols:>+9.1%}{r.r3:>+9.1%}"
          + ("   <- nocivo com a=3" if r.r3 > 0 else ""))
k0, k3 = int((P.ols > 0).sum()), int((P.r3 > 0).sum())
c2 = k0 >= 3 and k3 <= 0.10 * len(P)
print(f"\n  nocivas: OLS {k0} de {len(P)}, alpha=3 {k3} de {len(P)} "
      f"(barras >=3 e <= 10%) -> {'ok' if c2 else 'FALHOU'}")

section("CHECAGEM 3 E 4 — A CURVA EM N=40 AJUSTANDO NO GSE132203")
rows = []
s = D[NEW]
for r in range(30):
    idx = stratified_draw(s["chrono"], 40, np.random.default_rng(5000 + 11 * r))
    af, Cf = s["chrono"][idx], s["C12"][idx]
    Cp = Cf[np.random.default_rng(70_000 + r).permutation(40)]
    for dst in COHORTS:
        for c in clocks_for(NEW, dst):
            for kind, C_use in (("real", Cf), ("perm", Cp)):
                b = ridge_coefs(af, C_use, s["y"][c][idx], [0.0])[0.0]
                rows.append(dict(rep=r, dst=dst, clock=c, kind=kind,
                                 delta=net(dst, c, b, C_use.mean(axis=0))))
C = pd.DataFrame(rows)
C.to_csv(RES / "new_cohort_curve40.csv", index=False)
age = C[(C.kind == "real") & C.clock.isin(AGE_CLOCKS)].delta
dp = C[(C.kind == "real") & (C.clock == "DunedinPACE")].delta
perm = C[C.kind == "perm"].delta
c3 = (age > 0).mean() > 0.60
c4 = perm.median() > 0
print(f"  relogios de idade: mediana {age.median():+.1%}, nocivo em {(age > 0).mean():.0%}"
      f" (barra > 60%) -> {'ok' if c3 else 'FALHOU'}")
print(f"  DunedinPACE (sem barra): mediana {dp.median():+.1%}, nocivo em {(dp > 0).mean():.0%}")
print(f"  embaralhado: mediana {perm.median():+.1%} (barra > 0) -> {'ok' if c4 else 'FALHOU'}")
for c in AGE_CLOCKS + ["DunedinPACE"]:
    g = C[(C.kind == "real") & (C.clock == c)].delta
    if len(g):
        print(f"    {c:<13} mediana {g.median():+.1%}  nocivo {(g > 0).mean():.0%}")

section("FECHAMENTO")
print(f"  2. penalidade: {'ok' if c2 else 'FALHOU'}")
print(f"  3. dano em n=40 (idade), ajuste no GSE132203: {'ok' if c3 else 'FALHOU'}")
print(f"  4. referencia embaralhada positiva: {'ok' if c4 else 'FALHOU'}")
