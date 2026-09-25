#!/usr/bin/env python3
"""
Stage 32b — does the transport failure generalise to a pace-of-ageing clock?

Every result so far comes from three clocks that estimate AGE. DunedinPACE (stage
32a) estimates the PACE of ageing, was trained on a cohort not on GEO, and is
known to track monocyte subsets — a different kind of clock with a different
relation to blood composition. If the transported correction fails for it too,
the failure is a property of the correction and not of age clocks.

Same protocol as stages 23-31: fit with twelve types, measure with six, null
subtracted, counted per cell; y is DunedinPACE in place of clock age.

PRE-REGISTERED (fixed in the loop prompt before 32a finished)
  1. the stage 32a validation passed (scale and smoking positive control). Hard
     stop otherwise.
  2. n-curve: fitted on 40 age-stratified samples of GSE40279 and transported to
     the other three cohorts, the correction is harmful in more than 60% of draws.
  3. per cell at matched n (12 directed pairs): at least 3 harmful without a
     penalty, at most 1 at alpha = 3.
  4. controls -> cases inside GSE42861: paired excess above +1 point.
  Reported: the permuted reference at n = 40.

Usage:  .venv/bin/python analysis/32b_dunedinpace_transport.py
"""
import sys, gzip
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


CK = "DunedinPACE"
for t in COHORTS:
    a = pd.read_csv(CACHE / f"{t}_ages.csv", index_col=0)
    ix = a.index[a.chrono.notna()]
    dp = pd.read_csv(CACHE / f"{t}_dpace.csv", index_col=0)[CK].reindex(ix)
    assert dp.notna().all(), f"DunedinPACE faltando em {t}"
    D[t]["y"][CK] = dp.to_numpy()

NUL, BEF = {}, {}
for t in COHORTS:
    NUL[t] = null_mean(D[t]["y"][CK], D[t]["chrono"], D[t]["C6"])
    ib, base = inc_and_base(D[t]["y"][CK], D[t]["chrono"], D[t]["C6"])
    BEF[t] = ((ib - NUL[t]) / (1 - base), base)


def net(g, null, bef, b, cbar):
    yc = g["y"][CK] - (g["C12"][:, :-1] - cbar[:-1]) @ b
    ia, _ = inc_and_base(yc, g["chrono"], g["C6"])
    return (ia - null) / (1 - bef[1]) - bef[0]


section("CHECAGEM 2 — A CURVA EM N=40")
rows = []
fit = D["GSE40279"]
for r in range(30):
    idx = stratified_draw(fit["chrono"], 40, np.random.default_rng(4000 + 11 * r))
    af, Cf = fit["chrono"][idx], fit["C12"][idx]
    for kind in ("real", "perm"):
        C_use = Cf if kind == "real" else Cf[np.random.default_rng(60_000 + r).permutation(40)]
        b = ridge_coefs(af, C_use, fit["y"][CK][idx], [0.0])[0.0]
        for dst in ["GSE61151", "GSE50660", "GSE42861"]:
            rows.append(dict(rep=r, kind=kind, dst=dst,
                             delta=net(D[dst], NUL[dst], BEF[dst], b, C_use.mean(axis=0))))
C = pd.DataFrame(rows)
real, perm = C[C.kind == "real"].delta, C[C.kind == "perm"].delta
c2 = (real > 0).mean() > 0.60
print(f"  real: mediana {real.median():+.1%}, nocivo em {(real > 0).mean():.0%} (barra > 60%)"
      f" -> {'ok' if c2 else 'FALHOU'}")
print(f"  embaralhado: mediana {perm.median():+.1%}")

section("CHECAGEM 3 — CELULAS EM N CASADO, OLS vs ALPHA=3")
cells = []
for src, dst in permutations(COHORTS, 2):
    n = min(N[src], N[dst]); s = D[src]; out = {0.0: [], 3.0: []}
    for r in range(20):
        idx = (np.arange(N[src]) if n >= N[src]
               else stratified_draw(s["chrono"], n, np.random.default_rng(700 + 7 * r + n)))
        coefs = ridge_coefs(s["chrono"][idx], s["C12"][idx], s["y"][CK][idx], [0.0, 3.0])
        cbar = s["C12"][idx].mean(axis=0)
        for al, b in coefs.items():
            out[al].append(net(D[dst], NUL[dst], BEF[dst], b, cbar))
    cells.append(dict(src=src, dst=dst, n=n, ols=np.median(out[0.0]), r3=np.median(out[3.0])))
    print(f"  {src[3:]} -> {dst[3:]}  n={n:<4} OLS {cells[-1]['ols']:+.1%}  a=3 {cells[-1]['r3']:+.1%}")
P = pd.DataFrame(cells)
k0, k3 = int((P.ols > 0).sum()), int((P.r3 > 0).sum())
c3 = k0 >= 3 and k3 <= 1
print(f"\n  nocivas: OLS {k0} de 12, alpha=3 {k3} de 12 (barras >=3 e <=1) -> {'ok' if c3 else 'FALHOU'}")

section("CHECAGEM 4 — CONTROLES -> CASOS NO GSE42861")
samples, dis = None, None
with gzip.open(ROOT / "reference/data/GSE42861_series_matrix.txt.gz", "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("!series_matrix_table_begin"):
            break
        p = [v.strip('"') for v in line.rstrip("\n").split("\t")]
        if p[0] == "!Sample_geo_accession":
            samples = p[1:]
        elif p[0] == "!Sample_characteristics_ch1" and p[1].startswith("disease state"):
            dis = [x.partition(":")[2].strip() for x in p[1:]]
a42 = pd.read_csv(CACHE / "GSE42861_ages.csv", index_col=0)
ra = pd.Series(dis, index=samples).reindex(a42.index[a42.chrono.notna()]).str.contains(
    "rheumatoid", case=False).to_numpy()
d = D["GSE42861"]; ictl, icase = np.where(~ra)[0], np.where(ra)[0]


def sub(ix):
    return dict(chrono=d["chrono"][ix], C12=d["C12"][ix], C6=d["C6"][ix], y={CK: d["y"][CK][ix]})


gC = sub(icase)
nC = null_mean(gC["y"][CK], gC["chrono"], gC["C6"])
ibC, bC = inc_and_base(gC["y"][CK], gC["chrono"], gC["C6"]); befC = ((ibC - nC) / (1 - bC), bC)
ex = []
rng = np.random.default_rng(20260926)
for r in range(30):
    pm = rng.permutation(ictl); A, B = pm[: len(pm) // 2], pm[len(pm) // 2:]
    gA, gB = sub(A), sub(B)
    b = ridge_coefs(gA["chrono"], gA["C12"], gA["y"][CK], [0.0])[0.0]; cbar = gA["C12"].mean(axis=0)
    nB = null_mean(gB["y"][CK], gB["chrono"], gB["C6"], n_perm=300)
    ibB, bB = inc_and_base(gB["y"][CK], gB["chrono"], gB["C6"]); befB = ((ibB - nB) / (1 - bB), bB)
    ex.append(net(gC, nC, befC, b, cbar) - net(gB, nB, befB, b, cbar))
ex = np.array(ex)
c4 = np.median(ex) > 0.01
print(f"  excesso pareado: mediana {np.median(ex):+.1%}, IQR {np.quantile(ex,.25):+.1%} a "
      f"{np.quantile(ex,.75):+.1%} (barra > +1 ponto) -> {'ok' if c4 else 'FALHOU'}")

section("FECHAMENTO")
print(f"  2. curva n=40: {'ok' if c2 else 'FALHOU'}")
print(f"  3. celulas OLS vs alpha=3: {'ok' if c3 else 'FALHOU'}")
print(f"  4. controles -> casos: {'ok' if c4 else 'FALHOU'}")
