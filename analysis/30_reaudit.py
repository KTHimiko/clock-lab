#!/usr/bin/env python3
"""
Stage 30 — re-auditing the two manuscript claims computed before the audit.

Stages 23 to 29 changed the protocol: Horvath 2013 is out of every pair touching
GSE40279 (it was trained there), and damage is counted per clock, never as a
median over clocks (stage 26 showed pooling hides harm). Two claims the
manuscript still makes were computed under the old protocol:

  stage 19  cross-validation on the fitting cohort does not size the penalty —
            on the cohort that breaks the transport it picks alpha = 0.3 and
            leaves +2.9% of the damage standing
  stage 22  with the two cohorts deconvolved from different reference panels,
            the penalty still holds every pair

Both were medians over clocks, and both included Horvath 2013 on GSE40279. They
are recomputed here under the current protocol before the manuscript quotes
either again.

PRE-REGISTERED
  1. the cache is the stage 18 cache. Hard stop.
  2. CROSS-VALIDATION. For every (directed pair x clean clock) cell at matched n,
     the penalty is chosen by leave-one-out on the fitting subsample over the
     grid {0, 0.1, 0.3, 1, 3, 10} and the transport's net damage measured. The
     claim "cross-validation does not size the penalty" SURVIVES if the share of
     harmful cells under the cross-validated penalty exceeds the share under
     fixed alpha = 3 by at least ten points.
  3. PANEL MISMATCH. The fitting cohort's composition estimated from the Salas
     reference collapsed to six labels, the test cohort's from the Reinius
     reference, measured on the twelve-type panel of the test cohort (the
     stage 22 design). Per (pair x clean clock) cell: the share harmful under
     OLS and under alpha = 3. The claim "the penalty holds under panel mismatch"
     SURVIVES if alpha = 3 leaves no more than one harmful cell.

Usage:  .venv/bin/python analysis/30_reaudit.py
"""
import sys
from itertools import permutations
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_extended import TYPES12, COLLAPSE
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


GRIDA = [0.0, 0.1, 0.3, 1.0, 3.0, 10.0]


def ridge_loo(age, C, y, alphas):
    """Coefficients over the grid and the alpha leave-one-out would pick —
    the same computation as stage 18's, on the fitting subsample only."""
    Ct = partial_out(age, C[:, :-1])
    yt = partial_out(age, y.reshape(-1, 1)).ravel()
    mu, sd = Ct.mean(axis=0), Ct.std(axis=0) + 1e-12
    Z = (Ct - mu) / sd
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    scale = float((s ** 2).mean()); Uty = U.T @ yt
    out, loo = {}, {}
    for a in alphas:
        lam = a * scale
        out[a] = (Vt.T @ (s * Uty / (s ** 2 + lam))) / sd
        h = (U * (s ** 2 / (s ** 2 + lam))) @ U.T
        lev = np.clip(np.diag(h), 0, 1 - 1e-9)
        loo[a] = float(np.mean(((yt - h @ yt) / (1 - lev)) ** 2))
    return out, min(loo, key=loo.get)


def collapse6(C12):
    df = pd.DataFrame(C12, columns=TYPES12)
    c6 = pd.DataFrame({t: df[m].sum(axis=1) for t, m in COLLAPSE.items()})[TYPES6]
    return c6.div(c6.sum(axis=1), axis=0).to_numpy()


def cells(fit_panel, apply_panel, meas_panel, alphas_mode):
    rows = []
    for src, dst in permutations(COHORTS, 2):
        n = min(N[src], N[dst]); s, d = D[src], D[dst]
        for r in range(20):
            idx = (np.arange(N[src]) if n >= N[src]
                   else stratified_draw(s["chrono"], n, np.random.default_rng(950 + 7 * r + n)))
            af = s["chrono"][idx]; Cf = fit_panel(s)[idx]; cbar = Cf.mean(axis=0)
            for c in CLOCKS3:
                if not allowed(c, src, dst):
                    continue
                coefs, a_cv = ridge_loo(af, Cf, s["y"][c][idx], GRIDA)
                y = d["y"][c]; a = d["chrono"]; Mt = meas_panel(d)
                nb = NULLM[(dst, c, meas_panel.__name__)]
                ib, base = inc_and_base(y, a, Mt)
                before = (ib - nb) / (1 - base)
                picks = {"ols": 0.0, "a3": 3.0, "cv": a_cv} if alphas_mode == "cv" \
                    else {"ols": 0.0, "a3": 3.0}
                for label, al in picks.items():
                    yc = y - (apply_panel(d)[:, :-1] - cbar[:-1]) @ coefs[al]
                    ia, _ = inc_and_base(yc, a, Mt)
                    rows.append(dict(src=src, dst=dst, clock=c, pick=label, alpha=al,
                                     delta=(ia - nb) / (1 - base) - before))
        print(f"  {src[3:]} -> {dst[3:]}", flush=True)
    return pd.DataFrame(rows)


def six(d): return d["C6"]
def twelve(d): return d["C12"]
def salas6(d): return collapse6(d["C12"])


print("  nulos ...", flush=True)
NULLM = {}
for t in COHORTS:
    for c in CLOCKS3:
        for f in (six, twelve):
            NULLM[(t, c, f.__name__)] = null_mean(D[t]["y"][c], D[t]["chrono"], f(D[t]))

section("CHECAGEM 2 — A VALIDACAO CRUZADA DIMENSIONA A PENALIDADE?")
CV = cells(twelve, twelve, six, "cv")
CV.to_csv(RES / "reaudit_cv.csv", index=False)
cell = CV.groupby(["src", "dst", "clock", "pick"]).delta.median().unstack("pick").reset_index()
alpha_cv = CV[CV.pick == "cv"].groupby(["src", "dst", "clock"]).alpha.median()
h = {k: (cell[k] > 0).mean() for k in ("ols", "a3", "cv")}
print(f"  {len(cell)} celulas | nocivas: OLS {h['ols']:.0%}  alpha=3 {h['a3']:.0%}  "
      f"alpha por validacao cruzada {h['cv']:.0%}")
print(f"  alpha escolhido pela validacao cruzada, mediana das celulas: {alpha_cv.median():g}"
      f" (distribuicao {alpha_cv.value_counts().sort_index().to_dict()})")
bad = cell[cell.cv > 0].sort_values("cv", ascending=False)
for r in bad.head(6).itertuples():
    print(f"    {r.src[3:]+' -> '+r.dst[3:]:<16}{r.clock:<13} cv {r.cv:+.1%}  a3 {r.a3:+.1%}"
          f"  alpha_cv {alpha_cv.loc[(r.src, r.dst, r.clock)]:g}")
c2 = h["cv"] - h["a3"] >= 0.10
print(f"\n  diferenca {100*(h['cv']-h['a3']):+.0f} pontos (barra +10) -> afirmacao "
      f"{'SOBREVIVE' if c2 else 'NAO SOBREVIVE'}")

section("CHECAGEM 3 — PAINEIS DIFERENTES NAS DUAS PONTAS")
PM = cells(salas6, six, twelve, "fixed")
PM.to_csv(RES / "reaudit_panel.csv", index=False)
cellp = PM.groupby(["src", "dst", "clock", "pick"]).delta.median().unstack("pick").reset_index()
hp = {k: int((cellp[k] > 0).sum()) for k in ("ols", "a3")}
print(f"  {len(cellp)} celulas | nocivas: OLS {hp['ols']}  alpha=3 {hp['a3']}")
for r in cellp[cellp.a3 > 0].itertuples():
    print(f"    {r.src[3:]+' -> '+r.dst[3:]:<16}{r.clock:<13} OLS {r.ols:+.1%}  a3 {r.a3:+.1%}")
c3 = hp["a3"] <= 1
print(f"\n  -> afirmacao {'SOBREVIVE' if c3 else 'NAO SOBREVIVE'}")

section("FECHAMENTO")
print(f"  2. validacao cruzada nao dimensiona: {'sobrevive' if c2 else 'nao sobrevive'}")
print(f"  3. penalidade sob paineis diferentes: {'sobrevive' if c3 else 'nao sobrevive'}")
