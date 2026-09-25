#!/usr/bin/env python3
"""
Stage 35b — does the model-shift floor survive a fitting cohort four times larger?

Stage 24 split the transported correction's error into estimation noise (steep in
n, reproduced by permuted coefficients) and an n-independent floor, and read the
floor as model shift — the composition coefficients differ between cohorts. But
the largest fitting cohort was 656 (GSE40279), so "n-independent" rested on one
cohort's top end. With GSE55763 (2,639 unrelated adults after 35a dropped the
technical replicates) the fitting n can be pushed four times further. If the floor
were slow-decaying noise it would keep falling; if it is model shift it stays.

Fitting cohort GSE55763; test cohorts the four 450k cohorts and GSE132203 (EPIC).
Clocks: the age clocks that cleared 35a's checks, each out of any pair containing a
cohort it trained on (Hannum and Horvath 2013 out of GSE40279); counted per
(test cohort x clock) cell. DunedinPACE reported alongside, never pooled.

"error left" = what the six-type panel still sees of composition in the test
cohort after the correction, net of its own permutation null, as a share of
age-acceleration variance (stage 24's "after"). Permuted = the same fit with the
composition rows shuffled whole.

PRE-REGISTERED (written before any GSE55763 transport was computed)
  1. the 35a cache exists, with >= 2 usable age clocks. Hard stop.
  2. THE FLOOR PERSISTS: over the age-clock cells, median error left at the full n
     is at least half its median at n = 656. Pure 1/n noise predicts a quarter.
  3. THE FLOOR IS NOT NOISE: at the full n, median error left (real) exceeds the
     median permuted delta, in at least 2/3 of age-clock cells.
  4. the permuted reference still falls as ~1/n over 40..full: log-log slope of
     its median (sizes where positive) in [-1.3, -0.7].
  5. SMALL-N HARM from a third fitting cohort: at n = 40 the age clocks'
     correction is harmful in more than 60% of draws.
  6. THE PENALTY: at the full n and at n = 40, alpha = 3 leaves at most 10% of
     age-clock cells harmful (median over draws).

Usage:  .venv/bin/python analysis/35b_large_fit_floor.py
"""
import sys
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




BIG, NEW = "GSE55763", "GSE132203"
if not all((CACHE / f"{t}_{w}.csv").exists() for t in (BIG, NEW)
           for w in ("comp12", "comp6", "ages", "dpace", "usable_clocks")):
    sys.exit("parando: cache da 33a/35a ausente.")
USABLE = {t: pd.read_csv(CACHE / f"{t}_usable_clocks.csv").iloc[:, 0].tolist() for t in (BIG, NEW)}
for t in (BIG, NEW):
    a = pd.read_csv(CACHE / f"{t}_ages.csv", index_col=0)
    k = a.chrono.notna().to_numpy(); ix = a.index[k]
    D[t] = dict(chrono=a.chrono[k].to_numpy(),
                y={c: a[c][k].to_numpy() for c in CLOCKS},
                C12=pd.read_csv(CACHE / f"{t}_comp12.csv", index_col=0).loc[ix, TYPES12].to_numpy(),
                C6=pd.read_csv(CACHE / f"{t}_comp6.csv", index_col=0).loc[ix, TYPES6].to_numpy())
    N[t] = int(k.sum())
for t in COHORTS + [BIG, NEW]:
    aa = pd.read_csv(CACHE / f"{t}_ages.csv", index_col=0)
    D[t]["y"]["DunedinPACE"] = (pd.read_csv(CACHE / f"{t}_dpace.csv", index_col=0)["DunedinPACE"]
                                .reindex(aa.index[aa.chrono.notna()]).to_numpy())
AGE = [c for c in CLOCKS3 if c in USABLE[BIG]]
c1 = len(AGE) >= 2
print(f"  {BIG}: n = {N[BIG]}, relogios de idade utilizaveis {AGE} -> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("parando.")
DSTS = COHORTS + [NEW]


def clocks_for(dst):
    return [c for c in AGE if allowed(c, BIG, dst) and (dst != NEW or c in USABLE[NEW])]


print("  nulos ...", flush=True)
NUL, BEF = {}, {}
for t in DSTS:
    for c in AGE + ["DunedinPACE"]:
        NUL[(t, c)] = null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["C6"], n_perm=500)
        ib, base = inc_and_base(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        BEF[(t, c)] = ((ib - NUL[(t, c)]) / (1 - base), base)


def after(dst, c, b, cbar):
    d = D[dst]
    yc = d["y"][c] - (d["C12"][:, :-1] - cbar[:-1]) @ b
    ia, _ = inc_and_base(yc, d["chrono"], d["C6"])
    return (ia - NUL[(dst, c)]) / (1 - BEF[(dst, c)][1])


section("A CURVA ATE O N CHEIO DO GSE55763")
s = D[BIG]
SIZES = [40, 80, 160, 320, 656, 1300, N[BIG]]
rows = []
for n in SIZES:
    for r in range(30):
        idx = (np.arange(N[BIG]) if n >= N[BIG]
               else stratified_draw(s["chrono"], n, np.random.default_rng(35_000 + 13 * r + n)))
        af, Cf = s["chrono"][idx], s["C12"][idx]
        Cp = Cf[np.random.default_rng(90_000 + 17 * r + n).permutation(len(idx))]
        for dst in DSTS:
            for c in clocks_for(dst) + ["DunedinPACE"]:
                yf = s["y"][c][idx]
                if n < N[BIG] or r == 0:
                    for al, b in ridge_coefs(af, Cf, yf, [0.0, 3.0]).items():
                        aft = after(dst, c, b, Cf.mean(axis=0))
                        rows.append(dict(n=n, rep=r, dst=dst, clock=c,
                                         kind="real" if al == 0 else "ridge3",
                                         after=aft, delta=aft - BEF[(dst, c)][0]))
                b = ridge_coefs(af, Cp, yf, [0.0])[0.0]
                aft = after(dst, c, b, Cp.mean(axis=0))
                rows.append(dict(n=n, rep=r, dst=dst, clock=c, kind="perm",
                                 after=aft, delta=aft - BEF[(dst, c)][0]))
    print(f"  n = {n}: pronto", flush=True)
R = pd.DataFrame(rows)
R.to_csv(RES / "large_fit_curve.csv", index=False)

cell = (R.groupby(["n", "dst", "clock", "kind"])[["after", "delta"]].median().reset_index())
ag = cell[cell.clock.isin(AGE)]
print(f"\n  relogios de idade, mediana entre celulas:")
print(f"  {'n':>6}{'delta real':>12}{'erro restante':>15}{'embaralhado':>13}{'ridge a=3':>11}")
for n in SIZES:
    g = ag[ag.n == n].pivot_table(index=["dst", "clock"], columns="kind", values=["after", "delta"])
    print(f"  {n:>6}{g[('delta', 'real')].median():>+12.1%}{g[('after', 'real')].median():>+15.1%}"
          f"{g[('delta', 'perm')].median():>+13.1%}{g[('delta', 'ridge3')].median():>+11.1%}")
dpc = cell[(cell.clock == "DunedinPACE")]
print("\n  DunedinPACE (sem barra), mediana entre destinos:")
for n in SIZES:
    g = dpc[dpc.n == n].pivot_table(index="dst", columns="kind", values="delta")
    print(f"  {n:>6}  real {g.real.median():+.1%}  embaralhado {g.perm.median():+.1%}  a=3 {g.ridge3.median():+.1%}")

section("CHECAGENS")
full = N[BIG]
w = ag.pivot_table(index=["n", "dst", "clock"], columns="kind", values=["after", "delta"])
a656 = w.loc[656][("after", "real")].median(); afull = w.loc[full][("after", "real")].median()
c2 = afull >= 0.5 * a656
print(f"  2. erro restante: n=656 {a656:+.1%}, n={full} {afull:+.1%}, razao {afull / a656:.2f}"
      f" (barra >= 0,5; ruido puro ~0,25) -> {'ok' if c2 else 'FALHOU'}")
wf = w.loc[full]
share = float((wf[("after", "real")] > wf[("delta", "perm")]).mean())
c3 = share >= 2 / 3
print(f"  3. erro restante > embaralhado no n cheio em {share:.0%} das celulas (barra >= 67%)"
      f" -> {'ok' if c3 else 'FALHOU'}")
for (dst, c), v in wf.iterrows():
    print(f"     {dst:<10}{c:<13} restante {v[('after', 'real')]:+.1%}  embaralhado {v[('delta', 'perm')]:+.1%}")
pm = ag[ag.kind == "perm"].groupby("n").delta.median()
pos = pm[pm > 0]
slope = float(np.polyfit(np.log(pos.index.to_numpy(float)), np.log(pos.to_numpy()), 1)[0]) if len(pos) >= 3 else np.nan
c4 = -1.3 <= slope <= -0.7
print(f"  4. inclinacao log-log do embaralhado: {slope:.2f} em {len(pos)} tamanhos (barra [-1,3; -0,7])"
      f" -> {'ok' if c4 else 'FALHOU'}")
d40 = R[(R.n == 40) & (R.kind == "real") & R.clock.isin(AGE)].delta
c5 = (d40 > 0).mean() > 0.60
print(f"  5. n=40, relogios de idade: nocivo em {(d40 > 0).mean():.0%} dos sorteios, mediana {d40.median():+.1%}"
      f" (barra > 60%) -> {'ok' if c5 else 'FALHOU'}")
k_full = int((w.loc[full][("delta", "ridge3")] > 0).sum()); k_40 = int((w.loc[40][("delta", "ridge3")] > 0).sum())
nc = len(w.loc[full])
k_full0 = int((w.loc[full][("delta", "real")] > 0).sum()); k_400 = int((w.loc[40][("delta", "real")] > 0).sum())
c6 = k_full <= 0.10 * nc and k_40 <= 0.10 * nc
print(f"  6. celulas nocivas, OLS -> a=3: n cheio {k_full0} -> {k_full} de {nc}; n=40 {k_400} -> {k_40} de {nc}"
      f" (barra <= 10%) -> {'ok' if c6 else 'FALHOU'}")

section("FECHAMENTO")
for i, (lab, ok) in enumerate([("piso persiste", c2), ("piso nao e ruido", c3), ("embaralhado ~1/n", c4),
                                ("dano em n=40", c5), ("penalidade", c6)], start=2):
    print(f"  {i}. {lab}: {'ok' if ok else 'FALHOU'}")
