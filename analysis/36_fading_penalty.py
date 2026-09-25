#!/usr/bin/env python3
"""
Stage 36 — a penalty that fades with n, against one that does not.

The ridge penalty used since stage 21 is alpha x (mean eigenvalue of the fitting
cohort's standardised composition cross-product). That eigenvalue grows with n,
so alpha = 3 shrinks by the same proportion at every fitting size. Stage 35 found
the cost: with 2,639 fitting samples the unpenalised correction was better on
median (-4.0% against -3.0%). A textbook ridge has a fixed lambda, whose relative
weight fades as 1/n. Written here as alpha_n = 3 x 40 / n, it equals alpha = 3 at
n = 40 and is 0.045 at n = 2,639.

Patil, Du & Tibshirani (2024) separate covariate shift from regression shift (our
model shift) and show the optimal ridge level depends on which is present. For
estimation noise the fading penalty is the natural one. For model shift the error
does not shrink with n, so a penalty that vanishes with n cannot bound it. That is
the prediction under test.

Cohorts: the six cached (GSE40279, GSE61151, GSE50660, GSE42861, GSE132203,
GSE55763), every directed pair at matched n = min(n_A, n_B), 20 stratified draws.
Clocks: the three age clocks, each out of any pair containing a cohort it trained
on, and on GSE132203 only those that cleared 33a. Counted per (pair x clock) cell,
median over draws. DunedinPACE reported alongside, no bar.

PRE-REGISTERED (written before either penalty was computed on these pairs)
  1. the 33a and 35a caches exist. Hard stop.
  2. MODEL SHIFT NEEDS THE PROPORTIONAL PENALTY: at matched n, the fading penalty
     leaves at least twice as many harmful cells as alpha = 3, and at least 3 more.
  3. THE FADING PENALTY RECOVERS THE LARGE-N BENEFIT: fitted on all of GSE55763,
     its median over the 13 age-clock cells is within 0.3 points of the
     unpenalised median (alpha = 3 was 1.0 point away in stage 35).
  4. AGE EXTRAPOLATION (the hypothesis stage 35 wrote after its result): fitted on
     all of GSE55763 (ages 24-75), the error left in GSE40279 restricted to ages
     24-75 is less than half the error left in all of GSE40279, for both clocks
     clean on GSE40279.

Usage:  .venv/bin/python analysis/36_fading_penalty.py
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
ALL = COHORTS + [NEW, BIG]
for t in ALL:
    aa = pd.read_csv(CACHE / f"{t}_ages.csv", index_col=0)
    D[t]["y"]["DunedinPACE"] = (pd.read_csv(CACHE / f"{t}_dpace.csv", index_col=0)["DunedinPACE"]
                                .reindex(aa.index[aa.chrono.notna()]).to_numpy())
# GSE40279 restricted to the fitting cohort's age range, for check 4
R40 = "GSE40279r"
keep = (D["GSE40279"]["chrono"] >= 24) & (D["GSE40279"]["chrono"] <= 75)
D[R40] = dict(chrono=D["GSE40279"]["chrono"][keep], C12=D["GSE40279"]["C12"][keep],
              C6=D["GSE40279"]["C6"][keep],
              y={c: v[keep] for c, v in D["GSE40279"]["y"].items()})
N[R40] = int(keep.sum())
AGE = CLOCKS3
print(f"  coortes: {', '.join(f'{t} {N[t]}' for t in ALL)}; GSE40279 com 24-75 anos: {N[R40]}")


def clocks_for(src, dst):
    out = [c for c in AGE if allowed(c, src, dst)
           and (NEW not in (src, dst) or c in USABLE[NEW])
           and (BIG not in (src, dst) or c in USABLE[BIG])]
    return out + ["DunedinPACE"]


print("  nulos ...", flush=True)
NUL, BEF = {}, {}
for t in ALL + [R40]:
    for c in AGE + ["DunedinPACE"]:
        NUL[(t, c)] = null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["C6"], n_perm=500)
        ib, base = inc_and_base(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        BEF[(t, c)] = ((ib - NUL[(t, c)]) / (1 - base), base)


def after(dst, c, b, cbar):
    d = D[dst]
    yc = d["y"][c] - (d["C12"][:, :-1] - cbar[:-1]) @ b
    ia, _ = inc_and_base(yc, d["chrono"], d["C6"])
    return (ia - NUL[(dst, c)]) / (1 - BEF[(dst, c)][1])


def fade(n):
    return 3.0 * 40 / n


section("CHECAGEM 2 — PARES DIRECIONADOS EM N PAREADO")
cells = []
for src, dst in permutations(ALL, 2):
    n = min(N[src], N[dst]); s = D[src]; af = fade(n)
    res = {}
    for c in clocks_for(src, dst):
        res = {"ols": [], "a3": [], "fade": []}
        for r in range(20):
            idx = (np.arange(N[src]) if n >= N[src]
                   else stratified_draw(s["chrono"], n, np.random.default_rng(3600 + 7 * r + n)))
            co = ridge_coefs(s["chrono"][idx], s["C12"][idx], s["y"][c][idx], [0.0, 3.0, af])
            cbar = s["C12"][idx].mean(axis=0)
            for key, al in (("ols", 0.0), ("a3", 3.0), ("fade", af)):
                res[key].append(after(dst, c, co[al], cbar) - BEF[(dst, c)][0])
            if n >= N[src]:
                break
        cells.append(dict(src=src, dst=dst, n=n, clock=c, alpha_fade=af,
                          **{k: float(np.median(v)) for k, v in res.items()}))
    print(f"  {src[3:]} -> {dst[3:]} (n={n}, alpha_fade={af:.3f}): pronto", flush=True)
P = pd.DataFrame(cells)
P.to_csv(RES / "fading_penalty_cells.csv", index=False)
A = P[P.clock != "DunedinPACE"]
k0, k3, kf = int((A.ols > 0).sum()), int((A.a3 > 0).sum()), int((A.fade > 0).sum())
c2 = kf >= 2 * k3 and kf >= k3 + 3
print(f"\n  celulas de idade: {len(A)}; nocivas OLS {k0}, alpha=3 {k3}, penalidade que some {kf}")
print(f"  2. barra: que some >= 2x alpha=3 e >= +3 -> {'ok' if c2 else 'FALHOU'}")
for r in A[(A.fade > 0) | (A.a3 > 0)].sort_values("fade", ascending=False).itertuples():
    print(f"     {r.src[3:]+' -> '+r.dst[3:]:<16}{r.clock:<13} n={r.n:<5} OLS {r.ols:+.1%}  a=3 {r.a3:+.1%}  some {r.fade:+.1%}")
helped = A[A.ols < 0]
print(f"  onde OLS ajudou ({len(helped)} celulas), mediana: OLS {helped.ols.median():+.1%}, "
      f"a=3 {helped.a3.median():+.1%}, que some {helped.fade.median():+.1%}")
dp = P[P.clock == "DunedinPACE"]
print(f"  DunedinPACE (sem barra): nocivas OLS {(dp.ols > 0).sum()}, a=3 {(dp.a3 > 0).sum()}, "
      f"que some {(dp.fade > 0).sum()} de {len(dp)}")

section("CHECAGEM 3 E 4 — AJUSTE NO GSE55763 INTEIRO")
s = D[BIG]; idx = np.arange(N[BIG]); af = fade(N[BIG])
rows = []
for dst in COHORTS + [NEW, R40]:
    for c in clocks_for(BIG, "GSE40279" if dst == R40 else dst):
        if c == "DunedinPACE":
            continue
        co = ridge_coefs(s["chrono"], s["C12"], s["y"][c], [0.0, 3.0, af])
        cbar = s["C12"].mean(axis=0)
        rows.append(dict(dst=dst, clock=c,
                         **{k: after(dst, c, co[al], cbar) - BEF[(dst, c)][0]
                            for k, al in (("ols", 0.0), ("a3", 3.0), ("fade", af))},
                         left=after(dst, c, co[0.0], cbar)))
F = pd.DataFrame(rows)
F.to_csv(RES / "fading_penalty_full.csv", index=False)
M = F[F.dst != R40]
gap = abs(M.fade.median() - M.ols.median())
c3 = gap <= 0.003
print(f"  3. mediana em {len(M)} celulas: OLS {M.ols.median():+.1%}, a=3 {M.a3.median():+.1%}, "
      f"que some {M.fade.median():+.1%} (distancia {gap*100:.2f} ponto, barra <= 0,3) -> {'ok' if c3 else 'FALHOU'}")
full40 = F[F.dst == "GSE40279"].set_index("clock").left
rest40 = F[F.dst == R40].set_index("clock").left
c4 = all(rest40[c] < 0.5 * full40[c] for c in full40.index)
for c in full40.index:
    print(f"  4. {c:<13} erro restante GSE40279 inteiro {full40[c]:+.1%}, so 24-75 anos {rest40[c]:+.1%}")
print(f"     barra: restrito < metade do inteiro nos dois -> {'ok' if c4 else 'FALHOU'}")

section("FECHAMENTO")
for i, (lab, ok) in enumerate([("model shift pede a penalidade proporcional", c2),
                                ("a penalidade que some recupera o beneficio em n grande", c3),
                                ("extrapolacao de idade explica o piso no GSE40279", c4)], start=2):
    print(f"  {i}. {lab}: {'ok' if ok else 'FALHOU'}")
