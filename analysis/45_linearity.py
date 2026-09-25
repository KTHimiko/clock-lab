#!/usr/bin/env python3
"""
Stage 45 — is the composition effect linear, and does it matter if it is not?

Every correction here is linear in the proportions, and the decomposition into
estimation noise and model shift assumes it. If the true effect is curved, part
of what this project calls model shift could be a shared nonlinearity that every
cohort has, fitted at a different point of the curve in each. That would be a
different diagnosis with a different remedy.

Two tests, on the six blood cohorts, Levine 2018 and Horvath 2018:

  A. Is there curvature at all? Within each cohort, compare the linear fit
     y ~ age + C with one that adds squared terms for the four largest
     components (Neu, CD4mem, CD8mem, Mono), by F test, and report the extra
     variance explained.
  B. Does a quadratic correction transport better? Fit and transport both the
     linear and the quadratic correction at matched n over all directed pairs,
     and count harmful cells. If curvature were driving the harm, allowing for
     it should reduce it.

PRE-REGISTERED (written before either was computed)
  1. the cache loads. Hard stop.
  2. CURVATURE IS SMALL: the median extra R^2 from the squared terms is below
     0.01 within cohorts. If it fails, the linear assumption is doing real work
     and the decomposition needs the caveat stated, not assumed away.
  3. IT IS NOT THE MECHANISM: at matched n, the quadratic correction leaves at
     least as many harmful cells as the linear one. If it leaves fewer, part of
     the harm attributed to model shift is unmodelled curvature.

Usage:  .venv/bin/python analysis/45_linearity.py"""
import sys
from itertools import permutations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"; CACHE = RES / "cache"
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_extended import TYPES12
from model.clocks import CLOCKS
from model.deconvolution import TYPES as TYPES6
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




BLOOD = ["GSE40279", "GSE61151", "GSE50660", "GSE42861", "GSE132203", "GSE55763"]
AGE = ["Levine2018", "Horvath2018"]
BIG4 = ["Neu", "CD4mem", "CD8mem", "Mono"]
USABLE = {t: pd.read_csv(CACHE / f"{t}_usable_clocks.csv").iloc[:, 0].tolist()
          for t in ("GSE132203", "GSE55763")}
section("CHECAGEM 1 — CACHE")
D, N = {}, {}
for t in BLOOD:
    a = pd.read_csv(CACHE / f"{t}_ages.csv", index_col=0); a = a[a.chrono.notna()]
    D[t] = dict(chrono=a.chrono.to_numpy(float), y={c: a[c].to_numpy(float) for c in AGE},
                C12=pd.read_csv(CACHE / f"{t}_comp12.csv", index_col=0).loc[a.index, TYPES12].to_numpy(),
                C6=pd.read_csv(CACHE / f"{t}_comp6.csv", index_col=0).loc[a.index, TYPES6].to_numpy())
    N[t] = len(a)
print("  " + ", ".join(f"{t} {N[t]}" for t in BLOOD))
J = [TYPES12.index(x) for x in BIG4]


def quad(C):
    """composition with squared terms for the four largest types appended"""
    return np.column_stack([C[:, :-1], C[:, J] ** 2])


section("A — CURVATURA DENTRO DA COORTE")
rows = []
for t in BLOOD:
    d = D[t]; n = N[t]
    X0 = np.column_stack([np.ones(n), d["chrono"], d["C12"][:, :-1]])
    X1 = np.column_stack([X0, d["C12"][:, J] ** 2])
    for c in AGE:
        if t in USABLE and c not in USABLE[t]:
            continue
        y = d["y"][c]
        r0 = 1 - ((y - X0 @ np.linalg.lstsq(X0, y, rcond=None)[0]) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        r1 = 1 - ((y - X1 @ np.linalg.lstsq(X1, y, rcond=None)[0]) ** 2).sum() / ((y - y.mean()) ** 2).sum()
        df1, df2 = X1.shape[1] - X0.shape[1], n - X1.shape[1]
        F = ((r1 - r0) / df1) / ((1 - r1) / df2)
        rows.append(dict(cohort=t, clock=c, extra_r2=r1 - r0, F=F,
                         p=float(stats.f.sf(F, df1, df2)), n=n))
A = pd.DataFrame(rows); A.to_csv(RES / "linearity.csv", index=False)
print(f"  {'coorte':<11}{'relogio':<13}{'R2 extra':>10}{'F':>8}{'p':>11}")
for r in A.itertuples():
    print(f"  {r.cohort:<11}{r.clock:<13}{r.extra_r2:>10.4f}{r.F:>8.2f}{r.p:>11.2g}")
c2 = A.extra_r2.median() < 0.01
print(f"\n  2. R2 extra mediano {A.extra_r2.median():.4f} (barra < 0,01) -> {'ok' if c2 else 'FALHOU'}"
      f" | significativo em {(A.p < 0.05).sum()} de {len(A)} celulas")

section("B — A CORRECAO QUADRATICA TRANSPORTA MELHOR?")


def fit_apply(src, dst, c, n, rng, kind):
    s, d = D[src], D[dst]
    idx = (np.arange(N[src]) if n >= N[src] else stratified_draw(s["chrono"], n, rng))
    Cf = s["C12"][idx] if kind == "lin" else quad(s["C12"][idx])
    Cd = d["C12"][:, :-1] if kind == "lin" else quad(d["C12"])
    Cf_use = Cf if kind == "quad" else Cf[:, :-1]
    age = s["chrono"][idx]
    Z = partial_out(age, Cf_use); yt = partial_out(age, s["y"][c][idx].reshape(-1, 1)).ravel()
    mu, sd = Z.mean(0), Z.std(0) + 1e-12
    Zs = (Z - mu) / sd
    U, sv, Vt = np.linalg.svd(Zs, full_matrices=False)
    keep = sv > np.finfo(float).eps * max(len(idx), Zs.shape[1]) * sv.max()
    b = (Vt.T @ np.where(keep, (U.T @ yt) / sv, 0.0)) / sd
    yc = d["y"][c] - (Cd - Cf_use.mean(0)) @ b
    nul = NUL[(dst, c)]
    ia, _ = inc_and_base(yc, d["chrono"], d["C6"])
    return (ia - nul) / (1 - BEF[(dst, c)][1]) - BEF[(dst, c)][0]


NUL, BEF = {}, {}
for t in BLOOD:
    for c in AGE:
        NUL[(t, c)] = null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["C6"], n_perm=400)
        ib, base = inc_and_base(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        BEF[(t, c)] = ((ib - NUL[(t, c)]) / (1 - base), base)
cells = []
for src, dst in permutations(BLOOD, 2):
    n = min(N[src], N[dst])
    for c in AGE:
        if any(t in USABLE and c not in USABLE[t] for t in (src, dst)):
            continue
        res = {}
        for kind in ("lin", "quad"):
            v = [fit_apply(src, dst, c, n, np.random.default_rng(4500 + 7 * r + n), kind)
                 for r in range(1 if n >= N[src] else 20)]
            res[kind] = float(np.median(v))
        cells.append(dict(src=src, dst=dst, n=n, clock=c, **res))
B = pd.DataFrame(cells); B.to_csv(RES / "linearity_transport.csv", index=False)
kl, kq = int((B.lin > 0).sum()), int((B.quad > 0).sum())
c3 = kq >= kl
print(f"  celulas: {len(B)}; nocivas linear {kl}, quadratica {kq}")
print(f"  mediana: linear {B.lin.median():+.1%}, quadratica {B.quad.median():+.1%}")
print(f"  3. quadratica nao melhora (barra: nocivas >= {kl}) -> {'ok' if c3 else 'FALHOU'}")
section("FECHAMENTO")
print(f"  2. curvatura pequena: {'ok' if c2 else 'FALHOU'}")
print(f"  3. nao e o mecanismo: {'ok' if c3 else 'FALHOU'}")
