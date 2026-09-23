#!/usr/bin/env python3
"""
Stage 19 — the second moment, and whether the stage 18 fix reaches the anchor.

Stage 17 killed sample size as a sufficient explanation. Stage 18 killed
conditioning outright — GSE61151, the cohort stage 15 fitted on, is BETTER
conditioned (kappa 40.8) than every GSE40279 subsample including the full 656.
Two mechanisms down, and stage 15's reversed direction still unexplained.

WHAT THE ALGEBRA RULES OUT BEFORE ANY DATA IS TOUCHED, and this is why this
stage is not the one that was planned. The obvious next candidate was covariate
shift: the correction is a linear term extrapolated outside the composition
range it was fitted in, which is exactly what stage 12 demonstrated on purified
cells when neutrophils went from 50-75% to 100% and Horvath 2013's displacement
went from 5.52 to 13.67 years.

But the correction applied here is

    y_corrected = y - (C_test[:, :-1] - c_fit_mean[:-1]) @ b

and c_fit_mean is a constant vector. It contributes a constant to every sample,
and a constant changes no R-squared. **The difference in composition MEANS
between fitting and test cohort is therefore invisible to this metric.** Stage
12's purified-cell result was a displacement measurement, which is sensitive to
the mean; everything from stage 15 onward is a variance-share measurement, which
is not. Carrying the mechanism across without checking that would have been an
error of the same family as the stage 7 denominator.

What is left that can matter is the SECOND moment. The coefficient vector b is
the solution to a least-squares problem shaped by the fitting cohort's
composition covariance, and it is applied to a cohort with a different one. If
the two covariance structures disagree, b is the right answer to the wrong
geometry, and the damage should track that disagreement rather than any mean
shift.

THE DESIGN — three questions, all off the stage 18 cache, no series reloaded.

  Q1  Does stage 15's reversed direction reproduce under this machinery? The
      whole stage depends on it: stage 15 measured +4.9 points fitting on
      GSE61151 and testing on GSE40279, with a six-type panel for measurement
      and a twelve-type panel for fitting. If the current pipeline does not
      reproduce that, nothing below is about stage 15.

  Q2  Does the stage 18 fix reach it? Ridge at the penalty stage 18 found, and
      at the penalty cross-validation picks on GSE61151 itself. This is the
      question with practical consequences: a fix that works on subsamples of
      one cohort but not on the real cohort that failed is not a fix.

  Q3  Does covariance disagreement predict the damage? Measured two ways, on
      the composition block the correction actually uses:
        - Frobenius distance between the two correlation matrices
        - largest principal angle between their leading three eigenvector
          subspaces, which asks whether the dominant directions of composition
          variation are even the same directions
      Tested across the n-grid WITHIN each n x clock x cohort cell, because
      stage 18 was bitten by pooling clocks and reading a between-clock pattern
      as a within-stratum one.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. the cache must be the stage 18 cache: the four cohorts present, the panel
     numbers matching what stages 15, 17 and 18 all recorded (r = 0.789,
     MAE = 0.027). Hard stop, because a silently rebuilt cache would make every
     number here incomparable with the three stages it is being read against.
  2. REPRODUCTION OF STAGE 15. Fitting on GSE61151 and testing on GSE40279
     with OLS must come out positive — the correction must make things worse —
     for the three clocks external to GSE40279. Stage 15 reported +5.9, +4.9
     and +1.6 points for Horvath 2013, Horvath 2018 and Levine 2018 on its own
     scale. The scales are not identical (stage 15 quoted shares of the age
     residual before the null subtraction that stages 17-18 use), so the check
     is on SIGN and ORDER OF MAGNITUDE, not on the digits:
        median delta over the three clocks > 0, and at least two of three
        positive.
     If this fails, stage 15's reversed direction is not reproducible and that
     is the stage's result — the transport finding would then rest on a
     configuration nothing else can recover. Hard stop on interpretation, not
     on the run.
  3. the in-sample floor, now in the OTHER direction: a correction fitted on
     GSE61151 must drive GSE61151's own composition term down. If it does not,
     the fit is broken there and the transport says nothing.
  4. THE DEGENERATE CORNER AGAIN. Any ridge penalty reported as fixing the
     reversed direction must not be the shrink-to-nothing corner. Same
     two-sided rule as stage 18: at the penalty in question, the correction
     must beat OLS on GSE40279 AND still remove composition signal in a
     direction where OLS works — here, GSE40279 -> GSE61151 at full n, keeping
     at least 70% of the OLS benefit stage 18 measured there.
  5. Q3'S DECOUPLING. Covariance distance and n are correlated through the
     subsample, so the correlation is taken within each n x clock x cohort
     cell, never pooled. Pre-specified: covariance disagreement is the
     mechanism only if the median within-cell rho is positive and positive in
     more than 60% of cells. Stage 18's conditioning failed this at 48%, which
     is the coin-flip line this has to clear.
  6. the declared monocyte bias from stage 13 carries.

Usage:  .venv/bin/python analysis/19_covariance.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.linalg import subspace_angles

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_extended import TYPES12
from model.clocks import CLOCKS
from model.deconvolution import TYPES as TYPES6

OUT = ROOT / "results"; CACHE = OUT / "cache"
RNG = np.random.default_rng(20260923)
N_REPS = 30
GRID = [40, 60, 80, 120, 160, 184, 240, 320, 480]
ALPHAS = [0.0, 0.1, 0.3, 1.0, 3.0, 10.0]
VERDICT = ["Horvath2013", "Levine2018", "Horvath2018"]
COHORTS = ["GSE40279", "GSE61151", "GSE50660", "GSE42861"]
# clocks that have seen a cohort, so no verdict is taken from an in-sample cell
TRAINED_ON = {"Hannum2013": "GSE40279"}


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


def r2(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return 1 - (y - X @ beta).var() / y.var()


def inc_and_base(y, a, C):
    Xa = np.column_stack([np.ones(len(a)), a])
    base = r2(Xa, y)
    return r2(np.column_stack([Xa, C[:, :-1]]), y) - base, base


def null_mean(y, a, C, n_perm=1000, rng=RNG):
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


def ridge_path(age, C, y, alphas):
    """Composition coefficients over a scale-free penalty grid, plus the
    penalty leave-one-out would pick. Identical to stage 18's, which check 3
    there verified reproduces stage 17's OLS exactly at alpha = 0."""
    Ct = partial_out(age, C[:, :-1])
    yt = partial_out(age, y.reshape(-1, 1)).ravel()
    mu, sd = Ct.mean(axis=0), Ct.std(axis=0) + 1e-12
    Z = (Ct - mu) / sd
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    scale = float((s ** 2).mean())
    Uty = U.T @ yt
    out, loo = {}, {}
    for a in alphas:
        lam = a * scale
        out[a] = (Vt.T @ (s * Uty / (s ** 2 + lam))) / sd
        h = (U * (s ** 2 / (s ** 2 + lam))) @ U.T
        lev = np.clip(np.diag(h), 0, 1 - 1e-9)
        loo[a] = float(np.mean(((yt - h @ yt) / (1 - lev)) ** 2))
    return out, min(loo, key=loo.get)


def cov_distance(age_a, Ca, age_b, Cb, k=3):
    """
    How far apart are the two composition geometries?

    Both blocks are residualised on intercept and age first, because that is
    the space the coefficients live in, and correlation rather than covariance
    is compared so that a pure per-type scale difference does not masquerade as
    a structural one.
    """
    A = partial_out(age_a, Ca[:, :-1]); B = partial_out(age_b, Cb[:, :-1])
    Ra = np.corrcoef(A, rowvar=False); Rb = np.corrcoef(B, rowvar=False)
    Ra = np.nan_to_num(Ra); Rb = np.nan_to_num(Rb)
    frob = float(np.linalg.norm(Ra - Rb, "fro"))
    va = np.linalg.eigh(Ra)[1][:, -k:]
    vb = np.linalg.eigh(Rb)[1][:, -k:]
    ang = float(np.max(subspace_angles(va, vb)) * 180 / np.pi)
    return frob, ang


# ------------------------------------------------------------------ data ----
need = [CACHE / f"{t}_{w}.csv" for t in COHORTS
        for w in ("comp12", "comp6", "ages")] + [CACHE / "panel.csv"]
missing = [f.name for f in need if not f.exists()]
if missing:
    sys.exit(f"cache ausente ({missing[:3]}...). Rode analysis/18_conditioning.py antes.")

panel = pd.read_csv(CACHE / "panel.csv", index_col=0).iloc[:, 0]
D = {}
for tag in COHORTS:
    a = pd.read_csv(CACHE / f"{tag}_ages.csv", index_col=0)
    k = a.chrono.notna().to_numpy()
    D[tag] = dict(chrono=a.chrono[k].to_numpy(),
                  y={c: a[c][k].to_numpy() for c in CLOCKS},
                  C12=pd.read_csv(CACHE / f"{tag}_comp12.csv", index_col=0)
                        .loc[a.index[k], TYPES12].to_numpy(),
                  C6=pd.read_csv(CACHE / f"{tag}_comp6.csv", index_col=0)
                       .loc[a.index[k], TYPES6].to_numpy())


section("CHECAGEM 1 — O CACHE E O MESMO DAS ETAPAS 15, 17 E 18")
c1 = abs(panel.panel_r - 0.789) < 5e-3 and abs(panel.panel_mae - 0.027) < 5e-3
print(f"  painel: r = {panel.panel_r:.3f} (esperado 0.789) | "
      f"MAE = {panel.panel_mae:.3f} (esperado 0.027)")
print(f"  coortes: " + ", ".join(f"{t} n={len(D[t]['chrono'])}" for t in COHORTS))
print(f"  -> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("\n  parando: cache diferente torna tudo incomparavel com as etapas anteriores.")
print("  (vies de monocito da etapa 13 carregado; nenhum coeficiente por tipo e lido)")


def transport(src, dst, idx=None, alphas=ALPHAS, clocks=None):
    """Fit on src (optionally a subsample), measure on dst. Returns one row per
    clock x alpha. Measurement is always the cross-panel cell: fit with twelve
    types, measure the residual with six."""
    s, d = D[src], D[dst]
    idx = np.arange(len(s["chrono"])) if idx is None else idx
    a_fit, C_fit = s["chrono"][idx], s["C12"][idx]
    cbar = C_fit.mean(axis=0)
    rows = []
    for c in (clocks or CLOCKS):
        if TRAINED_ON.get(c) == dst:
            continue
        coefs, a_loo = ridge_path(a_fit, C_fit, s["y"][c][idx], alphas)
        ib, base = inc_and_base(d["y"][c], d["chrono"], d["C6"])
        nb = NULL[(dst, c)]
        before = (ib - nb) / (1 - base)
        for alpha, b in coefs.items():
            yc = d["y"][c] - (d["C12"][:, :-1] - cbar[:-1]) @ b
            ia, _ = inc_and_base(yc, d["chrono"], d["C6"])
            after = (ia - nb) / (1 - base)
            rows.append(dict(src=src, dst=dst, clock=c, alpha=alpha,
                             alpha_loo=a_loo, n=len(idx), before=before,
                             after=after, delta=after - before,
                             norm=float(np.linalg.norm(b))))
    return rows


print("\n  calculando os nulos de permutacao (uma vez por coorte x relogio) ...",
      flush=True)
NULL = {(t, c): null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        for t in COHORTS for c in CLOCKS}


section("CHECAGEM 2 — A DIRECAO INVERTIDA DA ETAPA 15 REPRODUZ?")
rev = pd.DataFrame(transport("GSE61151", "GSE40279"))
ro = rev[rev.alpha == 0.0]
print(f"  ajustado em GSE61151 (n={len(D['GSE61151']['chrono'])}), "
      f"testado em GSE40279 (n={len(D['GSE40279']['chrono'])}), OLS\n")
print(f"  {'relogio':<16}{'antes':>10}{'depois':>10}{'delta':>10}")
for _, r in ro.iterrows():
    print(f"  {r.clock:<16}{r.before:>10.1%}{r.after:>10.1%}{r.delta:>+10.1%}"
          + ("  pior" if r.delta > 0 else ""))
med = ro.delta.median()
npos = int((ro.delta > 0).sum())
c2 = (med > 0) and (npos >= 2)
print(f"\n  mediana {med:+.1%}, positivos {npos} de {len(ro)}")
print(f"  regra escrita antes: mediana > 0 e >= 2 de 3 positivos "
      f"-> {'ok' if c2 else 'FALHOU'}")
if not c2:
    print("  -> a direcao invertida da etapa 15 NAO reproduz sob esta maquinaria.\n"
          "     Isso e o resultado da etapa: o achado de transporte dependia de\n"
          "     uma configuracao que nada mais recupera. Nao remendar.")


section("CHECAGEM 3 — PISO EM CASA NA OUTRA DIRECAO")
s = D["GSE61151"]
c3 = True
for c in VERDICT:
    coefs, _ = ridge_path(s["chrono"], s["C12"], s["y"][c], [0.0])
    cbar = s["C12"].mean(axis=0)
    yc = s["y"][c] - (s["C12"][:, :-1] - cbar[:-1]) @ coefs[0.0]
    b0, _ = inc_and_base(s["y"][c], s["chrono"], s["C12"])
    a0, _ = inc_and_base(yc, s["chrono"], s["C12"])
    ok = a0 < b0
    c3 &= ok
    print(f"  {c:<16}{b0:>9.4f} -> {a0:>9.4f}  {'ok' if ok else 'FALHOU'}")
print(f"  -> {'ok' if c3 else 'FALHOU'}")


section("Q2 — O CONSERTO DA ETAPA 18 ALCANCA A DIRECAO INVERTIDA?")
piv = rev.pivot_table(index="alpha", columns="clock", values="delta")
print(f"  {'alpha':>7}" + "".join(f"{c:>15}" for c in piv.columns)
      + f"{'mediana':>11}")
for a in ALPHAS:
    print(f"  {a:>7}" + "".join(f"{piv.loc[a, c]:>+15.1%}" for c in piv.columns)
          + f"{piv.loc[a].median():>+11.1%}")
a_loo_rev = float(rev.alpha_loo.median())
print(f"\n  alpha que a validacao cruzada escolhe no GSE61151: {a_loo_rev:g}")
print(f"  norma L2 mediana por alpha: "
      + "  ".join(f"a={a}: {rev[rev.alpha == a]['norm'].median():.0f}" for a in ALPHAS))

# check 4 — the fix must not be the shrink-to-nothing corner
fwd = pd.DataFrame(transport("GSE40279", "GSE61151"))
fwd_ols = fwd[(fwd.alpha == 0.0) & fwd.clock.isin(VERDICT)].delta.median()
print(f"\n  CHECAGEM 4 — esquina degenerada")
print(f"  referencia: GSE40279 -> GSE61151 com OLS e n cheio = {fwd_ols:+.1%}")
print(f"  um alpha so conta se corrigir a direcao invertida E guardar >= 70% "
      f"disso ({0.7 * fwd_ols:+.1%})\n")
fixes = []
for a in ALPHAS:
    if a == 0.0:
        continue
    rv = piv.loc[a].median()
    fv = fwd[(fwd.alpha == a) & fwd.clock.isin(VERDICT)].delta.median()
    i = rv <= 0
    ii = fv <= 0.7 * fwd_ols
    if i and ii:
        fixes.append(a)
    print(f"    alpha={a:<5} invertida {rv:>+8.1%} {'ok ' if i else 'nao'}"
          f"   direta {fv:>+8.1%} {'ok ' if ii else 'nao'}"
          f"   {'<- CONSERTA AS DUAS' if (i and ii) else ''}")
c4 = len(fixes) > 0
print(f"\n  -> {'ok' if c4 else 'FALHOU'}: "
      + (f"alphas: {fixes}" if c4 else
         "nenhum alpha conserta a direcao invertida sem virar a esquina degenerada"))


section("Q3 — A GEOMETRIA DA COMPOSICAO DISCORDA ENTRE AS COORTES?")
print("  distancia entre as matrizes de CORRELACAO de composicao, depois de")
print("  tirar intercepto e idade — o espaco onde os coeficientes vivem.\n")
print(f"  {'par':<26}{'Frobenius':>12}{'maior angulo':>15}")
pairs = [("GSE61151", "GSE40279"), ("GSE40279", "GSE61151"),
         ("GSE40279", "GSE50660"), ("GSE40279", "GSE42861")]
seen, dist = set(), {}
for a, b in pairs:
    if (a, b) in seen or (b, a) in seen:
        continue
    seen.add((a, b))
    f, ang = cov_distance(D[a]["chrono"], D[a]["C12"], D[b]["chrono"], D[b]["C12"])
    dist[(a, b)] = (f, ang)
    print(f"  {a[3:]} <-> {b[3:]:<16}{f:>12.2f}{ang:>14.1f}°")
d15 = dist[("GSE61151", "GSE40279")][0]
others = [v[0] for k, v in dist.items() if k != ("GSE61151", "GSE40279")]
print(f"\n  o par da etapa 15 e o mais distante? "
      f"{d15:.2f} contra {min(others):.2f}-{max(others):.2f} "
      f"-> {'SIM' if d15 > max(others) else 'NAO'}")


section("Q3 — E ISSO PREVE O DANO, COM N FIXO?")
rows = []
for n in GRID:
    for r in range(N_REPS):
        idx = stratified_draw(D["GSE40279"]["chrono"], n,
                              np.random.default_rng(2000 + 13 * r + n))
        for dst in ["GSE61151", "GSE50660", "GSE42861"]:
            f, ang = cov_distance(D["GSE40279"]["chrono"][idx],
                                  D["GSE40279"]["C12"][idx],
                                  D[dst]["chrono"], D[dst]["C12"])
            for row in transport("GSE40279", dst, idx=idx, alphas=[0.0],
                                 clocks=VERDICT):
                row.update(rep=r, frob=f, angle=ang)
                rows.append(row)
    print(f"  n = {n:>4}: pronto", flush=True)

q3 = pd.DataFrame(rows)
q3.to_csv(OUT / "covariance.csv", index=False)

res = []
for (n, c, co), g in q3.groupby(["n", "clock", "dst"]):
    if g.rep.nunique() < 10:
        continue
    rf, _ = spearmanr(g.frob, g.delta)
    ra, _ = spearmanr(g.angle, g.delta)
    res.append(dict(n=n, rho_frob=rf, rho_angle=ra))
R = pd.DataFrame(res)
print(f"\n  correlacao DENTRO de cada celula n x relogio x coorte "
      f"({len(R)} celulas)\n")
print(f"  {'n':>5}{'rho(Frobenius)':>18}{'rho(angulo)':>16}")
for n, g in R.groupby("n"):
    print(f"  {n:>5}{g.rho_frob.median():>18.3f}{g.rho_angle.median():>16.3f}")
share = float((R.rho_frob > 0).mean())
c5 = (R.rho_frob.median() > 0) and (share > 0.60)
print(f"\n  rho(Frobenius) mediano {R.rho_frob.median():+.3f}, "
      f"positivo em {share:.0%} das celulas")
print(f"  regra escrita antes: mediana > 0 e mais de 60% positivas "
      f"-> {'ok' if c5 else 'FALHOU'}")
print(f"  (a linha do cara-ou-coroa e 50%; o condicionamento da etapa 18 deu 48%)")


section("CHECAGENS, FECHAMENTO")
for i, (name, ok) in enumerate([
        ("cache e o mesmo das etapas 15/17/18", c1),
        ("direcao invertida da etapa 15 reproduz", c2),
        ("piso em casa no GSE61151", c3),
        ("existe alpha que conserta as duas direcoes", c4),
        ("geometria da composicao preve o dano com n fixo", c5)], 1):
    print(f"  {i}. {name}: {'ok' if ok else 'FALHOU'}")
print("  6. vies de monocito da etapa 13: carregado")
print(f"\n  saida: results/covariance.csv ({len(q3)} linhas)")
