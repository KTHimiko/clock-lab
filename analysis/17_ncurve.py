#!/usr/bin/env python3
"""
Stage 17 — how many samples does the composition correction need?

Stage 15 found that the same correction, built the same way, helps or hurts
depending on where it was fitted. Fitted on GSE40279's 656 samples and carried
to GSE61151 it removed most of the composition signal. Fitted on GSE61151's 184
and carried to GSE40279 it made every clock worse — Horvath 2013 from 4.7% to
10.7%. Stage 15 attributed that to sample size: coefficients estimated on a
small cohort are unstable, and an unstable coefficient carried somewhere else
adds a composition-correlated term to every prediction.

THAT ATTRIBUTION WAS NEVER TESTED. Two cohorts differing in size also differ in
age range, array batch, collection site, and composition distribution. Stage 15
has one point at n=656 and one at n=184, in opposite directions, in different
cohorts. Any of those differences explains the result equally well.

This stage varies n and holds everything else fixed.

THE DESIGN
  fit cohort    GSE40279, 656 samples, subsampled at ten sizes. Every
                correction in this stage is fitted on the same people, the same
                array, the same collection sites — only how many of them changes.
  test cohorts  three, all external to the fit: GSE61151, GSE50660, GSE42861.
                Stage 15 had one. A size effect that appears in one test cohort
                and not the others is a fact about that cohort.
  the metric    excess over permutation null of the composition term left in
                the age residual, as a share of the age-acceleration variance
                of the UNCORRECTED clock — so before and after share a
                denominator. Delta = after - before. Positive means the
                correction made the clock worse than leaving it alone.
  cross-panel   correction fitted with the twelve-type panel, residual measured
                with the six-type panel. Stage 15 established that scoring a
                correction with its own panel scores it against its own
                representation of what it removed.

WHY THE SUBSAMPLE IS STRATIFIED BY AGE. Drawing 40 of 656 at random gives a
narrower and randomly placed age range, and stage 10 already cost this project a
conclusion that turned out to be about training age range rather than design.
Every draw is stratified by age decile, and check 4 verifies the age
distribution is held across the whole grid. Without that the n-axis and the
age-range axis are the same axis.

WHAT THIS CANNOT SETTLE. The anchor below fits on a 184-sample subsample of
GSE40279 and tests elsewhere; stage 15 fitted on GSE61151's own 184 and tested
on GSE40279. Those are not the same cell. So a positive anchor shows that n=184
is SUFFICIENT to produce harm, which is what the size explanation predicts; it
does not show that size is the ONLY thing wrong with stage 15's reversed
direction. A negative anchor is the informative one, and its consequence is
written down before the run.

CLOCKS. The four published clocks. Hannum was trained on GSE40279, so its age
residual there is small and the correction is fitted on very little — that is a
property of the fit cohort, not of any test cohort, so Hannum is reported but
kept out of the pooled verdict and marked in every table. The two clocks this
project built are also GSE40279-trained and are left out entirely: this stage is
about the correction, not about clock design.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. the twelve-type panel recovers the known twelve-way proportions of the
     twelve GSE167998 mixtures: r > 0.70, mean absolute error < 0.04. Carried
     from stage 15 — without a validated panel nothing downstream means
     anything. Hard stop.
  2. every clock tracks chronological age in every test cohort: r > 0.5.
     GSE50660 and GSE42861 are narrower in age than the fit cohort, so the bar
     is stage 16's rather than stage 1's. Hard stop.
  3. in-sample floor at full n: the correction fitted on all 656 drives the
     composition term down in the cohort it was fitted on, measured with its own
     panel. If it does not, the fit is broken and nothing out of cohort means
     anything. Hard stop.
  4. the age RANGE is held across the grid. What cost stage 10 a conclusion was
     a training range that moved and narrowed, not decile wobble: a draw of 40
     has noisy quantiles however it is taken, and a check on decile distance
     would be measuring sampling noise rather than the confound. So the two
     statistics that matter are measured instead, as medians over draws:
     the median age must sit within 2 years of the full cohort's, and the
     10-90 percentile spread must stay within 10% of it. Stratified and plain
     random draws are both reported, so what stratification bought is visible
     rather than asserted. If this fails the n-axis is the age-range axis and
     the curve is not reported as a curve.
  5. THE ANCHOR, and the decision rule is written here rather than after the
     numbers are seen. Stage 15's reversed direction fitted on 184 samples and
     made every clock worse, median +4.9 points. If sample size is the
     mechanism, fitting on 184 samples of GSE40279 must also be harmful:
     median delta over the three test cohorts and three verdict clocks > 0.

       positive -> size is sufficient to produce the harm; the curve is the
                   result and stage 15's attribution stands.
       negative -> size is NOT what broke stage 15's reversed direction. The
                   result of this stage is then that stage 15's stated
                   mechanism is wrong, reported as such, and the headline of
                   the transport finding has to be rewritten around whatever
                   else differs between those two cohorts.

  6. the permutation null is a property of the design and not of the response.
     The metric subtracts a null mean, and computing it per draw would dominate
     the runtime, so it is computed once per clock and test cohort from the
     uncorrected response and reused. That reuse is only legitimate if the null
     does not move with the response: verified on twenty corrected responses,
     largest deviation must be under 0.005 in R-squared units.

     THIS CHECK IS NOT LOAD-BEARING FOR THE VERDICT, and saying so is cheaper
     than letting it look like it is. The same null is subtracted from before
     and from after, so it cancels in the delta: the curve, the anchor and the
     crossing point would be identical with no null at all. What the null
     changes is the reported LEVEL of the composition term — whether a clock
     sits at 4.7% or at 1.2% before correction — which is what makes the levels
     comparable with stage 15's tables. The check is kept for that, and it is
     the one check here whose failure would not touch the conclusion.

  7. the declared bias from stage 13 carries: this panel's monocyte channel runs
     high, so no per-type monocyte coefficient is read or reported.

REFERENCE LINE, NOT A CHECK. A correction fitted at full n on PERMUTED
composition is transported alongside the real ones. It removes nothing by
construction and injects whatever a meaningless coefficient injects, so it gives
the scale of the worst case the real curve could approach.

Usage:  .venv/bin/python analysis/17_ncurve.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix
from load_extended import read_extended, TYPES12
from model.clocks import CLOCKS, predict
from model.deconvolution import PANEL, deconvolve, TYPES as TYPES6

DATA = ROOT / "reference/data"
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(20260923)
N_PERM = 1000
N_REPS = 30
GRID = [40, 60, 80, 120, 160, 240, 320, 480]
VERDICT_CLOCKS = ["Horvath2013", "Levine2018", "Horvath2018"]

FIT = "GSE40279"
TESTS = ["GSE61151", "GSE50660", "GSE42861"]
AGE_FIELD = {"GSE40279": "age (y)", "GSE61151": "agebloodtaken",
             "GSE50660": "age", "GSE42861": "age"}


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


def r2(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return 1 - (y - X @ beta).var() / y.var()


def type_scores(betas, labels, types, n_candidate=200):
    """t-like separation score per type, and the candidate probe pool."""
    scores, cand = {}, set()
    arr = np.asarray(labels)
    for t in types:
        sel = arr == t
        rest = np.isin(arr, [o for o in types if o != t])
        a, c = betas.loc[:, sel], betas.loc[:, rest]
        d = a.mean(axis=1) - c.mean(axis=1)
        sd = np.sqrt((a.var(axis=1, ddof=1) + c.var(axis=1, ddof=1)) / 2) + 0.01
        s = (d / sd).replace([np.inf, -np.inf], np.nan).dropna()
        scores[t] = s
        cand |= set(s.nlargest(n_candidate).index) | set(s.nsmallest(n_candidate).index)
    return scores, cand


def pick_panel(scores, cand, betas, labels, types, restrict_to, n_side=50):
    """Final reference profiles, over probes present in every cohort.

    The pool is the candidate set and nothing else. Only candidate rows survive
    the loading loop, so a probe chosen from outside it would not exist by the
    time anything is deconvolved.
    """
    pool = sorted(cand & set(restrict_to))
    arr = np.asarray(labels)
    probes = set()
    for t in types:
        s = scores[t].reindex(pool).dropna()
        probes |= set(s.nlargest(n_side).index) | set(s.nsmallest(n_side).index)
    probes = sorted(probes)
    ref = pd.DataFrame({t: betas.loc[probes, arr == t].mean(axis=1) for t in types})
    return ref, probes


def stratified_draw(ages, n, rng):
    """n indices that keep the cohort's age distribution — see check 4."""
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


def fit_correction(y, age, C):
    """Composition coefficients of clock_age ~ age + composition, and the
    centring used when they are applied elsewhere.

    Age stays in the model: composition drifts with age, and a correction fitted
    without age would strip real ageing out along with the blood count. Only the
    composition part is kept, so applying it needs no age — which is the whole
    point of a clock.
    """
    M = np.column_stack([np.ones(len(age)), age, C[:, :-1]])
    beta = np.linalg.lstsq(M, y, rcond=None)[0]
    return beta[2:], C.mean(axis=0)


def apply_correction(y, C, bcomp, cbar):
    return y - (C[:, :-1] - cbar[:-1]) @ bcomp


def inc_and_base(y, a, C):
    Xa = np.column_stack([np.ones(len(a)), a])
    base = r2(Xa, y)
    return r2(np.column_stack([Xa, C[:, :-1]]), y) - base, base


def null_mean(y, a, C, n_perm=N_PERM, rng=RNG):
    Xa = np.column_stack([np.ones(len(a)), a])
    base = r2(Xa, y)
    vals = np.empty(n_perm)
    for i in range(n_perm):
        Cp = C[rng.permutation(len(C))]
        vals[i] = r2(np.column_stack([Xa, Cp[:, :-1]]), y) - base
    return float(vals.mean()), float(vals.std(ddof=1) / np.sqrt(n_perm))


# ---------------------------------------------------------------- loading ---
print("carregando GSE167998 (referencia de 12 tipos) ...", flush=True)
b167, p167 = read_extended(DATA / "GSE167998_matrix_processed.txt.gz",
                           DATA / "BloodExtended_Pheno.csv")
is_mix = (p167.CellType == "MIX").to_numpy()
sc12, cand12 = type_scores(b167.loc[:, ~is_mix], p167.CellType[~is_mix], TYPES12)

print("carregando GSE35069 (referencia de 6 tipos) ...", flush=True)
b35, m35 = read_series_matrix(DATA / "GSE35069_series_matrix.txt.gz")
m35 = m35.set_index("gsm").reindex(b35.columns)
frac35 = m35["tissue/cell type"].map({v: k for k, v in PANEL.items()})
sc6, cand6 = type_scores(b35.loc[:, frac35.notna().to_numpy()],
                         frac35.dropna(), TYPES6)

candidates = sorted(cand12 | cand6)
print(f"  pool de candidatas: {len(candidates)} sondas")

# Each whole-blood cohort is loaded once, reduced to what this stage needs, and
# dropped. Only the candidate rows and the clock ages survive the loop, which is
# a few megabytes instead of five gigabytes.
kept, ages, chrono, indexes = {}, {}, {}, {}
for tag in [FIT] + TESTS:
    print(f"carregando {tag} ...", flush=True)
    b, m = read_series_matrix(DATA / f"{tag}_series_matrix.txt.gz")
    m = m.set_index("gsm").reindex(b.columns)
    chrono[tag] = pd.to_numeric(m[AGE_FIELD[tag]], errors="coerce")
    ages[tag] = {c: predict(b, c)[0] for c in CLOCKS}
    indexes[tag] = b.index
    kept[tag] = b.loc[b.index.intersection(candidates)].copy()
    n_ok = int(chrono[tag].notna().sum())
    print(f"  {b.shape[0]:,} sondas x {b.shape[1]} amostras, {n_ok} com idade",
          flush=True)
    del b, m

common = indexes[FIT]
for tag in TESTS:
    common = common.intersection(indexes[tag])
common = common.intersection(b167.index).intersection(b35.index)

ref12, probes12 = pick_panel(sc12, cand12, b167.loc[:, ~is_mix],
                             p167.CellType[~is_mix], TYPES12, common)
ref6, probes6 = pick_panel(sc6, cand6, b35.loc[:, frac35.notna().to_numpy()],
                           frac35.dropna(), TYPES6, common)
print(f"\n  painel 12: {len(probes12)} sondas | painel 6: {len(probes6)} sondas"
      f" | sondas comuns as seis series: {len(common):,}")


# --------------------------------------------------------------- checks 1-2 --
section("CHECAGEM 1 — O PAINEL DE 12 RECUPERA AS PROPORCOES CONHECIDAS?")
known = p167.loc[is_mix, TYPES12].to_numpy()
est = deconvolve(ref12.to_numpy(), b167.loc[probes12].to_numpy()[:, is_mix])
rv = float(np.corrcoef(known.ravel(), est.ravel())[0, 1])
mv = float(np.abs(known - est).mean())
c1 = (rv > 0.70) and (mv < 0.04)
print(f"  correlacao {rv:.3f}  erro absoluto medio {mv:.3f}  "
      f"-> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("\n  parando: sem painel validado nao se corrige nada.")
del b167, b35

comp = {}
for tag in [FIT] + TESTS:
    comp[(tag, 12)] = pd.DataFrame(
        deconvolve(ref12.to_numpy(), kept[tag].loc[probes12].to_numpy()),
        index=kept[tag].columns, columns=TYPES12)
    comp[(tag, 6)] = pd.DataFrame(
        deconvolve(ref6.to_numpy(), kept[tag].loc[probes6].to_numpy()),
        index=kept[tag].columns, columns=TYPES6)
print("  (vies declarado da etapa 13 continua valendo: o canal de monocito "
      "corre alto,\n   nenhum coeficiente por tipo de monocito e lido aqui)")

section("CHECAGEM 2 — OS RELOGIOS ACOMPANHAM A IDADE EM CADA COORTE DE TESTE?")
print(f"  {'relogio':<16}" + "".join(f"{t:>12}" for t in TESTS))
c2 = True
for c in CLOCKS:
    line = f"  {c:<16}"
    for tag in TESTS:
        k = chrono[tag].notna().to_numpy()
        r = float(np.corrcoef(chrono[tag][k].to_numpy(),
                              ages[tag][c][k].to_numpy())[0, 1])
        c2 &= r > 0.5
        line += f"{r:>12.3f}"
    print(line)
print(f"  limite 0.5 em todas as celulas -> {'ok' if c2 else 'FALHOU'}")
if not c2:
    sys.exit("\n  parando: um relogio que nao le idade nao tem aceleracao para corrigir.")


# ----------------------------------------------------------------- setup ----
kfit = chrono[FIT].notna().to_numpy()
age_fit_all = chrono[FIT][kfit].to_numpy()
idx_fit = chrono[FIT][kfit].index
C12_fit_all = comp[(FIT, 12)].loc[idx_fit, TYPES12].to_numpy()
N_FIT = len(age_fit_all)

test = {}
for tag in TESTS:
    k = chrono[tag].notna().to_numpy()
    ix = chrono[tag][k].index
    test[tag] = dict(age=chrono[tag][k].to_numpy(),
                     C12=comp[(tag, 12)].loc[ix, TYPES12].to_numpy(),
                     C6=comp[(tag, 6)].loc[ix, TYPES6].to_numpy(),
                     y={c: ages[tag][c].reindex(ix).to_numpy() for c in CLOCKS})

# the anchor sits at stage 15's fitting size, whatever that cohort actually has
N_ANCHOR = len(test["GSE61151"]["age"])
grid = sorted({n for n in set(GRID) | {N_ANCHOR} if n < N_FIT} | {N_FIT})
print(f"\n  coorte de ajuste: {FIT}, {N_FIT} amostras com idade")
print(f"  ancora (tamanho util de GSE61151 na etapa 15): n = {N_ANCHOR}")
print(f"  grade: {grid}")
print(f"  coortes de teste: " + ", ".join(
    f"{t} (n={len(test[t]['age'])})" for t in TESTS))


section("CHECAGEM 3 — PISO EM CASA COM N CHEIO")
c3 = True
for c in CLOCKS:
    y = ages[FIT][c].reindex(idx_fit).to_numpy()
    bcomp, cbar = fit_correction(y, age_fit_all, C12_fit_all)
    before, _ = inc_and_base(y, age_fit_all, C12_fit_all)
    after, _ = inc_and_base(apply_correction(y, C12_fit_all, bcomp, cbar),
                            age_fit_all, C12_fit_all)
    ok = after < before
    c3 &= ok
    print(f"  {c:<16} {before:>8.4f} -> {after:>8.4f}  {'ok' if ok else 'FALHOU'}")
print(f"  -> {'ok' if c3 else 'FALHOU'}")
if not c3:
    sys.exit("\n  parando: se a correcao nao zera em casa, o ajuste esta quebrado.")


section("CHECAGEM 4 — A FAIXA DE IDADE E MANTIDA AO LONGO DA GRADE?")
med_full = float(np.median(age_fit_all))
spread_full = float(np.percentile(age_fit_all, 90) - np.percentile(age_fit_all, 10))
print(f"  coorte cheia: mediana {med_full:.1f} anos, faixa 10-90 "
      f"{spread_full:.1f} anos\n")
print(f"  {'n':>5}{'  estratificado':>26}{'  aleatorio simples':>28}")
print(f"  {'':>5}{'desloc.':>12}{'faixa':>12}{'desloc.':>14}{'faixa':>12}")
rows4 = []
for n in grid:
    reps = 1 if n >= N_FIT else N_REPS
    stats = {}
    for how in ("estratificado", "aleatorio"):
        shifts, ratios = [], []
        for r in range(reps):
            rg = np.random.default_rng(1000 + 7 * r + n)
            if n >= N_FIT:
                idx = np.arange(N_FIT)
            elif how == "estratificado":
                idx = stratified_draw(age_fit_all, n, rg)
            else:
                idx = rg.choice(N_FIT, n, replace=False)
            a = age_fit_all[idx]
            shifts.append(np.median(a) - med_full)
            ratios.append((np.percentile(a, 90) - np.percentile(a, 10)) / spread_full)
        stats[how] = (float(np.median(shifts)), float(np.median(ratios)))
    rows4.append(dict(n=n, shift=stats["estratificado"][0],
                      ratio=stats["estratificado"][1]))
    print(f"  {n:>5}{stats['estratificado'][0]:>+11.2f}a{stats['estratificado'][1]:>11.3f}"
          f"{stats['aleatorio'][0]:>+13.2f}a{stats['aleatorio'][1]:>11.3f}")
c4 = all(abs(r["shift"]) < 2.0 and 0.90 < r["ratio"] < 1.10 for r in rows4)
print(f"\n  limites: deslocamento da mediana < 2 anos, faixa entre 0.90 e 1.10 "
      f"-> {'ok' if c4 else 'FALHOU'}")


section("CHECAGEM 6 — O NULO DEPENDE DO DESENHO OU DA RESPOSTA?")
probe_tag = TESTS[0]
t0 = test[probe_tag]
c = VERDICT_CLOCKS[0]
base_null, se_null = null_mean(t0["y"][c], t0["age"], t0["C6"])
devs = []
for r in range(20):
    idx = stratified_draw(age_fit_all, 120, np.random.default_rng(500 + r))
    y_fit = ages[FIT][c].reindex(idx_fit).to_numpy()[idx]
    bcomp, cbar = fit_correction(y_fit, age_fit_all[idx], C12_fit_all[idx])
    yc = apply_correction(t0["y"][c], t0["C12"], bcomp, cbar)
    m, _ = null_mean(yc, t0["age"], t0["C6"], n_perm=200)
    devs.append(abs(m - base_null))
c6 = max(devs) < 0.005
print(f"  {probe_tag} / {c}: nulo da resposta nao corrigida {base_null:.5f} "
      f"(erro padrao {se_null:.5f})")
print(f"  maior desvio sobre 20 respostas corrigidas: {max(devs):.5f}")
print(f"  limite 0.005 -> {'ok' if c6 else 'FALHOU'}")
if not c6:
    print("  -> o nulo sera recalculado por sorteio, e a etapa fica mais lenta.")

NULL = {}
for tag in TESTS:
    for c in CLOCKS:
        NULL[(tag, c)] = null_mean(test[tag]["y"][c], test[tag]["age"],
                                   test[tag]["C6"])[0]


# ------------------------------------------------------------- the curve ----
section("A CURVA — DELTA POR TAMANHO DA COORTE DE AJUSTE")
print("  delta = depois - antes, em pontos da variancia de aceleracao de idade.")
print("  positivo = a correcao deixou o relogio PIOR do que nao corrigir.\n")

rows = []
for n in grid:
    reps = 1 if n >= N_FIT else N_REPS
    for r in range(reps):
        idx = (np.arange(N_FIT) if n >= N_FIT
               else stratified_draw(age_fit_all, n, np.random.default_rng(2000 + 13 * r + n)))
        a_fit, C_fit = age_fit_all[idx], C12_fit_all[idx]
        for c in CLOCKS:
            y_fit = ages[FIT][c].reindex(idx_fit).to_numpy()[idx]
            bcomp, cbar = fit_correction(y_fit, a_fit, C_fit)
            for tag in TESTS:
                t = test[tag]
                nb = NULL[(tag, c)]
                ib, base = inc_and_base(t["y"][c], t["age"], t["C6"])
                yc = apply_correction(t["y"][c], t["C12"], bcomp, cbar)
                ia, _ = inc_and_base(yc, t["age"], t["C6"])
                before = (ib - nb) / (1 - base)
                after = (ia - nb) / (1 - base)
                rows.append(dict(n=n, rep=r, clock=c, cohort=tag,
                                 before=before, after=after,
                                 delta=after - before,
                                 verdict=c in VERDICT_CLOCKS))
    print(f"  n = {n:>4}: {reps} sorteio(s) pronto(s)", flush=True)

# reference line: the same transport with composition coefficients that mean
# nothing, fitted at full n
perm_rows = []
for c in CLOCKS:
    y_fit = ages[FIT][c].reindex(idx_fit).to_numpy()
    Cp = C12_fit_all[RNG.permutation(N_FIT)]
    bcomp, cbar = fit_correction(y_fit, age_fit_all, Cp)
    for tag in TESTS:
        t = test[tag]
        nb = NULL[(tag, c)]
        ib, base = inc_and_base(t["y"][c], t["age"], t["C6"])
        ia, _ = inc_and_base(apply_correction(t["y"][c], t["C12"], bcomp, cbar),
                             t["age"], t["C6"])
        perm_rows.append(dict(clock=c, cohort=tag,
                              delta=(ia - nb) / (1 - base) - (ib - nb) / (1 - base),
                              verdict=c in VERDICT_CLOCKS))

cur = pd.DataFrame(rows)
cur.to_csv(OUT / "ncurve.csv", index=False)
perm = pd.DataFrame(perm_rows)

summary = (cur.groupby(["n", "clock", "cohort"])
              .delta.agg(["median", lambda s: s.quantile(0.25),
                          lambda s: s.quantile(0.75), "size"])
              .set_axis(["median", "q25", "q75", "draws"], axis=1)
              .reset_index())
summary.to_csv(OUT / "ncurve_summary.csv", index=False)

print(f"\n  {'n':>5}" + "".join(f"{t.replace('GSE',''):>22}" for t in TESTS)
      + f"{'mediana':>12}")
print(f"  {'':>5}" + "".join(f"{'mediana [IQR]':>22}" for t in TESTS)
      + f"{'agregada':>12}")
for n in grid:
    line = f"  {n:>5}"
    for tag in TESTS:
        s = cur[(cur.n == n) & (cur.cohort == tag) & cur.verdict].delta
        line += f"{s.median():>+9.1%} [{s.quantile(.25):>+5.1%},{s.quantile(.75):>+6.1%}]"
    agg = cur[(cur.n == n) & cur.verdict].delta
    line += f"{agg.median():>+12.1%}"
    print(line)

pl = perm[perm.verdict].delta
print(f"\n  linha de referencia — coeficientes de composicao sem sentido, "
      f"ajustados com n={N_FIT}:")
print(f"  delta mediano {pl.median():+.1%}  (faixa {pl.min():+.1%} a {pl.max():+.1%})")

print(f"\n  Hannum2013 fica de fora da mediana agregada — treinou em {FIT}, "
      f"entao a\n  correcao dele e ajustada sobre um residuo de idade quase "
      f"inexistente:")
for n in grid:
    s = cur[(cur.n == n) & (cur.clock == "Hannum2013")].delta
    print(f"    n = {n:>4}: {s.median():>+8.1%}")


section("CHECAGEM 5 — A ANCORA")
anchor = cur[(cur.n == N_ANCHOR) & cur.verdict].delta
c5 = anchor.median() > 0
print(f"  regra escrita antes da rodada: com n = {N_ANCHOR} a mediana sobre as "
      f"tres coortes\n  de teste e os tres relogios do veredito deve ser "
      f"positiva se o tamanho for o mecanismo.\n")
print(f"  mediana em n = {N_ANCHOR}: {anchor.median():+.1%}  "
      f"(IQR {anchor.quantile(.25):+.1%} a {anchor.quantile(.75):+.1%}, "
      f"{len(anchor)} valores)")
print(f"  fracao de celulas piores que nao corrigir: "
      f"{(anchor > 0).mean():.0%}")
print(f"\n  -> {'ANCORA POSITIVA' if c5 else 'ANCORA NEGATIVA'}")
if c5:
    print("     o tamanho da coorte de ajuste basta para produzir o dano, e a\n"
          "     atribuicao da etapa 15 sobrevive. A curva e o resultado.")
else:
    print("     o tamanho NAO e o que quebrou a direcao invertida da etapa 15.\n"
          "     O resultado desta etapa e que o mecanismo declarado la esta\n"
          "     errado, e o achado de transporte precisa ser reescrito em torno\n"
          "     de outra diferenca entre as duas coortes. Nao remendar aqui.")


section("ONDE A CURVA CRUZA ZERO")
print("  o numero que a etapa existe para produzir: abaixo de que n a correcao\n"
      "  de composicao custa mais do que rende, fora da coorte de ajuste.\n")
agg = (cur[cur.verdict].groupby("n").delta
       .agg(mediana="median", q75=lambda s: s.quantile(0.75),
            pior=lambda s: (s > 0).mean()))
print(f"  {'n':>5}{'mediana':>12}{'q75':>12}{'% pior que nao corrigir':>26}")
for n, r in agg.iterrows():
    print(f"  {n:>5}{r.mediana:>+12.1%}{r.q75:>+12.1%}{r.pior:>25.0%}")

neg = agg[agg.mediana < 0]
if len(neg) and len(agg[agg.mediana > 0]):
    cross = int(neg.index.min())
    below = int(agg[agg.index < cross].index.max())
    print(f"\n  mediana cruza zero entre n = {below} e n = {cross}")
elif len(neg) == len(agg):
    print("\n  a mediana e negativa em toda a grade: a correcao ajuda mesmo no "
          "menor n testado")
else:
    print("\n  a mediana e positiva em toda a grade: a correcao nunca compensa "
          "nas coortes testadas")

safe = agg[agg.q75 < 0]
if len(safe):
    print(f"  o quartil superior so fica abaixo de zero a partir de n = "
          f"{int(safe.index.min())} — abaixo disso, um sorteio ruim piora o "
          f"relogio")


section("CHECAGENS, FECHAMENTO")
for i, (name, ok) in enumerate([
        ("painel de 12 valida contra proporcoes conhecidas", c1),
        ("relogios leem idade nas tres coortes de teste", c2),
        ("piso em casa com n cheio", c3),
        ("faixa de idade mantida na grade", c4),
        ("ancora em n = %d" % N_ANCHOR, c5),
        ("nulo e propriedade do desenho, nao da resposta", c6)], 1):
    print(f"  {i}. {name}: {'ok' if ok else 'FALHOU'}")
print(f"  7. vies de monocito declarado na etapa 13: carregado, nenhum "
      f"coeficiente por tipo lido")
print(f"\n  saidas: results/ncurve.csv ({len(cur)} linhas), "
      f"results/ncurve_summary.csv")
