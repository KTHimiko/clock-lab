#!/usr/bin/env python3
"""
Stage 23 — one of the verdict clocks was trained on the fitting cohort.

A reviewer of the manuscript draft asked whether Horvath 2013 or Horvath 2018
had been trained on any of the four cohorts. Nobody had checked. Hannum was
excluded from the start because it was trained on GSE40279; the other three
were assumed clean.

The primary source settles it. Horvath 2013, Additional file 1 (the 82-dataset
table with its 'Data Use' column), fetched through Europe PMC:

    data set 3   Blood WB, 450K, Training, 656 samples, Hannum 2012, GSE40279
    data set 44  Blood WB, 450K, Test,     689 samples, Liu 2013,    GSE42861

GSE40279 is a TRAINING set of Horvath 2013. GSE42861 was used only as a test
set, which does not make the clock in-sample there. GSE61151 and GSE50660 do not
appear. Horvath 2018 was trained on GSE80261, GSE50759, GSE104471, GSE77136,
GSE52026, E-MTAB-4385, GSE79056 and unpublished data — none of the four. Levine
2018 was trained on InCHIANTI (n = 456) — none of the four.

So Horvath 2013 falls under the rule Hannum was under all along: out of every
pair in which GSE40279 appears, as fitting cohort or as test cohort. Its share
of the training data was 656 of roughly four thousand, so the in-sample effect
should be milder than Hannum's — but milder is not absent, and stage 15 already
showed what a clock measured in the cohort it learned does to a share-of-residual
denominator.

This touches the project's headline. Stage 17's curve was fitted on GSE40279
with Horvath 2013 as one of three verdict clocks. Stages 18-22 carried it too.

WHAT THIS STAGE DOES. It does not rerun anything expensive. The saved per-clock
rows of stages 17 and 21 are filtered by the corrected rule, and the matched-n
comparisons that were never saved are recomputed off the stage 18 cache with the
rule applied. Then each headline is compared, old against corrected.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. the cache is the stage 18 cache: panel r = 0.789, MAE = 0.027. Hard stop.
  2. is the in-sample status visible in the data at all? Horvath 2013's age
     residual on GSE40279, relative to its residual on the other three cohorts,
     against the same ratio for Horvath 2018 and Levine. Reported, not
     pass/fail: it says how much the contamination should matter, not whether
     the rule applies. The rule applies on the documentary evidence.
  3. THE HEADLINES, each with its survival criterion written here:
       a. the n-curve: median delta at n = 40 still above +5%, and the median
          still crossing zero somewhere between n = 120 and n = 240
       b. the index ranks configurations: Spearman rho still above 0.70
       c. the asymmetry: still at least five of six pairs called
       d. the fix: ridge alpha = 3 still holds all twelve pairs at or below zero
     A headline that fails its criterion is reported as changed, in the
     synthesis and in both manuscripts. None is re-worded to survive.

Usage:  .venv/bin/python analysis/23_training_audit.py
"""
import sys
from itertools import permutations, combinations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

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


section("CHECAGEM 2 — O STATUS IN-SAMPLE APARECE NOS DADOS? (REPORTADO)")
print("  desvio-padrao do residuo de idade (relogio - idade), por coorte\n")
print(f"  {'relogio':<14}" + "".join(f"{t[3:]:>10}" for t in COHORTS)
      + f"{'40279 / outras':>17}")
for c in ["Hannum2013", "Horvath2013", "Horvath2018", "Levine2018"]:
    sds = []
    for t in COHORTS:
        d = D[t]; Xa = np.column_stack([np.ones(N[t]), d["chrono"]])
        res = d["y"][c] - Xa @ np.linalg.lstsq(Xa, d["y"][c], rcond=None)[0]
        sds.append(res.std())
    ratio = sds[0] / np.mean(sds[1:])
    mark = "   <- treinou no 40279" if c in TRAINED_ON else ""
    print(f"  {c:<14}" + "".join(f"{v:>10.2f}" for v in sds) + f"{ratio:>17.2f}{mark}")
print("\n  razao < 1 so para os relogios que treinaram ali = o treino aparece no residuo")


print("\n  nulos de permutacao ...", flush=True)
NULL, BEFORE = {}, {}
for t in COHORTS:
    for c in CLOCKS3:
        NULL[(t, c)] = null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        ib, base = inc_and_base(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        BEFORE[(t, c)] = ((ib - NULL[(t, c)]) / (1 - base), base)


def transport(src, dst, n, reps, alphas=(0.0,)):
    s, d = D[src], D[dst]
    rows = []
    for r in range(reps):
        idx = (np.arange(N[src]) if n >= N[src]
               else stratified_draw(s["chrono"], n, np.random.default_rng(700 + 7 * r + n)))
        af, Cf = s["chrono"][idx], s["C12"][idx]
        cbar = Cf.mean(axis=0)
        ti = transport_index(af, Cf, d["chrono"], d["C12"])
        for c in CLOCKS3:
            if not allowed(c, src, dst):
                continue
            coefs = ridge_coefs(af, Cf, s["y"][c][idx], list(alphas))
            before, base = BEFORE[(dst, c)]
            for alpha, b in coefs.items():
                yc = d["y"][c] - (d["C12"][:, :-1] - cbar[:-1]) @ b
                ia, _ = inc_and_base(yc, d["chrono"], d["C6"])
                rows.append(dict(src=src, dst=dst, n=len(idx), clock=c, alpha=alpha,
                                 tindex=ti,
                                 delta=(ia - NULL[(dst, c)]) / (1 - base) - before))
    return rows


section("3a — A CURVA DE N (ETAPA 17), ANTES E DEPOIS DA REGRA")
cur = pd.read_csv(RES / "ncurve.csv")
old = cur[cur.verdict]                                # Horvath13, Levine, Horvath18
new = old[old.clock != "Horvath2013"]                 # fitting cohort is GSE40279
print(f"  {'n':>5}{'antes (3 relogios)':>22}{'depois (2 relogios)':>22}{'H2013 sozinho':>16}")
for n in sorted(old.n.unique()):
    o = old[old.n == n].delta.median(); w = new[new.n == n].delta.median()
    h = cur[(cur.n == n) & (cur.clock == "Horvath2013")].delta.median()
    print(f"  {n:>5}{o:>+22.1%}{w:>+22.1%}{h:>+16.1%}")
g = new.groupby("n").delta.median()
cross_hi = g[g < 0].index.min(); cross_lo = g[g.index < cross_hi].index.max()
at40 = g.loc[40]
c3a = at40 > 0.05 and 120 <= cross_lo and cross_hi <= 240
worse40 = (new[new.n == 40].delta > 0).mean()
print(f"\n  corrigido: n=40 {at40:+.1%} (nocivo em {worse40:.0%}), cruza entre "
      f"{cross_lo} e {cross_hi}")
print(f"  criterio: n=40 > +5% e cruzamento entre 120 e 240 -> "
      f"{'SOBREVIVE' if c3a else 'MUDOU'}")


section("3b — O INDICE ORDENA AS CONFIGURACOES (ETAPA 21), COM A REGRA")
A = pd.read_csv(RES / "all_pairs.csv")
A0 = A[A.alpha == 0.0]
keep = A0.apply(lambda r: allowed(r.clock, r.src, r.dst), axis=1)
for label, frame in (("antes", A0), ("depois", A0[keep])):
    cfg = (frame.groupby(["src", "dst", "n"])
                .agg(indice=("tindex", "median"), dano=("delta", "median")).reset_index())
    rho, _ = spearmanr(cfg.indice, cfg.dano)
    if label == "depois":
        rho_new = rho; cfg_new = cfg
    print(f"  {label:<8} {len(cfg):>3} configuracoes, {len(frame):>5} linhas, rho = {rho:.3f}")
c3b = rho_new > 0.70
print(f"  criterio: rho > 0.70 -> {'SOBREVIVE' if c3b else 'MUDOU'}")
safe = cfg_new[cfg_new.indice <= 0.05]
print(f"  piso: abaixo de indice 0.05, {int((safe.dano > 0).sum())} de {len(safe)} "
      f"nocivas; acima, {(cfg_new[cfg_new.indice > 0.05].dano > 0).mean():.0%}")


section("3c — ASSIMETRIA EM N CASADO, COM A REGRA")
print(f"  {'par':<16}{'n':>6}{'indice A->B':>13}{'B->A':>9}{'dano A->B':>12}{'B->A':>9}"
      f"{'relogios':>10}{'acertou':>9}")
hits = []
for a, b in combinations(COHORTS, 2):
    n = min(N[a], N[b])
    ga = pd.DataFrame(transport(a, b, n, N_REPS)); gb = pd.DataFrame(transport(b, a, n, N_REPS))
    ia, ib = ga.tindex.median(), gb.tindex.median()
    da, db = ga.delta.median(), gb.delta.median()
    hit = (ia > ib) == (da > db); hits.append(hit)
    print(f"  {a[3:]+'/'+b[3:]:<16}{n:>6}{ia:>13.4f}{ib:>9.4f}{da:>+12.1%}{db:>+9.1%}"
          f"{ga.clock.nunique():>10}{'sim' if hit else 'NAO':>9}")
c3c = sum(hits) >= 5
print(f"\n  {sum(hits)} de {len(hits)} -> {'SOBREVIVE' if c3c else 'MUDOU'}")


section("3d — O CONSERTO NOS DOZE PARES, COM A REGRA")
print(f"  {'par':<20}{'n':>6}{'OLS':>9}{'ridge a=3':>12}{'relogios':>10}")
fix = []
for src, dst in permutations(COHORTS, 2):
    n = min(N[src], N[dst])
    g = pd.DataFrame(transport(src, dst, n, N_REPS, alphas=(0.0, ALPHA_FIX)))
    d0 = g[g.alpha == 0.0].delta.median(); d3 = g[g.alpha == ALPHA_FIX].delta.median()
    fix.append(d3 <= 0)
    print(f"  {src[3:]+' -> '+dst[3:]:<20}{n:>6}{d0:>+9.1%}{d3:>+12.1%}"
          f"{g.clock.nunique():>10}{'' if d3 <= 0 else '   <- nao zera'}")
c3d = all(fix)
print(f"\n  {sum(fix)} de 12 em zero ou abaixo -> {'SOBREVIVE' if c3d else 'MUDOU'}")


section("FECHAMENTO")
for k, (name, ok) in enumerate([("a curva de n", c3a), ("o indice ordena", c3b),
                                 ("a assimetria", c3c), ("o conserto", c3d)]):
    print(f"  3{'abcd'[k]}. {name}: {'SOBREVIVE' if ok else 'MUDOU'}")
