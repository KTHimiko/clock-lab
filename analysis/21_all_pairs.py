#!/usr/bin/env python3
"""
Stage 21 — all twelve directed pairs, because stage 20 rested on one.

Stage 20 found the quantity that governs the damage and tested it on
configurations built from a single fitting cohort. Two weaknesses were written
into its own closing section:

  - rho = 0.897 over 27 configurations, every one of them fitted on GSE40279.
    A ranking that holds inside one cohort's subsamples is a ranking of n, dressed
    up, until it is shown somewhere else.
  - the asymmetry — the property that separates the transport index from stage
    19's symmetric distance and the reason it explains the anchor — rests on ONE
    directed pair, GSE61151 -> GSE40279 against its own reverse.

Four cohorts give twelve directed pairs and six reversible ones. Everything runs
off the stage 18 cache, so this costs no reloading.

THE MATCHED-N RULE, which is what makes the asymmetry testable at all. A -> B at
n = 656 against B -> A at n = 184 confounds direction with size, and size is the
one thing already known to matter. Every asymmetry comparison therefore fits BOTH
directions at n = min(n_A, n_B). The index's prediction is then purely about
which covariance is doing the fitting.

CLOCKS. Three: Horvath 2013, Horvath 2018, Levine 2018. Hannum is dropped
entirely rather than included where it is safe — it was trained on GSE40279,
which appears here as a fitting cohort in six pairs and a test cohort in three,
and a clock that is in-sample in half the design is worse than a clock that is
absent from all of it.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. the cache is the stage 18 cache: panel r = 0.789, MAE = 0.027. Hard stop.
  2. in-sample floor for EVERY fitting cohort, not just GSE40279: a correction
     fitted on cohort A must drive A's own composition term down. Four cohorts,
     three clocks, twelve cells, all must pass. Hard stop — a cohort whose fit
     is broken cannot say anything about transport out of it.
  3. THE RANKING, at the bar stage 20's check 2b set and passed on one fitting
     cohort: Spearman rho > 0.70 between the transport index and the median
     damage, over every (directed pair x n) configuration. Now on four fitting
     cohorts instead of one. If it falls below 0.70 here, stage 20's rho = 0.897
     was a property of GSE40279's subsample ladder and not of the index.
  4. THE ASYMMETRY, which is the check this stage exists for. For each of the
     six reversible pairs, fitted at matched n, the index names one direction as
     worse. Does the damage agree?
        bar: correct in at least 5 of 6, with the exact binomial p reported.
     Five of six is p = 0.109 against a coin, which is not significance and is
     said so out loud — with six pairs there is no bar that both clears 0.05 and
     tolerates a single miss. The number reported is the count and its p, not a
     verdict dressed as one.
  5. THE FIX, generalised. Ridge at alpha = 3 — found by trial in stage 18,
     confirmed in stage 19, explained in stage 20 — must leave median damage at
     or below zero across all twelve directed pairs at matched n. This is the
     first test of it outside GSE40279-fitted corrections.
  6. the threshold is reported, not checked: stage 20 found damage crossing zero
     near an index of 0.05 on one fitting cohort. Where does it cross here?
  7. the declared monocyte bias from stage 13 carries.

Usage:  .venv/bin/python analysis/21_all_pairs.py
"""
import sys
from itertools import permutations, combinations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, binomtest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_extended import TYPES12
from model.clocks import CLOCKS
from model.deconvolution import TYPES as TYPES6

OUT = ROOT / "results"; CACHE = OUT / "cache"
RNG = np.random.default_rng(20260925)
N_REPS = 20
GRID = [40, 80, 160, 320, 480]
CLOCKS3 = ["Horvath2013", "Levine2018", "Horvath2018"]
COHORTS = ["GSE40279", "GSE61151", "GSE50660", "GSE42861"]
ALPHA_FIX = 3.0


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
    Ct = partial_out(age, C[:, :-1])
    Ct = Ct - Ct.mean(axis=0)
    return Ct.T @ Ct / len(Ct)


def transport_index(age_f, C_f, age_t, C_t, lam=0.0):
    Sf, St = sigma(age_f, C_f), sigma(age_t, C_t)
    M = np.linalg.inv(Sf + lam * np.eye(len(Sf)) * np.trace(Sf) / len(Sf))
    return float(np.trace(M @ Sf @ M @ St) / len(age_f))


def ridge_coefs(age, C, y, alphas):
    Ct = partial_out(age, C[:, :-1])
    yt = partial_out(age, y.reshape(-1, 1)).ravel()
    mu, sd = Ct.mean(axis=0), Ct.std(axis=0) + 1e-12
    Z = (Ct - mu) / sd
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    scale = float((s ** 2).mean())
    Uty = U.T @ yt
    return {a: (Vt.T @ (s * Uty / (s ** 2 + a * scale))) / sd for a in alphas}


# ------------------------------------------------------------------ data ----
need = [CACHE / f"{t}_{w}.csv" for t in COHORTS
        for w in ("comp12", "comp6", "ages")] + [CACHE / "panel.csv"]
if not all(f.exists() for f in need):
    sys.exit("cache ausente. Rode analysis/18_conditioning.py antes.")

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
N = {t: len(D[t]["chrono"]) for t in COHORTS}

section("CHECAGEM 1 — CACHE")
c1 = abs(panel.panel_r - 0.789) < 5e-3 and abs(panel.panel_mae - 0.027) < 5e-3
print(f"  painel r = {panel.panel_r:.3f} | MAE = {panel.panel_mae:.3f} "
      f"-> {'ok' if c1 else 'FALHOU'}")
print("  " + " | ".join(f"{t[3:]} n={N[t]}" for t in COHORTS))
if not c1:
    sys.exit("\n  parando: cache diferente.")
print("  (vies de monocito da etapa 13 carregado)")
print("  Hannum fora do desenho inteiro: treinou no GSE40279, que aqui e coorte\n"
      "  de ajuste em seis pares e de teste em tres.")

print("\n  nulos de permutacao ...", flush=True)
NULL = {(t, c): null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        for t in COHORTS for c in CLOCKS3}
BEFORE = {}
for t in COHORTS:
    for c in CLOCKS3:
        ib, base = inc_and_base(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        BEFORE[(t, c)] = ((ib - NULL[(t, c)]) / (1 - base), base)


section("CHECAGEM 2 — PISO EM CASA, NAS QUATRO COORTES DE AJUSTE")
c2 = True
print(f"  {'coorte':<12}" + "".join(f"{c:>16}" for c in CLOCKS3))
for t in COHORTS:
    s = D[t]; line = f"  {t[3:]:<12}"
    for c in CLOCKS3:
        b = ridge_coefs(s["chrono"], s["C12"], s["y"][c], [0.0])[0.0]
        cbar = s["C12"].mean(axis=0)
        yc = s["y"][c] - (s["C12"][:, :-1] - cbar[:-1]) @ b
        b0, _ = inc_and_base(s["y"][c], s["chrono"], s["C12"])
        a0, _ = inc_and_base(yc, s["chrono"], s["C12"])
        ok = a0 < b0; c2 &= ok
        line += f"{b0:>8.4f}->{a0:>6.4f}" if not ok else f"{b0:>9.4f}->{a0:.2f}"
    print(line)
print(f"  -> {'ok' if c2 else 'FALHOU'}")
if not c2:
    sys.exit("\n  parando: ajuste quebrado em alguma coorte.")


# ------------------------------------------------------- all directed pairs --
def run_pair(src, dst, n, reps, alphas=(0.0, ALPHA_FIX)):
    s, d = D[src], D[dst]
    rows = []
    for r in range(reps):
        idx = (np.arange(N[src]) if n >= N[src]
               else stratified_draw(s["chrono"], n, np.random.default_rng(900 + 7 * r + n)))
        a_f, C_f = s["chrono"][idx], s["C12"][idx]
        cbar = C_f.mean(axis=0)
        ti = transport_index(a_f, C_f, d["chrono"], d["C12"])
        for c in CLOCKS3:
            coefs = ridge_coefs(a_f, C_f, s["y"][c][idx], list(alphas))
            before, base = BEFORE[(dst, c)]
            for alpha, b in coefs.items():
                yc = d["y"][c] - (d["C12"][:, :-1] - cbar[:-1]) @ b
                ia, _ = inc_and_base(yc, d["chrono"], d["C6"])
                after = (ia - NULL[(dst, c)]) / (1 - base)
                rows.append(dict(src=src, dst=dst, n=len(idx), rep=r, clock=c,
                                 alpha=alpha, tindex=ti, before=before,
                                 delta=after - before))
    return rows


section("OS DOZE PARES DIRECIONADOS")
rows = []
for src, dst in permutations(COHORTS, 2):
    grid = sorted({g for g in GRID if g < N[src]} | {N[src]})
    for n in grid:
        rows += run_pair(src, dst, n, 1 if n >= N[src] else N_REPS)
    print(f"  {src[3:]} -> {dst[3:]}: n em {grid}", flush=True)

A = pd.DataFrame(rows)
A.to_csv(OUT / "all_pairs.csv", index=False)
ols = A[A.alpha == 0.0]


section("CHECAGEM 3 — O INDICE ORDENA AS CONFIGURACOES, EM QUATRO COORTES?")
cfg = (ols.groupby(["src", "dst", "n"])
          .agg(indice=("tindex", "median"), dano=("delta", "median"))
          .reset_index())
rho, pv = spearmanr(cfg.indice, cfg.dano)
c3 = rho > 0.70
print(f"  {len(cfg)} configuracoes, {ols.src.nunique()} coortes de ajuste")
print(f"  Spearman rho = {rho:.3f} (p = {pv:.2g})")
print(f"  barra: rho > 0.70 (a etapa 20 deu 0.897 com uma coorte so) "
      f"-> {'ok' if c3 else 'FALHOU'}")
for src in COHORTS:
    g = cfg[cfg.src == src]
    if len(g) > 3:
        rs, _ = spearmanr(g.indice, g.dano)
        print(f"    ajustando em {src[3:]:<8} ({len(g):>2} config.): rho = {rs:+.3f}")


section("CHECAGEM 4 — A ASSIMETRIA, EM N CASADO")
print("  cada par nas duas direcoes com n = min(n_A, n_B), para que so a\n"
      "  covariancia de ajuste mude.\n")
print(f"  {'par':<20}{'n':>6}{'indice A->B':>14}{'B->A':>10}"
      f"{'dano A->B':>12}{'B->A':>10}{'acertou?':>11}")
hits = []
for a, b in combinations(COHORTS, 2):
    n = min(N[a], N[b])
    ra = pd.DataFrame(run_pair(a, b, n, N_REPS, alphas=(0.0,)))
    rb = pd.DataFrame(run_pair(b, a, n, N_REPS, alphas=(0.0,)))
    ia, ib = ra.tindex.median(), rb.tindex.median()
    da, db = ra.delta.median(), rb.delta.median()
    pred_worse = a + "->" + b if ia > ib else b + "->" + a
    obs_worse = a + "->" + b if da > db else b + "->" + a
    hit = pred_worse == obs_worse
    hits.append(hit)
    print(f"  {a[3:]+'/'+b[3:]:<20}{n:>6}{ia:>14.4f}{ib:>10.4f}"
          f"{da:>+12.1%}{db:>+10.1%}{'sim' if hit else 'NAO':>11}")
k = sum(hits)
bt = binomtest(k, len(hits), 0.5, alternative="greater")
c4 = k >= 5
print(f"\n  acertou em {k} de {len(hits)} pares, binomial p = {bt.pvalue:.3f}")
print(f"  barra: >= 5 de 6 -> {'ok' if c4 else 'FALHOU'}")
print("  (5 de 6 da p = 0.109 contra uma moeda: com seis pares nao existe barra\n"
      "   que passe de 0.05 e tolere um erro. O numero e a contagem, nao veredito.)")


section("CHECAGEM 5 — O CONSERTO (ALPHA=3) FORA DO GSE40279")
print(f"  {'par':<22}{'n':>6}{'OLS':>10}{'ridge a=3':>12}")
c5 = True
for src, dst in permutations(COHORTS, 2):
    n = min(N[src], N[dst])
    r = pd.DataFrame(run_pair(src, dst, n, N_REPS))
    d0 = r[r.alpha == 0.0].delta.median()
    d3 = r[r.alpha == ALPHA_FIX].delta.median()
    ok = d3 <= 0
    c5 &= ok
    print(f"  {src[3:]+' -> '+dst[3:]:<22}{n:>6}{d0:>+10.1%}{d3:>+12.1%}"
          + ("" if ok else "   <- nao zera"))
print(f"\n  barra: mediana <= 0 nos doze pares -> {'ok' if c5 else 'FALHOU'}")


section("ONDE O DANO CRUZA ZERO, EM INDICE (REPORTADO, NAO CHECADO)")
cfg2 = cfg.sort_values("indice")
neg = cfg2[cfg2.dano <= 0]
pos = cfg2[cfg2.dano > 0]
if len(neg) and len(pos):
    hi_neg = neg.indice.max(); lo_pos = pos.indice.min()
    print(f"  maior indice ainda benefico: {hi_neg:.4f}")
    print(f"  menor indice ja nocivo:      {lo_pos:.4f}")
    print(f"  mediana dos dois:            {(hi_neg + lo_pos) / 2:.4f}")
print(f"  a etapa 20, com uma coorte de ajuste, achou a travessia perto de 0.05")
print(f"  fracao de configuracoes nocivas acima de indice 0.05: "
      f"{(cfg[cfg.indice > 0.05].dano > 0).mean():.0%}")
print(f"  fracao de configuracoes nocivas abaixo de 0.05:       "
      f"{(cfg[cfg.indice <= 0.05].dano > 0).mean():.0%}")


section("CHECAGENS, FECHAMENTO")
for i, (name, ok) in enumerate([
        ("cache identico", c1),
        ("piso em casa nas quatro coortes", c2),
        ("indice ordena as configuracoes (rho > 0.70)", c3),
        ("assimetria acerta >= 5 de 6 pares", c4),
        ("ridge alpha=3 nao prejudica em nenhum dos doze pares", c5)], 1):
    print(f"  {i}. {name}: {'ok' if ok else 'FALHOU'}")
print("  6. travessia em indice: reportada acima, nao checada")
print("  7. vies de monocito da etapa 13: carregado")
print(f"\n  saida: results/all_pairs.csv ({len(A)} linhas)")
