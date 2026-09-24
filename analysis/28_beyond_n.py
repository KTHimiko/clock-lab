#!/usr/bin/env python3
"""
Stage 28 — does the index know anything that the fitting size does not?

The manuscript quotes the transport index ranking 63 configurations at
p = 2e-24. A reviewer objected that the configurations share cohorts, so they
are not independent and the p-value is not what it claims. That is right, and
there is a deeper problem underneath it.

Within a directed pair the index falls with n by construction — it carries 1/n —
and the damage falls with n too. So a large part of any correlation pooled over
fitting sizes is guaranteed by n alone, and would appear for any quantity
proportional to 1/n whatever it knew about the cohorts. The question the index
has to answer is not whether it correlates with damage. It is whether it says
anything about WHICH pairs are dangerous that the fitting size does not.

THE DESIGN, off stage 26's saved rows (per clock, per draw, net damage).
Restricted to the fitting sizes every directed pair has — 40, 80 and 160 — so
that every pair contributes a complete profile: 12 pairs x 2 clocks x 3 sizes =
72 configurations. The index depends only on composition, so it is shared by the
two clocks of a pair; the pair is the block.

  - within-n: at each fitting size separately, Spearman between index and net
    damage over the 24 (pair x clock) cells. Here n is constant and cannot
    masquerade as anything.
  - block permutation (Winkler et al., NeuroImage 2015): whole index profiles
    are swapped between pairs, each pair keeping its own fitting sizes. The null
    distribution therefore contains all the correlation n produces, and the
    p-value measures only what the index adds across pairs.
  - leave-one-cohort-out: the pooled rho with every configuration touching one
    cohort removed, for each cohort in turn.

PRE-REGISTERED
  1. within-n: the median over the three fitting sizes of the within-n rho
     exceeds 0.3.
  2. block permutation: p < 0.05 over 20,000 permutations.
  3. leave-one-cohort-out: reported as a range, no bar.
  The pooled rho and its naive p are printed beside them for comparison with
  what the manuscript says.

Usage:  .venv/bin/python analysis/28_beyond_n.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


ROOT = Path(__file__).resolve().parent.parent
R = pd.read_csv(ROOT / "results" / "net_damage.csv")
NS = [40, 80, 160]
cfg = (R[R.n.isin(NS)].groupby(["src", "dst", "n", "clock"])
         .agg(index=("tindex", "median"), obs=("obs", "median")).reset_index())
cfg["pair"] = cfg.src + "->" + cfg.dst
pairs = sorted(cfg.pair.unique())
print(f"{len(cfg)} configuracoes, {len(pairs)} pares, tamanhos {NS}")
assert len(cfg) == len(pairs) * 2 * len(NS), "perfil incompleto em algum par"

section("REFERENCIA — O RHO POOLED E O P INGENUO (O QUE O MANUSCRITO DIZ)")
rho0, p0 = spearmanr(cfg["index"], cfg.obs)
print(f"  rho = {rho0:.3f}, p ingenuo = {p0:.1e}")

section("CHECAGEM 1 — DENTRO DE CADA N")
within = []
for n in NS:
    g = cfg[cfg.n == n]
    r, p = spearmanr(g["index"], g.obs)
    within.append(r)
    print(f"  n = {n:>3}: rho = {r:+.3f}  (p = {p:.2g}, {len(g)} celulas)")
c1 = float(np.median(within)) > 0.3
print(f"\n  mediana {np.median(within):+.3f} (barra > 0.3) -> "
      f"{'o indice sabe algo alem de n' if c1 else 'o indice NAO acrescenta alem de n'}")

section("CHECAGEM 2 — PERMUTACAO EM BLOCOS (PERFIS DE INDICE TROCADOS ENTRE PARES)")
prof = cfg.pivot_table(index="pair", columns="n", values="index", aggfunc="first")
rng = np.random.default_rng(20260928)
null = np.empty(20_000)
for i in range(len(null)):
    perm = dict(zip(pairs, rng.permutation(pairs)))
    shuffled = [prof.loc[perm[p], n] for p, n in zip(cfg.pair, cfg.n)]
    null[i] = spearmanr(shuffled, cfg.obs)[0]
p_block = float((np.sum(null >= rho0) + 1) / (len(null) + 1))
c2 = p_block < 0.05
print(f"  rho observado {rho0:.3f} | nulo em blocos: mediana {np.median(null):.3f}, "
      f"95o percentil {np.percentile(null, 95):.3f}")
print(f"  p em blocos = {p_block:.4f} (barra < 0.05) -> {'ok' if c2 else 'FALHOU'}")
print(f"  (o nulo ja carrega a correlacao que n produz sozinho: {np.median(null):.3f})")

section("CHECAGEM 3 — DEIXANDO UMA COORTE DE FORA (REPORTADO)")
for k in sorted(set(cfg.src)):
    g = cfg[(cfg.src != k) & (cfg.dst != k)]
    print(f"  sem {k}: rho = {spearmanr(g['index'], g.obs)[0]:+.3f} ({len(g)} config.)")

section("FECHAMENTO")
print(f"  1. dentro de n: {'ok' if c1 else 'FALHOU'} (mediana {np.median(within):+.3f})")
print(f"  2. permutacao em blocos: {'ok' if c2 else 'FALHOU'} (p = {p_block:.4f})")
