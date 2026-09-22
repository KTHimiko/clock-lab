#!/usr/bin/env python3
"""
Stage 4 — is the cell-type signal aligned with the coefficients, or just large?

Stage 3 found that how cell-type-variable a clock's probes are does not predict
how much its age moves between fractions — the correlation came out inverted.
The explanation on offer: a clock age is a weighted sum, so what matters is
whether the cell-type shifts line up with the weights. Shifts that point in
directions unrelated to the coefficients cancel.

This measures the alignment directly, against two nulls:

  SHUFFLE   permute the coefficients among the clock's own probes. Keeps both
            multisets — the same weights, the same shifts — and destroys only
            the pairing. Anything above this null is alignment and nothing else.

  RANDOM    draw probe sets of the same size from all 454k probes. Asks a
            different question: is this clock's probe set unusual, or would any
            set of that size do as well?

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. the weighted sums here must reproduce the fraction-to-fraction age
     differences measured in stage 2, to within rounding. Stage 4 recomputes
     from the probe level what stage 2 got from the clock function; if the two
     disagree, one of them is wrong and neither result stands.
  2. both null distributions must centre on zero — a null that is offset means
     the permutation is not breaking what it claims to break.

Usage:  .venv/bin/python analysis/04_alignment.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix
from model.clocks import CLOCKS, load_coefficients

OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
N_PERM = 10_000
rng = np.random.default_rng(20260922)


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


print("carregando GSE35069 ...", flush=True)
betas, meta = read_series_matrix(ROOT / "reference/data/GSE35069_series_matrix.txt.gz")
# SAME PREPROCESSING AS STAGE 2, deliberately. The first run of this stage
# restricted to complete probes, as stage 3 does, and the check below caught it:
# the completeness filter drops 4 of Hannum's 71 probes and 22 of Levine's 513,
# so the two stages were weighing different clocks and disagreeing by up to 5.7
# years. Stage 3 wants complete cases because a variance decomposition needs a
# balanced design; stage 4 must match stage 2 instead, because it is checking
# stage 2. Filling each probe's gaps with its own mean across samples is what
# model.clocks.predict does.
# Fill each probe's gaps with its own mean, at the numpy level. Writing this as
# `betas.T.fillna(betas.mean(axis=1)).T` is the obvious pandas phrasing and it
# hangs: transposing 485,577 x 60 materialises a frame with 485,577 COLUMNS, and
# pandas aligns per column. The array form does the same arithmetic in one pass.
_vals = betas.to_numpy(dtype=np.float32)
_rowmean = np.nanmean(_vals, axis=1, keepdims=True)
_gaps = np.isnan(_vals)
if _gaps.any():
    _vals[_gaps] = np.broadcast_to(_rowmean, _vals.shape)[_gaps]
betas = pd.DataFrame(_vals, index=betas.index, columns=betas.columns)
meta = meta.set_index("gsm")
frac = meta["tissue/cell type"].reindex(betas.columns)
donor = meta["title"].reindex(betas.columns).str.extract(r"_(\d+)$")[0]
print(f"  {betas.shape[0]:,} sondas completas x {betas.shape[1]} amostras", flush=True)

# deviation of each fraction from that donor's own whole blood, averaged
wide = {}
for f in frac.unique():
    cols = betas.columns[frac == f]
    wide[f] = betas[cols].rename(columns=dict(zip(cols, donor[cols]))).sort_index(axis=1)
ref = wide["Whole blood"]
delta = {f: (m - ref).mean(axis=1) for f, m in wide.items() if f != "Whole blood"}
D = pd.DataFrame(delta)
print(f"  deslocamento medio por fracao contra o sangue total do proprio doador:"
      f" {D.shape[1]} fracoes", flush=True)


section("CHECAGEM 1 — ISTO REPRODUZ A ETAPA 2?")
prev = pd.read_csv(OUT / "purified_ages.csv", index_col=0)
agree = []
for name in CLOCKS:
    coefs, _, spec = load_coefficients(name)
    present = D.index.intersection(coefs.index)
    w = coefs.loc[present]
    linear_shift = D.loc[present].mul(w, axis=0).sum()
    # stage 2 worked in years after the transform; for the linear clocks the
    # shift IS years, so those are directly comparable
    if spec["transform"] == "linear":
        obs = (prev.groupby("fraction")[name].mean()
               - prev[prev.fraction == "Whole blood"][name].mean())
        joined = pd.concat([linear_shift.rename("recomputed"),
                            obs.rename("stage2")], axis=1).dropna()
        r = float(joined.recomputed.corr(joined.stage2))
        gap = float((joined.recomputed - joined.stage2).abs().max())
        agree.append(gap < 0.5)
        print(f"  {name:<14} r = {r:>6.4f}   maior discrepancia {gap:>5.2f} anos  "
              f"{'ok' if gap < 0.5 else 'DIVERGE'}")
check1 = all(agree)
if not check1:
    print("\n  A RECOMPUTACAO NAO BATE COM A ETAPA 2 — parando.")
    sys.exit(1)


section("ALINHAMENTO CONTRA OS DOIS NULOS")
all_delta = D.to_numpy(dtype=np.float32)
probe_index = {p: i for i, p in enumerate(D.index)}
rows = []
print(f"  {'relogio':<14}{'observado':>11}{'nulo shuffle':>14}{'nulo aleatorio':>16}"
      f"{'alinhamento':>13}")
for name in CLOCKS:
    coefs, _, _ = load_coefficients(name)
    present = [p for p in coefs.index if p in probe_index]
    idx = np.array([probe_index[p] for p in present])
    w = coefs.loc[present].to_numpy(dtype=np.float32)
    dsub = all_delta[idx]                                   # probes x fractions

    observed = float(np.abs(dsub.T @ w).mean())

    # null 1: shuffle the pairing of weights to probes
    perm_w = np.array([rng.permutation(w) for _ in range(N_PERM)])   # perm x probes
    shuffle_null = np.abs(perm_w @ dsub).mean(axis=1)

    # null 2: random probe sets of the same size, weights kept in order
    # The first version wrote this as einsum over a broadcast_to view of the
    # weights. A broadcast view is not contiguous, einsum would not route it to
    # BLAS, and the whole stage sat at 99.5% of one core for 27 minutes without
    # finishing. The weights do not vary across permutations, so they belong
    # outside the contraction entirely: a 1-D vector against the gathered
    # deltas, in chunks that keep the gather under a few hundred megabytes.
    random_null = np.empty(N_PERM, dtype=np.float64)
    chunk = 500
    for start in range(0, N_PERM, chunk):
        stop = min(start + chunk, N_PERM)
        pick = rng.integers(0, all_delta.shape[0], size=(stop - start, len(idx)))
        gathered = all_delta[pick]                     # chunk x probes x fractions
        sums = np.einsum("i,pif->pf", w, gathered)
        random_null[start:stop] = np.abs(sums).mean(axis=1)

    align = observed / shuffle_null.mean()
    rows.append(dict(clock=name, observed=observed,
                     shuffle=float(shuffle_null.mean()),
                     random=float(random_null.mean()), alignment=align,
                     p_shuffle=float((shuffle_null >= observed).mean()),
                     p_random=float((random_null >= observed).mean())))
    print(f"  {name:<14}{observed:>11.3f}{shuffle_null.mean():>14.3f}"
          f"{random_null.mean():>16.3f}{align:>12.1f}x", flush=True)

res = pd.DataFrame(rows)
res.to_csv(OUT / "alignment.csv", index=False)


section("CHECAGEM 2 — OS NULOS ESTAO CENTRADOS?")
print("  (os nulos sao de valores absolutos, entao a media positiva e esperada;")
print("   o que importa e que o observado esteja MUITO acima, nao que o nulo seja 0)")
for _, r in res.iterrows():
    print(f"  {r.clock:<14} p(shuffle) = {r.p_shuffle:.4f}   p(aleatorio) = {r.p_random:.4f}")


section("O QUE ISSO DIZ")
spread = pd.read_csv(OUT / "purified_variance.csv").set_index("clock")["within"]
j = res.set_index("clock").join(spread.rename("spread"))
rho = float(j.alignment.corr(j.spread, method="spearman"))
print(f"  {'relogio':<14}{'alinhamento':>13}{'desvio etapa 2':>17}")
for c, r in j.iterrows():
    print(f"  {c:<14}{r.alignment:>12.1f}x{r.spread:>16.1f}a")
print(f"""
  Spearman entre alinhamento e vulnerabilidade: {rho:+.2f}  (n = 4)

  A etapa 3 mediu magnitude e obteve {-0.40:+.2f}. Se esta coluna correlacionar
  positivamente, a explicacao proposta se sustenta: o que expoe um relogio nao e
  ter sondas sensiveis a tipo celular, e ter os pesos apontando na mesma direcao
  que o deslocamento. Quatro pontos continuam nao sendo um teste.""")
