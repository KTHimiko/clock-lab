#!/usr/bin/env python3
"""
Cell composition from methylation — the stage 7 method, extracted.

A whole-blood profile is written as a non-negative mixture of purified
cell-type profiles that sums to one. The panel is built by picking, for each
type, the probes that separate it from the other five, and the proportions come
out of a constrained least squares.

Stage 7 validated this against the twelve GSE110554 mixtures, whose proportions
are known rather than estimated: r = 0.985 over the 72 cell-by-mixture values,
mean absolute error 0.021. Stage 7's copy stays inline in its own script, since
that is the run the validation describes; this module is the same method for
later stages to reuse, and every stage that uses it revalidates before trusting
it.
"""
import numpy as np
import pandas as pd
from scipy.optimize import nnls

PANEL = {"CD4T": "CD4+ T cells", "CD8T": "CD8+ T cells",
         "Bcell": "CD19+ B cells", "NK": "CD56+ NK cells",
         "Mono": "CD14+ Monocytes", "Neu": "Neutrophils"}
TYPES = list(PANEL)


def build_panel(betas, fraction, restrict_to=None, n_side=50, n_candidate=200):
    """Reference profiles over the probes that best separate the six types."""
    scores, cand = {}, set()
    for t, name in PANEL.items():
        sel = (fraction == name).to_numpy()
        rest = np.isin(fraction.to_numpy(), [PANEL[o] for o in TYPES if o != t])
        a, c = betas.loc[:, sel], betas.loc[:, rest]
        d = a.mean(axis=1) - c.mean(axis=1)
        sd = np.sqrt((a.var(axis=1, ddof=1) + c.var(axis=1, ddof=1)) / 2) + 0.01
        s = (d / sd).replace([np.inf, -np.inf], np.nan).dropna()
        scores[t] = s
        cand |= set(s.nlargest(n_candidate).index) | set(s.nsmallest(n_candidate).index)

    pool = sorted(cand if restrict_to is None else (cand & set(restrict_to)))
    probes = set()
    for t in TYPES:
        s = scores[t].reindex(pool).dropna()
        probes |= set(s.nlargest(n_side).index) | set(s.nsmallest(n_side).index)
    probes = sorted(probes)
    ref = pd.DataFrame({t: betas.loc[probes, (fraction == PANEL[t]).to_numpy()].mean(axis=1)
                        for t in TYPES})
    return ref


def deconvolve(ref, Y, weight=100.0):
    """
    Non-negative proportions summing to one, per column of Y.

    The sum-to-one constraint rides in as one heavily weighted extra equation,
    which keeps the whole thing inside plain NNLS.
    """
    R = np.vstack([ref, np.full((1, ref.shape[1]), weight)])
    out = np.empty((Y.shape[1], ref.shape[1]))
    for j in range(Y.shape[1]):
        y = np.concatenate([Y[:, j], [weight]])
        ok = ~np.isnan(y)
        w, _ = nnls(R[ok], y[ok])
        out[j] = w / w.sum() if w.sum() > 0 else np.nan
    return out
