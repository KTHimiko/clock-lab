"""
DunedinPACE (Belsky et al., eLife 2022), reimplemented from the published
model data of the DunedinPACE R package (GPL-3), which is kept in
reference/data/clocks/ and not redistributed here.

The score is not a plain weighted sum of raw betas. Each sample is first
quantile-normalised, over a background of 20,000 probes, to the Dunedin Study's
reference distribution, so that values sit on the scale of "years of biological
ageing per calendar year". Background probes absent from the array are filled
with their reference means before normalising; model probes missing after that
fall back to the model means.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
M = np.load(ROOT / "reference/data/clocks/dunedinpace.npz")


def predict(betas, min_coverage=0.8):
    gold, gmeans = M["gold_probes"], M["gold_means"]
    present = np.isin(gold, betas.index)
    if present.mean() < min_coverage or np.isin(M["model_probes"], betas.index).mean() < min_coverage:
        raise ValueError("DunedinPACE: cobertura de sondas insuficiente")
    X = np.tile(gmeans[:, None], (1, betas.shape[1]))
    X[present] = betas.loc[gold[present]].to_numpy(dtype=float)
    miss = np.isnan(X)
    X[miss] = np.broadcast_to(gmeans[:, None], X.shape)[miss]
    target = np.sort(gmeans)
    ranks = X.argsort(axis=0).argsort(axis=0)
    Xn = target[ranks]
    pos = {p: i for i, p in enumerate(gold)}
    rows = np.array([pos.get(p, -1) for p in M["model_probes"]])
    vals = np.where(rows[:, None] >= 0, Xn[np.clip(rows, 0, None)],
                    M["model_means"][:, None])
    score = M["intercept"] + M["weights"] @ vals
    return pd.Series(score, index=betas.columns, name="DunedinPACE"), float(present.mean())
