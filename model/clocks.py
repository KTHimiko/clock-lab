"""
Epigenetic clocks with published, open coefficients.

Each clock is a linear predictor over a small set of CpGs: a weighted sum of
beta values plus an intercept. Horvath's then passes the result through an
inverse transform, because it was fitted on a transformed age to handle the
fast methylation changes of childhood; the others predict age directly.

Coefficients come from the dnaMethyAge package's data files, which carry the
values published with each paper. Nothing here is refitted.
"""
import numpy as np
import pandas as pd
import pyreadr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CLOCK_DIR = ROOT / "reference/data/clocks"

CLOCKS = {
    "Horvath2013": dict(file="HorvathS2013.rda", transform="horvath", adult_age=20.0),
    "Hannum2013": dict(file="HannumG2013.rda", transform="linear"),
    "Levine2018": dict(file="LevineM2018.rda", transform="linear"),
    "Horvath2018": dict(file="HorvathS2018.rda", transform="horvath", adult_age=20.0),
}


def load_coefficients(name):
    spec = CLOCKS[name]
    obj = pyreadr.read_r(str(CLOCK_DIR / spec["file"]))
    coefs = next(iter(obj.values()))
    coefs.columns = ["probe", "coef"]
    intercept = 0.0
    mask = coefs.probe.str.lower().isin(["intercept", "(intercept)"])
    if mask.any():
        intercept = float(coefs.loc[mask, "coef"].iloc[0])
        coefs = coefs.loc[~mask]
    return coefs.set_index("probe")["coef"].astype(float), intercept, spec


def _inverse_horvath(x, adult_age=20.0):
    """Undo the transform Horvath fitted on: log(age+1)-log(adult+1) below
    adulthood, linear above. Applied to the linear predictor."""
    x = np.asarray(x, dtype=float)
    return np.where(x <= 0,
                    (1.0 + adult_age) * np.exp(x) - 1.0,
                    (1.0 + adult_age) * x + adult_age)


def predict(betas, name, min_coverage=0.8):
    """
    Predict epigenetic age for every column of `betas` (probes x samples).

    Missing probes are a fact of life across array versions, so coverage is
    reported rather than assumed. Below `min_coverage` the prediction is
    refused: a clock missing a fifth of its sites is not that clock any more.
    """
    coefs, intercept, spec = load_coefficients(name)
    present = coefs.index.intersection(betas.index)
    coverage = len(present) / len(coefs)
    if coverage < min_coverage:
        raise ValueError(f"{name}: only {coverage:.1%} of probes present")

    sub = betas.loc[present].astype(float)
    # a missing value at a present probe is filled with that probe's own mean
    # across samples, which is the least assuming choice available here
    sub = sub.T.fillna(sub.mean(axis=1)).T
    linear = intercept + sub.mul(coefs.loc[present], axis=0).sum(axis=0)

    if spec["transform"] == "horvath":
        age = _inverse_horvath(linear.values, spec.get("adult_age", 20.0))
    else:
        age = linear.values
    return pd.Series(age, index=betas.columns, name=name), coverage
