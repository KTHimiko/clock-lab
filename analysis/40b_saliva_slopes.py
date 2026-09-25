#!/usr/bin/env python3
"""
Stage 40b — the saliva immune-fraction slopes of stages 39-40, saved for the figure.

For each saliva cohort, clock age regressed on chronological age and the
immune fraction (EpiDISH, three-type step); the slope per +10 percentage points
of immune fraction with its standard error. GSE149747 uses its baseline samples.
Nothing new is estimated: these are the numbers stages 39 and 40 printed.

Usage:  .venv/bin/python analysis/40b_saliva_slopes.py
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "results/cache"


def slope(age, ic, y):
    X = np.column_stack([np.ones(len(y)), age, ic])
    b, *_ = np.linalg.lstsq(X, y, rcond=None); e = y - X @ b
    se = np.sqrt(e @ e / (len(y) - 3) * np.linalg.inv(X.T @ X)[2, 2])
    return b[2] / 10, se / 10


rows = []
for tag, array, label in [("GSE232891", "EPIC", "GSE232891"), ("GSE232332", "EPIC", "GSE232332"),
                          ("GSE78874", "450k", "GSE78874"), ("GSE78874qn", "450k", "GSE78874 (normalised)")]:
    a = pd.read_csv(CACHE / f"{tag}_ages.csv", index_col=0)
    ic = pd.read_csv(CACHE / f"{tag}_compmeas.csv", index_col=0).loc[a.index, "IC"].to_numpy()
    for c in ("Horvath2018", "Levine2018"):
        b, se = slope(a.chrono.to_numpy(float), ic, a[c].to_numpy(float))
        rows.append(dict(cohort=label, array=array, clock=c, n=len(a), slope=b, se=se))
d = pd.read_csv(CACHE / "GSE149747_saliva.csv", index_col=0)
d = d[d.tp == "Before"].dropna(subset=["age"])
for c in ("Horvath2018", "Levine2018"):
    b, se = slope(d.age.to_numpy(float), d.IC.to_numpy(float), d[c].to_numpy(float))
    rows.append(dict(cohort="GSE149747", array="EPIC", clock=c, n=len(d), slope=b, se=se))
out = pd.DataFrame(rows)
out.to_csv(ROOT / "results/saliva_slopes.csv", index=False)
print(out.round(2).to_string(index=False))
