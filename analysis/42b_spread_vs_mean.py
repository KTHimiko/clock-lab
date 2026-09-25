#!/usr/bin/env python3
"""
Stage 42b — between-cohort spread relative to the average effect.

Stage 42 read the blood/saliva difference as the sign of the composition effect.
A reader showed that too strong: with GSE78874 normalised, Horvath 2018's saliva
slope is no longer clearly of the opposite sign (−0.20, 95% CI −0.46 to +0.06)
yet 7 of 8 cells stay harmful; and Levine 2018 keeps its sign in all four saliva
cohorts while staying harmful at alpha = 3.

This restates it on a scale-free quantity. For each (tissue, axis, clock), a
random-effects fit of the per-cohort slopes gives the between-cohort standard
deviation tau (DerSimonian-Laird) and the inverse-variance mean; tau/|mean| says
how large the disagreement between cohorts is next to the effect itself. A shrunk
coefficient is near every cohort's own when that ratio is small.

Sources: results/sign_blood_saliva.csv (stage 42), normalised GSE78874 excluded.

Usage:  .venv/bin/python analysis/42b_spread_vs_mean.py
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
R = pd.read_csv(ROOT / "results/sign_blood_saliva.csv")
R = R[~R.cohort.str.contains("normalis")]
rows = []
for (tis, ax, c), g in R.groupby(["tissue", "axis", "clock"]):
    w = 1 / g.se ** 2
    m = float((w * g.slope).sum() / w.sum())
    Q = float((w * (g.slope - m) ** 2).sum()); df = len(g) - 1
    tau2 = max(0.0, (Q - df) / (w.sum() - (w ** 2).sum() / w.sum()))
    rows.append(dict(tissue=tis, axis=ax, clock=c, k=len(g), mean=m,
                     tau=np.sqrt(tau2), ratio=np.sqrt(tau2) / abs(m),
                     I2=max(0.0, (Q - df) / Q) if Q > 0 else 0.0,
                     both_signs=bool(((g.slope - 1.96 * g.se > 0).any()) and ((g.slope + 1.96 * g.se < 0).any()))))
D = pd.DataFrame(rows).sort_values(["tissue", "axis", "clock"])
D.to_csv(ROOT / "results/spread_vs_mean.csv", index=False)
print(f"  {'tecido':<8}{'eixo':<7}{'relogio':<13}{'k':>3}{'media':>9}{'tau':>8}{'tau/|media|':>13}{'I2':>7}  sinais opostos")
for r in D.itertuples():
    print(f"  {r.tissue:<8}{r.axis:<7}{r.clock:<13}{r.k:>3}{r.mean:>+9.2f}{r.tau:>8.2f}{r.ratio:>13.2f}{r.I2:>7.0%}"
          f"  {'sim' if r.both_signs else 'nao'}")
b = D[D.tissue == "blood"].ratio; s = D[D.tissue == "saliva"].ratio
print(f"\n  sangue {b.min():.2f} a {b.max():.2f} | saliva {s.min():.2f} a {s.max():.2f}"
      f" -> {'separam' if b.max() < s.min() else 'NAO separam'}")
