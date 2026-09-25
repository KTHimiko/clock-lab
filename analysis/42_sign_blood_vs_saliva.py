#!/usr/bin/env python3
"""
Stage 42 — why the penalty works in blood and not in saliva: does the sign hold?

Ridge shrinks a transported coefficient toward zero. If the target's true
coefficient has the same sign and a different size, shrinking toward zero moves
the correction toward any target on that side, and bounds the damage. If the
sign differs, shrinking can only reduce a wrong-direction correction, never make
it helpful. Stage 39 found a sign reversal in saliva (immune fraction). This
stage asks whether blood has one.

Per cohort, clock age regressed on chronological age and one composition
fraction; slope per +10 percentage points, with SE. Blood: the six blood cohorts,
twelve-type Salas composition, on the dominant myeloid-lymphoid axis (Neu,
neutrophils) and the naive/memory axis (CD8nv, naive CD8 T cells). Saliva: the
four cohorts of stages 38-40, immune fraction (results/saliva_slopes.csv, raw
GSE78874). Levine 2018 and Horvath 2018, usable in every cohort.

PRE-REGISTERED (written before any blood slope was computed)
  1. BLOOD KEEPS ITS SIGN on the dominant axis: for both clocks, the Neu slope
     has the same sign in all six blood cohorts.
  2. the same on the naive/memory axis (CD8nv).
  3. BETWEEN-COHORT HETEROGENEITY IS LARGER IN SALIVA: for both clocks, the
     Cochran I^2 of the saliva immune slopes exceeds the I^2 of the blood Neu
     slopes.

Usage:  .venv/bin/python analysis/42_sign_blood_vs_saliva.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))
from load_extended import TYPES12

CACHE = ROOT / "results/cache"; RES = ROOT / "results"
BLOOD = ["GSE40279", "GSE61151", "GSE50660", "GSE42861", "GSE132203", "GSE55763"]
TWO = ["Levine2018", "Horvath2018"]


def slope(age, x, y):
    X = np.column_stack([np.ones(len(y)), age, x])
    b, *_ = np.linalg.lstsq(X, y, rcond=None); e = y - X @ b
    return b[2] / 10, np.sqrt(e @ e / (len(y) - 3) * np.linalg.inv(X.T @ X)[2, 2]) / 10


def i2(b, se):
    w = 1 / se ** 2; m = (w * b).sum() / w.sum()
    Q = (w * (b - m) ** 2).sum(); df = len(b) - 1
    return max(0.0, (Q - df) / Q) if Q > 0 else 0.0, m


rows = []
for t in BLOOD:
    a = pd.read_csv(CACHE / f"{t}_ages.csv", index_col=0); a = a[a.chrono.notna()]
    C = pd.read_csv(CACHE / f"{t}_comp12.csv", index_col=0).loc[a.index, TYPES12]
    for c in TWO:
        for ax in ("Neu", "CD8nv"):
            b, se = slope(a.chrono.to_numpy(float), C[ax].to_numpy(float), a[c].to_numpy(float))
            rows.append(dict(tissue="blood", cohort=t, clock=c, axis=ax, slope=b, se=se))
S = pd.read_csv(RES / "saliva_slopes.csv")
S = S[~S.cohort.str.contains("normalised")]
for r in S.itertuples():
    rows.append(dict(tissue="saliva", cohort=r.cohort, clock=r.clock, axis="IC", slope=r.slope, se=r.se))
R = pd.DataFrame(rows); R.to_csv(RES / "sign_blood_saliva.csv", index=False)

print(f"  {'tecido':<8}{'eixo':<7}{'relogio':<13}" + "".join(f"{x[3:]:>14}" for x in BLOOD))
for (tis, ax, c), g in R[R.tissue == "blood"].groupby(["tissue", "axis", "clock"]):
    g = g.set_index("cohort").loc[BLOOD]
    print(f"  {tis:<8}{ax:<7}{c:<13}" + "".join(f"{b:>+8.2f} ({s:.2f})" for b, s in zip(g.slope, g.se)))
for c in TWO:
    g = R[(R.tissue == "saliva") & (R.clock == c)]
    print(f"  saliva  IC     {c:<13}" + "".join(f"{b:>+8.2f} ({s:.2f})" for b, s in zip(g.slope, g.se))
          + "   [" + ", ".join(g.cohort) + "]")

print()
ok1 = ok2 = ok3 = True
for c in TWO:
    for ax in ("Neu", "CD8nv"):
        g = R[(R.tissue == "blood") & (R.clock == c) & (R.axis == ax)]
        same = (np.sign(g.slope) == np.sign(g.slope.iloc[0])).all()
        I, m = i2(g.slope.to_numpy(), g.se.to_numpy())
        print(f"  sangue {ax:<6}{c:<13} mesmo sinal nas 6: {'sim' if same else 'NAO'}   I2 {I:.0%}  media ponderada {m:+.2f}")
        if ax == "Neu": ok1 &= same
        else: ok2 &= same
    gs = R[(R.tissue == "saliva") & (R.clock == c)]
    Is, ms = i2(gs.slope.to_numpy(), gs.se.to_numpy())
    Ib, _ = i2(*[R[(R.tissue == "blood") & (R.clock == c) & (R.axis == "Neu")][k].to_numpy() for k in ("slope", "se")])
    sgn = (np.sign(gs.slope) == np.sign(gs.slope.iloc[0])).all()
    print(f"  saliva IC    {c:<13} mesmo sinal nas 4: {'sim' if sgn else 'NAO'}   I2 {Is:.0%}  (sangue Neu {Ib:.0%})")
    ok3 &= Is > Ib
print(f"\n  1. sangue mantem o sinal no eixo Neu: {'ok' if ok1 else 'FALHOU'}")
print(f"  2. sangue mantem o sinal no eixo CD8nv: {'ok' if ok2 else 'FALHOU'}")
print(f"  3. heterogeneidade maior na saliva que no sangue (I2): {'ok' if ok3 else 'FALHOU'}")
