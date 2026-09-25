#!/usr/bin/env python3
"""
Stage 38c — control for 38b: does the saliva correction work at home?

38b found the transported correction harmful in 7 of 8 saliva cells at matched n,
and alpha = 3 did not rescue it. Before reading that as model shift between saliva
cohorts, the correction has to be shown to work inside a cohort, measured the same
way. Otherwise the failure is the correction or the measurement, not transport.

POST HOC (written after 38b, before this was run). Within each saliva cohort, 20
random half splits: fit on one half, apply to the other (the "home" transport,
same study, same array, same processing), and the full cohort applied to itself.
  H. the home transport is beneficial (median delta < 0) for both age clocks in
     all three cohorts. If it holds, 38b is transport. If it fails, 38b says
     nothing about transport.
  Also: condition number of the nine-type fitting composition, per cohort.

Usage:  .venv/bin/python analysis/38c_saliva_home.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
src = open(ROOT / "analysis/38b_saliva_transport.py").read()
exec(src.split('section("CHECAGEM 2 E 3')[0].split('"""', 2)[2])

section("CONTROLE — CORRECAO EM CASA")
okH = True
for t in SAL:
    d = D[t]; n = N[t]
    Z = partial_out(d["chrono"], d["C12"][:, :-1]); Z = (Z - Z.mean(0)) / Z.std(0)
    print(f"  {t}: numero de condicao da composicao de ajuste (padronizada) {np.linalg.cond(Z):.0f}")
    for c in AGE + ["DunedinPACE"]:
        dl = []
        for r in range(20):
            perm = np.random.default_rng(3900 + r).permutation(n)
            fi, te = np.sort(perm[: n // 2]), np.sort(perm[n // 2:])
            b = ridge_coefs(d["chrono"][fi], d["C12"][fi], d["y"][c][fi], [0.0])[0.0]
            cb = d["C12"][fi].mean(0)
            y, a, C12, C6 = d["y"][c][te], d["chrono"][te], d["C12"][te], d["C6"][te]
            nul = null_mean(y, a, C6, n_perm=200)
            ib, base = inc_and_base(y, a, C6)
            ia, _ = inc_and_base(y - (C12[:, :-1] - cb[:-1]) @ b, a, C6)
            dl.append((ia - ib) / (1 - base))
        b = ridge_coefs(d["chrono"], d["C12"], d["y"][c], [0.0])[0.0]
        l_full = after(t, c, b, d["C12"].mean(0))
        med = float(np.median(dl))
        if c in AGE:
            okH &= med < 0
        print(f"    {c:<13} metade->metade: mediana {med:+.1%} (nocivo {np.mean(np.array(dl) > 0):.0%})"
              f" | coorte inteira em si mesma: antes {BEF[(t, c)][0]:+.1%}, restante {l_full:+.1%}")
print(f"\n  H. correcao em casa benefica nos dois relogios de idade, nas tres coortes -> {'ok' if okH else 'FALHOU'}")
