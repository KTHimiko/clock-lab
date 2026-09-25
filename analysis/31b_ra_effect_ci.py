#!/usr/bin/env python3
"""
Stage 31b — a confidence interval for the arthritis effect under two adjustments.

Stage 31 reported that the Horvath 2018 arthritis effect on age acceleration is
-0.74 years with a whole-cohort composition adjustment and -1.39 years with a
controls-only one, "nearly a factor of two". A reader asked for intervals before
that phrase is used. Nonparametric bootstrap over people (1,000 resamples,
stratified by case status); both corrections are refitted inside every resample,
so the interval includes the uncertainty of the correction itself.

Reported for all three clocks: each effect and the difference between the two
adjustments, with percentile 95% intervals. No bar: this qualifies a sentence.

Usage:  .venv/bin/python analysis/31b_ra_effect_ci.py
"""
import numpy as np
src = open(__file__.replace("31b_ra_effect_ci.py", "31_controls_reference.py")).read()
exec(src.split('section("CHECAGEM 2')[0])


def effects(idx):
    dd_age = d["chrono"][idx]; C = d["C12"][idx]; r = ra[idx]
    ctl = np.where(~r)[0]
    Xa = np.column_stack([np.ones(len(idx)), dd_age, r.astype(float)])
    out = {}
    for c in CLOCKS3:
        y = d["y"][c][idx]
        b_all = ridge_coefs(dd_age, C, y, [0.0])[0.0]
        e1 = np.linalg.lstsq(Xa, y - (C[:, :-1] - C.mean(0)[:-1]) @ b_all, rcond=None)[0][2]
        b_ctl = ridge_coefs(dd_age[ctl], C[ctl], y[ctl], [0.0])[0.0]
        e2 = np.linalg.lstsq(Xa, y - (C[:, :-1] - C[ctl].mean(0)[:-1]) @ b_ctl, rcond=None)[0][2]
        out[c] = (e1, e2, e2 - e1)
    return out


n = len(ra); cases, ctls = np.where(ra)[0], np.where(~ra)[0]
point = effects(np.arange(n))
rng = np.random.default_rng(3131)
boot = {c: [] for c in CLOCKS3}
for _ in range(1000):
    idx = np.concatenate([rng.choice(cases, len(cases)), rng.choice(ctls, len(ctls))])
    e = effects(idx)
    for c in CLOCKS3:
        boot[c].append(e[c])
section("EFEITO DA AR, IC 95% POR BOOTSTRAP (1000)")
print(f"  {'relogio':<13}{'coorte inteira':>26}{'so controles':>26}{'diferenca':>26}")
for c in CLOCKS3:
    B = np.array(boot[c]); lo, hi = np.percentile(B, [2.5, 97.5], axis=0)
    p = point[c]
    print(f"  {c:<13}" + "".join(f"{p[k]:>+8.2f} [{lo[k]:+.2f}, {hi[k]:+.2f}]" for k in range(3)))
