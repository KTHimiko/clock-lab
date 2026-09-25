#!/usr/bin/env python3
"""
Stage 31 — model shift inside one study: controls as the reference for cases.

Stage 25 attributed the worst transports to model shift — the composition effect
itself differing between cohorts — and found, within GSE42861, that Levine 2018's
composition coefficients differ between rheumatoid arthritis cases and controls
(p = 0.013) while Horvath 2018's do not (p = 0.61). But between cohorts, model
shift is entangled with everything else that differs: array batch, laboratory,
collection, population.

Within GSE42861, cases and controls share the study, the array and the lab. A
correction fitted on controls and applied to cases — a reference-population design
some case-control analyses use — is a transport in which the only thing that
differs is disease. If it fails there beyond what estimation noise explains, the
failure cannot be batch.

THE DESIGN, repeated 30 times: controls split at random into halves A and B.
  control -> control   fit on A, apply to B: estimation noise, no model shift
  control -> case      fit on the SAME A, apply to the cases
Fitting set and size are identical in the two arms, so the paired difference in
net damage is the excess the disease adds. Scored as always (fit with twelve
types, residual measured with six, null subtracted, per clock).

CLOCKS. In GSE42861 all three are clean: GSE42861 was only a test set of Horvath
2013, and appears in no training of Levine 2018 or Horvath 2018.

The consequence for a study: the RA effect on age acceleration, estimated three
ways — unadjusted, adjusted within the whole cohort (the standard), and with a
correction fitted on controls and applied to all.

PRE-REGISTERED
  1. the cache and the metadata join are intact; controls and cases > 300. Hard stop.
  2. DIFFERENTIAL PREDICTION from stage 25:
       Levine 2018:   median paired excess (control->case minus control->control)
                      above +1 point
       Horvath 2018:  median paired excess below +1 point
     Both must hold for "model shift within a single study" to stand for the
     clock stage 25 flagged, and not for the one it cleared.
  3. Horvath 2013: reported, no prediction (stage 25 did not test it).
  4. the RA effect under the three adjustments: reported per clock.

Usage:  .venv/bin/python analysis/31_controls_reference.py
"""
import sys, gzip
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_extended import TYPES12
from model.clocks import CLOCKS
from model.deconvolution import TYPES as TYPES6

RES = ROOT / "results"; CACHE = RES / "cache"
RNG = np.random.default_rng(20260924)
N_REPS = 20
COHORTS = ["GSE40279", "GSE61151", "GSE50660", "GSE42861"]
CLOCKS3 = ["Horvath2013", "Levine2018", "Horvath2018"]
ALPHA_FIX = 3.0

# the corrected rule: a clock is out of any pair in which a cohort it trained on
# appears, whichever end of the transport that cohort is on
TRAINED_ON = {"Hannum2013": {"GSE40279"}, "Horvath2013": {"GSE40279"}}


def allowed(clock, *cohorts):
    return not (TRAINED_ON.get(clock, set()) & set(cohorts))


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


def r2(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return 1 - (y - X @ beta).var() / y.var()


def inc_and_base(y, a, C):
    Xa = np.column_stack([np.ones(len(a)), a])
    base = r2(Xa, y)
    return r2(np.column_stack([Xa, C[:, :-1]]), y) - base, base


def null_mean(y, a, C, n_perm=800, rng=RNG):
    Xa = np.column_stack([np.ones(len(a)), a])
    base = r2(Xa, y)
    return float(np.mean([r2(np.column_stack([Xa, C[rng.permutation(len(C))][:, :-1]]), y)
                          - base for _ in range(n_perm)]))


def stratified_draw(ages, n, rng):
    q = pd.qcut(pd.Series(ages), 10, labels=False, duplicates="drop").to_numpy()
    idx = []
    for b in np.unique(q):
        pool = np.where(q == b)[0]
        take = min(max(int(round(n * len(pool) / len(ages))), 1), len(pool))
        idx.extend(rng.choice(pool, take, replace=False))
    idx = np.array(sorted(idx))
    if len(idx) > n:
        idx = np.sort(rng.choice(idx, n, replace=False))
    elif len(idx) < n:
        rest = np.setdiff1d(np.arange(len(ages)), idx)
        idx = np.sort(np.concatenate([idx, rng.choice(rest, n - len(idx), replace=False)]))
    return idx


def partial_out(age, Y):
    Xa = np.column_stack([np.ones(len(age)), age])
    return Y - Xa @ np.linalg.pinv(Xa) @ Y


def sigma(age, C):
    Ct = partial_out(age, C[:, :-1]); Ct = Ct - Ct.mean(axis=0)
    return Ct.T @ Ct / len(Ct)


def transport_index(af, Cf, at, Ct):
    return float(np.trace(np.linalg.inv(sigma(af, Cf)) @ sigma(at, Ct)) / len(af))


def ridge_coefs(age, C, y, alphas):
    Ct = partial_out(age, C[:, :-1])
    yt = partial_out(age, y.reshape(-1, 1)).ravel()
    mu, sd = Ct.mean(axis=0), Ct.std(axis=0) + 1e-12
    Z = (Ct - mu) / sd
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    scale = float((s ** 2).mean()); Uty = U.T @ yt
    return {a: (Vt.T @ (s * Uty / (s ** 2 + a * scale))) / sd for a in alphas}


# ------------------------------------------------------------------ data ----
panel = pd.read_csv(CACHE / "panel.csv", index_col=0).iloc[:, 0]
D = {}
for tag in COHORTS:
    a = pd.read_csv(CACHE / f"{tag}_ages.csv", index_col=0)
    k = a.chrono.notna().to_numpy(); ix = a.index[k]
    D[tag] = dict(chrono=a.chrono[k].to_numpy(),
                  y={c: a[c][k].to_numpy() for c in CLOCKS},
                  C12=pd.read_csv(CACHE / f"{tag}_comp12.csv", index_col=0)
                        .loc[ix, TYPES12].to_numpy(),
                  C6=pd.read_csv(CACHE / f"{tag}_comp6.csv", index_col=0)
                       .loc[ix, TYPES6].to_numpy())
N = {t: len(D[t]["chrono"]) for t in COHORTS}

section("CHECAGEM 1 — CACHE")
c1 = abs(panel.panel_r - 0.789) < 5e-3 and abs(panel.panel_mae - 0.027) < 5e-3
print(f"  painel r = {panel.panel_r:.3f} | MAE = {panel.panel_mae:.3f} -> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("  parando: cache diferente.")


def header_disease(gse="GSE42861"):
    samples, field = None, None
    with gzip.open(ROOT / "reference/data" / f"{gse}_series_matrix.txt.gz", "rt",
                   errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            parts = [v.strip('"') for v in line.rstrip("\n").split("\t")]
            if parts[0] == "!Sample_geo_accession":
                samples = parts[1:]
            elif parts[0] == "!Sample_characteristics_ch1" and parts[1].startswith("disease state"):
                field = [p.partition(":")[2].strip() for p in parts[1:]]
    return pd.Series(field, index=samples)


a42 = pd.read_csv(CACHE / "GSE42861_ages.csv", index_col=0)
dis = header_disease().reindex(a42.index[a42.chrono.notna()])
ra = dis.str.contains("rheumatoid", case=False).to_numpy()
d = D["GSE42861"]
ictl, icase = np.where(~ra)[0], np.where(ra)[0]
c1 = dis.notna().all() and len(ictl) > 300 and len(icase) > 300
print(f"  GSE42861: controles {len(ictl)}, casos {len(icase)} -> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("  parando.")


def sub(idx):
    return dict(chrono=d["chrono"][idx], C12=d["C12"][idx], C6=d["C6"][idx],
                y={c: d["y"][c][idx] for c in CLOCKS})


print("  nulos por grupo ...", flush=True)
NUL, BEF = {}, {}
for lab, idx in (("case", icase),):
    g = sub(idx)
    for c in CLOCKS3:
        NUL[(lab, c)] = null_mean(g["y"][c], g["chrono"], g["C6"])
        ib, base = inc_and_base(g["y"][c], g["chrono"], g["C6"])
        BEF[(lab, c)] = ((ib - NUL[(lab, c)]) / (1 - base), base)


def score(g, key_null, key_bef, c, b, cbar):
    yc = g["y"][c] - (g["C12"][:, :-1] - cbar[:-1]) @ b
    ia, _ = inc_and_base(yc, g["chrono"], g["C6"])
    before, base = key_bef
    return (ia - key_null) / (1 - base) - before


rows = []
rng = np.random.default_rng(20260925)
for r in range(30):
    perm = rng.permutation(ictl)
    A, B = perm[: len(perm) // 2], perm[len(perm) // 2:]
    gA, gB, gC = sub(A), sub(B), sub(icase)
    for c in CLOCKS3:
        b = ridge_coefs(gA["chrono"], gA["C12"], gA["y"][c], [0.0])[0.0]
        cbar = gA["C12"].mean(axis=0)
        nB = null_mean(gB["y"][c], gB["chrono"], gB["C6"], n_perm=300)
        ibB, baseB = inc_and_base(gB["y"][c], gB["chrono"], gB["C6"])
        befB = ((ibB - nB) / (1 - baseB), baseB)
        rows.append(dict(rep=r, clock=c, arm="ctl->ctl", delta=score(gB, nB, befB, c, b, cbar)))
        rows.append(dict(rep=r, clock=c, arm="ctl->case",
                         delta=score(gC, NUL[("case", c)], BEF[("case", c)], c, b, cbar)))
R = pd.DataFrame(rows)
R.to_csv(RES / "controls_reference.csv", index=False)

section("CHECAGEM 2 — O EXCESSO DA DOENCA, PAREADO")
piv = R.pivot_table(index=["rep", "clock"], columns="arm", values="delta").reset_index()
piv["excess"] = piv["ctl->case"] - piv["ctl->ctl"]
print(f"  {'relogio':<13}{'ctl->ctl':>11}{'ctl->caso':>12}{'excesso':>10}{'IQR excesso':>22}")
ex = {}
for c in CLOCKS3:
    g = piv[piv.clock == c]
    ex[c] = g.excess.median()
    print(f"  {c:<13}{g['ctl->ctl'].median():>+11.1%}{g['ctl->case'].median():>+12.1%}"
          f"{ex[c]:>+10.1%}   [{g.excess.quantile(.25):+.1%}, {g.excess.quantile(.75):+.1%}]")
c2 = ex["Levine2018"] > 0.01 and ex["Horvath2018"] < 0.01
print(f"\n  previsao: Levine > +1 ponto e Horvath 2018 < +1 ponto -> "
      f"{'CONFIRMADA' if c2 else 'REFUTADA'}")

section("CHECAGEM 4 — O EFEITO DA AR SOB TRES AJUSTES (REPORTADO)")
print("  efeito da AR sobre a aceleracao de idade, em anos (caso - controle, ajustado por idade)\n")
Xa = np.column_stack([np.ones(len(ra)), d["chrono"], ra.astype(float)])
print(f"  {'relogio':<13}{'sem ajuste':>12}{'coorte inteira':>16}{'so controles':>15}")
for c in CLOCKS3:
    y = d["y"][c]
    e0 = np.linalg.lstsq(Xa, y, rcond=None)[0][2]
    b_all = ridge_coefs(d["chrono"], d["C12"], y, [0.0])[0.0]
    y_all = y - (d["C12"][:, :-1] - d["C12"].mean(0)[:-1]) @ b_all
    e1 = np.linalg.lstsq(Xa, y_all, rcond=None)[0][2]
    b_ctl = ridge_coefs(d["chrono"][ictl], d["C12"][ictl], y[ictl], [0.0])[0.0]
    y_ctl = y - (d["C12"][:, :-1] - d["C12"][ictl].mean(0)[:-1]) @ b_ctl
    e2 = np.linalg.lstsq(Xa, y_ctl, rcond=None)[0][2]
    print(f"  {c:<13}{e0:>+12.2f}{e1:>+16.2f}{e2:>+15.2f}")

section("FECHAMENTO")
print(f"  2. previsao diferencial (Levine sim, Horvath 2018 nao): "
      f"{'CONFIRMADA' if c2 else 'REFUTADA'}")
