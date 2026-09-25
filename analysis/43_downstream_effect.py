#!/usr/bin/env python3
"""
Stage 43 — the question a reader actually asks: does transport change the answer?

Every stage so far scored the correction by how much composition it leaves in
age acceleration. Nobody reports that number. What gets reported is an
association: smoking accelerates ageing by X years, this disease by Y. Stage 31
showed, for one cohort and one clock, that the choice of reference population
moved the arthritis effect from -0.74 to -1.39 years. This stage asks it
systematically, and against the right comparator.

For a target cohort with an exposure E, three estimates of the same quantity —
the coefficient of E in (clock age) ~ age + E:
  within      the composition correction fitted inside the target. This is what
              the target's own data support, and the practice the field treats
              as correct; it is the reference, not a truth.
  none        no composition correction at all.
  transported coefficients fitted on another cohort at size n and carried in.
The error of interest is |transported - within| in years: how far a reported
association moves because the correction came from elsewhere. A transport is
worse than useless when that error exceeds |none - within|, the error of simply
not correcting.

Targets and exposures (all already cached):
  GSE50660   smoking, ever (former+current) versus never      blood, 464
  GSE42861   rheumatoid arthritis, case versus control        blood, 689
  GSE232891  inflammatory bowel disease, case versus control  saliva, 552
  GSE232332  oesophageal cancer, case versus control          saliva, 265
Sources: the other blood cohorts for blood targets, the other saliva cohorts for
saliva targets, never pairing GSE232891 with GSE232332 (possible shared people),
and never scoring a clock on a pair containing a cohort it trained on.

PRE-REGISTERED (written before any effect was recomputed)
  1. the caches load and every target's exposure has both groups. Hard stop.
  2. TRANSPORT IS WORSE THAN NO CORRECTION, in the reported answer: at n = 40,
     the median |transported - within| exceeds the median |none - within| across
     (target x source x clock) cells.
  3. IT CHANGES CONCLUSIONS: at n = 40, in at least one cell the transported
     estimate differs from the within-cohort one by more than 1 year, or flips
     its sign.
  4. THE PENALTY HELPS HERE TOO: at alpha = 3 the median |transported - within|
     at n = 40 is below its unpenalised value.
  Reported without a bar: the same at each target's full fitting size, and the
  tail of absurd estimates, which a near-singular fit at n = 40 can produce. That
  tail is solver-dependent, so it is reported for exact least squares and for
  least squares with small singular values dropped, as numpy's lstsq does.

Usage:  .venv/bin/python analysis/43_downstream_effect.py
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


def header_fields(gse):
    samples, fields = None, {}
    with gzip.open(ROOT / "reference/data" / f"{gse}_series_matrix.txt.gz", "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            parts = [v.strip('"') for v in line.rstrip("\n").split("\t")]
            if parts[0] == "!Sample_geo_accession":
                samples = parts[1:]
            elif parts[0] == "!Sample_characteristics_ch1":
                for j, cell in enumerate(parts[1:]):
                    k, _, v = cell.partition(":")
                    fields.setdefault(k.strip(), [None] * len(samples))[j] = v.strip()
    return pd.DataFrame(fields, index=samples)


BLOOD = ["GSE40279", "GSE61151", "GSE50660", "GSE42861", "GSE132203", "GSE55763"]
SAL = ["GSE232891", "GSE232332", "GSE78874", "GSE149747"]
FIT_S = ["Epi", "Fib", "B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino"]
MEAS_S = ["Epi", "Fib", "IC"]
AGE = ["Levine2018", "Horvath2018"]
TRAINED_ON = {"Hannum2013": {"GSE40279"}, "Horvath2013": {"GSE40279"}}
USABLE = {t: pd.read_csv(CACHE / f"{t}_usable_clocks.csv").iloc[:, 0].tolist()
          for t in ("GSE132203", "GSE55763", "GSE232891", "GSE232332", "GSE78874")}


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


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


def ridge_coefs(age, C, y, alphas, rcond=None):
    """alpha = 0 is exact least squares. With rcond set, singular values below
    rcond * s_max are dropped, which is what numpy's lstsq does by default; the
    two differ only on a near-singular design, and that difference is itself
    reported in this stage."""
    Ct = partial_out(age, C[:, :-1]); yt = partial_out(age, y.reshape(-1, 1)).ravel()
    mu, sd = Ct.mean(0), Ct.std(0) + 1e-12
    Z = (Ct - mu) / sd
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    scale = float((s ** 2).mean()); Uty = U.T @ yt
    keep = s > (rcond * s.max() if rcond else 0.0)
    out = {}
    for a in alphas:
        f = np.where(keep, s * Uty / (s ** 2 + a * scale), 0.0)
        out[a] = (Vt.T @ f) / sd
    return out


section("CHECAGEM 1 — COORTES E EXPOSICOES")
D = {}
for t in BLOOD:
    a = pd.read_csv(CACHE / f"{t}_ages.csv", index_col=0); a = a[a.chrono.notna()]
    D[t] = dict(tis="blood", chrono=a.chrono.to_numpy(float),
                y={c: a[c].to_numpy(float) for c in CLOCKS},
                C=pd.read_csv(CACHE / f"{t}_comp12.csv", index_col=0).loc[a.index, TYPES12].to_numpy(),
                M=pd.read_csv(CACHE / f"{t}_comp6.csv", index_col=0).loc[a.index, TYPES6].to_numpy(),
                idx=a.index)
for t in SAL[:3]:
    a = pd.read_csv(CACHE / f"{t}_ages.csv", index_col=0)
    D[t] = dict(tis="saliva", chrono=a.chrono.to_numpy(float),
                y={c: a[c].to_numpy(float) for c in AGE},
                C=pd.read_csv(CACHE / f"{t}_compfit.csv", index_col=0).loc[a.index, FIT_S].to_numpy(),
                M=pd.read_csv(CACHE / f"{t}_compmeas.csv", index_col=0).loc[a.index, MEAS_S].to_numpy(),
                idx=a.index, group=a.group.to_numpy())
g = pd.read_csv(CACHE / "GSE149747_saliva.csv", index_col=0)
g = g[g.tp == "Before"].dropna(subset=["age"])
D["GSE149747"] = dict(tis="saliva", chrono=g.age.to_numpy(float),
                      y={c: g[c].to_numpy(float) for c in AGE},
                      C=None, M=g[MEAS_S].to_numpy(), idx=g.index)

EXPO = {}
h = header_fields("GSE50660").reindex(D["GSE50660"]["idx"])
sm = pd.to_numeric(h["smoking (0, 1 and 2, which represent never, former and current smokers)"], errors="coerce")
EXPO["GSE50660"] = ("fumo (alguma vez vs nunca)", (sm > 0).to_numpy(), sm.notna().to_numpy())
h = header_fields("GSE42861").reindex(D["GSE42861"]["idx"])
ra = h["disease state"].str.contains("rheumatoid", case=False)
EXPO["GSE42861"] = ("artrite reumatoide", ra.to_numpy(), ra.notna().to_numpy())
for t, lab, pos in (("GSE232891", "doenca inflamatoria intestinal", "Patient"),
                    ("GSE232332", "cancer de esofago", "Cancer")):
    grp = pd.Series(D[t]["group"])
    EXPO[t] = (lab, grp.str.contains(pos).to_numpy(), grp.notna().to_numpy())
ok1 = True
for t, (lab, e, k) in EXPO.items():
    n1, n0 = int((e & k).sum()), int((~e & k).sum())
    good = n1 >= 20 and n0 >= 20
    ok1 &= good
    print(f"  {t:<11}{lab:<32}expostos {n1:>4}  nao {n0:>4} -> {'ok' if good else 'FALHOU'}")
if not ok1:
    sys.exit("  parando.")


def eff(y, age, e, k):
    """coeficiente da exposicao em y ~ idade + exposicao, em anos"""
    X = np.column_stack([np.ones(int(k.sum())), age[k], e[k].astype(float)])
    return float(np.linalg.lstsq(X, y[k], rcond=None)[0][2])


def corrected(d, c, b, cbar):
    return d["y"][c] - (d["C"][:, :-1] - cbar[:-1]) @ b


def clocks_for(src, dst):
    out = []
    for c in AGE:
        if TRAINED_ON.get(c, set()) & {src, dst}:
            continue
        if any(t in USABLE and c not in USABLE[t] for t in (src, dst)):
            continue
        out.append(c)
    return out


SOURCES = {t: [s for s in (BLOOD if D[t]["tis"] == "blood" else SAL)
               if s != t and D[s]["C"] is not None
               and not (t in ("GSE232891", "GSE232332") and s in ("GSE232891", "GSE232332"))]
           for t in EXPO}
section("CHECAGEM 2, 3 E 4 — O EFEITO ESTIMADO SOB CADA CORRECAO")
rows = []
for t, (lab, e, k) in EXPO.items():
    d = D[t]
    for c in AGE:
        if t in USABLE and c not in USABLE[t]:
            continue
        e_none = eff(d["y"][c], d["chrono"], e, k)
        b_in = ridge_coefs(d["chrono"], d["C"], d["y"][c], [0.0])[0.0]
        e_within = eff(corrected(d, c, b_in, d["C"].mean(0)), d["chrono"], e, k)
        rows.append(dict(target=t, clock=c, src="-", n=0, kind="none", est=e_none, within=e_within))
        for s in SOURCES[t]:
            if c not in clocks_for(s, t):
                continue
            sd = D[s]; NS = len(sd["chrono"])
            for n in (40, NS):
                reps = 1 if n >= NS else 30
                for r in range(reps):
                    idx = (np.arange(NS) if n >= NS
                           else stratified_draw(sd["chrono"], n, np.random.default_rng(4300 + 13 * r + n)))
                    co = ridge_coefs(sd["chrono"][idx], sd["C"][idx], sd["y"][c][idx], [0.0, 3.0])
                    ct = ridge_coefs(sd["chrono"][idx], sd["C"][idx], sd["y"][c][idx], [0.0],
                                     rcond=np.finfo(float).eps * max(len(idx), sd["C"].shape[1] - 1))
                    cb = sd["C"][idx].mean(0)
                    for key, bb in (("ols", co[0.0]), ("a3", co[3.0]), ("ols_trunc", ct[0.0])):
                        rows.append(dict(target=t, clock=c, src=s, n=("40" if n < NS else "cheio"),
                                         kind=key, rep=r,
                                         est=eff(corrected(d, c, bb, cb), d["chrono"], e, k),
                                         within=e_within))
    print(f"  {t}: pronto", flush=True)
R = pd.DataFrame(rows)
R["err"] = (R.est - R.within).abs()
R.to_csv(RES / "downstream_effect.csv", index=False)

none = R[R.kind == "none"].groupby(["target", "clock"]).err.median()
print(f"\n  {'alvo':<11}{'relogio':<13}{'sem correcao':>14}{'dentro':>9}{'transp. n=40':>14}{'a=3 n=40':>11}{'transp. cheio':>15}")
for (t, c), g in R[R.kind != "none"].groupby(["target", "clock"]):
    w = g.within.iloc[0]
    m40 = g[(g.n == "40") & (g.kind == "ols")].err.median()
    a40 = g[(g.n == "40") & (g.kind == "a3")].err.median()
    mf = g[(g.n == "cheio") & (g.kind == "ols")].err.median()
    print(f"  {t:<11}{c:<13}{none[(t, c)]:>+14.2f}{w:>+9.2f}{m40:>+14.2f}{a40:>+11.2f}{mf:>+15.2f}")

cells40 = R[(R.n == "40") & (R.kind == "ols")].groupby(["target", "src", "clock"]).err.median()
cells40a = R[(R.n == "40") & (R.kind == "a3")].groupby(["target", "src", "clock"]).err.median()
nn = R[R.kind == "none"].set_index(["target", "clock"]).err
none_cells = pd.Series({ix: nn[(ix[0], ix[2])] for ix in cells40.index})
c2 = cells40.median() > none_cells.median()
print(f"\n  2. n=40: |transportado - dentro| mediano {cells40.median():.2f} ano contra "
      f"|sem correcao - dentro| {none_cells.median():.2f} -> {'ok' if c2 else 'FALHOU'}")
big = R[(R.n == "40") & (R.kind == "ols")]
flip = big[(np.sign(big.est) != np.sign(big.within)) | (big.err > 1.0)]
c3 = len(flip) > 0
print(f"  3. celulas com erro > 1 ano ou sinal trocado: {flip.groupby(['target','src','clock']).ngroups} "
      f"de {big.groupby(['target','src','clock']).ngroups} -> {'ok' if c3 else 'FALHOU'}")
for (t, s, c), gg in flip.groupby(["target", "src", "clock"]):
    print(f"     {t} <- {s[3:]:<8}{c:<13} dentro {gg.within.iloc[0]:+.2f}  transportado {gg.est.median():+.2f} ano")
c4 = cells40a.median() < cells40.median()
print(f"  4. alpha=3 reduz o erro em n=40: {cells40a.median():.2f} contra {cells40.median():.2f}"
      f" -> {'ok' if c4 else 'FALHOU'}")
section("CAUDA: ESTIMATIVAS ABSURDAS EM N=40 (REPORTADO, SEM BARRA)")
tis = {"GSE50660": "sangue", "GSE42861": "sangue", "GSE232891": "saliva", "GSE232332": "saliva"}
T = R[R.n == "40"].assign(tis=lambda x: x.target.map(tis))
print(f"  {'solucao':<12}{'sorteios':>9}{'|est|>10 anos':>15}{'pior |est|':>14}{'erro mediano':>14}")
for kind, lab in (("ols", "MQ exato"), ("ols_trunc", "MQ truncado"), ("a3", "ridge a=3")):
    g = T[T.kind == kind]
    print(f"  {lab:<12}{len(g):>9}{(g.est.abs() > 10).mean():>14.1%}{g.est.abs().max():>14.3g}{g.err.median():>14.2f}")
for kind in ("ols", "ols_trunc", "a3"):
    g = T[T.kind == kind]
    print(f"    {kind:<10} por tecido, |est|>10: " + "  ".join(
        f"{t} {gg.est.abs().gt(10).mean():.0%}" for t, gg in g.groupby("tis")))

section("FECHAMENTO")
for i, (lab, ok) in enumerate([("transporte pior que nao corrigir", c2),
                                ("muda conclusoes", c3), ("penalidade ajuda", c4)], start=2):
    print(f"  {i}. {lab}: {'ok' if ok else 'FALHOU'}")
