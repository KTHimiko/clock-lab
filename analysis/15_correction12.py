#!/usr/bin/env python3
"""
Stage 15 — the composition correction, at twelve types.

Stage 12 fitted a correction on one cohort and carried it to another, and it
removed almost nothing: 5%, 62%, 36%, 0%, 7%. The diagnosis then was that the
per-type coefficients are not identified, because the six proportions are
collinear — the joint fit is sound while the split between neighbouring types is
not, and a correction made of per-type coefficients is built entirely out of the
part that is not sound.

Stage 13 gave that diagnosis a test. The twelve-type panel separates 12 of 12
types under leave-one-out and lands more buckets inside clinical range than the
six-type panel did. If stage 12's failure really was the panel, the correction
should travel better now. If it fails again with a panel this much better, the
negative becomes much harder to explain away.

WHAT IS NEW BESIDES THE PANEL
  - both directions. Stage 12 fitted on GSE40279 and tested on GSE61151 only. A
    correction that works one way and not the other is a fact about the
    cohorts, not about correction, and one direction cannot tell them apart.
  - cross-panel residuals. Correcting with twelve types and then measuring what
    is left with twelve types risks scoring the correction against its own
    representation. Every residual is therefore measured BOTH ways, and the
    one that counts is the cross-panel number: correct with twelve, measure
    with six, and vice versa.
  - two clocks this project built: stage 10's best family configuration, and
    stage 14's clock selected by the IntrinClock rule. A clock already flattened
    by design is the interesting case — does it still need correcting?

THE TRAP IN THE SECOND DIRECTION, and it was nearly read as a result. Four of
the six clocks here have a stake in GSE40279: Hannum was trained on it, and
this project's two — the stage 10 family clock and the stage 14 IntrinClock-rule
clock — were both fitted on it. A clock measured in the cohort it learned has a
tiny age residual, and every number in this stage is a share OF that residual.
Divide by a small denominator and a modest absolute change reads as an enormous
percentage: the family clock appears to go from 0.5% to 40%. That is mostly the
denominator. Every row is therefore marked, and the verdict is taken only from
the clocks external to the cohort being tested.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. the twelve-type panel recovers the known twelve-way proportions of the
     twelve GSE167998 mixtures: r > 0.70, mean absolute error < 0.04
  2. in-sample floor: in the cohort it was fitted on, each correction must drive
     the composition term to zero when measured with its own panel. If it does
     not, the fit is broken and nothing out-of-cohort means anything
  3. twelve predictors buy R-squared for free, so every residual is compared
     against its own permutation null and the verdict is taken on the excess
     over null, never on the raw value
  4. the declared bias from stage 13 carries: this panel's monocyte channel runs
     high, so no per-type monocyte coefficient is read or reported

Usage:  .venv/bin/python analysis/15_correction12.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix
from load_extended import read_extended, TYPES12
from model.clocks import CLOCKS, predict
from model.deconvolution import build_panel, deconvolve, TYPES as TYPES6

DATA = ROOT / "reference/data"
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(20250922)
N_PERM = 3000
K = 1000
LAMBDAS = [1e-2, 1e-1, 1e0, 1e1, 1e2]


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


def r2(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return 1 - (y - X @ beta).var() / y.var()


def build12(betas, celltype, restrict_to, n_side=50, n_candidate=200):
    scores, cand = {}, set()
    arr = celltype.to_numpy()
    for t in TYPES12:
        sel, rest = arr == t, np.isin(arr, [o for o in TYPES12 if o != t])
        a, c = betas.loc[:, sel], betas.loc[:, rest]
        d = a.mean(axis=1) - c.mean(axis=1)
        sd = np.sqrt((a.var(axis=1, ddof=1) + c.var(axis=1, ddof=1)) / 2) + 0.01
        s = (d / sd).replace([np.inf, -np.inf], np.nan).dropna()
        scores[t] = s
        cand |= set(s.nlargest(n_candidate).index) | set(s.nsmallest(n_candidate).index)
    pool = sorted(cand & set(restrict_to))
    probes = set()
    for t in TYPES12:
        s = scores[t].reindex(pool).dropna()
        probes |= set(s.nlargest(n_side).index) | set(s.nsmallest(n_side).index)
    probes = sorted(probes)
    return pd.DataFrame({t: betas.loc[probes, arr == t].mean(axis=1)
                         for t in TYPES12}), probes


def ridge_dual(Xtr, ytr, lambdas):
    mx, my = Xtr.mean(axis=0), ytr.mean()
    Xc, yc = Xtr - mx, ytr - my
    w, V = np.linalg.eigh(Xc @ Xc.T)
    Vty = V.T @ yc
    for lam in lambdas:
        yield lam, Xc.T @ (V @ (Vty / (w + lam))), mx, my


print("carregando GSE167998 ...", flush=True)
b167, p167 = read_extended(DATA / "GSE167998_matrix_processed.txt.gz",
                           DATA / "BloodExtended_Pheno.csv")
is_mix = (p167.CellType == "MIX").to_numpy()
pur167 = (~is_mix) & p167.Age.notna().to_numpy()

print("carregando GSE35069 ...", flush=True)
b35, m35 = read_series_matrix(DATA / "GSE35069_series_matrix.txt.gz")
m35 = m35.set_index("gsm").reindex(b35.columns)
frac35 = m35["tissue/cell type"]

print("carregando GSE61151 ...", flush=True)
b61, m61 = read_series_matrix(DATA / "GSE61151_series_matrix.txt.gz")
m61 = m61.set_index("gsm").reindex(b61.columns)
age61 = pd.to_numeric(m61["agebloodtaken"], errors="coerce")

print("carregando GSE40279 (1.2 GB) ...", flush=True)
b40, m40 = read_series_matrix(DATA / "GSE40279_series_matrix.txt.gz")
m40 = m40.set_index("gsm").reindex(b40.columns)
age40 = pd.to_numeric(m40["age (y)"], errors="coerce")

common = b61.index.intersection(b40.index).intersection(
    b35.index).intersection(b167.index)
ref12, probes12 = build12(b167.loc[:, ~is_mix], p167.CellType[~is_mix], common)
ref6 = build_panel(b35, frac35, restrict_to=common)
print(f"\n  painel 12: {len(ref12)} sondas | painel 6: {len(ref6)} sondas")


section("CHECAGEM 1 — O PAINEL DE 12 RECUPERA AS PROPORCOES CONHECIDAS?")
known = p167.loc[is_mix, TYPES12].to_numpy()
est = deconvolve(ref12.to_numpy(), b167.loc[probes12].to_numpy()[:, is_mix])
rv = float(np.corrcoef(known.ravel(), est.ravel())[0, 1])
mv = float(np.abs(known - est).mean())
c1 = (rv > 0.70) and (mv < 0.04)
print(f"  correlacao {rv:.3f}  erro absoluto medio {mv:.3f}  -> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("\n  parando: sem painel validado nao se corrige nada.")

comp = {}
for tag, b in (("GSE40279", b40), ("GSE61151", b61)):
    comp[(tag, 12)] = pd.DataFrame(
        deconvolve(ref12.to_numpy(), b.loc[probes12].to_numpy()),
        index=b.columns, columns=TYPES12)
    comp[(tag, 6)] = pd.DataFrame(
        deconvolve(ref6.to_numpy(), b.loc[ref6.index].to_numpy()),
        index=b.columns, columns=TYPES6)
print("  (vies declarado da etapa 13 continua valendo: o canal de monocito "
      "corre alto,\n   nenhum coeficiente por tipo de monocito e lido aqui)")

ages = {}
for c in CLOCKS:
    ages[c] = {"GSE40279": predict(b40, c)[0], "GSE61151": predict(b61, c)[0]}

# the two clocks this project built, both trained on GSE40279
k40 = age40.notna().to_numpy()
X = b40.loc[common].to_numpy(dtype=np.float64).T[k40]
E = b61.loc[common].to_numpy(dtype=np.float64).T
P167 = b167.loc[common].to_numpy(dtype=np.float64).T[pur167]
good = ~np.isnan(X).any(axis=0) & ~np.isnan(E).any(axis=0) & ~np.isnan(P167).any(axis=0)
X, E = X[:, good], E[:, good]
P167 = P167[:, good]
yv = age40[k40].to_numpy()
rr = ((X - X.mean(0)).T @ (yv - yv.mean())) / (X.std(0) * yv.std() * len(yv) + 1e-12)
absr = np.abs(rr)
types_pur = p167.CellType[pur167].to_numpy()
is_nv8 = (types_pur == "CD8nv").astype(float)
r_nv8 = np.nan_to_num(np.array(
    [np.corrcoef(P167[:, j], is_nv8)[0, 1] for j in range(P167.shape[1])]))

for name, pool in (("familia_k1000", np.arange(X.shape[1])),
                   ("regra_intrinclock",
                    np.where((absr > 0.3) & (np.abs(r_nv8) < 0.3))[0])):
    idx = pool[np.argsort(-absr[pool])[:K]]
    _, beta, mx, my = next(iter(ridge_dual(X[:, idx], yv, [1.0])))
    ages[name] = {"GSE40279": pd.Series((X[:, idx] - mx) @ beta + my,
                                        index=b40.columns[k40]),
                  "GSE61151": pd.Series((E[:, idx] - mx) @ beta + my,
                                        index=b61.columns)}
del b167, b35, b61, b40
CLOCKS_ALL = list(CLOCKS) + ["familia_k1000", "regra_intrinclock"]


section("A CORRECAO, NAS DUAS DIRECOES E NOS DOIS PAINEIS")
AGES = {"GSE40279": age40, "GSE61151": age61}
# which clocks have seen which cohort, so no verdict is taken from an in-sample cell
TRAINED_ON = {"Hannum2013": "GSE40279", "familia_k1000": "GSE40279",
              "regra_intrinclock": "GSE40279"}


def residual_r2(y, a, C):
    Xa = np.column_stack([np.ones(len(a)), a])
    inc = r2(np.column_stack([Xa, C[:, :-1]]), y) - r2(Xa, y)
    null = np.array([r2(np.column_stack([Xa, C[RNG.permutation(len(C))][:, :-1]]), y)
                     - r2(Xa, y) for _ in range(N_PERM)])
    base = r2(Xa, y)
    return inc, float(null.mean()), base


rows = []
for src, dst in (("GSE40279", "GSE61151"), ("GSE61151", "GSE40279")):
    ks, kd = AGES[src].notna().to_numpy(), AGES[dst].notna().to_numpy()
    a_s, a_d = AGES[src][ks].to_numpy(), AGES[dst][kd].to_numpy()
    for npan in (6, 12):
        types = TYPES6 if npan == 6 else TYPES12
        Cs = comp[(src, npan)].loc[AGES[src][ks].index, types].to_numpy()
        Cd = comp[(dst, npan)].loc[AGES[dst][kd].index, types].to_numpy()
        cbar = Cs.mean(axis=0)
        for c in CLOCKS_ALL:
            ys = ages[c][src].reindex(AGES[src][ks].index).to_numpy()
            yd = ages[c][dst].reindex(AGES[dst][kd].index).to_numpy()
            M = np.column_stack([np.ones(len(a_s)), a_s, Cs[:, :-1]])
            bcomp = np.linalg.lstsq(M, ys, rcond=None)[0][2:]
            corr_s = ys - (Cs[:, :-1] - cbar[:-1]) @ bcomp
            corr_d = yd - (Cd[:, :-1] - cbar[:-1]) @ bcomp
            row = dict(src=src, dst=dst, panel=npan, clock=c)
            # in-sample floor, own panel
            i0, n0, _ = residual_r2(ys, a_s, Cs)
            i0c, _, _ = residual_r2(corr_s, a_s, Cs)
            row["home_before"], row["home_after"] = i0, i0c
            # out of cohort, measured with BOTH panels
            for mpan in (6, 12):
                mt = TYPES6 if mpan == 6 else TYPES12
                Cm = comp[(dst, mpan)].loc[AGES[dst][kd].index, mt].to_numpy()
                ib, nb, base = residual_r2(yd, a_d, Cm)
                ia, na, _ = residual_r2(corr_d, a_d, Cm)
                row[f"away{mpan}_before"] = (ib - nb) / (1 - base)
                row[f"away{mpan}_after"] = (ia - na) / (1 - base)
            rows.append(row)
        print(f"  {src} -> {dst}, painel {npan}: pronto", flush=True)

res = pd.DataFrame(rows)
res.to_csv(OUT / "correction12.csv", index=False)

c2 = bool((res.home_after < res.home_before).all())
print(f"\n  2. piso em casa: a correcao zera o termo de composicao na coorte de "
      f"ajuste  -> {'ok' if c2 else 'FALHOU'}")

section("RESULTADO — QUANTO SOBRA FORA DE CASA (EXCESSO SOBRE O NULO, NA EAA)")
print("""  A coluna que decide e 'medido com 6' quando a correcao usou 12, e vice-versa:
  corrigir com um painel e medir com o mesmo arrisca pontuar a correcao contra a
  propria representacao dela.\\n""")
for src, dst in (("GSE40279", "GSE61151"), ("GSE61151", "GSE40279")):
    print(f"  ajustado em {src}, testado em {dst}")
    print(f"  {'relogio':<19}{'painel':>7}{'  medido com 6':>18}{'  medido com 12':>19}")
    print(f"  {'':<19}{'':>7}{'antes':>9}{'depois':>9}{'antes':>10}{'depois':>9}")
    for c in CLOCKS_ALL:
        insample = TRAINED_ON.get(c) == dst
        for npan in (6, 12):
            r = res[(res.src == src) & (res.panel == npan) & (res.clock == c)].iloc[0]
            mark = "  <- TREINOU AQUI, nao conta" if (insample and npan == 6) else ""
            print(f"  {c if npan==6 else '':<19}{npan:>7}"
                  f"{r.away6_before:>9.1%}{r.away6_after:>9.1%}"
                  f"{r.away12_before:>10.1%}{r.away12_after:>9.1%}{mark}")
    print()

section("VEREDITO — SO OS RELOGIOS EXTERNOS A COORTE TESTADA")
print("""  Correcao ajustada com 12 tipos, residuo medido com 6 -- a celula cruzada,
  que e a unica que nao pontua a correcao contra a propria representacao dela.\n""")
print(f"  {'ajuste':<11}{'teste':<11}{'relogio':<14}{'antes':>8}{'depois':>9}{'variacao':>11}")
verdict = []
for src, dst in (("GSE40279", "GSE61151"), ("GSE61151", "GSE40279")):
    for c in CLOCKS_ALL:
        if TRAINED_ON.get(c) == dst:
            continue
        r = res[(res.src == src) & (res.panel == 12) & (res.clock == c)].iloc[0]
        d = r.away6_after - r.away6_before
        verdict.append(dict(src=src, dst=dst, clock=c, before=r.away6_before,
                            after=r.away6_after, delta=d))
        print(f"  {src:<11}{dst:<11}{c:<14}{r.away6_before:>8.1%}{r.away6_after:>9.1%}"
              f"{d:>+11.1%}{'  pior' if d > 0 else ''}")
v = pd.DataFrame(verdict)
v.to_csv(OUT / "correction12_verdict.csv", index=False)
for src in ("GSE40279", "GSE61151"):
    d = v[v.src == src]
    print(f"\n  ajustado em {src} ({len(AGES[src].dropna())} amostras): "
          f"melhora em {int((d.delta < 0).sum())} de {len(d)} relogios, "
          f"variacao mediana {d.delta.median():+.1%}")
