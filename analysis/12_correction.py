#!/usr/bin/env python3
"""
Stage 12 — does correcting for cell composition actually work?

This is what practitioners do: estimate each donor's blood composition from the
methylation, and subtract the part of the epigenetic age that composition
explains. Stage 7 showed that part is small in real blood — 1 to 1.6 years of
standard deviation — and stage 11 showed it cannot be designed away. So the
remaining question is whether it can be subtracted away, and whether the
subtraction survives leaving the cohort it was fitted on.

THE CORRECTION. On GSE40279, fit  clock_age ~ chronological_age + composition,
and keep only the composition coefficients. Applying them needs no age, so the
correction is usable on a sample whose age is unknown — which is the whole point
of a clock. Fitting WITH age in the model matters: composition itself drifts with
age, and a correction fitted without age would strip out real ageing along with
the blood count.

THE TEST THAT MATTERS. Apply those coefficients, unchanged, to GSE61151 and ask
whether the composition-attributable variance there actually falls. A correction
that only works where it was fitted is not a correction.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. the panel is revalidated on GSE110554's twelve mixtures, whose proportions
     are known: correlation > 0.80 and mean absolute error < 0.08, the stage 7
     bar. Nothing is corrected with an unvalidated panel
  2. physiology in both whole-blood cohorts: median neutrophil-like fraction
     between 0.40 and 0.75
  3. the two cohorts must be comparable in composition, or transporting a
     linear term between them is extrapolation rather than correction. Median
     neutrophil fraction within 0.15 of each other
  4. in-sample floor: the correction must remove the composition variance in
     GSE40279 itself. If it fails there, the fit is broken and the
     out-of-cohort result means nothing

THE EXTRAPOLATION, DECLARED. The corrected clock is also measured on purified
cells, where composition is 100% of one type. The correction learned its
coefficients where neutrophils run 50-75% and B cells 1-4%. Applying it at 100%
is far outside that range, and a linear term extrapolated fifty percentage
points is not evidence about anything except itself. It is reported because the
displacement axis is how this whole project measures exposure, and refusing to
report it would be worse — but it is not the verdict.

Usage:  .venv/bin/python analysis/12_correction.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix
from model.clocks import CLOCKS, predict
from model.deconvolution import build_panel, deconvolve, TYPES, PANEL

DATA = ROOT / "reference/data"
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(20250922)
MIX_COL = {"CD4T": "cd4t", "CD8T": "cd8t", "Bcell": "bcell",
           "NK": "nk", "Mono": "mono", "Neu": "neu"}
N_PERM = 5000


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


def r2(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return 1 - (y - X @ beta).var() / y.var()


print("carregando GSE35069 (painel + medida) ...", flush=True)
b35, m35 = read_series_matrix(DATA / "GSE35069_series_matrix.txt.gz")
m35 = m35.set_index("gsm").reindex(b35.columns)
frac35 = m35["tissue/cell type"]
donor35 = m35["title"].str.extract(r"_(\d+)$")[0]

print("carregando GSE110554 (validacao do painel) ...", flush=True)
b110, m110 = read_series_matrix(DATA / "GSE110554_series_matrix.txt.gz")
m110 = m110.set_index("gsm").reindex(b110.columns)

print("carregando GSE61151 (teste) ...", flush=True)
b61, m61 = read_series_matrix(DATA / "GSE61151_series_matrix.txt.gz")
m61 = m61.set_index("gsm").reindex(b61.columns)
age61 = pd.to_numeric(m61["agebloodtaken"], errors="coerce")

print("carregando GSE40279 (onde a correcao e ajustada, 1.2 GB) ...", flush=True)
b40, m40 = read_series_matrix(DATA / "GSE40279_series_matrix.txt.gz")
m40 = m40.set_index("gsm").reindex(b40.columns)
age40 = pd.to_numeric(m40["age (y)"], errors="coerce")

common = b35.index.intersection(b110.index).intersection(
    b61.index).intersection(b40.index)
ref = build_panel(b35, frac35, restrict_to=common)
print(f"\n  painel: {len(ref)} sondas presentes nas quatro coortes")


section("CHECAGEM 1 — O PAINEL AINDA RECUPERA AS PROPORCOES CONHECIDAS?")
is_mix = (m110["cell type"] == "MIX").to_numpy()
known = m110.loc[is_mix, [MIX_COL[t] for t in TYPES]].apply(
    pd.to_numeric, errors="coerce").to_numpy() / 100.0
est = deconvolve(ref.to_numpy(), b110.loc[ref.index].to_numpy()[:, is_mix])
r_mix = float(np.corrcoef(known.ravel(), est.ravel())[0, 1])
mae_mix = float(np.abs(known - est).mean())
c1 = (r_mix > 0.80) and (mae_mix < 0.08)
print(f"  correlacao {r_mix:.3f}  erro absoluto medio {mae_mix:.3f}  "
      f"-> {'ok' if c1 else 'FALHOU'}")
del b110


section("CHECAGENS 2 E 3 — COMPOSICAO NAS DUAS COORTES DE SANGUE TOTAL")
comp40 = pd.DataFrame(deconvolve(ref.to_numpy(), b40.loc[ref.index].to_numpy()),
                      index=b40.columns, columns=TYPES)
comp61 = pd.DataFrame(deconvolve(ref.to_numpy(), b61.loc[ref.index].to_numpy()),
                      index=b61.columns, columns=TYPES)
print(f"  {'tipo':<8}{'GSE40279 mediana':>18}{'GSE61151 mediana':>18}{'diferenca':>12}")
for t in TYPES:
    a, b_ = comp40[t].median(), comp61[t].median()
    print(f"  {t:<8}{a:>18.3f}{b_:>18.3f}{b_-a:>+12.3f}")
c2 = all(0.40 <= c["Neu"].median() <= 0.75 for c in (comp40, comp61))
c3 = abs(comp40["Neu"].median() - comp61["Neu"].median()) < 0.15
print(f"\n  2. neutrofilos em faixa adulta nas duas  -> {'ok' if c2 else 'FALHOU'}")
print(f"  3. as duas coortes sao comparaveis em composicao  -> {'ok' if c3 else 'FALHOU'}")

ages = {}
for c in CLOCKS:
    ages[c] = (predict(b40, c)[0], predict(b61, c)[0], predict(b35, c)[0])

# the best family clock from stage 10: 1000 probes ranked by age correlation,
# ridge with lambda = 1, trained on GSE40279
k40 = age40.notna().to_numpy()
X = b40.loc[common].to_numpy(dtype=np.float64).T[k40]
E = b61.loc[common].to_numpy(dtype=np.float64).T
P = b35.loc[common].to_numpy(dtype=np.float64).T
good = ~np.isnan(X).any(axis=0) & ~np.isnan(E).any(axis=0) & ~np.isnan(P).any(axis=0)
X, E, P = X[:, good], E[:, good], P[:, good]
yv = age40[k40].to_numpy()
rr = ((X - X.mean(0)).T @ (yv - yv.mean())) / (X.std(0) * yv.std() * len(yv) + 1e-12)
idx = np.argsort(-np.abs(rr))[:1000]
mx, my = X[:, idx].mean(0), yv.mean()
Xc = X[:, idx] - mx
w_, V = np.linalg.eigh(Xc @ Xc.T)
beta_f = Xc.T @ (V @ ((V.T @ (yv - my)) / (w_ + 1.0)))
ages["familia_k1000"] = (
    pd.Series((X[:, idx] - mx) @ beta_f + my, index=b40.columns[k40]),
    pd.Series((E[:, idx] - mx) @ beta_f + my, index=b61.columns),
    pd.Series((P[:, idx] - mx) @ beta_f + my, index=b35.columns))
del b40, b61

if not (c1 and c2 and c3):
    sys.exit("\n  UMA CHECAGEM FALHOU — parando antes de corrigir nada.")


section("A CORRECAO, AJUSTADA NO GSE40279 E LEVADA PARA O GSE61151")
C40 = comp40.loc[age40[k40].index, TYPES].to_numpy()
k61 = age61.notna().to_numpy()
C61 = comp61.loc[age61[k61].index, TYPES].to_numpy()
a40, a61 = age40[k40].to_numpy(), age61[k61].to_numpy()
cbar = C40.mean(axis=0)

rows = []
for name, (p40, p61, p35) in ages.items():
    y40 = p40.reindex(age40[k40].index).to_numpy()
    y61 = p61.reindex(age61[k61].index).to_numpy()

    # fit with age in the model, keep only the composition coefficients
    M = np.column_stack([np.ones(len(a40)), a40, C40[:, :-1]])
    bta, *_ = np.linalg.lstsq(M, y40, rcond=None)
    bcomp = bta[2:]

    def corrected(y, C):
        return y - (C[:, :-1] - cbar[:-1]) @ bcomp

    c40_, c61_ = corrected(y40, C40), corrected(y61, C61)

    def increment(y, a, C):
        Xa = np.column_stack([np.ones(len(a)), a])
        return r2(np.column_stack([Xa, C[:, :-1]]), y) - r2(Xa, y)

    inc40, inc40c = increment(y40, a40, C40), increment(c40_, a40, C40)
    inc61, inc61c = increment(y61, a61, C61), increment(c61_, a61, C61)
    null = np.array([increment(c61_, a61, C61[RNG.permutation(len(C61))])
                     for _ in range(N_PERM)])
    rows.append(dict(clock=name,
                     inc40=inc40, inc40c=inc40c, inc61=inc61, inc61c=inc61c,
                     null61=float(null.mean()), p61c=float((null >= inc61c).mean()),
                     r61=float(np.corrcoef(a61, y61)[0, 1]),
                     r61c=float(np.corrcoef(a61, c61_)[0, 1]),
                     mae61=float(np.abs(y61 - a61).mean()),
                     mae61c=float(np.abs(c61_ - a61).mean())))

res = pd.DataFrame(rows)
res.to_csv(OUT / "correction.csv", index=False)
print(f"  {'relogio':<15}{'R2 comp em casa':>17}{'depois':>9}"
      f"{'R2 comp FORA':>15}{'depois':>9}{'nulo':>8}{'p':>8}")
for _, r in res.iterrows():
    print(f"  {r.clock:<15}{r.inc40:>17.4f}{r.inc40c:>9.4f}"
          f"{r.inc61:>15.4f}{r.inc61c:>9.4f}{r.null61:>8.4f}{r.p61c:>8.3f}")
c4 = bool((res.inc40c < res.inc40).all())
print(f"\n  4. em casa a correcao remove mesmo a variancia de composicao  "
      f"-> {'ok' if c4 else 'FALHOU'}")
print("""
  'R2 comp' e quanto a composicao explica da idade epigenetica ALEM da idade
  cronologica. 'em casa' e no GSE40279, onde os coeficientes foram ajustados --
  ali a queda e garantida e nao prova nada. A coluna que vale e 'FORA'.
  O nulo embaralha a composicao entre pessoas: seis preditores compram R2 de
  graca, e ele mede quanto.""")


section("O QUE A CORRECAO CUSTOU EM ACURACIA")
print(f"  {'relogio':<15}{'r antes':>10}{'r depois':>11}{'erro antes':>13}{'erro depois':>14}")
for _, r in res.iterrows():
    print(f"  {r.clock:<15}{r.r61:>+10.3f}{r.r61c:>+11.3f}"
          f"{r.mae61:>12.1f}a{r.mae61c:>13.1f}a")


section("EXTRAPOLACAO — O DESLOCAMENTO EM CELULAS PURAS")
print("""  A correcao aprendeu onde neutrofilo vai de 50 a 75% e linfocito B de 1 a 4%.
  Aqui ela e aplicada a amostras que sao 100% de um tipo -- cinquenta pontos
  percentuais fora da faixa onde os coeficientes foram estimados. Isto NAO e
  veredito sobre a correcao; e o que acontece quando se extrapola um termo
  linear, relatado porque o eixo de deslocamento e como este projeto mede
  exposicao o tempo todo.\\n""")
sel = frac35.isin(list(PANEL.values())).to_numpy()
onehot = np.zeros((int(sel.sum()), len(TYPES)))
fr = frac35[sel].to_numpy()
for j, t in enumerate(TYPES):
    onehot[fr == PANEL[t], j] = 1.0
dono = donor35.to_numpy()[sel]
print(f"  {'relogio':<15}{'desloc. antes':>15}{'desloc. depois':>16}")
drows = []
for name, (_, _, p35) in ages.items():
    v = p35.to_numpy()[sel]
    for tag, vals in (("antes", v),
                      ("depois", v - (onehot[:, :-1] - cbar[:-1]) @
                       np.linalg.lstsq(np.column_stack(
                           [np.ones(len(a40)), a40, C40[:, :-1]]),
                           ages[name][0].reindex(age40[k40].index).to_numpy(),
                           rcond=None)[0][2:])):
        w = pd.DataFrame({"d": dono, "f": fr, "v": vals}).pivot_table(
            index="d", columns="f", values="v")
        drows.append(dict(clock=name, when=tag, disp=float(w.std(axis=1).mean())))
d = pd.DataFrame(drows).pivot(index="clock", columns="when", values="disp")
d.to_csv(OUT / "correction_displacement.csv")
for name in ages:
    print(f"  {name:<15}{d.loc[name,'antes']:>14.2f}a{d.loc[name,'depois']:>15.2f}a")
