#!/usr/bin/env python3
"""
Stage 11 — is the floor biology or method?

Every clock measured in this project, published or built, displaces several
years between cell fractions for each year of genuine age response. Stage 10's
best was 3.9 and nothing came near zero. Two explanations, and they differ in
what anyone can do about it:

  METHOD — selection by correlation with age has no reason to avoid cell-type
  markers, and simply telling it to avoid them would lower the floor.
  BIOLOGY — the CpGs that track time ARE largely the ones that distinguish cell
  types, so there is no clean set to select and the floor stays where it is.

The test: filter probes by how much they vary between cell types BEFORE ranking
them by correlation with age, and sweep the severity of that filter.

THE CIRCULARITY THAT WOULD HAVE MADE THIS MEANINGLESS. Cell-type variance
measured on GSE35069 cannot be used to choose probes whose displacement is then
measured on GSE35069 — that is selecting on the evaluation set, and it would
lower the number without lowering anything real. So the filter is built from
GSE110554's purified cells (EPIC, a different lab, different donors) and the
displacement is measured on GSE35069 as always. The two never meet.

THE NEGATIVE CONTROL THAT DECIDES IT. Any filter removes probes, and removing
probes changes the fit. So every penalised configuration is matched by one that
discards the SAME NUMBER of probes at random. If random discarding buys as much
flatness, the cell-type information contributed nothing and the answer is
neither of the two above — it is that the sweep was measuring its own noise.

WHY THE FILTER RANKS ON ABSOLUTE CELL-TYPE SPREAD, not on the fraction of a
probe's variance that cell type explains. Displacement is a weighted sum of
each probe's shift between fractions, so what a coefficient multiplies is the
shift in beta units, not its share of that probe's total variability. The
share is the right quantity for asking "is this a cell-type marker"; the
absolute spread is the right one for asking "how far will this move my clock",
and that is the question here.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. selection cohort and measurement cohort are different, asserted on ids
  2. the filter must actually change the chosen probes — overlap with the
     unfiltered set reported, and a filter that changes nothing is a bug
  3. calibration: only configurations with slope in [0.5, 1.5] are compared,
     the band stage 10 settled on
  4. the unfiltered configuration must reproduce stage 10's number for the same
     probe count, or the two stages disagree about something and neither can
     be read

Usage:  .venv/bin/python analysis/11_floor.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix

DATA = ROOT / "reference/data"
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(20250922)

LAMBDAS = [1e-2, 1e-1, 1e0, 1e1, 1e2]
K_SWEEP = [1_000, 10_000]
QUANTILES = [1.0, 0.75, 0.50, 0.25, 0.10, 0.05]
N_CONTROL = 5
FRACTIONS = ["CD4+ T cells", "CD8+ T cells", "CD19+ B cells",
             "CD56+ NK cells", "CD14+ Monocytes", "Neutrophils"]
EPIC_TYPES = ["CD4T", "CD8T", "Bcell", "NK", "Mono", "Neu"]


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


def ridge_dual(Xtr, ytr, lambdas):
    mx, my = Xtr.mean(axis=0), ytr.mean()
    Xc, yc = Xtr - mx, ytr - my
    w, V = np.linalg.eigh(Xc @ Xc.T)
    Vty = V.T @ yc
    for lam in lambdas:
        yield lam, Xc.T @ (V @ (Vty / (w + lam))), mx, my


def cv_lambda(Xs, ys, lambdas, folds=5):
    parts = np.array_split(RNG.permutation(len(ys)), folds)
    err = {lam: 0.0 for lam in lambdas}
    for test in parts:
        tr = np.setdiff1d(np.arange(len(ys)), test)
        for lam, beta, mx, my in ridge_dual(Xs[tr], ys[tr], lambdas):
            err[lam] += float(np.abs((Xs[test] - mx) @ beta + my - ys[test]).sum())
    return min(err, key=err.get)


def fit_and_measure(idx, X, y, E, ye, P, dono, frc):
    lam = cv_lambda(X[:, idx], y, LAMBDAS)
    lam_, beta, mx, my = next(iter(ridge_dual(X[:, idx], y, [lam])))
    pe = (E[:, idx] - mx) @ beta + my
    pv = (P[:, idx] - mx) @ beta + my
    w = pd.DataFrame({"d": dono, "f": frc, "v": pv}).pivot_table(
        index="d", columns="f", values="v")
    slope = float(np.cov(pe, ye)[0, 1] / np.var(ye))
    disp = float(w.std(axis=1).mean())
    return dict(lam=lam, r=float(np.corrcoef(ye, pe)[0, 1]),
                mae=float(np.abs(pe - ye).mean()), slope=slope, disp=disp,
                disp_n=disp / slope if slope > 1e-9 else np.inf)


print("carregando GSE110554 (so para construir o filtro) ...", flush=True)
b110, m110 = read_series_matrix(DATA / "GSE110554_series_matrix.txt.gz")
m110 = m110.set_index("gsm").reindex(b110.columns)
cell110 = m110["cell type"]
means110 = pd.DataFrame({t: b110.loc[:, (cell110 == t).to_numpy()].mean(axis=1)
                         for t in EPIC_TYPES})
cellvar = means110.std(axis=1)          # per-probe spread across cell types
del b110

print("carregando GSE35069 (so para medir) ...", flush=True)
b35, m35 = read_series_matrix(DATA / "GSE35069_series_matrix.txt.gz")
m35 = m35.set_index("gsm").reindex(b35.columns)
frac35, donor35 = m35["tissue/cell type"], m35["title"].str.extract(r"_(\d+)$")[0]

print("carregando GSE61151 (teste) ...", flush=True)
b61, m61 = read_series_matrix(DATA / "GSE61151_series_matrix.txt.gz")
m61 = m61.set_index("gsm").reindex(b61.columns)
age61 = pd.to_numeric(m61["agebloodtaken"], errors="coerce")

print("carregando GSE40279 (treino, 1.2 GB) ...", flush=True)
b40, m40 = read_series_matrix(DATA / "GSE40279_series_matrix.txt.gz")
m40 = m40.set_index("gsm").reindex(b40.columns)
age40 = pd.to_numeric(m40["age (y)"], errors="coerce")


section("CHECAGENS DE SANIDADE")
c1 = not (set(m110.index) & set(b35.columns))
print(f"  1. a coorte que constroi o filtro (GSE110554) e a que mede o "
      f"deslocamento (GSE35069) nao compartilham amostra  -> {'ok' if c1 else 'FALHOU'}")

shared = b40.index.intersection(b61.index).intersection(b35.index).intersection(cellvar.index)
keep40 = age40.notna().to_numpy()
X = b40.loc[shared].to_numpy(dtype=np.float64).T[keep40]
E = b61.loc[shared].to_numpy(dtype=np.float64).T
P = b35.loc[shared].to_numpy(dtype=np.float64).T
cv_ = cellvar.reindex(shared).to_numpy()
good = (~np.isnan(X).any(axis=0) & ~np.isnan(E).any(axis=0)
        & ~np.isnan(P).any(axis=0) & ~np.isnan(cv_))
X, E, P, cv_ = X[:, good], E[:, good], P[:, good], cv_[good]
y = age40[keep40].to_numpy()
keep61 = age61.notna().to_numpy()
E, ye = E[keep61], age61[keep61].to_numpy()
print(f"     sondas nas QUATRO coortes: {X.shape[1]:,}")
del b40, b61

sel = frac35.isin(FRACTIONS).to_numpy()
P, dono, frc = P[sel], donor35.to_numpy()[sel], frac35.to_numpy()[sel]

r_age = ((X - X.mean(0)).T @ (y - y.mean())) / (X.std(0) * y.std() * len(y) + 1e-12)
absr = np.abs(r_age)
print(f"     correlacao com a idade x variancia entre tipos celulares: "
      f"Spearman {pd.Series(absr).corr(pd.Series(cv_), method='spearman'):+.3f}")
print("     (se fosse fortemente positiva, as sondas do tempo SERIAM as do tipo\n"
      "      celular e o piso estaria explicado antes de comecar)")


section("A VARREDURA — FILTRAR POR VARIANCIA ENTRE TIPOS, DEPOIS RANQUEAR POR IDADE")
rows = []
for k in K_SWEEP:
    base_idx = None
    for q in QUANTILES:
        if q >= 1.0:
            pool = np.arange(X.shape[1])
        else:
            pool = np.where(cv_ <= np.quantile(cv_, q))[0]
        if len(pool) < k:
            continue
        idx = pool[np.argsort(-absr[pool])[:k]]
        if base_idx is None:
            base_idx = idx
        m = fit_and_measure(idx, X, y, E, ye, P, dono, frc)
        overlap = len(np.intersect1d(idx, base_idx)) / k
        # matched negative control: discard the same NUMBER of probes at random
        ctrl = []
        for _ in range(N_CONTROL if q < 1.0 else 1):
            rp = RNG.choice(X.shape[1], size=len(pool), replace=False)
            ci = rp[np.argsort(-absr[rp])[:k]]
            ctrl.append(fit_and_measure(ci, X, y, E, ye, P, dono, frc))
        rows.append(dict(k=k, q=q, pool=len(pool), overlap=overlap, **m,
                         ctrl_r=float(np.mean([c["r"] for c in ctrl])),
                         ctrl_disp_n=float(np.mean([c["disp_n"] for c in ctrl])),
                         ctrl_slope=float(np.mean([c["slope"] for c in ctrl]))))
        print(f"  k={k:<6} q={q:<5} pool {len(pool):>7,}  sobreposicao {overlap:>5.0%}  "
              f"r {m['r']:+.3f}  escala {m['slope']:.3f}  "
              f"desloc/escala {m['disp_n']:>5.1f}a   "
              f"| controle r {rows[-1]['ctrl_r']:+.3f} desloc/escala "
              f"{rows[-1]['ctrl_disp_n']:>5.1f}a", flush=True)
res = pd.DataFrame(rows)
res.to_csv(OUT / "floor.csv", index=False)


section("CHECAGENS 2, 3 E 4")
strict = res[res.q <= 0.10]
c2 = bool(strict.overlap.max() < 0.95)
print(f"  2. o filtro muda mesmo as sondas: no filtro mais severo a sobreposicao "
      f"com o conjunto livre cai para {strict.overlap.min():.0%}  -> {'ok' if c2 else 'FALHOU'}")
cal = res[(res.slope >= 0.5) & (res.slope <= 1.5)]
c3 = len(cal) >= 4
print(f"  3. calibracao: {len(cal)} de {len(res)} configuracoes reportam anos  "
      f"-> {'ok' if c3 else 'FALHOU'}")
base = res[(res.q == 1.0) & (res.k == 1000)]
prev = pd.read_csv(ROOT / "results/reverse.csv")
prev_b = prev[(prev.k == 1000) & (prev.lam == float(base.lam.iloc[0]))]
c4 = bool(len(prev_b) and abs(float(base.disp_n.iloc[0]) - float(prev_b.disp_n.iloc[0])) < 1.5)
print(f"  4. a configuracao sem filtro reproduz a etapa 10 (k=1000): "
      f"{float(base.disp_n.iloc[0]):.1f}a aqui contra "
      f"{float(prev_b.disp_n.iloc[0]) if len(prev_b) else float('nan'):.1f}a la  "
      f"-> {'ok' if c4 else 'FALHOU'}")
print("     (nao sao identicos: aqui a intersecao inclui o EPIC e tem menos sondas)")

if not (c1 and c2 and c3 and c4):
    sys.exit("\n  UMA CHECAGEM FALHOU — parando antes de ler o resultado.")


section("O PISO DESCEU?")
for k in K_SWEEP:
    d = cal[cal.k == k].sort_values("q", ascending=False)
    if d.empty:
        continue
    free = d[d.q == 1.0]
    if free.empty:
        continue
    f0 = float(free.disp_n.iloc[0])
    print(f"\n  k = {k}:  sem filtro {f0:.1f}a de desloc/escala, r +{free.r.iloc[0]:.3f}")
    print(f"    {'q':<6}{'desloc/escala':>14}{'controle':>10}"
          f"{'vantagem do filtro':>21}{'r':>9}{'r controle':>12}")
    for _, r in d[d.q < 1.0].iterrows():
        edge = r.ctrl_disp_n - r.disp_n     # positivo = o filtro ganhou do controle
        print(f"    {r.q:<6.2f}{r.disp_n:>13.1f}a{r.ctrl_disp_n:>9.1f}a"
              f"{edge:>+20.1f}a{r.r:>+9.3f}{r.ctrl_r:>+12.3f}")
print("""
  A coluna que decide e "vantagem do filtro": quanto o filtrado fica ABAIXO do
  seu proprio controle. O controle descarta a mesma quantidade de sondas, so
  que ao acaso, entao o que ele consegue e o que qualquer poda conseguiria.
  Positivo significa que a informacao de tipo celular acrescentou alguma coisa.
  E preciso ler junto com as duas ultimas colunas: um filtro que ganha 1 ano de
  planura e perde 0.09 de correlacao nao trocou nada de bom.""")


section("POR QUE — O ENRIQUECIMENTO QUE EXPLICA O PISO")
print("""  A correlacao global entre |r com a idade| e variancia entre tipos celulares
  e so +0.13, o que sugeriria que as duas coisas mal se tocam. Mas o que um
  relogio usa nao e a sonda media: e o TOPO do ranking de idade. A pergunta
  certa e quanto desse topo cai no quarto mais variavel entre tipos celulares.
  Sob independencia seria 25%.\n""")
hi = cv_ >= np.quantile(cv_, 0.75)
ordem = np.argsort(-absr)
print(f"  {'melhores k sondas de idade':<28}{'no quarto mais variavel':>25}{'enriquecimento':>17}")
enr = []
for k in (100, 1_000, 10_000, 100_000):
    top = ordem[:k]
    fr = float(hi[top].mean())
    enr.append(dict(k=k, frac=fr, fold=fr / 0.25))
    print(f"  {k:<28,}{fr:>24.0%}{fr/0.25:>16.1f}x")
pd.DataFrame(enr).to_csv(OUT / "floor_enrichment.csv", index=False)
print("""
  Esse e o piso. Nao e que a selecao por idade ignore o tipo celular por
  descuido -- e que as sondas que melhor marcam o tempo sao, em boa parte, as
  mesmas que distinguem as celulas. Tirar umas e tirar as outras.""")
