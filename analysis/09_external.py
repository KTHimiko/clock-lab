#!/usr/bin/env python3
"""
Stage 9 — does dilution still work away from home?

Stage 8 found that cell-type displacement is governed by the L2 norm of a
clock's coefficients, so spreading the same total weight over more probes buys
flatness. It also carried a caveat that swallowed half the claim: accuracy there
was cross-validated INSIDE the cohort the family trained on. A diluted ridge
over 458,674 probes is exactly the kind of model that can memorise a cohort's
batch structure and lose everything when it travels.

GSE40279 settles it. 656 whole blood samples, ages 19 to 101, a different study
on a different population. A clock trained on GSE61151 has never seen one of
these people. Accuracy measured here is the real thing.

THE TRAP IN THIS DATASET, DECLARED BEFORE ANY NUMBER. **GSE40279 is the cohort
the Hannum clock was trained on.** Its accuracy here is a memory, not a
prediction. It is still reported — because leaving it out would hide the fact —
but it is marked, and it is never used as the bar the family has to clear.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. the three cohorts share no sample, asserted on GSM ids
  2. age tripwire on the new cohort: nothing below 18 or above 110. The loader
     bug that stage 5 fixed produced exactly this signature, and a new dataset
     is where it would show up again
  3. the published clocks must work here, or the file was parsed wrong rather
     than the clocks being wrong. Horvath 2013 must reach r > 0.75 — it has no
     stake in this cohort, so it is the honest probe
  4. at least 300,000 probes shared and complete across all three cohorts,
     otherwise the family cannot be rebuilt as stage 8 built it

THE TRAP THIS STAGE WALKED INTO, AND HOW BOTH AXES CHANGED BECAUSE OF IT.
Correlation is scale-invariant. At the degenerate end of a ridge sweep the
coefficients collapse towards zero, the predictions become the training mean
plus an epsilon multiple of a real age direction — and r stays high while every
prediction is wrong by thirteen years. Displacement, measured in years, collapses
by the same epsilon. So the degenerate corner reads as "accurate AND perfectly
flat" and dominates any frontier drawn on those two columns. It is not a clock;
it is a constant with a rumour of a direction.

Both axes are therefore normalised by the predictor's SCALE: the slope of its
prediction on real age in the external cohort. Displacement becomes years of
cell-type spread per year of genuine age response, and epsilon cancels out of
both. The raw columns are printed too, because the trap is worth seeing.

Usage:  .venv/bin/python analysis/09_external.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix
from model.clocks import CLOCKS, predict

DATA = ROOT / "reference/data"
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)

K_PROBES = [100, 1_000, 10_000, 100_000, None]
LAMBDAS = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e6]
FRACTIONS = ["CD4+ T cells", "CD8+ T cells", "CD19+ B cells",
             "CD56+ NK cells", "CD14+ Monocytes", "Neutrophils"]
HANNUM_HOME = "Hannum2013"      # trained on GSE40279


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


def ridge_dual(Xtr, ytr, lambdas):
    mx, my = Xtr.mean(axis=0), ytr.mean()
    Xc, yc = Xtr - mx, ytr - my
    w, V = np.linalg.eigh(Xc @ Xc.T)
    Vty = V.T @ yc
    for lam in lambdas:
        yield lam, Xc.T @ (V @ (Vty / (w + lam))), mx, my


print("carregando GSE35069 (medida de deslocamento) ...", flush=True)
b35, m35 = read_series_matrix(DATA / "GSE35069_series_matrix.txt.gz")
m35 = m35.set_index("gsm").reindex(b35.columns)
frac35 = m35["tissue/cell type"]
donor35 = m35["title"].str.extract(r"_(\d+)$")[0]

print("carregando GSE61151 (treino) ...", flush=True)
b61, m61 = read_series_matrix(DATA / "GSE61151_series_matrix.txt.gz")
m61 = m61.set_index("gsm").reindex(b61.columns)
age61 = pd.to_numeric(m61["agebloodtaken"], errors="coerce")

print("carregando GSE40279 (teste externo, 1.2 GB) ...", flush=True)
b40, m40 = read_series_matrix(DATA / "GSE40279_series_matrix.txt.gz")
m40 = m40.set_index("gsm").reindex(b40.columns)
age40 = pd.to_numeric(m40["age (y)"], errors="coerce")
print(f"  {b40.shape[0]:,} sondas x {b40.shape[1]} amostras", flush=True)


section("CHECAGENS DE SANIDADE")
ids = [set(b35.columns), set(b61.columns), set(b40.columns)]
c1 = not (ids[0] & ids[1] or ids[0] & ids[2] or ids[1] & ids[2])
print(f"  1. nenhuma amostra compartilhada entre as tres coortes  "
      f"-> {'ok' if c1 else 'FALHOU'}")

bad = ((age40 < 18) | (age40 > 110)).sum()
c2 = (bad == 0) and age40.notna().all()
print(f"  2. idades do GSE40279: {age40.min():.0f} a {age40.max():.0f}, "
      f"{int(age40.isna().sum())} faltando, {int(bad)} fora de 18-110  "
      f"-> {'ok' if c2 else 'FALHOU'}")

pub, predict_cache = [], {}
for c in CLOCKS:
    try:
        pv, cov = predict(b40, c)
        predict_cache[c] = pv
    except ValueError as e:
        print(f"     {c:<14} RECUSADO: {e}"); continue
    ok = age40.notna() & pv.notna()
    pub.append(dict(clock=c, coverage=cov,
                    r=float(np.corrcoef(age40[ok], pv[ok])[0, 1]),
                    mae=float(np.abs(pv[ok] - age40[ok]).mean())))
pub = pd.DataFrame(pub)
c3 = float(pub.set_index("clock").loc["Horvath2013", "r"]) > 0.75
print(f"  3. relogios publicados no GSE40279 (Horvath2013 e a sonda honesta):")
for _, r in pub.iterrows():
    tag = "  <- TREINOU AQUI, isto e memoria" if r.clock == HANNUM_HOME else ""
    print(f"     {r.clock:<14} r {r.r:+.3f}  erro {r.mae:>4.1f}a  "
          f"cobertura {r.coverage:.0%}{tag}")
print(f"     -> {'ok' if c3 else 'FALHOU'}")

shared = b61.index.intersection(b35.index).intersection(b40.index)
keep61 = age61.notna().to_numpy()
X = b61.loc[shared].to_numpy(dtype=np.float64).T[keep61]
P = b35.loc[shared].to_numpy(dtype=np.float64).T
E = b40.loc[shared].to_numpy(dtype=np.float64).T
good = ~np.isnan(X).any(axis=0) & ~np.isnan(P).any(axis=0) & ~np.isnan(E).any(axis=0)
X, P, E = X[:, good], P[:, good], E[:, good]
y = age61[keep61].to_numpy()
ye = age40.to_numpy()
c4 = X.shape[1] >= 300_000
print(f"  4. sondas comuns e completas nas tres: {X.shape[1]:,}  "
      f"-> {'ok' if c4 else 'FALHOU'}")
del b61, b40

if not (c1 and c2 and c3 and c4):
    sys.exit("\n  UMA CHECAGEM FALHOU — parando antes de ler o resultado.")

sel = frac35.isin(FRACTIONS).to_numpy()
P, dono, frc = P[sel], donor35.to_numpy()[sel], frac35.to_numpy()[sel]


section("A FAMILIA, TREINADA NO GSE61151 E TESTADA EM 656 ESTRANHOS")
r_full = ((X - X.mean(0)).T @ (y - y.mean())) / (X.std(0) * y.std() * len(y) + 1e-12)
order = np.argsort(-np.abs(r_full))
rows = []
for k in K_PROBES:
    idx = order if k is None else order[:k]
    for lam, beta, mx, my in ridge_dual(X[:, idx], y, LAMBDAS):
        pe = (E[:, idx] - mx) @ beta + my
        pv = (P[:, idx] - mx) @ beta + my
        w = pd.DataFrame({"d": dono, "f": frc, "v": pv}).pivot_table(
            index="d", columns="f", values="v")
        slope = float(np.cov(pe, ye)[0, 1] / np.var(ye))
        disp = float(w.std(axis=1).mean())
        rows.append(dict(k=(X.shape[1] if k is None else k), lam=lam,
                         r_ext=float(np.corrcoef(ye, pe)[0, 1]),
                         mae_ext=float(np.abs(pe - ye).mean()),
                         l2=float(np.sqrt((beta ** 2).sum())),
                         slope=slope, disp=disp,
                         disp_n=disp / slope if slope > 1e-9 else np.inf))
res = pd.DataFrame(rows)

# In-cohort accuracy is recomputed HERE rather than read from stage 8, because
# stage 8 intersected two cohorts and this stage intersects three: its full-probe
# model has 458,674 probes and this one has 447,564. They are different models,
# and putting one's cross-validated r beside the other's external r would be the
# cross-stage comparison this project keeps getting wrong.
print("\n  validacao cruzada refeita sobre AS MESMAS sondas ...", flush=True)
RNG = np.random.default_rng(20250922)
n = len(y)
folds = np.array_split(RNG.permutation(n), 5)
cvpred = {(k, lam): np.empty(n) for k in K_PROBES for lam in LAMBDAS}
for test in folds:
    tr = np.setdiff1d(np.arange(n), test)
    Xtr, ytr = X[tr], y[tr]
    rr = ((Xtr - Xtr.mean(0)).T @ (ytr - ytr.mean())) / (
        Xtr.std(0) * ytr.std() * len(tr) + 1e-12)
    o = np.argsort(-np.abs(rr))
    for k in K_PROBES:
        i = o if k is None else o[:k]
        for lam, beta, mx, my in ridge_dual(Xtr[:, i], ytr, LAMBDAS):
            cvpred[(k, lam)][test] = (X[np.ix_(test, i)] - mx) @ beta + my
res["r_cv"] = [float(np.corrcoef(y, cvpred[(None if r.k == X.shape[1] else int(r.k),
                                            r.lam)])[0, 1]) for _, r in res.iterrows()]
res["queda"] = res.r_cv - res.r_ext
res.to_csv(OUT / "external.csv", index=False)

print(f"\n  {'sondas':>8}{'lambda':>9}{'r casa':>9}{'r fora':>9}{'erro':>8}"
      f"{'escala':>9}{'desloc.':>9}{'desloc/escala':>15}")
for _, r in res.iterrows():
    print(f"  {r.k:>8.0f}{r.lam:>9.0e}{r.r_cv:>+9.3f}{r.r_ext:>+9.3f}"
          f"{r.mae_ext:>7.1f}a{r.slope:>9.3f}{r.disp:>8.2f}a{r.disp_n:>14.1f}a")
print("""
  Leia a coluna 'escala' junto com 'r fora'. Onde a escala cai para perto de
  zero, o relogio esta prevendo uma constante: o r continua alto porque
  correlacao nao ve escala, mas o erro vai para 13 anos e o deslocamento vai a
  zero pelo mesmo motivo. A ultima coluna e a que compara: anos de deslocamento
  por ano de resposta real a idade.""")


section("A DILUICAO SOBREVIVE A VIAGEM?")
q = res[(res.r_cv > 0.5) & (res.disp > 1e-3)].copy()
rho = float(np.log10(q.l2).corr(q.queda, method="spearman"))
print(f"  sobre as {len(q)} configuracoes que eram relogios em casa (r_cv > 0.5):")
print(f"    Spearman(log norma L2, queda de r ao sair de casa) = {rho:+.3f}")
print("""
  Positivo significa que o relogio concentrado perde MAIS ao viajar, e a
  diluicao defende duas vezes: contra composicao e contra a coorte. Negativo
  significaria que a diluicao so funcionava em casa, e a etapa 8 teria vendido
  uma defesa que nao existe.\\n""")
conc = q.nlargest(5, "l2")[["k", "lam", "r_cv", "r_ext", "queda", "disp"]]
dilu = q.nsmallest(5, "l2")[["k", "lam", "r_cv", "r_ext", "queda", "disp"]]
for tag, d in (("mais concentrados", conc), ("mais diluidos", dilu)):
    print(f"  os 5 {tag}:")
    for _, r in d.iterrows():
        print(f"    {r.k:>7.0f} sondas, lambda {r.lam:.0e}: casa {r.r_cv:+.3f} "
              f"-> fora {r.r_ext:+.3f}  (queda {r.queda:+.3f}), desloc {r.disp:.2f}a")


section("A FRONTEIRA, COM OS DOIS EIXOS INVARIANTES A ESCALA")
front = []
for _, r in res.sort_values("r_ext", ascending=False).iterrows():
    if not front or r.disp_n < min(f["disp_n"] for f in front):
        front.append(dict(k=r.k, lam=r.lam, r_ext=r.r_ext, mae_ext=r.mae_ext,
                          slope=r.slope, disp=r.disp, disp_n=r.disp_n))
pd.DataFrame(front).to_csv(OUT / "external_pareto.csv", index=False)
print(f"  {'sondas':>8}{'lambda':>9}{'r fora':>9}{'erro':>8}{'escala':>9}"
      f"{'desloc/escala':>15}")
for f in front:
    print(f"  {f['k']:>8.0f}{f['lam']:>9.0e}{f['r_ext']:>+9.3f}"
          f"{f['mae_ext']:>7.1f}a{f['slope']:>9.3f}{f['disp_n']:>14.1f}a")

pubd = pd.read_csv(ROOT / "results/frontier_published.csv").set_index("clock")
print("\n  os publicados, na mesma medida (escala tirada do GSE40279,")
print("  deslocamento do GSE35069):")
for _, r in pub.iterrows():
    pv, _ = predict(b35, r.clock)
    pe, _ = predict_cache[r.clock], None
    sl = float(np.cov(pe[age40.notna()].to_numpy(),
                      age40.dropna().to_numpy())[0, 1] / np.var(age40.dropna()))
    dn = float(pubd.loc[r.clock, "disp"]) / sl
    better = res[(res.disp_n < dn) & (res.r_ext >= r.r)]
    tag = "  [TREINOU NESTA COORTE]" if r.clock == HANNUM_HOME else ""
    print(f"    {r.clock:<14} r {r.r:+.3f}  escala {sl:.3f}  "
          f"desloc/escala {dn:>5.2f}a   {len(better)} configuracoes sao mais "
          f"planas E pelo menos tao precisas{tag}")
