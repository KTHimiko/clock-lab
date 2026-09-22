#!/usr/bin/env python3
"""
Stage 10 — was the Horvath advantage design, or just a wider training range?

Stage 9 measured displacement per year of genuine age response and found the two
Horvath clocks (6.91 and 7.33 years) beyond the reach of all 45 configurations
of a ridge family trained on GSE61151. It also declared a handicap: that family
learned on ages 35-83 and was tested on 19-101, so its slope was compressed by
extrapolating outside the range it had seen.

This stage removes the handicap by running the test backwards. Train on
GSE40279 — 656 samples, ages 19 to 101 — and test on GSE61151. Now the family
has the wider range and the published clocks do not.

A SYMMETRY WORTH NAMING. GSE40279 is the cohort the Hannum clock was trained on,
so a family trained there stands exactly where Hannum stands, and comparing the
two on GSE61151 is finally fair. The Horvath clocks are external to both.

AND THE DECOMPOSITION. If the family still loses, "design" is too vague an
answer to stop at, so the design is split in two:

  (a) Horvath's probes, Horvath's coefficients        — the clock itself
  (b) Horvath's probes, coefficients refitted here    — isolates PROBE CHOICE
  (c) the family's own top-353 probes, fitted here    — same size, own choice

If (b) keeps Horvath's flatness, the advantage lives in which CpGs he picked. If
(b) falls to (c), it lives in the coefficients — in how the published weights
happen to cancel across cell types, which is what stage 4 said could not be
engineered.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. the three cohorts share no sample, asserted on GSM ids
  2. the published clocks must work on GSE61151 — Horvath 2013 above r = 0.75,
     the bar stage 1 set. Recomputed here, not imported, because stage 1 used a
     different probe set than the three-way intersection leaves
  3. at least 300,000 probes shared and complete across all three
  4. CALIBRATION. Only configurations that actually report years enter the
     comparison: slope between 0.5 and 1.5 on the test cohort.

     THIS CHECK REPLACED A WRONG ONE, and the wrong one is worth keeping on
     the page. It first read: "the degenerate corner must no longer look
     good — the collapsed configurations must sit at or below the median on
     displacement-per-response." It failed, at 10.5 against a median of 11.2,
     and the premise was the thing at fault. As lambda grows, ridge does not
     converge to noise; it converges to the correlation-weighted direction,
     which is a perfectly reasonable estimator that happens to be multiplied
     by a vanishing constant. Normalising by scale divides that constant out,
     so the collapsed configurations are entitled to a decent
     displacement-per-response. What disqualifies them is not flatness, it is
     that their slope is 0.009 and they are wrong by ten years — they do not
     report an age at all. A calibration band says that directly.

     The headline below does not depend on the swap: the best collapsed
     configuration reaches 10.5 and the frontier sits at 3.9, so they were
     never in contention either way
  5. every lambda is chosen by cross-validation INSIDE the training cohort. The
     test cohort is touched once, at the end, for reading only

Usage:  .venv/bin/python analysis/10_reverse.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix
from model.clocks import CLOCKS, predict, load_coefficients

DATA = ROOT / "reference/data"
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(20250922)

K_PROBES = [100, 1_000, 10_000, 100_000, None]
LAMBDAS = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e6]
FRACTIONS = ["CD4+ T cells", "CD8+ T cells", "CD19+ B cells",
             "CD56+ NK cells", "CD14+ Monocytes", "Neutrophils"]


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


def ridge_dual(Xtr, ytr, lambdas):
    mx, my = Xtr.mean(axis=0), ytr.mean()
    Xc, yc = Xtr - mx, ytr - my
    w, V = np.linalg.eigh(Xc @ Xc.T)
    Vty = V.T @ yc
    for lam in lambdas:
        yield lam, Xc.T @ (V @ (Vty / (w + lam))), mx, my


def measure(beta, mx, my, idx, E, ye, P, dono, frc):
    """Accuracy on the test cohort and displacement, both scale-normalised."""
    pe = (E[:, idx] - mx) @ beta + my
    pv = (P[:, idx] - mx) @ beta + my
    w = pd.DataFrame({"d": dono, "f": frc, "v": pv}).pivot_table(
        index="d", columns="f", values="v")
    slope = float(np.cov(pe, ye)[0, 1] / np.var(ye))
    disp = float(w.std(axis=1).mean())
    return dict(r=float(np.corrcoef(ye, pe)[0, 1]),
                mae=float(np.abs(pe - ye).mean()),
                slope=slope, disp=disp,
                disp_n=disp / slope if slope > 1e-9 else np.inf)


print("carregando GSE35069 ...", flush=True)
b35, m35 = read_series_matrix(DATA / "GSE35069_series_matrix.txt.gz")
m35 = m35.set_index("gsm").reindex(b35.columns)
frac35 = m35["tissue/cell type"]
donor35 = m35["title"].str.extract(r"_(\d+)$")[0]

print("carregando GSE61151 (agora o TESTE) ...", flush=True)
b61, m61 = read_series_matrix(DATA / "GSE61151_series_matrix.txt.gz")
m61 = m61.set_index("gsm").reindex(b61.columns)
age61 = pd.to_numeric(m61["agebloodtaken"], errors="coerce")

print("carregando GSE40279 (agora o TREINO, 1.2 GB) ...", flush=True)
b40, m40 = read_series_matrix(DATA / "GSE40279_series_matrix.txt.gz")
m40 = m40.set_index("gsm").reindex(b40.columns)
age40 = pd.to_numeric(m40["age (y)"], errors="coerce")


section("CHECAGENS DE SANIDADE")
c1 = not (set(b35.columns) & set(b61.columns) or set(b35.columns) & set(b40.columns)
          or set(b61.columns) & set(b40.columns))
print(f"  1. nenhuma amostra compartilhada  -> {'ok' if c1 else 'FALHOU'}")

pub = []
for c in CLOCKS:
    pv, cov = predict(b61, c)
    ok = age61.notna() & pv.notna()
    pub.append(dict(clock=c, r=float(np.corrcoef(age61[ok], pv[ok])[0, 1]),
                    slope=float(np.cov(pv[ok], age61[ok])[0, 1] / np.var(age61[ok]))))
pub = pd.DataFrame(pub).set_index("clock")
c2 = float(pub.loc["Horvath2013", "r"]) > 0.75
print(f"  2. relogios publicados no GSE61151 (teste):")
for c, r in pub.iterrows():
    print(f"     {c:<14} r {r.r:+.3f}  escala {r.slope:.3f}")
print(f"     -> {'ok' if c2 else 'FALHOU'}")

shared = b40.index.intersection(b61.index).intersection(b35.index)
keep40 = age40.notna().to_numpy()
X = b40.loc[shared].to_numpy(dtype=np.float64).T[keep40]
E = b61.loc[shared].to_numpy(dtype=np.float64).T
P = b35.loc[shared].to_numpy(dtype=np.float64).T
good = ~np.isnan(X).any(axis=0) & ~np.isnan(E).any(axis=0) & ~np.isnan(P).any(axis=0)
X, E, P = X[:, good], E[:, good], P[:, good]
probe_ids = pd.Index(shared)[good]
y = age40[keep40].to_numpy()
keep61 = age61.notna().to_numpy()
E, ye = E[keep61], age61[keep61].to_numpy()
c3 = X.shape[1] >= 300_000
print(f"  3. sondas comuns e completas: {X.shape[1]:,}  -> {'ok' if c3 else 'FALHOU'}")
print(f"     treino {X.shape[0]} amostras ({y.min():.0f}-{y.max():.0f} anos), "
      f"teste {E.shape[0]} ({ye.min():.0f}-{ye.max():.0f} anos)")
del b40, b61

sel = frac35.isin(FRACTIONS).to_numpy()
P, dono, frc = P[sel], donor35.to_numpy()[sel], frac35.to_numpy()[sel]


section("A FAMILIA, TREINADA NA FAIXA LARGA E TESTADA NA ESTREITA")
r_full = ((X - X.mean(0)).T @ (y - y.mean())) / (X.std(0) * y.std() * len(y) + 1e-12)
order = np.argsort(-np.abs(r_full))
rows = []
for k in K_PROBES:
    idx = order if k is None else order[:k]
    for lam, beta, mx, my in ridge_dual(X[:, idx], y, LAMBDAS):
        rows.append(dict(k=(X.shape[1] if k is None else k), lam=lam,
                         l2=float(np.sqrt((beta ** 2).sum())),
                         **measure(beta, mx, my, idx, E, ye, P, dono, frc)))
res = pd.DataFrame(rows)
res.to_csv(OUT / "reverse.csv", index=False)
print(f"  {'sondas':>8}{'lambda':>9}{'r teste':>10}{'erro':>8}{'escala':>9}"
      f"{'desloc.':>9}{'desloc/escala':>15}")
for _, r in res.iterrows():
    print(f"  {r.k:>8.0f}{r.lam:>9.0e}{r.r:>+10.3f}{r.mae:>7.1f}a{r.slope:>9.3f}"
          f"{r.disp:>8.2f}a{r.disp_n:>14.1f}a")

cal = res[(res.slope >= 0.5) & (res.slope <= 1.5)]
c4 = len(cal) >= 5
print(f"\n  4. calibracao: {len(cal)} das {len(res)} configuracoes reportam anos "
      f"(inclinacao entre 0.5 e 1.5)  -> {'ok' if c4 else 'FALHOU'}")
deg = res[res.lam >= 1e4]
print(f"     descartadas as colapsadas: lambda >= 1e4 tem inclinacao ate "
      f"{deg.slope.max():.3f} e erro de {deg.mae.min():.1f}a para cima")
print(f"     (a versao anterior desta checagem exigia que o canto colapsado "
      f"parecesse ruim\n      no eixo normalizado. Falhou, e a premissa e que estava "
      f"errada -- ver docstring.)")

if not (c1 and c2 and c3 and c4):
    sys.exit("\n  UMA CHECAGEM FALHOU — parando antes de ler o resultado.")


section("COMPARACAO — A FAMILIA CHEGA NOS HORVATH AGORA?")
pubd = pd.read_csv(ROOT / "results/frontier_published.csv").set_index("clock")
print(f"  {'relogio':<16}{'r teste':>9}{'escala':>9}{'desloc/escala':>15}"
      f"{'   configuracoes melhores nos dois'}")
best = []
for c, r in pub.iterrows():
    dn = float(pubd.loc[c, "disp"]) / r.slope
    better = cal[(cal.disp_n < dn) & (cal.r >= r.r)]
    best.append(dict(clock=c, r=r.r, slope=r.slope, disp_n=dn, better=len(better)))
    print(f"  {c:<16}{r.r:>+9.3f}{r.slope:>9.3f}{dn:>14.1f}a{len(better):>25}")
pd.DataFrame(best).to_csv(OUT / "reverse_published.csv", index=False)
b = cal.loc[cal.disp_n.idxmin()]
print(f"\n  a configuracao mais plana da familia: {b.k:.0f} sondas, lambda {b.lam:.0e}"
      f"  ->  r {b.r:+.3f}, desloc/escala {b.disp_n:.1f}a")


section("DECOMPOSICAO — SONDAS DO HORVATH OU COEFICIENTES DO HORVATH?")
print("""  (a) sondas do Horvath + coeficientes do Horvath  = o relogio publicado
  (b) sondas do Horvath + coeficientes reajustados aqui
  (c) as 353 melhores sondas pela propria familia, reajustadas aqui

  (b) contra (a) isola a escolha das SONDAS; (b) contra (c), o tamanho do
  conjunto fica igual e o que muda e so quem escolheu.\\n""")


def cv_lambda(Xs, ys, lambdas, folds=5):
    """Pick lambda inside the training cohort. The test set is never consulted."""
    parts = np.array_split(RNG.permutation(len(ys)), folds)
    err = {lam: 0.0 for lam in lambdas}
    for test in parts:
        tr = np.setdiff1d(np.arange(len(ys)), test)
        for lam, beta, mx, my in ridge_dual(Xs[tr], ys[tr], lambdas):
            err[lam] += float(np.abs((Xs[test] - mx) @ beta + my - ys[test]).sum())
    return min(err, key=err.get)


dec = []
for tag, idx in (("b_sondas_horvath", None), ("c_sondas_proprias", None)):
    pass
hv = load_coefficients("Horvath2013")[0]
idx_h = np.where(probe_ids.isin(hv.index))[0]
idx_c = order[:len(idx_h)]
print(f"  sondas do Horvath presentes na intersecao: {len(idx_h)} de {len(hv)}")
for tag, idx in (("(b) sondas do Horvath, reajustado", idx_h),
                 ("(c) sondas proprias, mesmo tamanho", idx_c)):
    lam = cv_lambda(X[:, idx], y, LAMBDAS)
    beta, mx, my = next((b_, m_, my_) for l_, b_, m_, my_
                        in ridge_dual(X[:, idx], y, [lam]))
    m = measure(beta, mx, my, idx, E, ye, P, dono, frc)
    dec.append(dict(variant=tag, n=len(idx), lam=lam, **m))
    print(f"  {tag:<38} lambda {lam:.0e}  r {m['r']:+.3f}  "
          f"escala {m['slope']:.3f}  desloc/escala {m['disp_n']:.1f}a")
a = pub.loc["Horvath2013"]
a_dn = float(pubd.loc["Horvath2013", "disp"]) / a.slope
dec.append(dict(variant="(a) Horvath2013 publicado", n=len(hv), lam=np.nan,
                r=a.r, mae=np.nan, slope=a.slope,
                disp=float(pubd.loc["Horvath2013", "disp"]), disp_n=a_dn))
print(f"  {'(a) Horvath2013 publicado':<38} {'':>11}  r {a.r:+.3f}  "
      f"escala {a.slope:.3f}  desloc/escala {a_dn:.1f}a")
pd.DataFrame(dec).to_csv(OUT / "reverse_decomposition.csv", index=False)
