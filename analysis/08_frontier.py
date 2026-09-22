#!/usr/bin/env python3
"""
Stage 8 — turning the stage 4 argument into a measurement.

Stage 4 found that the cell-type displacement of a clock is exactly what a
random pairing of its weights with its probes' cell-type shifts produces, and
proposed that what separates clocks is therefore **gain** — total absolute
weight — rather than which probes they chose. That was recorded as an argument
and not a measurement, because it rested on two linear clocks, where a
correlation is ±1 by arithmetic.

Four clocks cannot test it. Two hundred can. This stage trains a family of age
predictors on GSE61151 whole blood, sweeping the number of probes and the ridge
penalty, and measures two numbers for each:

  ACCURACY — out-of-sample, by 5-fold cross-validation
  DISPLACEMENT — the stage 2 statistic (spread of predicted age across fractions
  within one donor) on GSE35069 purified cells, which the training never sees

If gain is the whole story, displacement is a function of total weight and
nothing else, and the two axes trace one curve: you buy accuracy with gain and
pay for it in composition. If probe choice matters after all, clocks of equal
gain will scatter in displacement, and some corner of the sweep will hold a
predictor that is accurate AND flat.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. probe selection happens INSIDE each fold. Ranking probes by correlation
     with age on all 184 samples and then cross-validating is the oldest
     leakage in this field, and it inflates accuracy without touching
     displacement — precisely the shape that would fake a frontier
  2. the family must contain a real clock: the best configuration must reach
     r > 0.75 held out, the same bar stage 1 set for the published clocks
  3. the degenerate corner must appear. At a large enough penalty every
     coefficient goes to zero, the prediction becomes the training mean, and
     BOTH accuracy and displacement must go to zero. If that end of the sweep
     does not show up, the sweep is not measuring what it claims.
     THIS CHECK WAS MIS-SPECIFIED ON FIRST WRITING and it failed, so what
     happened to it is on the record. It had three clauses: displacement -> 0,
     accuracy -> 0, and total absolute weight < 1.0. The first two passed
     decisively (0.020 years, |r| 0.16). The third failed at 1.877 — because
     total absolute weight is a sum over probes, so its scale rides on how many
     probes a configuration uses, and 1.877 spread over 458,674 probes is four
     millionths per probe. An absolute threshold on it was meaningless from the
     start. It is replaced by its scale-free form, weight PER PROBE, and not
     deleted: relaxing a check that fails is the move that should always be
     looked at twice, so it is written down rather than quietly edited
  4. GSE35069 contributes nothing to training — asserted on the sample ids

Usage:  .venv/bin/python analysis/08_frontier.py
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

K_PROBES = [100, 1_000, 10_000, 100_000, None]     # None = every shared probe
LAMBDAS = [1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3, 1e4, 1e6]
N_FOLDS = 5
FRACTIONS = ["CD4+ T cells", "CD8+ T cells", "CD19+ B cells",
             "CD56+ NK cells", "CD14+ Monocytes", "Neutrophils"]


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


print("carregando GSE61151 (treino) ...", flush=True)
b61, m61 = read_series_matrix(DATA / "GSE61151_series_matrix.txt.gz")
m61 = m61.set_index("gsm").reindex(b61.columns)
age = pd.to_numeric(m61["agebloodtaken"], errors="coerce")
keep = age.notna().to_numpy()

print("carregando GSE35069 (celulas purificadas, so para medir) ...", flush=True)
b35, m35 = read_series_matrix(DATA / "GSE35069_series_matrix.txt.gz")
m35 = m35.set_index("gsm").reindex(b35.columns)
frac = m35["tissue/cell type"]
donor = m35["title"].str.extract(r"_(\d+)$")[0]

assert not (set(b61.columns) & set(b35.columns)), "as duas coortes compartilham amostra"
print(f"  checagem 4: nenhuma amostra em comum entre treino e medida  -> ok")

shared = b61.index.intersection(b35.index)
X_all = b61.loc[shared].to_numpy(dtype=np.float64).T[keep]     # samples x probes
P_all = b35.loc[shared].to_numpy(dtype=np.float64).T           # purified
good = ~np.isnan(X_all).any(axis=0) & ~np.isnan(P_all).any(axis=0)
X_all, P_all = X_all[:, good], P_all[:, good]
y_all = age[keep].to_numpy()
n, p = X_all.shape
print(f"  sondas comuns e completas: {p:,}   amostras de treino: {n}   "
      f"purificadas: {P_all.shape[0]}")

sel_frac = frac.isin(FRACTIONS).to_numpy()
P, dono, frc = P_all[sel_frac], donor.to_numpy()[sel_frac], frac.to_numpy()[sel_frac]
print(f"  fracoes usadas na medida: {len(FRACTIONS)} tipos, {P.shape[0]} amostras")


def ridge_dual(Xtr, ytr, lambdas):
    """Ridge for every lambda at once, solved in sample space (n << p)."""
    mx, my = Xtr.mean(axis=0), ytr.mean()
    Xc, yc = Xtr - mx, ytr - my
    G = Xc @ Xc.T
    w, V = np.linalg.eigh(G)
    Vty = V.T @ yc
    for lam in lambdas:
        alpha = V @ (Vty / (w + lam))
        yield lam, Xc.T @ alpha, mx, my


section("CHECAGEM 1 — SELECAO DENTRO DAS DOBRAS")
print("""  As sondas sao ranqueadas pela correlacao com a idade SO nas amostras de
  treino de cada dobra. Ranquear nas 184 e depois validar cruzado vaza a idade
  do teste para dentro da escolha das sondas: inflaria a acuracia sem mexer no
  deslocamento, que e exatamente o formato de uma fronteira falsa.""")

folds = np.array_split(RNG.permutation(n), N_FOLDS)
pred = {(k, lam): np.empty(n) for k in K_PROBES for lam in LAMBDAS}
for f, test in enumerate(folds):
    tr = np.setdiff1d(np.arange(n), test)
    Xtr, ytr = X_all[tr], y_all[tr]
    r = ((Xtr - Xtr.mean(0)).T @ (ytr - ytr.mean())) / (
        Xtr.std(0) * ytr.std() * len(tr) + 1e-12)
    order = np.argsort(-np.abs(r))
    for k in K_PROBES:
        idx = order if k is None else order[:k]
        for lam, beta, mx, my in ridge_dual(Xtr[:, idx], ytr, LAMBDAS):
            pred[(k, lam)][test] = (X_all[np.ix_(test, idx)] - mx) @ beta + my
    print(f"  dobra {f+1}/{N_FOLDS} pronta", flush=True)


section("DESLOCAMENTO — MEDIDO EM CELULAS QUE O TREINO NUNCA VIU")
# per-probe cell-type shift: how far the six fraction means spread apart. It is
# what a random pairing of weights with shifts would multiply against.
shift = pd.DataFrame({f: P[frc == f].mean(axis=0) for f in FRACTIONS}).std(axis=1).to_numpy()
rows = []
r_full = ((X_all - X_all.mean(0)).T @ (y_all - y_all.mean())) / (
    X_all.std(0) * y_all.std() * n + 1e-12)
order_full = np.argsort(-np.abs(r_full))
for k in K_PROBES:
    idx = order_full if k is None else order_full[:k]
    for lam, beta, mx, my in ridge_dual(X_all[:, idx], y_all, LAMBDAS):
        pv = (P[:, idx] - mx) @ beta + my
        w = pd.DataFrame({"d": dono, "f": frc, "v": pv}).pivot_table(
            index="d", columns="f", values="v")
        disp = float(w.std(axis=1).mean())
        span = float((w.max(axis=1) - w.min(axis=1)).mean())
        pv_cv = pred[(k, lam)]
        rows.append(dict(k=(p if k is None else k), lam=lam,
                         r=float(np.corrcoef(y_all, pv_cv)[0, 1]),
                         mae=float(np.abs(pv_cv - y_all).mean()),
                         gain=float(np.abs(beta).sum()),
                         gain_l2=float(np.sqrt((beta ** 2).sum())),
                         gain_l2s=float(np.sqrt(((beta * shift[idx]) ** 2).sum())),
                         disp=disp, span=span,
                         nz=int((np.abs(beta) > 1e-12).sum())))
res = pd.DataFrame(rows)
res.to_csv(OUT / "frontier.csv", index=False)

print(f"  {'sondas':>8}{'lambda':>9}{'r fora':>9}{'erro':>8}{'ganho':>10}"
      f"{'desloc.':>10}{'ampl.':>8}")
for _, r in res.iterrows():
    print(f"  {r.k:>8.0f}{r.lam:>9.0e}{r.r:>+9.3f}{r.mae:>7.1f}a{r.gain:>10.1f}"
          f"{r.disp:>9.2f}a{r.span:>7.1f}a")


section("CHECAGEM 2 E 3 — A FAMILIA CONTEM UM RELOGIO, E CONTEM O DEGENERADO?")
best = res.loc[res.r.idxmax()]
c2 = best.r > 0.75
print(f"  melhor configuracao: {best.k:.0f} sondas, lambda {best.lam:.0e}  "
      f"-> r = {best.r:.3f}, erro {best.mae:.1f}a  {'ok' if c2 else 'FALHOU'}")
deg = res[res.lam == max(LAMBDAS)].copy()
deg["gpp"] = deg.gain / deg.k
c3 = bool((deg.gpp.max() < 1e-4) and (deg.disp.max() < 0.5) and (deg.r.abs().max() < 0.3))
print(f"  canto degenerado (lambda {max(LAMBDAS):.0e}): deslocamento ate "
      f"{deg.disp.max():.3f}a, |r| ate {deg.r.abs().max():.3f}, peso por sonda ate "
      f"{deg.gpp.max():.1e}  -> {'ok' if c3 else 'FALHOU'}")
print("  (a terceira clausula era peso TOTAL < 1.0 e falhou em 1.877; ver docstring)")
if not (c2 and c3):
    sys.exit("\n  UMA CHECAGEM FALHOU — parando antes de ler o resultado.")


section("A PERGUNTA DA ETAPA 4 — DESLOCAMENTO E FUNCAO DO GANHO?")
print("""  A etapa 4 disse: o deslocamento e o que um pareamento ALEATORIO dos pesos
  com os desvios de tipo celular produz, logo o que separa relogios e ganho.
  Mas "ganho" tem mais de uma leitura, e elas nao sao a mesma coisa:

    L1  = soma dos |coeficientes|          -- o que a etapa 4 usou
    L2  = raiz da soma dos coeficientes^2  -- o que um pareamento aleatorio
                                              realmente produz
    L2s = idem, pesado pelo desvio de cada sonda entre tipos celulares

  Espalhar o MESMO peso total sobre mais sondas mantem L1 e derruba L2. Se o
  deslocamento seguir L2 e nao L1, a etapa 4 estava certa no mecanismo e errada
  na grandeza -- e isso e uma diferenca pratica, porque diz que dilucao e uma
  defesa.\n""")
q = res[res.disp > 1e-3].copy()
ld = np.log10(q.disp)
print(f"  {'resumo de beta':<8}{'Spearman':>11}{'inclinacao':>13}{'dispersao residual':>21}")
for tag in ("L1", "L2", "L2s"):
    lg = np.log10(q[{"L1": "gain", "L2": "gain_l2", "L2s": "gain_l2s"}[tag]])
    rho = float(pd.Series(lg).corr(ld, method="spearman"))
    fit = np.polyfit(lg, ld, 1)
    resid = ld - np.polyval(fit, lg)
    print(f"  {tag:<8}{rho:>+11.3f}{fit[0]:>13.2f}{10**resid.std():>20.2f}x")
print("""
  Inclinacao 1.0 e dispersao 1.00x seriam proporcionalidade exata. A dispersao
  que sobra e o quanto a ESCOLHA das sondas importa depois de descontado o
  resumo -- ou seja, o quanto a etapa 4 deixou de fora.""")


section("A FRONTEIRA — DA PARA SER PRECISO SEM LER COMPOSICAO?")
front = []
for _, r in res.sort_values("r", ascending=False).iterrows():
    if not front or r.disp < min(f["disp"] for f in front):
        front.append(dict(k=r.k, lam=r.lam, r=r.r, mae=r.mae,
                          gain=r.gain, disp=r.disp))
print(f"  {'sondas':>8}{'lambda':>9}{'r fora':>9}{'erro':>8}{'ganho':>10}{'desloc.':>10}")
for f in front:
    print(f"  {f['k']:>8.0f}{f['lam']:>9.0e}{f['r']:>+9.3f}{f['mae']:>7.1f}a"
          f"{f['gain']:>10.1f}{f['disp']:>9.2f}a")
pd.DataFrame(front).to_csv(OUT / "frontier_pareto.csv", index=False)

# The published clocks, measured with THE SAME statistic on THE SAME six
# fractions. Stage 2's number is over ten fractions and is not comparable —
# reading it next to this table would be the cross-stage mistake this project
# keeps making, so it is recomputed here instead of imported.
from model.clocks import CLOCKS, predict
b35_sel = b35.loc[:, sel_frac]
print(f"\n  para comparar, os relogios publicados sobre AS MESMAS 6 fracoes,\n"
      f"  com a mesma estatistica (recalculado, nao importado da etapa 2):")
pub = []
for c in CLOCKS:
    pv, _ = predict(b35_sel, c)
    w = pd.DataFrame({"d": dono, "f": frc, "v": pv.to_numpy()}).pivot_table(
        index="d", columns="f", values="v")
    pub.append(dict(clock=c, disp=float(w.std(axis=1).mean()),
                    span=float((w.max(axis=1) - w.min(axis=1)).mean())))
pub = pd.DataFrame(pub)
pub.to_csv(OUT / "frontier_published.csv", index=False)
for _, r in pub.iterrows():
    better = res[(res.disp < r.disp) & (res.r > 0.75)]
    note = (f"  <- {len(better)} configuracoes da familia sao mais planas E "
            f"tem r>0.75" if len(better) else "")
    print(f"    {r.clock:<14} deslocamento {r.disp:.2f}a  amplitude {r.span:.1f}a{note}")
print("""
  UMA RESSALVA QUE MUDA A LEITURA. Os dois eixos nao sao igualmente justos.
  O deslocamento e medido do mesmo jeito para todos, no GSE35069, e nenhum dos
  dois lados treinou ali -- esse eixo e limpo. A acuracia nao e: a familia foi
  treinada nesta coorte e validada cruzado DENTRO dela, enquanto os relogios
  publicados chegam de fora, sem nunca terem visto estes 188 tubos. A familia
  joga em casa.

  Entao a afirmacao que sobrevive nao e "da para fazer melhor que o Hannum".
  E esta: o deslocamento cai pela metade sem perder acuracia DENTRO da coorte,
  e quem manda nele e a norma L2 dos coeficientes. Diluir o mesmo peso por mais
  sondas e uma defesa que funciona. Se a acuracia dessa versao diluida
  atravessa para outra coorte, esta etapa nao testou.""")
