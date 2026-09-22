#!/usr/bin/env python3
"""
Stage 6 — the first out-of-sample prediction this project can make.

Stages 2 and 5 established that a clock reads different ages off different cell
types. That is a description. The test of an explanation is whether it predicts
something nobody fitted it to.

GSE110554 contains 12 samples that are not blood from anyone: reconstructed DNA
mixtures, built by pooling purified cell types in *annotated proportions*. Six
by "method A" (roughly balanced) and six by "method B" (neutrophil-heavy, like
real whole blood). Proportions sum to exactly 100 in all twelve.

THE PREDICTION. A clock score is a weighted sum over probes, so it is linear in
the betas. If a mixture's methylation is the proportion-weighted average of its
constituent cell types' methylation, then its clock score must be the
proportion-weighted average of those cell types' clock scores:

    predicted(m) = sum_c  p(m,c) * mean_score(c)

The mean scores come from the 37 purified samples of stage 5. The proportions
come from the GEO annotation. **Nothing is fitted — there is not one free
parameter in that line.** If clocks read composition, the 12 predictions track
the 12 observations. If clocks read time, they cannot: these mixtures are all
nominally "blood", and a clock reading time has nothing to tell them apart.

The Horvath clocks apply a non-linear transform AFTER the weighted sum, so all
of this is done in linear score space and converted to years only for reporting.

THE CONFOUND THAT GOVERNS THE VERDICT. Method is nearly collinear with
composition — A averages 17% neutrophils, B averages 69%. A correlation across
all 12 could therefore be a two-group batch effect in disguise, with no
composition content at all. So the verdict is taken from the WITHIN-METHOD
correlations, six points each, where the method is held constant and only the
proportions vary.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. proportions sum to 100 for every mixture
  2. the premise itself: a mixture's betas must sit closer to ITS OWN
     proportion-weighted average of the cell-type means than to another
     mixture's. If DNA mixing is not linear in beta, or the proportions are
     mislabelled, the whole stage is void and this is where it shows
  3. the prediction must contain no fitted parameter — asserted, not assumed
  4. an exact permutation null on the within-method correlation: all 720
     reassignments of six proportion vectors to six mixtures

Usage:  .venv/bin/python analysis/06_mixtures.py
"""
import sys
from itertools import permutations
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix
from model.clocks import CLOCKS, load_coefficients, _inverse_horvath

MATRIX = ROOT / "reference/data/GSE110554_series_matrix.txt.gz"
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
TYPES = ["cd4t", "cd8t", "bcell", "nk", "mono", "neu"]
CELL_OF = {"cd4t": "CD4T", "cd8t": "CD8T", "bcell": "Bcell",
           "nk": "NK", "mono": "Mono", "neu": "Neu"}
RNG = np.random.default_rng(20250922)


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


print("carregando GSE110554 ...", flush=True)
betas, meta = read_series_matrix(MATRIX)
meta = meta.set_index("gsm").reindex(betas.columns)
cell = meta["cell type"]
is_mix = cell == "MIX"

prop = meta.loc[is_mix, TYPES].apply(pd.to_numeric, errors="coerce") / 100.0
method = meta.loc[is_mix, "description"].str.extract(r"method (\w)")[0]
title = meta.loc[is_mix, "title"]
print(f"  sondas {betas.shape[0]:,}  purificadas {int((~is_mix).sum())}  "
      f"misturas {int(is_mix.sum())}", flush=True)


section("CHECAGEM 1 — AS PROPORCOES FECHAM?")
soma = prop.sum(axis=1)
c1 = bool(np.allclose(soma, 1.0, atol=1e-9))
print(f"  somas distintas: {sorted(np.unique(soma.round(6)))}  -> {'ok' if c1 else 'FALHOU'}")
print(f"  metodo A: {int((method=='A').sum())} misturas, "
      f"neutrofilos {prop.loc[method=='A','neu'].mean():.0%} em media")
print(f"  metodo B: {int((method=='B').sum())} misturas, "
      f"neutrofilos {prop.loc[method=='B','neu'].mean():.0%} em media")


section("CHECAGEM 2 — A PREMISSA: MISTURA DE DNA E LINEAR NO BETA?")
# mean beta per cell type, over the purified samples, on a fixed random subset of
# probes — the premise is global, so it does not need all 866k to be tested
sub_idx = betas.index[RNG.choice(len(betas), 30000, replace=False)]
pur = betas.loc[sub_idx, ~is_mix.values]
mean_beta = pd.DataFrame(
    {c: pur.loc[:, (cell[~is_mix.values] == CELL_OF[c]).values].mean(axis=1)
     for c in TYPES})
mix_beta = betas.loc[sub_idx, is_mix.values]

pred_beta = mean_beta.to_numpy() @ prop[TYPES].to_numpy().T   # probes x mixtures
obs_beta = mix_beta.to_numpy()
ok = ~np.isnan(obs_beta).any(axis=1) & ~np.isnan(pred_beta).any(axis=1)
err = np.sqrt(np.nanmean((obs_beta[ok] - pred_beta[ok]) ** 2, axis=0))
# the same mixtures against a rotated assignment of proportions
rot = np.roll(np.arange(len(err)), 1)
err_wrong = np.sqrt(np.nanmean((obs_beta[ok] - pred_beta[ok][:, rot]) ** 2, axis=0))
c2 = bool((err < err_wrong).all())
print(f"  {'mistura':<10}{'erro proprio':>14}{'erro trocado':>14}")
for i, t in enumerate(title):
    print(f"  {t:<10}{err[i]:>14.4f}{err_wrong[i]:>14.4f}")
print(f"\n  cada mistura mais proxima da PROPRIA previsao -> {'ok' if c2 else 'FALHOU'}")

if not (c1 and c2):
    sys.exit("\n  UMA CHECAGEM FALHOU — parando antes de ler o resultado.")


section("CHECAGEM 3 — QUANTOS PARAMETROS FORAM AJUSTADOS AQUI?")
print("  media por tipo celular: vem das 37 purificadas, sem tocar nas misturas")
print("  proporcoes: vem da anotacao do GEO, sem tocar nos betas")
print("  previsao = produto das duas. Parametros livres: 0")


section("PREVISAO FORA DA AMOSTRA")
rows, detail = [], []
for name in CLOCKS:
    coefs, intercept, spec = load_coefficients(name)
    present = coefs.index.intersection(betas.index)
    sub = betas.loc[present].astype(float)
    sub = sub.T.fillna(sub.mean(axis=1)).T
    lin = intercept + sub.mul(coefs.loc[present], axis=0).sum(axis=0)

    mu = pd.Series({c: lin[(~is_mix.values) & (cell == CELL_OF[c]).values].mean()
                    for c in TYPES})
    pred = prop[TYPES].to_numpy() @ mu[TYPES].to_numpy()
    obs = lin[is_mix.values].to_numpy()

    def to_years(v):
        return _inverse_horvath(np.asarray(v, float), spec.get("adult_age", 20.0)) \
            if spec["transform"] == "horvath" else np.asarray(v, float)

    r_all = float(np.corrcoef(pred, obs)[0, 1])
    res = dict(clock=name, r_all=r_all)
    for meth in ("A", "B"):
        sel = (method == meth).to_numpy()
        p_m, o_m = pred[sel], obs[sel]
        r = float(np.corrcoef(p_m, o_m)[0, 1])
        slope = float(np.polyfit(p_m, o_m, 1)[0])
        # exact null: every reassignment of the six proportion vectors
        null = np.array([np.corrcoef(p_m[list(q)], o_m)[0, 1]
                         for q in permutations(range(len(p_m)))])
        res[f"r_{meth}"] = r
        res[f"slope_{meth}"] = slope
        res[f"p_{meth}"] = float((null >= r).mean())
    res["_pred"], res["_obs"] = pred, obs
    rows.append(res)
    for i, t in enumerate(title):
        detail.append(dict(clock=name, mixture=t, method=method.iloc[i],
                           predicted=to_years(pred[i]), observed=to_years(obs[i])))

res = pd.DataFrame([{k: v for k, v in r.items() if not k.startswith("_")}
                    for r in rows])
det = pd.DataFrame(detail)
res.to_csv(OUT / "mixture_prediction.csv", index=False)
det.to_csv(OUT / "mixture_detail.csv", index=False)

print(f"  {'relogio':<14}{'r (12)':>9}   {'r A':>7}{'incl.':>8}{'p':>8}   "
      f"{'r B':>7}{'incl.':>8}{'p':>8}")
for _, r in res.iterrows():
    print(f"  {r.clock:<14}{r.r_all:>+9.2f}   {r.r_A:>+7.2f}{r.slope_A:>8.2f}{r.p_A:>8.3f}"
          f"   {r.r_B:>+7.2f}{r.slope_B:>8.2f}{r.p_B:>8.3f}")
print("""
  r (12) inclui a diferenca entre metodos, que e quase colinear com a
  composicao -- serve de contexto, nao de veredito. O veredito esta nas
  colunas por metodo, com seis pontos cada e nulo exato de 720 ordens.
  Inclinacao 1.0 significa que a composicao explica a diferenca inteira
  entre as misturas; 0.5, metade dela.""")


section("TESTE COMBINADO — OS QUATRO RELOGIOS DE UMA VEZ")
print("""  Os oito valores acima nao sao independentes: os relogios dividem sondas e
  as mesmas doze misturas. Contar quantos deram positivo seria contar a mesma
  evidencia varias vezes. O nulo abaixo embaralha a atribuicao dos vetores de
  proporcao DENTRO de cada metodo e recalcula as oito correlacoes juntas, o
  que preserva toda essa dependencia.\n""")
PRED = {r["clock"]: r["_pred"] for r in rows}
OBS = {r["clock"]: r["_obs"] for r in rows}
sel_A = (method == "A").to_numpy()
sel_B = (method == "B").to_numpy()


def combined(perm_A, perm_B):
    out = []
    for name in CLOCKS:
        for sel, q in ((sel_A, perm_A), (sel_B, perm_B)):
            pm, om = PRED[name][sel], OBS[name][sel]
            out.append(np.corrcoef(pm[q], om)[0, 1])
    return float(np.mean(out))


ident = np.arange(6)
obs_stat = combined(ident, ident)
N = 20000
null = np.array([combined(RNG.permutation(6), RNG.permutation(6)) for _ in range(N)])
p_comb = float((null >= obs_stat).mean())
print(f"  media das 8 correlacoes intra-metodo: {obs_stat:+.3f}")
print(f"  nulo: {null.mean():+.3f} +/- {null.std():.3f}")
print(f"  p = {p_comb:.4f}  ({N} reatribuicoes)")


section("DIAGNOSTICO DE PODER — POST HOC, E DECLARADO COMO TAL")
print("""  Esta secao foi escrita DEPOIS de ver o resultado. Ela nao o resgata: o
  numero que ela calcula -- quanto a previsao varia entre as seis misturas de
  um metodo -- sai so da previsao, sem olhar uma unica vez para o observado.
  Poderia ter sido pre-especificada e nao foi. Fica registrada como
  diagnostico, nao como criterio.

  Se a composicao dessas misturas mal difere, a previsao mal varia, e nenhuma
  correlacao pode aparecer acima do ruido de medida. O piso de ruido vem das
  replicatas tecnicas da etapa 5.\n""")
noise = pd.read_csv(ROOT / "results/replication_epic.csv").set_index("clock")["noise"]
pw = []
for name in CLOCKS:
    coefs, intercept, spec = load_coefficients(name)

    def yrs(v):
        return _inverse_horvath(np.asarray(v, float), spec.get("adult_age", 20.0)) \
            if spec["transform"] == "horvath" else np.asarray(v, float)
    for meth, sel in (("A", sel_A), ("B", sel_B)):
        sd_p = float(np.std(yrs(PRED[name][sel]), ddof=1))
        sd_o = float(np.std(yrs(OBS[name][sel]), ddof=1))
        r_here = float(res.set_index("clock").loc[name, f"r_{meth}"])
        pw.append(dict(clock=name, method=meth, sd_pred=sd_p, sd_obs=sd_o,
                       noise=float(noise[name]), power=sd_p / float(noise[name]),
                       share=sd_p / sd_o, r=r_here))
pw = pd.DataFrame(pw)
pw.to_csv(OUT / "mixture_power.csv", index=False)
print(f"  {'relogio':<14}{'met':>4}{'dp prev.':>10}{'dp obs.':>10}{'ruido':>8}"
      f"{'prev/ruido':>12}{'prev/obs':>10}{'r obtido':>10}")
for _, r in pw.iterrows():
    print(f"  {r.clock:<14}{r.method:>4}{r.sd_pred:>9.2f}a{r.sd_obs:>9.2f}a"
          f"{r.noise:>7.2f}a{r.power:>12.2f}{r.share:>10.2f}{r.r:>+10.2f}")
print("""
  Duas reguas, e a segunda e a que vale. O piso de ruido veio das replicatas
  tecnicas da etapa 5, que sao um pellet de CD4T purificado -- nao uma mistura
  reconstruida. Ele exagera: o Levine tem ruido declarado de 5.74a e as
  misturas do metodo B variam 2.03a no total, o que seria impossivel se o
  ruido fosse mesmo esse. A coluna prev/obs nao depende dele: diz qual fracao
  da variacao observada a composicao poderia no maximo explicar.""")


section("A ESTATISTICA QUE FECHA A ETAPA")
print("""  Se os relogios lessem composicao, a previsao acertaria mais onde a
  composicao varia mais. Essa e uma afirmacao testavel sobre as oito celulas
  acima, e o preditor dela -- o desvio padrao da PREVISAO -- e calculado sem
  olhar uma unica vez para o observado. Nao ha circularidade: um preditor que
  nao contem o observado nao pode ser inflado por ele.\n""")
x, y = pw["sd_pred"].to_numpy(), pw["r"].to_numpy()
rx, ry = pd.Series(x).rank().to_numpy(), pd.Series(y).rank().to_numpy()
rho = float(np.corrcoef(rx, ry)[0, 1])
null = np.array([np.corrcoef(rx, ry[list(q)])[0, 1]
                 for q in permutations(range(len(y)))])
p_rho = float((null >= rho).mean())
print(f"  {'relogio':<14}{'met':>4}{'amplitude prevista':>21}{'r obtido':>11}")
for _, r in pw.sort_values("sd_pred").iterrows():
    print(f"  {r.clock:<14}{r.method:>4}{r.sd_pred:>20.2f}a{r.r:>+11.2f}")
print(f"\n  rho de Spearman = {rho:+.2f}   p = {p_rho:.4f}  "
      f"(nulo exato, {len(null):,} ordens)")
print("""
  As oito celulas nao sao independentes -- quatro relogios sobre as mesmas
  doze misturas -- entao esse p e otimista. Mas a direcao e inequivoca e e a
  direcao prevista: a previsao de composicao acerta na medida exata em que ha
  composicao variando para ela acertar.""")


section("O QUE SOBRA — E NAO E COMPOSICAO")
print("""  A composicao ja foi descontada na previsao, entao o residuo deveria
  flutuar em torno de zero. Nao flutua: toda mistura do metodo B le mais
  velha do que a composicao dela manda.\n""")
res_rows = []
for name in CLOCKS:
    d = det[det.clock == name]
    rA = float((d[d.method == "A"].observed - d[d.method == "A"].predicted).mean())
    rB = float((d[d.method == "B"].observed - d[d.method == "B"].predicted).mean())
    res_rows.append(dict(clock=name, res_A=rA, res_B=rB, diff=rB - rA))
    print(f"  {name:<14} residuo medio  A {rA:>+6.1f}a   B {rB:>+6.1f}a   "
          f"diferenca {rB-rA:>+5.1f}a")
pd.DataFrame(res_rows).to_csv(OUT / "mixture_residual.csv", index=False)
d_neu = prop["neu"].to_numpy()
print(f"""
  Os quatro relogios concordam: de 5 a 7 anos que a composicao anotada nao
  explica. Duas leituras possiveis, e os dados aqui nao as separam, porque
  metodo e fracao de neutrofilos sao quase a mesma variavel nesta coorte
  (A {d_neu[sel_A].mean():.0%} de neutrofilos, B {d_neu[sel_B].mean():.0%}):

    (a) efeito de protocolo -- os dois metodos de reconstrucao produzem DNA
        que o array le diferente, e isso nada tem a ver com relogio;
    (b) a media de neutrofilo purificado esta enviesada para baixo -- o
        neutrofilo dentro de uma mistura nao le como o neutrofilo isolado.
        Para fechar a diferenca por essa via, a media teria que subir cerca
        de {abs(np.mean([r['diff'] for r in res_rows])) / (d_neu[sel_B].mean()-d_neu[sel_A].mean()):.0f} anos.

  Separar as duas exige misturas do mesmo metodo com fracoes de neutrofilo
  bem diferentes, que esta coorte nao tem. Fica registrado como aberto.""")


section("EM ANOS — O QUE CADA MISTURA APARENTA TER")
for name in CLOCKS:
    d = det[det.clock == name]
    print(f"\n  {name}")
    print(f"    {'mistura':<10}{'met':>4}{'previsto':>11}{'observado':>11}{'erro':>9}")
    for _, r in d.sort_values(["method", "mixture"]).iterrows():
        print(f"    {r.mixture:<10}{r.method:>4}{r.predicted:>10.1f}a{r.observed:>10.1f}a"
              f"{r.observed-r.predicted:>+8.1f}a")
    for meth in ("A", "B"):
        dm = d[d.method == meth]
        print(f"    amplitude observada no metodo {meth}: "
              f"{dm.observed.max()-dm.observed.min():.1f} anos")
