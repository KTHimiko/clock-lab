#!/usr/bin/env python3
"""
Stage 7 — how much of a person's "biological age" is their blood count?

Every stage so far has worked on cells someone sorted in a lab. Nobody buys a
sorted-cell test; they send a tube of blood. So: estimate each donor's cell
composition from their methylation, then ask how much of what the clock calls
their biological age is that composition.

THE REFERENCE. Cell composition is recovered by writing a whole-blood profile as
a non-negative mixture of purified cell-type profiles that sums to one. The
panel is built from GSE35069's six sorted types — the SAME dataset stages 2-4
used, which makes this a place to be careful rather than a convenience.

THE VALIDATION, AND WHY IT IS HARSH. GSE110554's twelve reconstructed mixtures
have proportions that are known, not estimated. So the panel is tested by
recovering them — a different cohort, a different array, a different lab from
the panel it is built on. That is strictly harder than the job it will then do
(450K panel on 450K whole blood), so passing here means passing there.

ACCEPTANCE, FIXED BEFORE THE RESULT IS READ. The deconvolution is used on real
people only if, on those twelve mixtures:
  1. correlation between estimated and known proportion > 0.80 across all 72
     cell-by-mixture values
  2. mean absolute error < 0.08 — eight percentage points
  3. no cell type systematically off by more than 0.15 in mean
And on GSE61151, a physiology check that has nothing to do with the panel:
  4. median estimated neutrophil fraction between 0.40 and 0.75, which is what
     adult whole blood contains

A KNOWN LIMITATION, STATED UP FRONT. The panel has six types. Real blood also
holds eosinophils and basophils, which this panel cannot name and will fold into
whatever resembles them — granulocyte signal, so mostly neutrophils. Estimates
are therefore "neutrophil-like granulocytes", not neutrophils.

Usage:  .venv/bin/python analysis/07_deconvolution.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import nnls

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix
from model.clocks import CLOCKS, predict

DATA = ROOT / "reference/data"
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(20250922)

# the six types the two studies have in common, panel name -> GSE35069 fraction
PANEL = {"CD4T": "CD4+ T cells", "CD8T": "CD8+ T cells",
         "Bcell": "CD19+ B cells", "NK": "CD56+ NK cells",
         "Mono": "CD14+ Monocytes", "Neu": "Neutrophils"}
TYPES = list(PANEL)
MIX_COL = {"CD4T": "cd4t", "CD8T": "cd8t", "Bcell": "bcell",
           "NK": "nk", "Mono": "mono", "Neu": "neu"}
N_PER_SIDE = 50          # probes per cell type per direction
N_CANDIDATE = 200        # generous pool, trimmed after intersecting platforms


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


def deconvolve(ref, Y, weight=100.0):
    """
    Non-negative proportions summing to one, per column of Y.

    The sum-to-one constraint is imposed by appending a row of `weight` to the
    reference and `weight` to the target: a heavily weighted equation saying the
    proportions add up. That keeps the whole thing inside plain NNLS.
    """
    R = np.vstack([ref, np.full((1, ref.shape[1]), weight)])
    out = np.empty((Y.shape[1], ref.shape[1]))
    for j in range(Y.shape[1]):
        y = np.concatenate([Y[:, j], [weight]])
        ok = ~np.isnan(y)
        w, _ = nnls(R[ok], y[ok])
        out[j] = w / w.sum() if w.sum() > 0 else np.nan
    return out


section("PAINEL DE REFERENCIA — GSE35069")
print("carregando ...", flush=True)
b35, m35 = read_series_matrix(DATA / "GSE35069_series_matrix.txt.gz")
m35 = m35.set_index("gsm").reindex(b35.columns)
frac = m35["tissue/cell type"]

# candidate probes: for each type, the ones that separate it from the other five
cand = set()
means = {}
for t, name in PANEL.items():
    sel = (frac == name).to_numpy()
    rest = np.isin(frac.to_numpy(), [PANEL[o] for o in TYPES if o != t])
    a, c = b35.loc[:, sel], b35.loc[:, rest]
    d = a.mean(axis=1) - c.mean(axis=1)
    sd = np.sqrt((a.var(axis=1, ddof=1) + c.var(axis=1, ddof=1)) / 2) + 0.01
    tstat = (d / sd).replace([np.inf, -np.inf], np.nan).dropna()
    cand |= set(tstat.nlargest(N_CANDIDATE).index)
    cand |= set(tstat.nsmallest(N_CANDIDATE).index)
    means[t] = tstat
print(f"  pool de candidatas: {len(cand):,} sondas "
      f"({N_CANDIDATE} por tipo por direcao)")

print("carregando GSE110554 (EPIC) ...", flush=True)
b110, m110 = read_series_matrix(DATA / "GSE110554_series_matrix.txt.gz")
m110 = m110.set_index("gsm").reindex(b110.columns)
print("carregando GSE61151 (450K, 710 MB) ...", flush=True)
b61, m61 = read_series_matrix(DATA / "GSE61151_series_matrix.txt.gz")
m61 = m61.set_index("gsm").reindex(b61.columns)

shared = sorted(cand & set(b110.index) & set(b61.index))
print(f"\n  presentes nas tres plataformas: {len(shared):,} de {len(cand):,}")

# final selection, taken only among probes all three arrays carry, so nothing is
# chosen that the targets cannot supply
probes = set()
for t in TYPES:
    s = means[t].reindex(shared).dropna()
    probes |= set(s.nlargest(N_PER_SIDE).index) | set(s.nsmallest(N_PER_SIDE).index)
probes = sorted(probes)
print(f"  painel final: {len(probes):,} sondas "
      f"({N_PER_SIDE} por tipo por direcao, apos a intersecao)")

ref = pd.DataFrame({t: b35.loc[probes, (frac == PANEL[t]).to_numpy()].mean(axis=1)
                    for t in TYPES})
ref.to_csv(OUT / "deconv_reference.csv")
del b35


section("VALIDACAO — RECUPERAR AS PROPORCOES CONHECIDAS DAS 12 MISTURAS")
is_mix = (m110["cell type"] == "MIX").to_numpy()
known = m110.loc[is_mix, [MIX_COL[t] for t in TYPES]].apply(
    pd.to_numeric, errors="coerce").to_numpy() / 100.0
est = deconvolve(ref.to_numpy(), b110.loc[probes].to_numpy()[:, is_mix])

titles = m110.loc[is_mix, "title"].to_numpy()
print(f"  {'mistura':<10}" + "".join(f"{t:>14}" for t in TYPES))
for i, t in enumerate(titles):
    print(f"  {t:<10}" + "".join(
        f"{known[i,j]:>6.2f}/{est[i,j]:<7.2f}" for j in range(len(TYPES))))
print("  (conhecido/estimado)")

r_all = float(np.corrcoef(known.ravel(), est.ravel())[0, 1])
mae = float(np.abs(known - est).mean())
bias = est.mean(axis=0) - known.mean(axis=0)
print(f"\n  correlacao sobre os 72 valores : {r_all:.3f}  "
      f"-> {'ok' if r_all > 0.80 else 'FALHOU'}")
print(f"  erro absoluto medio            : {mae:.3f}  "
      f"-> {'ok' if mae < 0.08 else 'FALHOU'}")
print(f"  {'vies por tipo':<32}: " + "  ".join(
    f"{t} {bias[j]:+.2f}" for j, t in enumerate(TYPES)))
print(f"  {'maior vies':<32}: {np.abs(bias).max():.3f}  "
      f"-> {'ok' if np.abs(bias).max() < 0.15 else 'FALHOU'}")
print(f"  {'por tipo, correlacao':<32}: " + "  ".join(
    f"{t} {np.corrcoef(known[:,j], est[:,j])[0,1]:+.2f}" for j, t in enumerate(TYPES)))

ok_valid = (r_all > 0.80) and (mae < 0.08) and (np.abs(bias).max() < 0.15)
del b110


section("APLICACAO — 188 SANGUES TOTAIS DO GSE61151")
comp = pd.DataFrame(deconvolve(ref.to_numpy(), b61.loc[probes].to_numpy()),
                    index=b61.columns, columns=TYPES)
print(f"  {'tipo':<8}{'mediana':>10}{'p10':>8}{'p90':>8}")
for t in TYPES:
    print(f"  {t:<8}{comp[t].median():>10.3f}{comp[t].quantile(.1):>8.3f}"
          f"{comp[t].quantile(.9):>8.3f}")
med_neu = float(comp["Neu"].median())
ok_phys = 0.40 <= med_neu <= 0.75
print(f"\n  mediana de granulocitos tipo neutrofilo: {med_neu:.3f}  "
      f"-> {'ok (faixa adulta 0.40-0.75)' if ok_phys else 'FORA DA FISIOLOGIA'}")

# POST HOC, e declarado. A checagem 4 que eu pre-especifiquei olhava so o
# neutrofilo, que e 64% do sangue e por isso o mais facil de acertar. Os tipos
# menores nao tinham checagem nenhuma. Abaixo estao todos, contra a faixa
# clinica de adulto em fracao de leucocitos totais. Isto NAO e criterio de
# aceitacao -- foi escrito depois -- mas fica no registro porque muda como as
# secoes seguintes devem ser lidas.
CLIN = {"CD4T": (0.07, 0.15), "CD8T": (0.03, 0.08), "Bcell": (0.01, 0.04),
        "NK": (0.02, 0.06), "Mono": (0.02, 0.08), "Neu": (0.40, 0.75)}
print("\n  todos os tipos contra a faixa clinica (post hoc, nao e criterio):")
fora = []
for t in TYPES:
    lo, hi = CLIN[t]
    med = float(comp[t].median())
    mark = "ok" if lo <= med <= hi else ("ABAIXO" if med < lo else "ACIMA")
    if mark != "ok":
        fora.append(t)
    print(f"    {t:<8} mediana {med:.3f}   faixa {lo:.2f}-{hi:.2f}   {mark}")
if fora:
    print(f"""
  {', '.join(fora)} fora da faixa. As proporcoes sao correlacionadas entre si,
  entao o ajuste conjunto pode estar certo enquanto a divisao entre tipos
  vizinhos esta errada -- e o erro aqui tem a cara disso: sobra CD4T e NK,
  falta CD8T e Bcell, que sao exatamente os pares que competem pelas mesmas
  sondas linfoides. Consequencia pratica: o teste CONJUNTO adiante (o
  incremento de R2, que usa os seis preditores de uma vez) resiste a isso; os
  coeficientes POR TIPO nao, e nao devem ser lidos como medida de tipo.

  Vale notar tambem o que a validacao nas misturas nao cobre. Uma mistura
  reconstruida e, por construcao, uma combinacao linear de celulas separadas
  -- parecida demais com o painel. Acertar 0.985 nela nao garante acertar em
  sangue de verdade, que traz eosinofilo e basofilo que o painel nao tem e
  celulas que ninguem tirou do corpo. A validacao era necessaria; suficiente
  ela nao e.""")

if not (ok_valid and ok_phys):
    sys.exit("\n  ACEITACAO FALHOU — a deconvolucao nao e usada em pessoas reais.")

age = pd.to_numeric(m61["agebloodtaken"], errors="coerce")
ages = {}
for name in CLOCKS:
    ages[name], _ = predict(b61, name)
comp["age"] = age.values
comp.to_csv(OUT / "deconv_gse61151.csv")
del b61


section("QUANTO DA IDADE EPIGENETICA E O HEMOGRAMA?")
print("""  Tres regressoes aninhadas por relogio. A composicao entra depois da idade
  cronologica, entao o incremento de R2 e o que ela explica ALEM do que a idade
  ja explicava -- o que importa, porque a propria composicao muda com a idade e
  esse pedaco pertence a idade, nao ao hemograma.

  O nulo embaralha as linhas da matriz de composicao entre as pessoas: seis
  preditores em 188 amostras compram R2 de graca, e isso mede quanto.\\n""")

keep = age.notna().to_numpy()
X_age = np.column_stack([np.ones(keep.sum()), age[keep].to_numpy()])
C = comp.loc[keep, TYPES].to_numpy()[:, :-1]   # one type dropped: they sum to 1


def r2(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    return 1 - resid.var() / y.var()


rows, coefs_by_type = [], {}
N_PERM = 5000
for name in CLOCKS:
    y = ages[name][keep].to_numpy()
    r2_age = r2(X_age, y)
    r2_both = r2(np.column_stack([X_age, C]), y)
    r2_comp = r2(np.column_stack([np.ones(len(C)), C]), y)
    inc = r2_both - r2_age

    null = np.empty(N_PERM)
    for i in range(N_PERM):
        null[i] = r2(np.column_stack([X_age, C[RNG.permutation(len(C))]]), y) - r2_age
    p = float((null >= inc).mean())

    # the same thing in years: how far the composition term moves a prediction
    beta, *_ = np.linalg.lstsq(np.column_stack([X_age, C]), y, rcond=None)
    comp_term = C @ beta[2:]
    rows.append(dict(clock=name, r2_age=r2_age, r2_comp=r2_comp, r2_both=r2_both,
                     increment=inc, null=float(null.mean()), p=p,
                     sd_years=float(comp_term.std()),
                     span_years=float(comp_term.max() - comp_term.min())))

    # Per-type coefficients, identified: all six proportions, no intercept, so
    # each coefficient is what the clock reads off a sample made entirely of that
    # type. With an intercept they would be identified only up to a constant,
    # because the six proportions sum to one.
    Cf = comp.loc[keep, TYPES].to_numpy()
    bf, *_ = np.linalg.lstsq(np.column_stack([age[keep].to_numpy(), Cf]), y, rcond=None)
    coefs_by_type[name] = pd.Series(bf[1:], index=TYPES)

res = pd.DataFrame(rows)
res.to_csv(OUT / "deconv_variance.csv", index=False)
print(f"  {'relogio':<14}{'R2 idade':>10}{'R2 comp':>9}{'R2 ambos':>10}"
      f"{'incr.':>8}{'nulo':>8}{'p':>8}{'dp anos':>9}{'ampl.':>8}")
for _, r in res.iterrows():
    print(f"  {r.clock:<14}{r.r2_age:>10.3f}{r.r2_comp:>9.3f}{r.r2_both:>10.3f}"
          f"{r.increment:>8.3f}{r.null:>8.3f}{r.p:>8.4f}{r.sd_years:>9.2f}{r.span_years:>8.1f}")
print("""
  R2 comp e a composicao SOZINHA, e ela e alta so porque composicao anda com a
  idade -- nao leia essa coluna como efeito. O incremento e o efeito.
  dp anos e a amplitude sao o termo de composicao convertido em anos: quanto o
  relogio desloca a leitura de uma pessoa por causa do sangue dela, na coorte.""")


section("RECONCILIACAO — O SANGUE TOTAL CONHECE A ORDEM DAS CELULAS SEPARADAS?")
print("""  Este e o fecho, e e tambem a camada que ja produziu os piores erros deste
  projeto, entao vale dizer exatamente o que esta sendo comparado.

  A esquerda: o deslocamento por tipo celular medido em CELULAS SEPARADAS numa
  bancada -- etapa 2 no GSE35069, etapa 5 no GSE110554.
  A direita: o coeficiente por tipo de uma regressao sobre 188 tubos de SANGUE
  TOTAL, onde ninguem separou nada e a composicao foi inferida da metilacao.

  Sao medidas de coisas diferentes por caminhos diferentes. Se a ordem bater,
  o efeito atravessa do experimento de bancada para o exame que se vende.

  Uma ressalva de circularidade, declarada: o painel de deconvolucao foi
  construido do GSE35069, o mesmo da etapa 2. A comparacao com a etapa 5, que
  e de outra coorte e outra plataforma, nao tem esse problema.\n""")
coefs = pd.DataFrame(coefs_by_type)
coefs.to_csv(OUT / "deconv_coefficients.csv")
from itertools import permutations
s2 = pd.read_csv(ROOT / "results/purified_ages.csv", index_col=0)
S2_NAME = {t: PANEL[t] for t in TYPES}
off5 = pd.read_csv(ROOT / "results/replication_offsets.csv", index_col=0)
rec = []
for name in CLOCKS:
    w = s2.pivot_table(index="donor", columns="fraction", values=name)
    a2 = pd.Series({t: (w.sub(w.mean(axis=1), axis=0)).mean()[S2_NAME[t]] for t in TYPES})
    a5 = off5[name].reindex(TYPES)
    b = coefs[name].reindex(TYPES)
    b = b - b.mean()
    line = dict(clock=name)
    for tag, a in (("etapa2", a2 - a2.mean()), ("etapa5", a5 - a5.mean())):
        ra, rb = a.rank().to_numpy(), b.rank().to_numpy()
        rho = float(np.corrcoef(ra, rb)[0, 1])
        null = np.array([np.corrcoef(ra, rb[list(q)])[0, 1] for q in permutations(range(6))])
        line[f"rho_{tag}"] = rho
        line[f"p_{tag}"] = float((null >= rho).mean())
    rec.append(line)
rec = pd.DataFrame(rec)
rec.to_csv(OUT / "deconv_reconciliation.csv", index=False)
print(f"  {'relogio':<14}{'rho vs etapa2':>15}{'p':>8}{'rho vs etapa5':>16}{'p':>8}")
for _, r in rec.iterrows():
    print(f"  {r.clock:<14}{r.rho_etapa2:>+15.2f}{r.p_etapa2:>8.3f}"
          f"{r.rho_etapa5:>+16.2f}{r.p_etapa5:>8.3f}")
print("\n  coeficientes por tipo, centrados (anos que o tipo puro desloca a leitura):")
print("  " + f"{'relogio':<14}" + "".join(f"{t:>9}" for t in TYPES))
for name in CLOCKS:
    b = coefs[name] - coefs[name].mean()
    print("  " + f"{name:<14}" + "".join(f"{b[t]:>+9.1f}" for t in TYPES))
