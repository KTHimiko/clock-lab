#!/usr/bin/env python3
"""
Stage 5 — does the stage 2 finding survive in other people, on another array?

Stage 2 found epigenetic age moving up to 35 years between cell fractions drawn
from ONE person on ONE day. That is either biology or it is GSE35069: six men,
one lab, one 450K run, one sorting protocol. This stage takes the claim to
GSE110554 — Salas et al., Illumina EPIC, cells separated by magnetic beads in a
different lab.

THE DESIGN IS NOT THE SAME, and that governs everything below. GSE35069 is
paired: ten fractions per donor, so the within-person spread is measured
directly. GSE110554 is not — each purified sample is one cell type from one
DIFFERENT donor, 37 samples, ages 19 to 59. The paired statistic cannot be
computed here at all.

So the replication is of the claim, not of the arithmetic. If cell type shifts
epigenetic age, then in an unpaired cohort each cell type must carry its own
systematic offset from its donor's real age. The quantity is `agegap` (predicted
minus chronological), and the statistic is how far the six cell-type mean agegaps
spread apart. Under the null — clocks read time, and cell type is irrelevant —
that spread is what shuffling the cell-type labels produces.

A constant caveat, stated before the numbers: EPIC is missing 19 of Horvath2013's
353 probes and 6 of Hannum's 71, and `predict` drops a missing probe from the
sum. That shifts every prediction on this platform by the same constant, so
absolute agegap is not comparable to stage 2 or to the papers. Every statistic
here is a CONTRAST between cell types, where the constant cancels.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. clock probe coverage on EPIC >= 0.8 for every clock, the threshold
     `predict` already enforces — EPIC is not a superset of 450K
  2. the clocks must still track chronological age in these 37 samples: median
     r across the four clocks > 0.5. Below that, the cohort or the parsing is
     broken and no contrast within it means anything
  3. the ruler. Th2535-1 and Th2535-2 are the same donor, same cell type, same
     plate — technical replicates, visible only after the stage 5 loader fix.
     Their difference IS the measurement noise. If the cell-type spread is not
     clearly larger than it, there is no effect to report

  4. the confound. A clock regresses toward the cohort mean, so an older donor
     gets a negative agegap for free. The cell types here are NOT age-balanced —
     B cells average 40.5 years and neutrophils 26.8 — and those two sit at
     opposite ends of the raw result, which is what the confound would produce.
     Every headline number is therefore recomputed on agegap residualised on
     chronological age, and the verdict is taken from the residualised version

Usage:  .venv/bin/python analysis/05_replication.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix
from model.clocks import CLOCKS, predict

MATRIX = ROOT / "reference/data/GSE110554_series_matrix.txt.gz"
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(20250922)
N_PERM = 20000

# Name mapping between the two studies. This is exactly the cross-stage layer
# that has produced six wrong results in this project, so it is written out and
# printed rather than inferred from string similarity.
TO_STAGE2 = {"CD4T": "CD4+ T cells", "CD8T": "CD8+ T cells",
             "Bcell": "CD19+ B cells", "NK": "CD56+ NK cells",
             "Mono": "CD14+ Monocytes", "Neu": "Neutrophils"}


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


print("carregando GSE110554 (EPIC, 208 MB comprimidos) ...", flush=True)
betas, meta = read_series_matrix(MATRIX)
meta = meta.set_index("gsm").reindex(betas.columns)
cell = meta["cell type"]
age = pd.to_numeric(meta["age"], errors="coerce")
print(f"  sondas {betas.shape[0]:,}  amostras {betas.shape[1]}", flush=True)


section("CHECAGEM 1 — COBERTURA DAS SONDAS NO EPIC")
ages, cov = {}, {}
for name in CLOCKS:
    try:
        pred, c = predict(betas, name)
    except ValueError as e:
        print(f"  {name:<14} RECUSADO: {e}")
        continue
    ages[name], cov[name] = pred, c
    print(f"  {name:<14} cobertura {c:>6.1%}  {'ok' if c >= 0.8 else 'ABAIXO DO PISO'}")
if not ages:
    sys.exit("nenhum relogio tem cobertura suficiente no EPIC — parar")
CLK = list(ages)

A = pd.DataFrame(ages)
A["cell"], A["age"] = cell.values, age.values
A["replicate"] = meta["replicate"].values if "replicate" in meta.columns else None
for name in CLK:
    A[name + "_gap"] = A[name] - A["age"]

pure = A[A.cell != "MIX"].copy()
print(f"\n  amostras purificadas {len(pure)} (MIX {int((A.cell=='MIX').sum())} "
      f"excluidas: sao misturas reconstruidas, sem doador nem idade)")
print("  " + ", ".join(f"{k} {v}" for k, v in pure.cell.value_counts().items()))
print(f"  idades {pure.age.min():.0f} a {pure.age.max():.0f}, "
      f"{pure.age.nunique()} valores distintos em {len(pure)} amostras")


section("CHECAGEM 2 — OS RELOGIOS AINDA ACOMPANHAM A IDADE CRONOLOGICA?")
rs = {}
for name in CLK:
    rs[name] = float(np.corrcoef(pure["age"], pure[name])[0, 1])
    print(f"  {name:<14} r = {rs[name]:>6.3f}   erro medio "
          f"{float(np.abs(pure[name+'_gap']).mean()):>5.1f}a")
med_r = float(np.median(list(rs.values())))
check2 = med_r > 0.5
print(f"\n  mediana dos r = {med_r:.3f}  -> {'ok' if check2 else 'FALHOU'}")


section("CHECAGEM 3 — A REGUA: REPLICATAS TECNICAS DO MESMO DOADOR")
rep = pure[pure["replicate"].notna()]
if len(rep) == 2:
    for name in CLK:
        v = rep[name].values
        print(f"  {name:<14} {v[0]:>6.1f}a  vs {v[1]:>6.1f}a   "
              f"diferenca {abs(v[0]-v[1]):>5.2f}a")
    noise = {n: float(abs(rep[n].values[0] - rep[n].values[1])) for n in CLK}
    print(f"\n  mesma pessoa, mesmo tipo celular ({rep.cell.iloc[0]}), mesma placa.")
    print("  Essa diferenca e ruido de medida, e serve de piso para tudo abaixo.")
else:
    noise = {n: np.nan for n in CLK}
    print(f"  {len(rep)} replicatas encontradas — checagem indisponivel")

# the replicates are one donor measured twice; averaging them avoids counting
# that person twice in a cohort where every other point is a distinct person
if len(rep) == 2:
    keep = pure.drop(index=rep.index[1:])
    for name in CLK:
        keep.loc[rep.index[0], name] = rep[name].mean()
        keep.loc[rep.index[0], name + "_gap"] = rep[name + "_gap"].mean()
    pure = keep
    print(f"  as duas foram promediadas numa amostra so -> n = {len(pure)}")

if not check2:
    print("\n  CHECAGEM 2 FALHOU — parando antes de ler o resultado.")
    sys.exit(1)


section("RESULTADO — DESLOCAMENTO DE IDADE EPIGENETICA POR TIPO CELULAR")
print(f"  {'tipo':<8}{'n':>3}" + "".join(f"{c.replace('2013','').replace('2018',''):>12}" for c in CLK))
order = pure.groupby("cell")[CLK[0] + "_gap"].mean().sort_values().index
for f in order:
    sub = pure[pure.cell == f]
    print(f"  {f:<8}{len(sub):>3}" + "".join(f"{sub[c+'_gap'].mean():>+12.1f}" for c in CLK))
print("\n  (valores em anos: media de predito menos idade real, dentro do tipo)")

section("CHECAGEM 4 — O CONFUNDIMENTO COM A IDADE DOS DOADORES")
by_age = pure.groupby("cell")["age"].agg(["mean", "min", "max"])
print("  idade cronologica media por tipo (os tipos NAO sao balanceados):")
for f, r in by_age.sort_values("mean").iterrows():
    print(f"    {f:<8} {r['mean']:>5.1f} anos  ({r['min']:.0f}-{r['max']:.0f})")
print(f"  amplitude entre as medias: {by_age['mean'].max()-by_age['mean'].min():.1f} anos")
print("""
  Um relogio encolhe a escala: prediz alto para jovem e baixo para velho, entao
  agegap cai com a idade mesmo sem nenhum efeito de tipo celular. Com os tipos
  desbalanceados em 13.7 anos, parte do resultado bruto pode ser so isso. O
  residuo abaixo tira de cada agegap a reta que o liga a idade cronologica.""")

labels = pure.cell.values
slopes = {}
for name in CLK:
    g = pure[name + "_gap"].values
    sl, ic = np.polyfit(pure["age"].values, g, 1)
    slopes[name] = sl
    pure[name + "_res"] = g - (sl * pure["age"].values + ic)
    print(f"  {name:<14} inclinacao agegap~idade = {sl:>+6.3f} a/a")


def spread_test(values, labels, rng, n_perm):
    obs = float(pd.Series(values).groupby(labels).mean().std(ddof=0))
    null = np.empty(n_perm)
    for i in range(n_perm):
        null[i] = pd.Series(rng.permutation(values)).groupby(labels).mean().std(ddof=0)
    return obs, float(null.mean()), float((null >= obs).mean())


rows = []
for name in CLK:
    o_r, n_r, p_r = spread_test(pure[name + "_gap"].values, labels, RNG, N_PERM)
    o_s, n_s, p_s = spread_test(pure[name + "_res"].values, labels, RNG, N_PERM)
    rows.append(dict(clock=name, coverage=cov[name], r_age=rs[name],
                     slope=slopes[name],
                     spread_raw=o_r, null_raw=n_r, p_raw=p_r,
                     spread=o_s, null=n_s, ratio=o_s / n_s, p=p_s,
                     noise=noise[name]))
res = pd.DataFrame(rows)

section("TESTE — O ESPALHAMENTO ENTRE TIPOS E MAIOR QUE O ACASO?")
print("  bruto = agegap cru | residuo = agegap sem a reta da idade (o que vale)\n")
print(f"  {'relogio':<14}{'bruto':>8}{'p':>8}   {'residuo':>8}{'nulo':>7}{'razao':>7}{'p':>8}{'ruido':>8}")
for _, r in res.iterrows():
    print(f"  {r.clock:<14}{r.spread_raw:>7.2f}a{r.p_raw:>8.4f}   "
          f"{r.spread:>7.2f}a{r.null:>6.2f}a{r.ratio:>7.2f}{r.p:>8.4f}{r.noise:>7.2f}a")
res.to_csv(OUT / "replication_epic.csv", index=False)


section("COMPARACAO COM A ETAPA 2 — A MESMA ORDEM DE TIPOS CELULARES?")
s2 = pd.read_csv(ROOT / "results/purified_ages.csv", index_col=0)
print("  mapeamento de nomes usado (escrito a mao, nao inferido):")
for k, v in TO_STAGE2.items():
    ok = v in set(s2.fraction)
    print(f"    {k:<8} -> {v:<18} {'ok' if ok else 'NAO EXISTE NA ETAPA 2'}")
assert all(v in set(s2.fraction) for v in TO_STAGE2.values())

print()
rho_rows = []
for name in CLK:
    # stage 2: centre each donor on his own mean, so what is left is the
    # fraction's offset, free of that donor's age
    w = s2.pivot_table(index="donor", columns="fraction", values=name)
    off2 = (w.sub(w.mean(axis=1), axis=0)).mean()
    a = pd.Series({k: off2[v] for k, v in TO_STAGE2.items()})
    b = pure.groupby("cell")[name + "_res"].mean().reindex(a.index)
    ra, rb = a.rank().values, b.rank().values
    rho = float(np.corrcoef(ra, rb)[0, 1])
    # with six points a rho reads high by accident, so it gets its own null:
    # all 720 orderings of the six labels, exactly, not a sample of them
    from itertools import permutations
    null = np.array([np.corrcoef(ra, rb[list(q)])[0, 1]
                     for q in permutations(range(6))])
    p = float((null >= rho).mean())
    rho_rows.append(dict(clock=name, rho=rho, p=p))
    print(f"  {name:<14} rho de Spearman entre as 6 medias = {rho:>+6.2f}   "
          f"p = {p:.3f}  (nulo exato, 720 ordens)")
    print("    " + "  ".join(f"{k}:{a[k]:+.0f}/{b[k]:+.0f}" for k in a.index))
pd.DataFrame(rho_rows).to_csv(OUT / "replication_rank.csv", index=False)
print("\n  (etapa2/etapa5, em anos; a etapa 2 centrada por doador, a 5 em residuo\n"
      "   de idade. Sao duas coortes, duas plataformas e dois protocolos de\n"
      "   separacao: a concordancia de ordem e o teste mais exigente aqui.)")
