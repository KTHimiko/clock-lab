#!/usr/bin/env python3
"""
Stage 14 — the floor, re-tested on the axis stage 11 could not see.

Stage 11 concluded that the exposure floor is biology: filtering probes by how
much they vary between cell types never beat an unfiltered clock. Tomusiak et
al. (2024, Communications Biology) report the opposite — IntrinClock, a clock
built to be invariant across ten immune cell types.

The two are not necessarily in conflict, and the reason is stage 13's finding.
Stage 11's filter was built from six cell types. It could not express the
naive-versus-memory axis, which is where the largest cell-type age differences
live. IntrinClock's published criterion is specific about exactly that axis:
keep CpGs with |r| > 0.3 against chronological age and |r| < 0.3 against a
sample being naive CD8. Stage 11 may simply have filtered on the wrong axis
because its reference panel did not contain the right one.

This stage puts three selection rules side by side, all trained the same way on
GSE40279, all tested on GSE61151, all measured against the same twelve-type
displacement:

  (a) no filter — stage 10's configuration
  (b) stage 11's filter, cell-type spread over six types
  (c) the IntrinClock rule, correlation with being naive CD8
  (d) a matched random control for each filter, discarding the same COUNT

DISPLACEMENT IS MEASURED ON TWELVE TYPES NOW. GSE35069, which every earlier
stage used, has six. GSE167998 has twelve but is not paired — each purified
sample is a different donor — so the statistic is stage 5's: age gap
residualised on chronological age, then spread across the twelve type means.
Both axes are normalised by scale, per stage 9.

SANITY CHECK THAT DECIDES WHETHER ANY OF THIS IS REAL
  1. the naive/memory gap must reproduce the published one. Tomusiak reports
     naive CD8 reading 15-20 years younger than effector memory CD8 from the
     same individual. If that does not appear here, the twelve-type labelling
     or the clocks are wrong and nothing below can be read. Required: naive CD8
     at least 8 years below memory CD8 on at least two of the four clocks
  2. the three cohorts share no sample
  3. calibration: slope on the test cohort within [0.5, 1.5]
  4. each filter is judged against a matched random control that randomises
     ONLY the cell-type criterion.

     The first version of this control was wrong and the fix matters. It drew
     the replacement pool from all 441,010 probes, so for rule (c) — whose pool
     is only 4,575 — it handed the control four thousand random probes and then
     asked it to build a clock. The control then lost on age accuracy, and the
     filter's apparent advantage included "I kept probes that correlate with
     age", which is not what is being tested. The control now draws from the
     probes that already pass |r| > 0.3 against age, so the only thing being
     randomised is the cell-type criterion.

     For rule (c) the decisive comparison needs no control at all: both (a) and
     (c) take the top 1,000 probes by correlation with age, and in a cohort of
     656 those all clear |r| > 0.3 anyway, so the age criterion is not binding
     and the only difference between the two rules is the naive-CD8 exclusion

Usage:  .venv/bin/python analysis/14_floor_naive.py
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

DATA = ROOT / "reference/data"
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(20250922)
LAMBDAS = [1e-2, 1e-1, 1e0, 1e1, 1e2]
K = 1000
N_CONTROL = 5


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


def ridge_dual(Xtr, ytr, lambdas):
    mx, my = Xtr.mean(axis=0), ytr.mean()
    Xc, yc = Xtr - mx, ytr - my
    w, V = np.linalg.eigh(Xc @ Xc.T)
    Vty = V.T @ yc
    for lam in lambdas:
        yield lam, Xc.T @ (V @ (Vty / (w + lam))), mx, my


def cv_lambda(Xs, ys, folds=5):
    parts = np.array_split(RNG.permutation(len(ys)), folds)
    err = {lam: 0.0 for lam in LAMBDAS}
    for test in parts:
        tr = np.setdiff1d(np.arange(len(ys)), test)
        for lam, beta, mx, my in ridge_dual(Xs[tr], ys[tr], LAMBDAS):
            err[lam] += float(np.abs((Xs[test] - mx) @ beta + my - ys[test]).sum())
    return min(err, key=err.get)


def displacement12(pred, ages_pur, types_pur):
    """Stage 5's statistic: age gap residualised on age, spread across types."""
    gap = pred - ages_pur
    sl, ic = np.polyfit(ages_pur, gap, 1)
    resid = gap - (sl * ages_pur + ic)
    return float(pd.Series(resid).groupby(types_pur).mean().std(ddof=0))


print("carregando GSE167998 (12 tipos, medida de deslocamento) ...", flush=True)
b167, p167 = read_extended(DATA / "GSE167998_matrix_processed.txt.gz",
                           DATA / "BloodExtended_Pheno.csv")
pur = (p167.CellType != "MIX").to_numpy() & p167.Age.notna().to_numpy()
b167p, p167p = b167.loc[:, pur], p167[pur]
ages_pur = p167p.Age.to_numpy(dtype=float)
types_pur = p167p.CellType.to_numpy()
print(f"  {int(pur.sum())} purificadas com idade, {len(set(types_pur))} tipos, "
      f"{p167p['Subject.ID'].nunique()} doadores", flush=True)

print("carregando GSE61151 (teste) ...", flush=True)
b61, m61 = read_series_matrix(DATA / "GSE61151_series_matrix.txt.gz")
m61 = m61.set_index("gsm").reindex(b61.columns)
age61 = pd.to_numeric(m61["agebloodtaken"], errors="coerce")

print("carregando GSE40279 (treino, 1.2 GB) ...", flush=True)
b40, m40 = read_series_matrix(DATA / "GSE40279_series_matrix.txt.gz")
m40 = m40.set_index("gsm").reindex(b40.columns)
age40 = pd.to_numeric(m40["age (y)"], errors="coerce")


section("CHECAGEM 1 — O VAO VIRGEM/MEMORIA APARECE COMO NA LITERATURA?")
print("""  Tomusiak et al. relatam CD8 virgem lendo 15 a 20 anos mais novo que CD8 de
  memoria efetora da mesma pessoa. Aqui os doadores sao diferentes, entao o que
  se compara e o desvio de idade de cada tipo depois de descontar a idade real
  de quem doou.\\n""")
pub_off = {}
for c in CLOCKS:
    pv, _ = predict(b167p, c)
    gap = pv.to_numpy() - ages_pur
    sl, ic = np.polyfit(ages_pur, gap, 1)
    pub_off[c] = pd.Series(gap - (sl * ages_pur + ic)).groupby(types_pur).mean()
off = pd.DataFrame(pub_off)
off.to_csv(OUT / "offsets12.csv")
print(f"  {'tipo':<8}" + "".join(f"{c.replace('2013','').replace('2018',''):>13}" for c in CLOCKS))
for t in off.index:
    print(f"  {t:<8}" + "".join(f"{off.loc[t,c]:>+13.1f}" for c in CLOCKS))
gaps = {c: float(off.loc["CD8mem", c] - off.loc["CD8nv", c]) for c in CLOCKS}
print(f"\n  CD8 memoria menos CD8 virgem:  " +
      "  ".join(f"{c.replace('2013','').replace('2018','')} {g:+.1f}a" for c, g in gaps.items()))
c1 = sum(g >= 8 for g in gaps.values()) >= 2
print(f"  relogios com vao >= 8 anos: {sum(g >= 8 for g in gaps.values())}/4  "
      f"-> {'ok' if c1 else 'FALHOU'}")

c2 = not (set(b167.columns) & set(b61.columns) or set(b40.columns) & set(b61.columns))
print(f"\n  2. nenhuma amostra compartilhada  -> {'ok' if c2 else 'FALHOU'}")

common = b40.index.intersection(b61.index).intersection(b167.index)
k40 = age40.notna().to_numpy()
X = b40.loc[common].to_numpy(dtype=np.float64).T[k40]
E = b61.loc[common].to_numpy(dtype=np.float64).T
P = b167p.loc[common].to_numpy(dtype=np.float64).T
good = ~np.isnan(X).any(axis=0) & ~np.isnan(E).any(axis=0) & ~np.isnan(P).any(axis=0)
X, E, P = X[:, good], E[:, good], P[:, good]
ids = pd.Index(common)[good]
y = age40[k40].to_numpy()
k61 = age61.notna().to_numpy()
E, ye = E[k61], age61[k61].to_numpy()
print(f"  sondas comuns e completas: {X.shape[1]:,}")

if not (c1 and c2):
    sys.exit("\n  UMA CHECAGEM FALHOU — parando antes de ler o resultado.")


section("OS TRES CRITERIOS DE SELECAO")
r_age = ((X - X.mean(0)).T @ (y - y.mean())) / (X.std(0) * y.std() * len(y) + 1e-12)
absr = np.abs(r_age)

# (b) stage 11's filter: spread across the SIX collapsed types, the information
# that stage 11 actually had
SIX = {"CD4T": ["CD4nv", "CD4mem", "Treg"], "CD8T": ["CD8nv", "CD8mem"],
       "Bcell": ["Bnv", "Bmem"], "NK": ["NK"], "Mono": ["Mono"],
       "Neu": ["Neu", "Eos", "Bas"]}
m6 = pd.DataFrame({k: P[:, :][np.isin(types_pur, v)].mean(axis=0)
                   for k, v in SIX.items()})
spread6 = m6.std(axis=1).to_numpy()

# (c) the IntrinClock rule, on the axis it names
is_nv8 = (types_pur == "CD8nv").astype(float)
r_nv8 = np.array([np.corrcoef(P[:, j], is_nv8)[0, 1] for j in range(P.shape[1])])
r_nv8 = np.nan_to_num(r_nv8)

rules = {
    "(a) sem filtro": np.arange(X.shape[1]),
    "(b) etapa 11: 6 tipos": np.where(spread6 <= np.quantile(spread6, 0.75))[0],
    "(c) IntrinClock: CD8 virgem": np.where((absr > 0.3) & (np.abs(r_nv8) < 0.3))[0],
}
for name, pool in rules.items():
    print(f"  {name:<30} elegiveis {len(pool):>8,}")


section("RESULTADO — ACURACIA E DESLOCAMENTO EM 12 TIPOS")
def run(pool, k=K):
    if len(pool) < k:
        return None
    idx = pool[np.argsort(-absr[pool])[:k]]
    lam = cv_lambda(X[:, idx], y)
    _, beta, mx, my = next(iter(ridge_dual(X[:, idx], y, [lam])))
    pe = (E[:, idx] - mx) @ beta + my
    pv = (P[:, idx] - mx) @ beta + my
    slope = float(np.cov(pe, ye)[0, 1] / np.var(ye))
    disp = displacement12(pv, ages_pur, types_pur)
    return dict(lam=lam, r=float(np.corrcoef(ye, pe)[0, 1]),
                mae=float(np.abs(pe - ye).mean()), slope=slope, disp=disp,
                disp_n=disp / slope if slope > 1e-9 else np.inf)

universe = np.where(absr > 0.3)[0]
print(f"  universo do controle (|r| com a idade > 0.3): {len(universe):,} sondas")
print(f"  as 1000 melhores por idade tem |r| de {absr[np.argsort(-absr)[:K]].min():.2f} "
      f"para cima, entao o criterio de idade da regra (c) nao morde a selecao "
      f"final")

rows = []
for name, pool in rules.items():
    m = run(pool)
    if m is None:
        print(f"  {name:<30} POUCAS SONDAS ({len(pool)})")
        continue
    ctrl = None
    if name != "(a) sem filtro":
        cs = [run(RNG.choice(universe, size=min(len(pool), len(universe)),
                             replace=False))
              for _ in range(N_CONTROL)]
        ctrl = {k: float(np.mean([c[k] for c in cs])) for k in ("r", "disp_n", "slope")}
    rows.append(dict(rule=name, n_pool=len(pool), **m,
                     ctrl_r=ctrl["r"] if ctrl else np.nan,
                     ctrl_disp_n=ctrl["disp_n"] if ctrl else np.nan))
res = pd.DataFrame(rows)
res.to_csv(OUT / "floor_naive.csv", index=False)

cal = res[(res.slope >= 0.5) & (res.slope <= 1.5)]
print(f"  3. calibracao: {len(cal)}/{len(res)} dentro de [0.5, 1.5]")
print(f"\n  {'criterio':<30}{'r teste':>9}{'escala':>8}{'desloc.':>9}"
      f"{'desloc/escala':>15}{'controle':>10}{'vantagem':>10}")
for _, r in res.iterrows():
    edge = (r.ctrl_disp_n - r.disp_n) if np.isfinite(r.ctrl_disp_n) else np.nan
    print(f"  {r.rule:<30}{r.r:>+9.3f}{r.slope:>8.3f}{r.disp:>8.2f}a"
          f"{r.disp_n:>14.1f}a"
          f"{'' if not np.isfinite(r.ctrl_disp_n) else f'{r.ctrl_disp_n:>9.1f}a'}"
          f"{'' if not np.isfinite(edge) else f'{edge:>+10.1f}'}")

base = res[res.rule == "(a) sem filtro"].iloc[0]
print(f"\n  os publicados, no mesmo eixo de 12 tipos:")
for c in CLOCKS:
    pv, _ = predict(b61, c)
    ok = age61.notna() & pv.notna()
    sl = float(np.cov(pv[ok], age61[ok])[0, 1] / np.var(age61[ok]))
    pvp, _ = predict(b167p, c)
    d = displacement12(pvp.to_numpy(), ages_pur, types_pur)
    print(f"    {c:<14} r {float(np.corrcoef(age61[ok], pv[ok])[0,1]):+.3f}  "
          f"escala {sl:.3f}  desloc/escala {d/sl:>5.1f}a")
print(f"""
  A coluna 'vantagem' e o quanto cada filtro fica ABAIXO do seu proprio
  controle aleatorio de mesmo tamanho. A comparacao com '(a) sem filtro'
  responde a pergunta da etapa 11 no eixo que ela nao tinha: da para ficar mais
  plano do que nao filtrar nada, agora que o eixo virgem/memoria esta visivel?""")
