#!/usr/bin/env python3
"""
Stage 22 — what happens when the two cohorts were deconvolved by different teams.

Stage 21 wrote the limitation into its own closing section: every cohort's
composition in this project comes from the same deconvolution with the same
panel, and that is not the situation the result is about. The real case is two
independent groups — one publishes composition coefficients estimated with their
reference panel, another applies them to proportions estimated with theirs.

WHAT THE MISMATCH ACTUALLY IS, because the obvious framing is impossible. "Fit
with six types and apply to twelve" cannot be done: the coefficient vectors do
not have the same length. Coefficients only transfer between analyses that use
the SAME cell-type labels. So the mismatch that matters is the same nominal
types estimated from different reference data:

  direct-6     the six-type panel built from Reinius (GSE35069) — the panel
               stages 7 to 12 used
  collapsed-6  the twelve-type panel built from Salas (GSE167998), collapsed
               onto the same six labels by the map in scripts/load_extended.py

Same six names, different donors, different probes, different array generation.
The literature says they differ where this design needs them to: the Reinius
reference was purified from Swedish male donors only and discriminates NK and
granulocytes least well, while Salas is the later and finer library.

THE TWO ARMS, and the measurement is deliberately the same in both

  matched      fit on A with direct-6(A)      -> apply to B with direct-6(B)
  mismatched   fit on A with collapsed-6(A)   -> apply to B with direct-6(B)

Both arms APPLY with direct-6(B), and both are measured with the twelve-type
composition of B, which comes from the other reference entirely. So the
measurement's relationship to the applied panel is identical in the two arms,
and the only thing that differs is which panel estimated the proportions the
coefficients were fitted on. Any gap between the arms is the mismatch and
nothing else.

The mismatched arm is the two-teams case written out: group A published
coefficients from the Salas panel, group B applies them to Reinius-derived
proportions, both of them calling the columns CD4T, CD8T, Bcell, NK, Mono, Neu.

THE INDEX NEEDS NO SPECIAL HANDLING. Sigma_fit is taken from the fitting panel
on A and Sigma_test from the applying panel on B, which in the mismatched arm
are different panels. If two references disagree about the same cell type, that
disagreement IS a covariance mismatch, and the index should already see it. That
is the interesting prediction and it is check 5.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. the cache is the stage 18 cache: panel r = 0.789, MAE = 0.027. Hard stop.
  2. THE TWO PANELS MUST BE MEASURING THE SAME THING. If direct-6 and
     collapsed-6 disagree wholesale, this stage is about a broken panel rather
     than about panel mismatch, and the result would mean nothing. Per cell type
     and per cohort, the Pearson correlation between the two estimates: the
     median over the six types and four cohorts must exceed 0.50. Reported per
     type as well, because the literature predicts specific weak spots and if
     the weak spots land where it says they should, that is a validation rather
     than a coincidence. Hard stop.
  3. in-sample floor in both arms, on all four fitting cohorts.
  4. THE QUESTION. Does the mismatch add damage? Paired over the twelve directed
     pairs at matched n, the median of delta(mismatched) - delta(matched).
     Pre-specified as reported rather than as a pass/fail: it is a magnitude,
     and the honest output is the number and its sign, not a verdict. What IS
     pass/fail is check 5.
  5. DOES THE INDEX STILL PREDICT WHEN THE PANELS DISAGREE? Spearman rho over
     the configurations of the mismatched arm, at the bar stage 21 used and
     passed on matched panels: rho > 0.70.
        passes -> the diagnostic survives the case it will actually be used in,
                  which is two groups who never spoke
        fails  -> the index needs both cohorts deconvolved the same way, and
                  that is a limitation that has to travel with every claim made
                  about it since stage 20
  6. DOES THE FIX STILL FIX? Ridge at alpha = 3 must hold median damage at or
     below zero across the twelve pairs of the mismatched arm.
  7. the declared monocyte bias from stage 13 carries, and matters more here:
     the monocyte channel of the twelve-type panel runs high, so monocyte is one
     of the six types where the two panels are expected to disagree most.

Usage:  .venv/bin/python analysis/22_panel_mismatch.py
"""
import sys
from itertools import permutations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_extended import TYPES12, COLLAPSE
from model.clocks import CLOCKS
from model.deconvolution import TYPES as TYPES6

OUT = ROOT / "results"; CACHE = OUT / "cache"
RNG = np.random.default_rng(20260926)
N_REPS = 20
GRID = [40, 80, 160, 320, 480]
CLOCKS3 = ["Horvath2013", "Levine2018", "Horvath2018"]
COHORTS = ["GSE40279", "GSE61151", "GSE50660", "GSE42861"]
ALPHA_FIX = 3.0


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


def r2(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return 1 - (y - X @ beta).var() / y.var()


def inc_and_base(y, a, C):
    Xa = np.column_stack([np.ones(len(a)), a])
    base = r2(Xa, y)
    return r2(np.column_stack([Xa, C[:, :-1]]), y) - base, base


def null_mean(y, a, C, n_perm=800, rng=RNG):
    Xa = np.column_stack([np.ones(len(a)), a])
    base = r2(Xa, y)
    return float(np.mean([r2(np.column_stack([Xa, C[rng.permutation(len(C))][:, :-1]]), y)
                          - base for _ in range(n_perm)]))


def stratified_draw(ages, n, rng):
    q = pd.qcut(pd.Series(ages), 10, labels=False, duplicates="drop").to_numpy()
    idx = []
    for b in np.unique(q):
        pool = np.where(q == b)[0]
        take = min(max(int(round(n * len(pool) / len(ages))), 1), len(pool))
        idx.extend(rng.choice(pool, take, replace=False))
    idx = np.array(sorted(idx))
    if len(idx) > n:
        idx = np.sort(rng.choice(idx, n, replace=False))
    elif len(idx) < n:
        rest = np.setdiff1d(np.arange(len(ages)), idx)
        idx = np.sort(np.concatenate([idx, rng.choice(rest, n - len(idx), replace=False)]))
    return idx


def partial_out(age, Y):
    Xa = np.column_stack([np.ones(len(age)), age])
    return Y - Xa @ np.linalg.pinv(Xa) @ Y


def sigma(age, C):
    Ct = partial_out(age, C[:, :-1])
    Ct = Ct - Ct.mean(axis=0)
    return Ct.T @ Ct / len(Ct)


def transport_index(age_f, C_f, age_t, C_t, lam=0.0):
    Sf, St = sigma(age_f, C_f), sigma(age_t, C_t)
    M = np.linalg.inv(Sf + lam * np.eye(len(Sf)) * np.trace(Sf) / len(Sf))
    return float(np.trace(M @ Sf @ M @ St) / len(age_f))


def ridge_coefs(age, C, y, alphas):
    Ct = partial_out(age, C[:, :-1])
    yt = partial_out(age, y.reshape(-1, 1)).ravel()
    mu, sd = Ct.mean(axis=0), Ct.std(axis=0) + 1e-12
    Z = (Ct - mu) / sd
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    scale = float((s ** 2).mean())
    Uty = U.T @ yt
    return {a: (Vt.T @ (s * Uty / (s ** 2 + a * scale))) / sd for a in alphas}


# ------------------------------------------------------------------ data ----
need = [CACHE / f"{t}_{w}.csv" for t in COHORTS
        for w in ("comp12", "comp6", "ages")] + [CACHE / "panel.csv"]
if not all(f.exists() for f in need):
    sys.exit("cache ausente. Rode analysis/18_conditioning.py antes.")

panel = pd.read_csv(CACHE / "panel.csv", index_col=0).iloc[:, 0]
D = {}
for tag in COHORTS:
    a = pd.read_csv(CACHE / f"{tag}_ages.csv", index_col=0)
    k = a.chrono.notna().to_numpy()
    ix = a.index[k]
    c12 = pd.read_csv(CACHE / f"{tag}_comp12.csv", index_col=0).loc[ix, TYPES12]
    # the same six labels, from the other reference, via the published map
    col6 = pd.DataFrame({t: c12[members].sum(axis=1) for t, members in COLLAPSE.items()})
    col6 = col6[TYPES6].div(col6[TYPES6].sum(axis=1), axis=0)
    D[tag] = dict(chrono=a.chrono[k].to_numpy(),
                  y={c: a[c][k].to_numpy() for c in CLOCKS},
                  C12=c12.to_numpy(),
                  direct6=pd.read_csv(CACHE / f"{tag}_comp6.csv", index_col=0)
                            .loc[ix, TYPES6].to_numpy(),
                  collapsed6=col6.to_numpy())
N = {t: len(D[t]["chrono"]) for t in COHORTS}

section("CHECAGEM 1 — CACHE")
c1 = abs(panel.panel_r - 0.789) < 5e-3 and abs(panel.panel_mae - 0.027) < 5e-3
print(f"  painel r = {panel.panel_r:.3f} | MAE = {panel.panel_mae:.3f} "
      f"-> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("\n  parando: cache diferente.")


section("CHECAGEM 2 — OS DOIS PAINEIS MEDEM A MESMA COISA?")
print("  correlacao entre direct-6 (Reinius) e collapsed-6 (Salas), por tipo\n")
print(f"  {'tipo':<10}" + "".join(f"{t[3:]:>10}" for t in COHORTS) + f"{'mediana':>10}")
rs = {}
for j, ty in enumerate(TYPES6):
    row, line = [], f"  {ty:<10}"
    for tag in COHORTS:
        v = float(np.corrcoef(D[tag]["direct6"][:, j], D[tag]["collapsed6"][:, j])[0, 1])
        row.append(v); line += f"{v:>10.3f}"
    rs[ty] = float(np.median(row))
    print(line + f"{rs[ty]:>10.3f}")
med_r = float(np.median(list(rs.values())))
c2 = med_r > 0.50
print(f"\n  mediana sobre os seis tipos: {med_r:.3f} (barra 0.50) "
      f"-> {'ok' if c2 else 'FALHOU'}")
worst = min(rs, key=rs.get)
print(f"  pior tipo: {worst} ({rs[worst]:.3f}) — a literatura prevê NK e "
      f"granulocitos\n  como os pontos fracos do Reinius, e o vies de monocito "
      f"da etapa 13 continua\n  declarado do lado do painel de doze.")
if not c2:
    sys.exit("\n  parando: os paineis discordam demais; a etapa seria sobre um "
             "painel quebrado.")

print("\n  nulos de permutacao ...", flush=True)
NULL, BEFORE = {}, {}
for t in COHORTS:
    for c in CLOCKS3:
        NULL[(t, c)] = null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["C12"])
        ib, base = inc_and_base(D[t]["y"][c], D[t]["chrono"], D[t]["C12"])
        BEFORE[(t, c)] = ((ib - NULL[(t, c)]) / (1 - base), base)


def run(src, dst, n, reps, fit_panel, alphas=(0.0, ALPHA_FIX)):
    """Fit on src with `fit_panel`, apply to dst with direct6, measure with C12."""
    s, d = D[src], D[dst]
    rows = []
    for r in range(reps):
        idx = (np.arange(N[src]) if n >= N[src]
               else stratified_draw(s["chrono"], n, np.random.default_rng(700 + 7 * r + n)))
        a_f = s["chrono"][idx]; C_f = s[fit_panel][idx]
        cbar = C_f.mean(axis=0)
        ti = transport_index(a_f, C_f, d["chrono"], d["direct6"])
        for c in CLOCKS3:
            coefs = ridge_coefs(a_f, C_f, s["y"][c][idx], list(alphas))
            before, base = BEFORE[(dst, c)]
            for alpha, b in coefs.items():
                yc = d["y"][c] - (d["direct6"][:, :-1] - cbar[:-1]) @ b
                ia, _ = inc_and_base(yc, d["chrono"], d["C12"])
                rows.append(dict(src=src, dst=dst, n=len(idx), rep=r, clock=c,
                                 alpha=alpha, arm=fit_panel, tindex=ti,
                                 delta=(ia - NULL[(dst, c)]) / (1 - base) - before))
    return rows


section("CHECAGEM 3 — PISO EM CASA NOS DOIS BRACOS")
c3 = True
for arm in ("direct6", "collapsed6"):
    for t in COHORTS:
        s = D[t]
        for c in CLOCKS3:
            b = ridge_coefs(s["chrono"], s[arm], s["y"][c], [0.0])[0.0]
            cbar = s[arm].mean(axis=0)
            yc = s["y"][c] - (s[arm][:, :-1] - cbar[:-1]) @ b
            b0, _ = inc_and_base(s["y"][c], s["chrono"], s[arm])
            a0, _ = inc_and_base(yc, s["chrono"], s[arm])
            c3 &= a0 < b0
    print(f"  braco {arm}: {'ok' if c3 else 'FALHOU'}")
print(f"  -> {'ok' if c3 else 'FALHOU'}")
if not c3:
    sys.exit("\n  parando.")


section("OS DOZE PARES, NOS DOIS BRACOS")
rows = []
for src, dst in permutations(COHORTS, 2):
    grid = sorted({g for g in GRID if g < N[src]} | {N[src]})
    for n in grid:
        reps = 1 if n >= N[src] else N_REPS
        rows += run(src, dst, n, reps, "direct6")
        rows += run(src, dst, n, reps, "collapsed6")
    print(f"  {src[3:]} -> {dst[3:]}", flush=True)
A = pd.DataFrame(rows)
A.to_csv(OUT / "panel_mismatch.csv", index=False)
ols = A[A.alpha == 0.0]


section("CHECAGEM 4 — QUANTO O DESENCONTRO DE PAINEL CUSTA (REPORTADO)")
print(f"  {'par':<22}{'n':>6}{'mesmo painel':>15}{'painel diferente':>19}{'custo':>10}")
# the matched n is almost never on the fitting cohort's grid — 184 is not in
# GSE40279's — so these configurations are run explicitly rather than selected
# out of the grid, which is what stage 21 did and what the first version of this
# stage failed to do (it selected empty frames and printed NaN).
MATCH = {}
for a, b in permutations(COHORTS, 2):
    n = min(N[a], N[b])
    MATCH[(a, b)] = pd.DataFrame(run(a, b, n, N_REPS, "direct6")
                                 + run(a, b, n, N_REPS, "collapsed6"))
pairs, costs = [], []
for src, dst in permutations(COHORTS, 2):
    n = min(N[src], N[dst])
    g = MATCH[(src, dst)][MATCH[(src, dst)].alpha == 0.0]
    dm = g[g.arm == "direct6"].delta.median()
    dx = g[g.arm == "collapsed6"].delta.median()
    pairs.append((src, dst)); costs.append(dx - dm)
    print(f"  {src[3:]+' -> '+dst[3:]:<22}{n:>6}{dm:>+15.1%}{dx:>+19.1%}"
          f"{dx - dm:>+10.1%}")
costs = np.array(costs)
print(f"\n  custo mediano do desencontro: {np.median(costs):+.1%}")
print(f"  pares onde o desencontro piora: {int((costs > 0).sum())} de {len(costs)}")
print(f"  pares nocivos (delta > 0): mesmo painel "
      f"{int((ols[ols.arm=='direct6'].groupby(['src','dst']).delta.median() > 0).sum())}"
      f" de 12, painel diferente "
      f"{int((ols[ols.arm=='collapsed6'].groupby(['src','dst']).delta.median() > 0).sum())}"
      f" de 12")


section("CHECAGEM 5 — O INDICE AINDA PREVE COM OS PAINEIS DESENCONTRADOS?")
for arm, label in (("direct6", "mesmo painel"), ("collapsed6", "painel diferente")):
    cfg = (ols[ols.arm == arm].groupby(["src", "dst", "n"])
             .agg(indice=("tindex", "median"), dano=("delta", "median")).reset_index())
    rho, pv = spearmanr(cfg.indice, cfg.dano)
    if arm == "collapsed6":
        c5 = rho > 0.70; rho_x = rho
    print(f"  {label:<20} {len(cfg):>3} config.  rho = {rho:+.3f}  (p = {pv:.2g})")
print(f"\n  barra no braco desencontrado: rho > 0.70 (a etapa 21 deu 0.907 com\n"
      f"  paineis casados) -> {'ok' if c5 else 'FALHOU'}")
if not c5:
    print("  -> o indice precisa das duas coortes deconvoluidas do mesmo jeito.\n"
          "     Essa limitacao passa a viajar com tudo que foi dito desde a etapa 20.")


section("CHECAGEM 6 — O CONSERTO AGUENTA O DESENCONTRO?")
print(f"  {'par':<22}{'OLS':>10}{'ridge a=3':>12}")
c6 = True
for src, dst in permutations(COHORTS, 2):
    n = min(N[src], N[dst])
    g = MATCH[(src, dst)]
    g = g[g.arm == "collapsed6"]
    d0 = g[g.alpha == 0.0].delta.median()
    d3 = g[g.alpha == ALPHA_FIX].delta.median()
    ok = d3 <= 0; c6 &= ok
    print(f"  {src[3:]+' -> '+dst[3:]:<22}{d0:>+10.1%}{d3:>+12.1%}"
          + ("" if ok else "   <- nao zera"))
print(f"\n  barra: mediana <= 0 nos doze -> {'ok' if c6 else 'FALHOU'}")


section("CHECAGENS, FECHAMENTO")
for i, (name, ok) in enumerate([
        ("cache identico", c1),
        ("os dois paineis medem a mesma coisa", c2),
        ("piso em casa nos dois bracos", c3),
        ("custo do desencontro: reportado, nao checado", True),
        ("indice preve com paineis desencontrados (rho > 0.70)", c5),
        ("ridge alpha=3 aguenta o desencontro", c6)], 1):
    print(f"  {i}. {name}: {'ok' if ok else 'FALHOU'}")
print("  7. vies de monocito da etapa 13: carregado, e relevante aqui")
print(f"\n  saida: results/panel_mismatch.csv ({len(A)} linhas)")
