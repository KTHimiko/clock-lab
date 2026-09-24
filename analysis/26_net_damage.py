#!/usr/bin/env python3
"""
Stage 26 — predicting the SIGN of the damage, not only its rank.

The transport index predicts the estimation error a transported correction
leaves, never the net effect: net = what is left minus what was there. That is
why its safety floor was one-sided and why the zone between 0.05 and 0.16 held
both outcomes. The reviewer proposed closing it in advance. With b the fitted
coefficients and S_B the target's composition covariance,

    net damage = b' S_B b - 2 b' S_B beta_B

and if beta is shared between cohorts, E[b' S_B b] = beta' S_B beta +
sigma^2 * index, which gives the pre-transport estimate

    E[net] ~ 2 sigma_A^2 * index - b' S_B b                       (P_est)

computable from the fitting cohort (b, sigma_A^2) and the target's proportions
(S_B) alone — no clock values in the target, no outcome.

Stage 25 showed beta is NOT always shared, and that the worst transport in the
project is a specification failure. So P_est is expected to fail exactly there.
An oracle that knows the target's own full-cohort beta,

    oracle = b' S_B b - 2 b' S_B beta_B_hat                        (P_oracle)

gives the ceiling: how well any calculation in twelve-type space can track a
measurement made with the six-type panel. The gap between P_est and P_oracle is
what the unseeable component costs.

Both are divided by the target's age-residual variance so they sit on the same
scale as the measured share. Configurations as in stage 21: every directed pair
of the four cohorts, fitting sizes 40 to the full cohort, 20 draws, with the two
clocks clean everywhere (Levine 2018, Horvath 2018).

PRE-REGISTERED
  1. the cache is the stage 18 cache. Hard stop.
  2. the ceiling: Spearman(P_oracle, observed net) over configurations is
     reported first. If it is below 0.7, the twelve-versus-six measurement caps
     every predictor, and the bars below are read against it.
  3. P_est calls the SIGN of the configuration's median net damage correctly in
     at least 75% of configurations.
  4. P_est ranks the observed net damage at Spearman rho >= 0.7.
  5. its sign errors concentrate where stage 25's specification term is large:
     the median specification term of the misses exceeds that of the hits.
     Reported with the numbers; this is the stage 25 mechanism predicting where
     the stage 26 predictor breaks.

Usage:  .venv/bin/python analysis/26_net_damage.py
"""
import sys
from itertools import permutations, combinations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_extended import TYPES12
from model.clocks import CLOCKS
from model.deconvolution import TYPES as TYPES6

RES = ROOT / "results"; CACHE = RES / "cache"
RNG = np.random.default_rng(20260924)
N_REPS = 20
COHORTS = ["GSE40279", "GSE61151", "GSE50660", "GSE42861"]
CLOCKS3 = ["Horvath2013", "Levine2018", "Horvath2018"]
ALPHA_FIX = 3.0

# the corrected rule: a clock is out of any pair in which a cohort it trained on
# appears, whichever end of the transport that cohort is on
TRAINED_ON = {"Hannum2013": {"GSE40279"}, "Horvath2013": {"GSE40279"}}


def allowed(clock, *cohorts):
    return not (TRAINED_ON.get(clock, set()) & set(cohorts))


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
    Ct = partial_out(age, C[:, :-1]); Ct = Ct - Ct.mean(axis=0)
    return Ct.T @ Ct / len(Ct)


def transport_index(af, Cf, at, Ct):
    return float(np.trace(np.linalg.inv(sigma(af, Cf)) @ sigma(at, Ct)) / len(af))


def ridge_coefs(age, C, y, alphas):
    Ct = partial_out(age, C[:, :-1])
    yt = partial_out(age, y.reshape(-1, 1)).ravel()
    mu, sd = Ct.mean(axis=0), Ct.std(axis=0) + 1e-12
    Z = (Ct - mu) / sd
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    scale = float((s ** 2).mean()); Uty = U.T @ yt
    return {a: (Vt.T @ (s * Uty / (s ** 2 + a * scale))) / sd for a in alphas}


# ------------------------------------------------------------------ data ----
panel = pd.read_csv(CACHE / "panel.csv", index_col=0).iloc[:, 0]
D = {}
for tag in COHORTS:
    a = pd.read_csv(CACHE / f"{tag}_ages.csv", index_col=0)
    k = a.chrono.notna().to_numpy(); ix = a.index[k]
    D[tag] = dict(chrono=a.chrono[k].to_numpy(),
                  y={c: a[c][k].to_numpy() for c in CLOCKS},
                  C12=pd.read_csv(CACHE / f"{tag}_comp12.csv", index_col=0)
                        .loc[ix, TYPES12].to_numpy(),
                  C6=pd.read_csv(CACHE / f"{tag}_comp6.csv", index_col=0)
                       .loc[ix, TYPES6].to_numpy())
N = {t: len(D[t]["chrono"]) for t in COHORTS}

section("CHECAGEM 1 — CACHE")
c1 = abs(panel.panel_r - 0.789) < 5e-3 and abs(panel.panel_mae - 0.027) < 5e-3
print(f"  painel r = {panel.panel_r:.3f} | MAE = {panel.panel_mae:.3f} -> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("  parando: cache diferente.")


def fwl(age, C, y):
    X0 = np.column_stack([np.ones(len(age)), age])
    P = X0 @ np.linalg.pinv(X0)
    Ct = C[:, :-1] - P @ C[:, :-1]; yt = y - P @ y
    b = np.linalg.lstsq(Ct, yt, rcond=None)[0]
    s2 = float(((yt - Ct @ b) ** 2).sum() / (len(y) - 2 - Ct.shape[1]))
    return b, s2


section("CHECAGEM 1 — CACHE")
c1 = abs(panel.panel_r - 0.789) < 5e-3 and abs(panel.panel_mae - 0.027) < 5e-3
print(f"  -> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("  parando.")

CL = ["Levine2018", "Horvath2018"]
print("  nulos ...", flush=True)
NULL, BEFORE, RV, BFULL = {}, {}, {}, {}
for t in COHORTS:
    Xa = np.column_stack([np.ones(N[t]), D[t]["chrono"]])
    for c in CL:
        NULL[(t, c)] = null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        ib, base = inc_and_base(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        BEFORE[(t, c)] = ((ib - NULL[(t, c)]) / (1 - base), base)
        RV[(t, c)] = float((D[t]["y"][c] - Xa @ np.linalg.lstsq(Xa, D[t]["y"][c],
                                                                 rcond=None)[0]).var())
        BFULL[(t, c)] = fwl(D[t]["chrono"], D[t]["C12"], D[t]["y"][c])[0]
SIG = {t: sigma(D[t]["chrono"], D[t]["C12"]) for t in COHORTS}

spec_tab = pd.read_csv(RES / "beta_heterogeneity.csv")
SPEC = {(r.src, r.dst, r.clock): r.spec for r in spec_tab.itertuples()}

GRID = [40, 80, 160, 320, 480]
rows = []
for src, dst in permutations(COHORTS, 2):
    grid = sorted({g for g in GRID if g < N[src]} | {N[src]})
    S_B = SIG[dst]
    for n in grid:
        reps = 1 if n >= N[src] else 20
        for r in range(reps):
            idx = (np.arange(N[src]) if n >= N[src]
                   else stratified_draw(D[src]["chrono"], n, np.random.default_rng(800 + 7 * r + n)))
            af, Cf = D[src]["chrono"][idx], D[src]["C12"][idx]
            ti = transport_index(af, Cf, D[dst]["chrono"], D[dst]["C12"])
            cbar = Cf.mean(axis=0)
            for c in CL:
                b, s2 = fwl(af, Cf, D[src]["y"][c][idx])
                p_est = (2 * s2 * ti - b @ S_B @ b) / RV[(dst, c)]
                p_orc = (b @ S_B @ b - 2 * b @ S_B @ BFULL[(dst, c)]) / RV[(dst, c)]
                yc = D[dst]["y"][c] - (D[dst]["C12"][:, :-1] - cbar[:-1]) @ b
                ia, _ = inc_and_base(yc, D[dst]["chrono"], D[dst]["C6"])
                before, base = BEFORE[(dst, c)]
                obs = (ia - NULL[(dst, c)]) / (1 - base) - before
                rows.append(dict(src=src, dst=dst, n=len(idx), rep=r, clock=c,
                                 p_est=p_est, p_orc=p_orc, obs=obs, tindex=ti,
                                 spec=SPEC[(src, dst, c)]))
    print(f"  {src[3:]} -> {dst[3:]}: pronto", flush=True)
R = pd.DataFrame(rows)
R.to_csv(RES / "net_damage.csv", index=False)
cfg = (R.groupby(["src", "dst", "n", "clock"])
         .agg(p_est=("p_est", "median"), p_orc=("p_orc", "median"),
              tindex=("tindex", "median"),
              obs=("obs", "median"), spec=("spec", "first")).reset_index())

section("CHECAGEM 2 — O TETO: O ORACULO ACOMPANHA A MEDICAO?")
rho_o, _ = spearmanr(cfg.p_orc, cfg.obs)
sign_o = ((cfg.p_orc > 0) == (cfg.obs > 0)).mean()
print(f"  {len(cfg)} configuracoes | Spearman {rho_o:.3f} | sinal certo {sign_o:.0%}")
cap = rho_o < 0.7
print(f"  {'teto BAIXO: as barras abaixo sao lidas contra ele' if cap else 'teto ok'}")

section("CHECAGEM 3 E 4 — O PREVISOR DE ANTEMAO (P_est)")
rho_e, pe = spearmanr(cfg.p_est, cfg.obs)
hit = (cfg.p_est > 0) == (cfg.obs > 0)
sign_e = hit.mean()
c3 = sign_e >= 0.75; c4 = rho_e >= 0.7
print(f"  sinal certo em {sign_e:.0%} das configuracoes (barra 75%) -> {'ok' if c3 else 'FALHOU'}")
print(f"  Spearman com o dano liquido observado: {rho_e:.3f} (barra 0.7) -> {'ok' if c4 else 'FALHOU'}")
print(f"\n  por coorte de ajuste:")
for s in COHORTS:
    g = cfg[cfg.src == s]
    print(f"    {s[3:]:<8} sinal {((g.p_est > 0) == (g.obs > 0)).mean():>4.0%}  "
          f"rho {spearmanr(g.p_est, g.obs)[0]:+.3f}  ({len(g)} config.)")
print(f"\n  matriz de sinal (previsto x observado):")
tab = pd.crosstab(np.where(cfg.p_est > 0, "prev nocivo", "prev benefico"),
                  np.where(cfg.obs > 0, "obs nocivo", "obs benefico"))
print("   " + tab.to_string().replace("\n", "\n   "))

section("POSTERIOR — O PREVISOR SOFISTICADO BATE O INDICE SIMPLES?")
print("  escrito depois das checagens 3 e 4 falharem; nao e checagem, e comparacao\n")
rho_i, _ = spearmanr(cfg.tindex, cfg.obs)
print(f"  mesmas {len(cfg)} configuracoes, dano liquido observado:")
print(f"    indice de transporte    rho {rho_i:.3f}   (so ordena; nao da sinal)")
print(f"    P_est (parecer)         rho {rho_e:.3f}   sinal {sign_e:.0%}")
print(f"    oraculo (beta do alvo)  rho {rho_o:.3f}   sinal {sign_o:.0%}")
safe = cfg[cfg.tindex <= 0.05]
print(f"\n  piso do indice nestas configuracoes: abaixo de 0.05, "
      f"{int((safe.obs > 0).sum())} de {len(safe)} nocivas")
fb = cfg[cfg.p_est <= 0]
print(f"  'seguro' segundo P_est: {int((fb.obs > 0).sum())} de {len(fb)} nocivas")


section("CHECAGEM 5 — OS ERROS DE SINAL MORAM ONDE O BETA DIFERE?")
m_miss, m_hit = cfg[~hit].spec.median(), cfg[hit].spec.median()
c5 = m_miss > m_hit
print(f"  termo de especificacao mediano: erros {m_miss:+.1%} | acertos {m_hit:+.1%}"
      f" -> {'sim' if c5 else 'nao'}")
worst = cfg[~hit].sort_values("spec", ascending=False).head(6)
print(f"\n  erros com maior termo de especificacao:")
for r in worst.itertuples():
    print(f"    {r.src[3:]+' -> '+r.dst[3:]:<16} n={r.n:<4} {r.clock:<12} "
          f"prev {r.p_est:+.1%}  obs {r.obs:+.1%}  espec {r.spec:+.1%}")

section("FECHAMENTO")
print(f"  2. teto do oraculo: rho {rho_o:.3f}, sinal {sign_o:.0%}")
print(f"  3. sinal de P_est >= 75%: {'ok' if c3 else 'FALHOU'} ({sign_e:.0%})")
print(f"  4. rho de P_est >= 0.7: {'ok' if c4 else 'FALHOU'} ({rho_e:.3f})")
print(f"  5. erros onde o beta difere: {'sim' if c5 else 'nao'}")
