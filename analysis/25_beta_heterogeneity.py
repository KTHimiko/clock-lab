#!/usr/bin/env python3
"""
Stage 25 — does the composition effect itself differ between cohorts?

Stage 24 found that a correction fitted on all 656 samples of GSE40279 — where
estimation noise is negligible — still leaves +4.5% of composition in the age
residual of other cohorts, and in one pair doubles it. At home, measured with the
other panel, the same correction removes 85-95%. So the floor is transport, and
not measurement.

The reviewer's fourth point named the candidate. GSE42861 is a rheumatoid
arthritis case-control study and GSE50660 a smoking cohort; both conditions
reshape blood, and may also reshape how composition maps onto a clock. The
transport index assumes the true coefficient beta is the same in both cohorts.
If it is not, the transported correction carries a specification error the index
cannot see, and no fitting size removes.

That has a name. When the input distribution shifts but the input-output map is
shared, it is covariate shift; when the map itself differs, it is model shift
(also regression shift, or posterior drift). The linear-regression literature
treats the two together (Lei et al., ICML 2021; Patil, Du & Tibshirani, ICML
2024). Stages 17-24 measured only the first.

THE DESIGN — full-cohort fits, off the stage 18 cache, plus the disease and
smoking fields read from the series headers (the table itself is not reloaded).

  For each cohort k and clock c, the composition coefficients beta_k are fitted
  on the whole cohort by Frisch-Waugh-Lovell (age partialled out), with the usual
  covariance V_k = s_k^2 (C'C)^-1.

  (a) heterogeneity: for every pair, a Wald test of beta_A = beta_B,
      (b_A - b_B)' (V_A + V_B)^-1 (b_A - b_B) on 11 degrees of freedom
  (b) the specification term: if estimation were perfect, the correction fitted
      on A and applied in B leaves (beta_A - beta_B)' S_B (beta_A - beta_B). Its
      plug-in estimate is biased upward by the estimation noise of both fits, so
      the bias-corrected version subtracts tr(S_B V_A) + tr(S_B V_B)
  (c) within-cohort: beta fitted separately in RA cases and controls (GSE42861),
      and in ever and never smokers (GSE50660)
  (d) whether carrying disease and smoking as covariates in the fit reduces the
      between-cohort heterogeneity

CLOCKS. Levine 2018 and Horvath 2018, clean in all four cohorts after stage 23.

SANITY CHECKS AND PRE-REGISTERED CRITERIA
  1. the cache is the stage 18 cache. Hard stop.
  2. the metadata join: every cohort sample in the cache must find its disease or
     smoking field, and each within-cohort group must be large enough to fit
     eleven coefficients — over 150. Hard stop.

     THIS CHECK WAS WRITTEN WRONG THE FIRST TIME, the seventh in this project
     revisited after failing. It demanded over 100 current smokers in GSE50660,
     a number assumed rather than looked up. GSE50660 has 22 current, 263 former
     and 179 never smokers: the join was intact, and the check stopped the run on
     a false belief about the cohort. Two things had been folded into one check —
     join integrity and group size — and are now separate. With 22 samples,
     eleven coefficients cannot be estimated, so the current-versus-never
     contrast is infeasible and the smoking contrast becomes EVER versus NEVER
     (285 against 179). It is the weaker contrast: a former smoker's blood
     partly recovers, so any difference it finds understates the current one.
  3. HETEROGENEITY IS ESTABLISHED if at least half of the (pair x clock) Wald
     tests reject at a Bonferroni-corrected 0.05 over the twelve tests.
  4. THE SPECIFICATION TERM EXPLAINS THE FLOOR if its bias-corrected estimate
     ranks the observed full-n leftover across the (directed pair x clock)
     configurations at Spearman rho >= 0.5.
  5. within-cohort heterogeneity (RA; smoking): reported with p-values, no bar —
     it is the mechanism question, not a pass/fail.
  6. adjustment: reported as the median change in the Wald statistic for pairs
     involving GSE42861 or GSE50660.

Usage:  .venv/bin/python analysis/25_beta_heterogeneity.py
"""
import sys, gzip
from itertools import permutations, combinations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, chi2

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


def header_fields(gse):
    """Sample characteristics from the series header only; the beta table below
    it is never read."""
    samples, fields = None, {}
    with gzip.open(ROOT / "reference/data" / f"{gse}_series_matrix.txt.gz", "rt",
                   errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            parts = [v.strip('"') for v in line.rstrip("\n").split("\t")]
            if parts[0] == "!Sample_geo_accession":
                samples = parts[1:]
            elif parts[0] == "!Sample_characteristics_ch1":
                for j, cell in enumerate(parts[1:]):
                    k, _, v = cell.partition(":")
                    fields.setdefault(k.strip(), [None] * len(samples))[j] = v.strip()
    return pd.DataFrame(fields, index=samples)


def fwl_fit(age, C, y, Z=None):
    """Composition coefficients with age (and optional covariates Z) partialled
    out, and their covariance."""
    X0 = np.column_stack([np.ones(len(age)), age] + ([] if Z is None else [Z]))
    P = X0 @ np.linalg.pinv(X0)
    Ct = C[:, :-1] - P @ C[:, :-1]
    yt = y - P @ y
    b = np.linalg.lstsq(Ct, yt, rcond=None)[0]
    dof = len(y) - X0.shape[1] - Ct.shape[1]
    s2 = float(((yt - Ct @ b) ** 2).sum() / dof)
    V = s2 * np.linalg.pinv(Ct.T @ Ct)
    return b, V


def wald(bA, VA, bB, VB):
    d = bA - bB
    stat = float(d @ np.linalg.pinv(VA + VB) @ d)
    return stat, float(chi2.sf(stat, len(d)))


section("CHECAGEM 1 — CACHE")
c1 = abs(panel.panel_r - 0.789) < 5e-3 and abs(panel.panel_mae - 0.027) < 5e-3
print(f"  -> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("  parando.")

section("CHECAGEM 2 — METADADOS DE DOENCA E TABAGISMO")
ages42 = pd.read_csv(CACHE / "GSE42861_ages.csv", index_col=0)
ages50 = pd.read_csv(CACHE / "GSE50660_ages.csv", index_col=0)
h42 = header_fields("GSE42861").reindex(ages42.index[ages42.chrono.notna()])
h50 = header_fields("GSE50660").reindex(ages50.index[ages50.chrono.notna()])
ra = (h42["disease state"].str.contains("rheumatoid", case=False)).astype(float).to_numpy()
smk42 = h42["smoking status"].map({"never": 0, "ex": 1, "occasional": 1,
                                   "current": 2}).to_numpy(dtype=float)
kcol = [c for c in h50.columns if c.startswith("smoking")][0]
smk50 = pd.to_numeric(h50[kcol], errors="coerce").to_numpy()
n_ra, n_ctl = int(ra.sum()), int((1 - ra).sum())
n_cur, n_nev = int((smk50 == 2).sum()), int((smk50 == 0).sum())
miss = int(h42["disease state"].isna().sum() + np.isnan(smk50).sum())
print(f"  GSE42861: AR {n_ra}, controles {n_ctl}; tabagismo sem codigo: "
      f"{int(np.isnan(smk42).sum())} ({sorted(h42['smoking status'].dropna().unique())})")
print(f"  GSE50660: atuais {n_cur}, nunca {n_nev}, ex {int((smk50 == 1).sum())}")
n_ever = int((smk50 >= 1).sum())
join_ok = miss == 0
size_ok = min(n_ra, n_ctl, n_ever, n_nev) > 150
c2 = join_ok and size_ok
print(f"  juncao sem faltantes: {'ok' if join_ok else 'FALHOU'}")
print(f"  grupos > 150 (AR {n_ra}, ctl {n_ctl}, alguma vez {n_ever}, nunca {n_nev}): "
      f"{'ok' if size_ok else 'FALHOU'}")
print(f"  (atual vs nunca e inviavel: {n_cur} fumantes atuais para 11 coeficientes)")
if not c2:
    sys.exit("  parando: juncao de metadados quebrada.")

CL = ["Levine2018", "Horvath2018"]
FIT = {(t, c): fwl_fit(D[t]["chrono"], D[t]["C12"], D[t]["y"][c]) for t in COHORTS for c in CL}
SIG = {t: sigma(D[t]["chrono"], D[t]["C12"]) for t in COHORTS}


section("CHECAGEM 3 — O BETA DIFERE ENTRE COORTES?")
alpha_b = 0.05 / (len(list(combinations(COHORTS, 2))) * len(CL))
print(f"  Wald de beta_A = beta_B, 11 g.l., Bonferroni a {alpha_b:.4f}\n")
print(f"  {'par':<16}" + "".join(f"{c:>24}" for c in CL))
W, sig = {}, []
for a, b in combinations(COHORTS, 2):
    line = f"  {a[3:]+'/'+b[3:]:<16}"
    for c in CL:
        st, p = wald(*FIT[(a, c)], *FIT[(b, c)])
        W[(a, b, c)] = st; sig.append(p < alpha_b)
        line += f"{st:>12.1f} p={p:<9.1e}{'*' if p < alpha_b else ' '}"
    print(line)
c3 = np.mean(sig) >= 0.5
print(f"\n  {int(np.sum(sig))} de {len(sig)} rejeitam -> heterogeneidade "
      f"{'ESTABELECIDA' if c3 else 'NAO ESTABELECIDA'}")


section("CHECAGEM 4 — O TERMO DE ESPECIFICACAO EXPLICA O PISO?")
print("  observado: composicao que sobra no resíduo de B com a correcao ajustada\n"
      "  em A cheio (fracao da variancia de aceleracao, medida no painel de 6, menos\n"
      "  o nulo). previsto: (bA-bB)' S_B (bA-bB) - tr(S_B VA) - tr(S_B VB), sobre a\n"
      "  variancia residual de idade em B.\n")
print("  nulos ...", flush=True)
NULL, BEFORE = {}, {}
for t in COHORTS:
    for c in CL:
        NULL[(t, c)] = null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        ib, base = inc_and_base(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        BEFORE[(t, c)] = ((ib - NULL[(t, c)]) / (1 - base), base)
rows = []
for a, b in permutations(COHORTS, 2):
    for c in CL:
        bA, VA = FIT[(a, c)]; bB, VB = FIT[(b, c)]; S = SIG[b]
        d = bA - bB
        spec = float(d @ S @ d - np.trace(S @ VA) - np.trace(S @ VB))
        Xa = np.column_stack([np.ones(N[b]), D[b]["chrono"]])
        rv = float((D[b]["y"][c] - Xa @ np.linalg.lstsq(Xa, D[b]["y"][c], rcond=None)[0]).var())
        cbar = D[a]["C12"].mean(axis=0)
        yc = D[b]["y"][c] - (D[b]["C12"][:, :-1] - cbar[:-1]) @ bA
        ia, _ = inc_and_base(yc, D[b]["chrono"], D[b]["C6"])
        before, base = BEFORE[(b, c)]
        after = (ia - NULL[(b, c)]) / (1 - base)
        rows.append(dict(src=a, dst=b, clock=c, spec=spec / rv, after=after,
                         delta=after - before, before=before))
S4 = pd.DataFrame(rows)
S4.to_csv(RES / "beta_heterogeneity.csv", index=False)
print(f"  {'par':<20}{'relogio':<13}{'previsto':>10}{'sobra obs.':>12}{'delta obs.':>12}")
for r in S4.sort_values("spec", ascending=False).itertuples():
    print(f"  {r.src[3:]+' -> '+r.dst[3:]:<20}{r.clock:<13}{r.spec:>+10.1%}"
          f"{r.after:>+12.1%}{r.delta:>+12.1%}")
rho4, p4 = spearmanr(S4.spec, S4.after)
c4 = rho4 >= 0.5
print(f"\n  Spearman(previsto, sobra) = {rho4:.3f} (p = {p4:.2g}), {len(S4)} configuracoes"
      f" -> {'EXPLICA' if c4 else 'NAO EXPLICA'}")


section("CHECAGEM 5 — DENTRO DA COORTE: DOENCA E TABAGISMO MUDAM O BETA? (REPORTADO)")
for label, t, g1, g0 in (("AR vs controle", "GSE42861", ra == 1, ra == 0),
                         ("fumou alguma vez vs nunca", "GSE50660", smk50 >= 1, smk50 == 0)):
    d = D[t]
    for c in CL:
        f1 = fwl_fit(d["chrono"][g1], d["C12"][g1], d["y"][c][g1])
        f0 = fwl_fit(d["chrono"][g0], d["C12"][g0], d["y"][c][g0])
        st, p = wald(*f1, *f0)
        print(f"  {t} {label:<24}{c:<13} Wald {st:>7.1f}  p = {p:.2g}")


section("CHECAGEM 6 — AJUSTAR POR DOENCA/TABAGISMO REDUZ A HETEROGENEIDADE? (REPORTADO)")
Zs = {"GSE42861": np.column_stack([ra, np.nan_to_num(smk42, nan=1.0)]),
      "GSE50660": np.column_stack([(smk50 == 1).astype(float), (smk50 == 2).astype(float)])}
FITZ = {(t, c): (fwl_fit(D[t]["chrono"], D[t]["C12"], D[t]["y"][c], Zs.get(t))
                 if t in Zs else FIT[(t, c)]) for t in COHORTS for c in CL}
ch = []
for a, b in combinations(COHORTS, 2):
    if not ({a, b} & set(Zs)):
        continue
    for c in CL:
        st0 = W[(a, b, c)]; st1, _ = wald(*FITZ[(a, c)], *FITZ[(b, c)])
        ch.append(st1 - st0)
        print(f"  {a[3:]+'/'+b[3:]:<16}{c:<13} Wald {st0:>7.1f} -> {st1:>7.1f}")
print(f"\n  mudanca mediana no Wald: {np.median(ch):+.1f}")


section("FECHAMENTO")
print(f"  1. cache: {'ok' if c1 else 'FALHOU'}")
print(f"  2. metadados: {'ok' if c2 else 'FALHOU'}")
print(f"  3. heterogeneidade: {'ESTABELECIDA' if c3 else 'NAO ESTABELECIDA'}")
print(f"  4. termo de especificacao explica o piso: {'SIM' if c4 else 'NAO'} (rho {rho4:.3f})")
