#!/usr/bin/env python3
"""
Stage 18 — is it the row count, or the conditioning? And does ridge fix it?

Stage 17 built the curve and broke its own explanation. Harm rises steeply as
the fitting cohort shrinks — +12.2 points at n=40, worse than doing nothing in
89% of draws — but 184 samples of GSE40279 do NOT reproduce the +4.9 that stage
15 measured fitting on GSE61151's own 184. Size is a strong cause and not a
sufficient one, so it is not the cause stage 15 named.

The candidate left over is in stage 12's own numbers, unconnected until now:
median CD8T is 0.036 in GSE40279 and **0.000** in GSE61151. A proportion pinned
at the non-negativity boundary of NNLS carries no information about its own
coefficient at any sample size. If what governs the damage is the CONDITIONING
of the composition matrix rather than the number of rows in it, then n matters
only because more rows usually buy better conditioning — and a badly conditioned
cohort of 184 behaves like a well-conditioned cohort of 40.

That is question (a), and it is worth answering because the diagnostic it
implies is computable by a practitioner on their own single cohort, with no
second cohort to validate against. Sample size is not such a diagnostic: stage
17 shows n is not sufficient.

Question (b) is whether the damage can be removed rather than only predicted.
Meredith et al. (2019) give the mechanism — collinear cell-type predictors,
variance inflation above 100, 83% of coefficients flipping sign in simulation.
Inflated coefficients are exactly what a ridge penalty is for, and the
eight-method benchmark the field uses (Genome Biology 2016) contains no
penalised adjustment at all. So it is both an obvious fix and an untried one.

THE TRAP IN QUESTION (B), and this project has fallen into its twin before.
Ridge shrinks the composition coefficients toward zero, and a correction shrunk
to zero IS no correction. So "ridge removes the harm" is guaranteed at large
enough lambda and means nothing — it is stage 8's degenerate corner again, where
"dilution is a defence" turned out to be a model shrinking itself rather than
resisting anything. A lambda only counts here if it removes the harm at small n
AND keeps the benefit at large n. That two-sided requirement is check 4, written
before the run.

THE DESIGN
  Everything from stage 17 is held: fit cohort GSE40279 subsampled on the same
  grid with the same age-decile stratification, the same three external test
  cohorts, the same cross-panel metric (fit with twelve types, measure the
  residual with six).

  What is new per draw:
    - the conditioning of the composition block, measured AFTER partialling out
      the intercept and chronological age, since that is the system whose
      solution the composition coefficients are: condition number, smallest
      singular value, and the largest variance inflation factor, which is the
      quantity Meredith reports as 113.7
    - the L2 norm of the fitted composition coefficients, which stage 17 did not
      record and which is the thing the harm should track if inflation is the
      mechanism
    - ridge over a scale-free lambda grid, plus the lambda a practitioner would
      actually pick, by leave-one-out cross-validation on the fitting subsample

  Ridge is fitted by Frisch-Waugh-Lovell: partial [1, age] out of both the clock
  age and the composition, then penalise only the composition block. At lambda=0
  this reproduces stage 17's coefficients exactly, which check 3 verifies.
  Columns are standardised before penalising and the coefficients are returned to
  the original scale, because the proportions run from ~0.6 for neutrophils to
  ~0.005 for basophils and an unstandardised penalty would fall almost entirely
  on the rare types.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. the twelve-type panel recovers the known twelve-way proportions of the
     twelve GSE167998 mixtures: r > 0.70, mean absolute error < 0.04. Carried
     from stages 15 and 17. Hard stop.
  2. in-sample floor at full n: the OLS correction drives the composition term
     down in the cohort it was fitted on. Carried. Hard stop.
  3. ridge at lambda = 0 reproduces stage 17's OLS correction to 1e-8 in the
     coefficients and to 1e-6 in the resulting delta. If the two paths disagree,
     the ridge implementation is wrong and question (b) cannot be asked. Hard
     stop.
  4. THE DEGENERATE-CORNER CHECK, which is what keeps question (b) honest. At
     the largest lambda the correction must collapse toward doing nothing —
     delta -> 0 and coefficient norm -> 0 — and that corner is NOT a fix. A
     lambda is reported as a fix only if, at that single lambda:
         (i)  at n = 40 the median delta is <= 0, i.e. no worse than not
              correcting, and strictly better than OLS at the same n, and
         (ii) at n = 656 it retains at least 70% of the OLS benefit.
     If no lambda satisfies both, the answer to question (b) is no, and it is
     reported as no. Shrinking a correction into silence is not a method.
  5. THE DECOUPLING TEST for question (a). Sample size and conditioning are
     correlated by construction, so a raw association between conditioning and
     harm proves nothing. The test is WITHIN each n stratum, where n is
     constant: Spearman correlation between the draw's condition number and its
     delta, computed separately at each n, over the 30 draws x 3 clocks x 3
     cohorts in that stratum.
         conditioning is the mechanism -> positive rho in most strata
         it is not                     -> rho scattered around zero, and the
                                          stage's answer to (a) is that
                                          conditioning adds nothing over n
  6. THE PLACEMENT, and this is the number the stage exists for. Stage 15's
     fitting cohort, GSE61151, has its own composition matrix. Where does its
     conditioning sit against GSE40279's subsamples?
         if it lands beyond the median of GSE40279 draws at n <= 120, then a
         cohort of 184 that is conditioned like a cohort of 40 is exactly what
         the failed anchor needs, and the conditioning hypothesis carries the
         explanation stage 15 got wrong
         if it lands near GSE40279's full-cohort value, conditioning does not
         explain the anchor either, and both mechanisms are out
  7. the declared bias from stage 13 carries: this panel's monocyte channel runs
     high, so no per-type monocyte coefficient is read or reported.

A NOTE ON THE CACHE. Loading five series matrices costs about fifteen minutes
and produces two small things: clock ages and composition estimates. Those are
written to results/cache_*.csv on the first run and reused after. The cache is
derived, gitignored with the rest of results/, and rebuilt by deleting it.

Usage:  .venv/bin/python analysis/18_conditioning.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_extended import TYPES12
from model.clocks import CLOCKS
from model.deconvolution import TYPES as TYPES6

OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
CACHE = OUT / "cache"
RNG = np.random.default_rng(20260923)
N_REPS = 30
GRID = [40, 60, 80, 120, 160, 240, 320, 480]
VERDICT_CLOCKS = ["Horvath2013", "Levine2018", "Horvath2018"]
ALPHAS = [0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]

FIT = "GSE40279"
TESTS = ["GSE61151", "GSE50660", "GSE42861"]


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


def r2(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return 1 - (y - X @ beta).var() / y.var()


def inc_and_base(y, a, C):
    Xa = np.column_stack([np.ones(len(a)), a])
    base = r2(Xa, y)
    return r2(np.column_stack([Xa, C[:, :-1]]), y) - base, base


def null_mean(y, a, C, n_perm=1000, rng=RNG):
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
    """Residualise on [1, age] — the Frisch-Waugh-Lovell half of the fit."""
    Xa = np.column_stack([np.ones(len(age)), age])
    P = Xa @ np.linalg.pinv(Xa)
    return Y - P @ Y


def conditioning(age, C):
    """
    Conditioning of the composition block as the fit actually sees it.

    Measured after partialling out intercept and age, and on the same columns
    the correction uses (the last type is dropped, since proportions sum to one
    and the full set is singular by construction rather than by cohort).
    """
    Ct = partial_out(age, C[:, :-1])
    Ct = Ct - Ct.mean(axis=0)
    sv = np.linalg.svd(Ct, compute_uv=False)
    sv = sv[sv > 0]
    kappa = float(sv[0] / sv[-1]) if len(sv) == Ct.shape[1] else np.inf
    # largest variance inflation factor, the quantity Meredith reports as 113.7
    vif = []
    for j in range(Ct.shape[1]):
        others = np.delete(Ct, j, axis=1)
        rj = r2(np.column_stack([np.ones(len(Ct)), others]), Ct[:, j])
        vif.append(1.0 / max(1e-12, 1.0 - rj))
    return kappa, float(sv[-1] / np.sqrt(len(Ct))), float(np.max(vif))


def ridge_path(age, C, y, alphas):
    """
    Composition coefficients over a scale-free lambda grid.

    lambda = alpha * mean eigenvalue of the standardised cross-product, so the
    grid means the same thing at every n and every cohort. alpha = 0 is OLS.
    """
    Ct = partial_out(age, C[:, :-1])
    yt = partial_out(age, y.reshape(-1, 1)).ravel()
    mu, sd = Ct.mean(axis=0), Ct.std(axis=0) + 1e-12
    Z = (Ct - mu) / sd
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    scale = float((s ** 2).mean())
    Uty = U.T @ yt
    out = {}
    for a in alphas:
        lam = a * scale
        b = Vt.T @ (s * Uty / (s ** 2 + lam))
        out[a] = b / sd
    # the lambda a practitioner would pick: leave-one-out over the grid
    loo = {}
    for a in alphas:
        lam = a * scale
        h = (U * (s ** 2 / (s ** 2 + lam))) @ U.T
        resid = yt - h @ yt
        lev = np.clip(np.diag(h), 0, 1 - 1e-9)
        loo[a] = float(np.mean((resid / (1 - lev)) ** 2))
    return out, min(loo, key=loo.get)


# ------------------------------------------------------------------ data ----
def build_cache():
    from load_geo import read_series_matrix
    from load_extended import read_extended
    from model.deconvolution import PANEL, build_panel, deconvolve
    DATA = ROOT / "reference/data"
    AGE_FIELD = {"GSE40279": "age (y)", "GSE61151": "agebloodtaken",
                 "GSE50660": "age", "GSE42861": "age"}
    from model.clocks import predict

    print("cache ausente — carregando as series (uma vez) ...", flush=True)
    b167, p167 = read_extended(DATA / "GSE167998_matrix_processed.txt.gz",
                               DATA / "BloodExtended_Pheno.csv")
    is_mix = (p167.CellType == "MIX").to_numpy()
    b35, m35 = read_series_matrix(DATA / "GSE35069_series_matrix.txt.gz")
    m35 = m35.set_index("gsm").reindex(b35.columns)
    frac35 = m35["tissue/cell type"]

    def scores_and_cand(betas, labels, types, n_candidate=200):
        sc, cand = {}, set()
        arr = np.asarray(labels)
        for t in types:
            a = betas.loc[:, arr == t]
            c = betas.loc[:, np.isin(arr, [o for o in types if o != t])]
            d = a.mean(axis=1) - c.mean(axis=1)
            sd = np.sqrt((a.var(axis=1, ddof=1) + c.var(axis=1, ddof=1)) / 2) + 0.01
            s = (d / sd).replace([np.inf, -np.inf], np.nan).dropna()
            sc[t] = s
            cand |= set(s.nlargest(n_candidate).index) | set(s.nsmallest(n_candidate).index)
        return sc, cand

    sc12, cand12 = scores_and_cand(b167.loc[:, ~is_mix], p167.CellType[~is_mix], TYPES12)
    lab6 = frac35.map({v: k for k, v in PANEL.items()})
    sc6, cand6 = scores_and_cand(b35.loc[:, lab6.notna().to_numpy()],
                                 lab6.dropna(), TYPES6)
    candidates = sorted(cand12 | cand6)

    kept, ages_all, chrono, indexes = {}, {}, {}, {}
    for tag in [FIT] + TESTS:
        print(f"  {tag} ...", flush=True)
        b, m = read_series_matrix(DATA / f"{tag}_series_matrix.txt.gz")
        m = m.set_index("gsm").reindex(b.columns)
        chrono[tag] = pd.to_numeric(m[AGE_FIELD[tag]], errors="coerce")
        ages_all[tag] = pd.DataFrame({c: predict(b, c)[0] for c in CLOCKS})
        indexes[tag] = b.index
        kept[tag] = b.loc[b.index.intersection(candidates)].copy()
        del b, m

    common = indexes[FIT]
    for tag in TESTS:
        common = common.intersection(indexes[tag])
    common = common.intersection(b167.index).intersection(b35.index)

    def final_panel(sc, cand, betas, labels, types, n_side=50):
        pool = sorted(cand & set(common))
        arr = np.asarray(labels)
        probes = set()
        for t in types:
            s = sc[t].reindex(pool).dropna()
            probes |= set(s.nlargest(n_side).index) | set(s.nsmallest(n_side).index)
        probes = sorted(probes)
        return pd.DataFrame({t: betas.loc[probes, arr == t].mean(axis=1)
                             for t in types}), probes

    ref12, probes12 = final_panel(sc12, cand12, b167.loc[:, ~is_mix],
                                  p167.CellType[~is_mix], TYPES12)
    ref6, probes6 = final_panel(sc6, cand6, b35.loc[:, lab6.notna().to_numpy()],
                                lab6.dropna(), TYPES6)

    known = p167.loc[is_mix, TYPES12].to_numpy()
    est = deconvolve(ref12.to_numpy(), b167.loc[probes12].to_numpy()[:, is_mix])
    panel_r = float(np.corrcoef(known.ravel(), est.ravel())[0, 1])
    panel_mae = float(np.abs(known - est).mean())

    CACHE.mkdir(exist_ok=True)
    for tag in [FIT] + TESTS:
        pd.DataFrame(deconvolve(ref12.to_numpy(), kept[tag].loc[probes12].to_numpy()),
                     index=kept[tag].columns, columns=TYPES12
                     ).to_csv(CACHE / f"{tag}_comp12.csv")
        pd.DataFrame(deconvolve(ref6.to_numpy(), kept[tag].loc[probes6].to_numpy()),
                     index=kept[tag].columns, columns=TYPES6
                     ).to_csv(CACHE / f"{tag}_comp6.csv")
        ages_all[tag].assign(chrono=chrono[tag]).to_csv(CACHE / f"{tag}_ages.csv")
    pd.Series({"panel_r": panel_r, "panel_mae": panel_mae}).to_csv(CACHE / "panel.csv")
    print("  cache escrito em results/cache/", flush=True)


need = [CACHE / f"{t}_{w}.csv" for t in [FIT] + TESTS
        for w in ("comp12", "comp6", "ages")] + [CACHE / "panel.csv"]
if not all(f.exists() for f in need):
    build_cache()
else:
    print("usando results/cache/ (apague o diretorio para reconstruir)")

panel = pd.read_csv(CACHE / "panel.csv", index_col=0).iloc[:, 0]
data = {}
for tag in [FIT] + TESTS:
    a = pd.read_csv(CACHE / f"{tag}_ages.csv", index_col=0)
    k = a.chrono.notna().to_numpy()
    data[tag] = dict(
        chrono=a.chrono[k].to_numpy(),
        y={c: a[c][k].to_numpy() for c in CLOCKS},
        C12=pd.read_csv(CACHE / f"{tag}_comp12.csv", index_col=0)
              .loc[a.index[k], TYPES12].to_numpy(),
        C6=pd.read_csv(CACHE / f"{tag}_comp6.csv", index_col=0)
             .loc[a.index[k], TYPES6].to_numpy())


section("CHECAGEM 1 — O PAINEL DE 12 (CARREGADA DAS ETAPAS 15 E 17)")
c1 = (panel.panel_r > 0.70) and (panel.panel_mae < 0.04)
print(f"  correlacao {panel.panel_r:.3f}  erro absoluto medio {panel.panel_mae:.3f}"
      f"  -> {'ok' if c1 else 'FALHOU'}")
if not c1:
    sys.exit("\n  parando: sem painel validado nao se corrige nada.")
print("  (vies de monocito da etapa 13 continua declarado; nenhum coeficiente "
      "por tipo e lido)")

fit = data[FIT]
N_FIT = len(fit["chrono"])
N_ANCHOR = len(data["GSE61151"]["chrono"])
grid = sorted({n for n in set(GRID) | {N_ANCHOR} if n < N_FIT} | {N_FIT})
print(f"\n  ajuste {FIT} n={N_FIT} | teste " +
      ", ".join(f"{t} n={len(data[t]['chrono'])}" for t in TESTS))
print(f"  grade: {grid}")


section("CHECAGEM 2 — PISO EM CASA COM N CHEIO (OLS)")
c2 = True
for c in CLOCKS:
    y = fit["y"][c]
    coefs, _ = ridge_path(fit["chrono"], fit["C12"], y, [0.0])
    cbar = fit["C12"].mean(axis=0)
    yc = y - (fit["C12"][:, :-1] - cbar[:-1]) @ coefs[0.0]
    before, _ = inc_and_base(y, fit["chrono"], fit["C12"])
    after, _ = inc_and_base(yc, fit["chrono"], fit["C12"])
    ok = after < before
    c2 &= ok
    print(f"  {c:<16}{before:>9.4f} -> {after:>9.4f}  {'ok' if ok else 'FALHOU'}")
print(f"  -> {'ok' if c2 else 'FALHOU'}")
if not c2:
    sys.exit("\n  parando: ajuste quebrado em casa.")


section("CHECAGEM 3 — RIDGE EM LAMBDA=0 REPRODUZ O OLS DA ETAPA 17")


def ols_stage17(y, age, C):
    M = np.column_stack([np.ones(len(age)), age, C[:, :-1]])
    return np.linalg.lstsq(M, y, rcond=None)[0][2:]


dmax_b, dmax_d = 0.0, 0.0
for c in VERDICT_CLOCKS:
    for n in (60, 184, N_FIT):
        idx = (np.arange(N_FIT) if n >= N_FIT
               else stratified_draw(fit["chrono"], n, np.random.default_rng(9)))
        b17 = ols_stage17(fit["y"][c][idx], fit["chrono"][idx], fit["C12"][idx])
        b18, _ = ridge_path(fit["chrono"][idx], fit["C12"][idx], fit["y"][c][idx], [0.0])
        dmax_b = max(dmax_b, float(np.abs(b17 - b18[0.0]).max()))
        t = data["GSE61151"]; cbar = fit["C12"][idx].mean(axis=0)
        d17 = inc_and_base(t["y"][c] - (t["C12"][:, :-1] - cbar[:-1]) @ b17,
                           t["chrono"], t["C6"])[0]
        d18 = inc_and_base(t["y"][c] - (t["C12"][:, :-1] - cbar[:-1]) @ b18[0.0],
                           t["chrono"], t["C6"])[0]
        dmax_d = max(dmax_d, abs(d17 - d18))
c3 = (dmax_b < 1e-8) and (dmax_d < 1e-6)
print(f"  maior diferenca de coeficiente {dmax_b:.2e} (limite 1e-8)")
print(f"  maior diferenca de delta       {dmax_d:.2e} (limite 1e-6)")
print(f"  -> {'ok' if c3 else 'FALHOU'}")
if not c3:
    sys.exit("\n  parando: as duas implementacoes discordam, a pergunta (b) nao pode ser feita.")


# ----------------------------------------------------------------- run ------
section("RODANDO A GRADE — CONDICIONAMENTO, NORMA E RIDGE POR SORTEIO")
NULL = {(tag, c): null_mean(data[tag]["y"][c], data[tag]["chrono"], data[tag]["C6"])
        for tag in TESTS for c in CLOCKS}
BEFORE = {}
for tag in TESTS:
    for c in CLOCKS:
        ib, base = inc_and_base(data[tag]["y"][c], data[tag]["chrono"], data[tag]["C6"])
        BEFORE[(tag, c)] = ((ib - NULL[(tag, c)]) / (1 - base), base)

rows = []
for n in grid:
    reps = 1 if n >= N_FIT else N_REPS
    for r in range(reps):
        idx = (np.arange(N_FIT) if n >= N_FIT
               else stratified_draw(fit["chrono"], n, np.random.default_rng(2000 + 13 * r + n)))
        a_fit, C_fit = fit["chrono"][idx], fit["C12"][idx]
        kappa, min_sv, vif = conditioning(a_fit, C_fit)
        cbar = C_fit.mean(axis=0)
        for c in CLOCKS:
            coefs, a_loo = ridge_path(a_fit, C_fit, fit["y"][c][idx], ALPHAS)
            for alpha, b in coefs.items():
                for tag in TESTS:
                    t = data[tag]
                    before, base = BEFORE[(tag, c)]
                    yc = t["y"][c] - (t["C12"][:, :-1] - cbar[:-1]) @ b
                    ia, _ = inc_and_base(yc, t["chrono"], t["C6"])
                    after = (ia - NULL[(tag, c)]) / (1 - base)
                    rows.append(dict(n=n, rep=r, clock=c, cohort=tag, alpha=alpha,
                                     kappa=kappa, min_sv=min_sv, vif=vif,
                                     norm=float(np.linalg.norm(b)),
                                     alpha_loo=a_loo, before=before, after=after,
                                     delta=after - before,
                                     verdict=c in VERDICT_CLOCKS))
    print(f"  n = {n:>4}: {reps} sorteio(s)", flush=True)

df = pd.DataFrame(rows)
df.to_csv(OUT / "conditioning.csv", index=False)
ols = df[df.alpha == 0.0]
V = ols[ols.verdict]


section("(A) CONDICIONAMENTO CONTRA CONTAGEM DE LINHAS")
print("  condicionamento do bloco de composicao, depois de tirar intercepto e idade\n")
cond = (ols.groupby("n")[["kappa", "vif", "min_sv", "norm"]].median())
print(f"  {'n':>5}{'kappa':>12}{'VIF max':>12}{'menor v.s.':>13}{'norma L2':>12}")
for n, r in cond.iterrows():
    print(f"  {n:>5}{r.kappa:>12.1f}{r.vif:>12.1f}{r.min_sv:>13.5f}{r.norm:>12.1f}")

print("\n  CHECAGEM 5 — teste de desacoplamento, dentro de cada estrato de n")
print("  (n constante dentro da linha, entao a correlacao nao pode ser n disfarcado)\n")
print(f"  {'n':>5}{'rho(kappa, delta)':>22}{'p':>10}{'rho(norma, delta)':>22}{'p':>10}")
rhos = []
for n in grid:
    s = V[V.n == n]
    if s.rep.nunique() < 3:
        print(f"  {n:>5}{'— sorteio unico —':>22}")
        continue
    rk, pk = spearmanr(s.kappa, s.delta)
    rn, pn = spearmanr(s["norm"], s.delta)
    rhos.append(rk)
    print(f"  {n:>5}{rk:>22.3f}{pk:>10.3g}{rn:>22.3f}{pn:>10.3g}")
c5 = (np.array(rhos) > 0).mean() > 0.5
print(f"\n  rho positivo em {int((np.array(rhos) > 0).sum())} de {len(rhos)} estratos"
      f" -> {'ok' if c5 else 'FALHOU'}")
print("  ok aqui significa: condicionamento explica dano ALEM do que n explica.")


section("CHECAGEM 6 — A PLACAGEM DO GSE61151, A COORTE DE AJUSTE DA ETAPA 15")
k61, sv61, vif61 = conditioning(data["GSE61151"]["chrono"], data["GSE61151"]["C12"])
print(f"  GSE61151 (n={N_ANCHOR}, a coorte que a etapa 15 ajustou):")
print(f"    kappa {k61:.1f} | VIF max {vif61:.1f} | menor valor singular {sv61:.5f}")
kap_by_n = {n: float(ols[ols.n == n].kappa.median()) for n in grid}
# the smallest grid n whose typical draw is ALREADY better conditioned than
# GSE61151 — i.e. GSE61151 sits at or below that rung of the ladder
better = [n for n in grid if kap_by_n[n] <= k61]
equiv_n = min(better) if better else None
print(f"\n  onde isso cai entre os sorteios do GSE40279:")
for n in grid:
    mark = "  <- GSE61151 fica aqui" if n == equiv_n else ""
    print(f"    n = {n:>4}: kappa mediano {kap_by_n[n]:>10.1f}{mark}")
c6 = equiv_n is not None and equiv_n <= 120
print(f"\n  GSE61151 e tao mal condicionado quanto um sorteio de n <= "
      f"{equiv_n if equiv_n else 'nenhum na grade'}")
print(f"  regra escrita antes: explica a ancora se esse n for <= 120 "
      f"-> {'ok' if c6 else 'FALHOU'}")
if not c6:
    print("  -> condicionamento tambem NAO explica a ancora da etapa 17.\n"
          "     Os dois mecanismos candidatos caem, e o que resta e diferenca\n"
          "     de coorte que nenhuma das duas grandezas captura. Reportar assim.")


section("(B) RIDGE — CHECAGEM 4, A ESQUINA DEGENERADA")
piv = (df[df.verdict].groupby(["alpha", "n"]).delta.median().unstack())
print(f"  delta mediano por alpha e n (positivo = pior que nao corrigir)\n")
print(f"  {'alpha':>7}" + "".join(f"{n:>9}" for n in grid))
for a in ALPHAS:
    print(f"  {a:>7}" + "".join(f"{piv.loc[a, n]:>+9.1%}" for n in grid))

nrm = df[df.verdict].groupby("alpha")["norm"].median()
print(f"\n  norma L2 mediana dos coeficientes, por alpha:")
print("  " + "  ".join(f"a={a}: {nrm[a]:.1f}" for a in ALPHAS))

ols40 = piv.loc[0.0, 40]
ols_full = piv.loc[0.0, N_FIT]
print(f"\n  referencias do OLS: n=40 {ols40:+.1%} | n={N_FIT} {ols_full:+.1%}")
print(f"  um alpha so conta como conserto se, nele:")
print(f"    (i)  n=40 tiver delta <= 0 e melhor que o OLS ({ols40:+.1%})")
print(f"    (ii) n={N_FIT} guardar >= 70% do ganho do OLS "
      f"(<= {0.7 * ols_full:+.1%})\n")
fixes = []
for a in ALPHAS:
    if a == 0.0:
        continue
    d40, dfull = piv.loc[a, 40], piv.loc[a, N_FIT]
    i = (d40 <= 0) and (d40 < ols40)
    ii = dfull <= 0.7 * ols_full
    if i and ii:
        fixes.append(a)
    print(f"    alpha={a:<5} n=40 {d40:>+8.1%} {'ok ' if i else 'nao'}"
          f"   n={N_FIT} {dfull:>+8.1%} {'ok ' if ii else 'nao'}"
          f"   {'<- CONSERTO' if (i and ii) else ''}")
c4 = len(fixes) > 0
print(f"\n  -> {'ok' if c4 else 'FALHOU'}: "
      + (f"alphas que consertam: {fixes}" if c4
         else "nenhum alpha remove o dano sem matar o beneficio. "
              "A resposta a (b) e NAO."))

loo = df[df.verdict].groupby("n").alpha_loo.median()
print(f"\n  o alpha que a validacao cruzada escolheria, por n:")
print("  " + "  ".join(f"n={n}: {loo[n]:g}" for n in grid))


section("CHECAGENS, FECHAMENTO")
for i, (name, ok) in enumerate([
        ("painel de 12 (carregada)", c1),
        ("piso em casa com n cheio", c2),
        ("ridge em lambda=0 reproduz o OLS da etapa 17", c3),
        ("esquina degenerada: existe alpha que conserta sem matar", c4),
        ("desacoplamento: condicionamento explica alem de n", c5),
        ("placagem do GSE61151 explica a ancora", c6)], 1):
    print(f"  {i}. {name}: {'ok' if ok else 'FALHOU'}")
print("  7. vies de monocito da etapa 13: carregado")
print(f"\n  saida: results/conditioning.csv ({len(df)} linhas)")
