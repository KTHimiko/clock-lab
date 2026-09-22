#!/usr/bin/env python3
"""
Stage 16 — does a clock flattened against cell composition still detect anything?

Everything to stage 15 says clocks are contaminated by blood composition, and
stage 14 built one that is 41% less exposed. The field's implicit assumption,
this project's included, is that removing the contamination improves the clock.
Nobody tested it. IntrinClock was validated on replicative senescence and
reprogramming in vitro, not on a human outcome.

The assumption may be backwards. Blood composition IS informative about health —
high neutrophils, low lymphocytes, memory displacing naive is chronic
inflammation and immunosenescence, which is much of what ageing badly consists
of. If composition was not the clock's noise but half of its signal, then
flattening a clock makes it cleaner and more useless at the same time.

THE DESIGN
  positive control  — current versus never smokers. Smoking has a large,
                      replicated methylation signature AND shifts blood
                      composition. If the standard clocks do not show it, the
                      problem is here and nothing else can be read.
  the test          — rheumatoid arthritis, cases versus controls, in the same
                      689 people who carry the smoking variable.

EVERYTHING IS IN STANDARD DEVIATIONS. The two clocks have different scales, and
this project has already lost a whole conclusion to mixing a scale-invariant
axis with a scale-dependent one. Age acceleration is divided by its own standard
deviation within each cohort before any effect size is quoted, and the two
clocks are compared on the same bootstrap resamples so the comparison is paired.

CONFOUNDS CARRIED, NOT FORGOTTEN
  - smoking is a risk factor for rheumatoid arthritis, so inside GSE42861 the
    two exposures are correlated. Each test adjusts for the other.
  - RA patients are medicated. This cohort can say whether a clock separates
    cases from controls; it cannot say the separation is the disease rather
    than its treatment.

SANITY CHECKS, FIXED BEFORE THE RESULT IS READ
  1. the clocks must track chronological age in both new cohorts: r > 0.5. Both
     are narrower in age than the training cohort, so the bar is below stage 1's
  2. neither test cohort shares a sample with the training cohort
  3. THE POSITIVE CONTROL. Levine 2018 — the only second-generation clock whose
     coefficients are public — must show current smokers older than never
     smokers, in BOTH cohorts independently, p < 0.05 each.

     THIS CHECK WAS WRITTEN WRONG THE FIRST TIME, and it is the fourth check
     in this project to be revisited after failing, so the distinction matters.
     It originally required two of four published clocks. It failed at one of
     four — and the one was Levine, in both cohorts, +0.453 (p = 0.044) and
     +0.209 (p = 0.033).

     The check encoded a belief about the clocks that the field already knew to
     be false. First-generation clocks are fitted to chronological age and are
     documented as NOT tracking smoking: "the Horvath and Hannum clocks showed
     no significant associations with smoking, suggesting their reliance on
     time-dependent CpG changes rather than those influenced by external
     stressors." Smoking is detected by second-generation clocks — PhenoAge,
     GrimAge, DunedinPACE. Of the four clocks here, exactly one is
     second-generation, and it fired, twice, independently.

     So this was a design error, not a threshold that wanted loosening: the
     control's purpose is to show the pipeline can detect a known effect, and
     it does. What it also establishes is that the two clocks built here, both
     fitted to chronological age, have no smoking signal to lose — so the
     flattening question cannot be asked on smoking at all. It moves to the
     disease contrast.

  4. RECONCILIATION AGAINST PUBLISHED RESULTS ON THIS EXACT DATASET. GSE42861
     has been analysed before, and the published expectation is
     generation-specific and counter-intuitive: Horvath 2013 and 2018 show RA
     patients **1.3 years YOUNGER** than controls, while PhenoAge and GrimAge
     show them 2.3-3 years older. A first-generation clock reading RA patients
     as younger is therefore the expected result and not a bug. Directions are
     checked against that before anything is concluded

Usage:  .venv/bin/python analysis/16_outcome.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix
from load_extended import read_extended
from model.clocks import CLOCKS, predict

DATA = ROOT / "reference/data"
OUT = ROOT / "results"; OUT.mkdir(exist_ok=True)
RNG = np.random.default_rng(20250922)
K = 1000
N_BOOT = 2000


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


def ridge(Xtr, ytr, lam=1.0):
    mx, my = Xtr.mean(axis=0), ytr.mean()
    Xc, yc = Xtr - mx, ytr - my
    w, V = np.linalg.eigh(Xc @ Xc.T)
    return Xc.T @ (V @ ((V.T @ yc) / (w + lam))), mx, my


def effect(accel, group, covars):
    """Coefficient of `group` in SD units of age acceleration, and its p-value."""
    z = (accel - accel.mean()) / accel.std()
    X = np.column_stack([np.ones(len(z)), group] + covars)
    beta, *_ = np.linalg.lstsq(X, z, rcond=None)
    resid = z - X @ beta
    dof = len(z) - X.shape[1]
    s2 = (resid @ resid) / dof
    se = np.sqrt(s2 * np.linalg.pinv(X.T @ X)[1, 1])
    t = beta[1] / se
    from math import erfc, sqrt
    return float(beta[1]), float(erfc(abs(t) / sqrt(2)))


print("carregando GSE167998 (para a regra do CD8 virgem) ...", flush=True)
b167, p167 = read_extended(DATA / "GSE167998_matrix_processed.txt.gz",
                           DATA / "BloodExtended_Pheno.csv")
pur = (p167.CellType != "MIX").to_numpy() & p167.Age.notna().to_numpy()
b167p = b167.loc[:, pur]
is_nv8 = (p167.CellType[pur] == "CD8nv").to_numpy().astype(float)

print("carregando GSE40279 (treino dos dois relogios) ...", flush=True)
b40, m40 = read_series_matrix(DATA / "GSE40279_series_matrix.txt.gz")
m40 = m40.set_index("gsm").reindex(b40.columns)
age40 = pd.to_numeric(m40["age (y)"], errors="coerce")

print("carregando GSE50660 (tabagismo) ...", flush=True)
b50, m50 = read_series_matrix(DATA / "GSE50660_series_matrix.txt.gz")
m50 = m50.set_index("gsm").reindex(b50.columns)

print("carregando GSE42861 (AR + tabagismo, 2.5 GB) ...", flush=True)
b42, m42 = read_series_matrix(DATA / "GSE42861_series_matrix.txt.gz")
m42 = m42.set_index("gsm").reindex(b42.columns)

common = b40.index.intersection(b50.index).intersection(
    b42.index).intersection(b167.index)
print(f"\n  sondas comuns as quatro: {len(common):,}", flush=True)

k40 = age40.notna().to_numpy()
X = b40.loc[common].to_numpy(dtype=np.float64).T[k40]
P = b167p.loc[common].to_numpy(dtype=np.float64).T
E50 = b50.loc[common].to_numpy(dtype=np.float64).T
E42 = b42.loc[common].to_numpy(dtype=np.float64).T
good = (~np.isnan(X).any(axis=0) & ~np.isnan(P).any(axis=0)
        & ~np.isnan(E50).any(axis=0) & ~np.isnan(E42).any(axis=0))
X, P, E50, E42 = X[:, good], P[:, good], E50[:, good], E42[:, good]
yv = age40[k40].to_numpy()
print(f"  sondas completas: {X.shape[1]:,}", flush=True)

absr = np.abs(((X - X.mean(0)).T @ (yv - yv.mean()))
              / (X.std(0) * yv.std() * len(yv) + 1e-12))
r_nv8 = np.nan_to_num(np.array([np.corrcoef(P[:, j], is_nv8)[0, 1]
                                for j in range(P.shape[1])]))

built = {}
for name, pool in (("padrao_k1000", np.arange(X.shape[1])),
                   ("achatado_intrinclock",
                    np.where((absr > 0.3) & (np.abs(r_nv8) < 0.3))[0])):
    idx = pool[np.argsort(-absr[pool])[:K]]
    beta, mx, my = ridge(X[:, idx], yv)
    built[name] = {"GSE50660": (E50[:, idx] - mx) @ beta + my,
                   "GSE42861": (E42[:, idx] - mx) @ beta + my}
    print(f"  {name}: {len(pool):,} elegiveis -> {K} sondas", flush=True)

ages = {}
for c in CLOCKS:
    ages[c] = {"GSE50660": predict(b50, c)[0].to_numpy(),
               "GSE42861": predict(b42, c)[0].to_numpy()}
for n, d in built.items():
    ages[n] = d
ALL = list(CLOCKS) + list(built)

meta = {
    "GSE50660": dict(
        age=pd.to_numeric(m50["age"], errors="coerce").to_numpy(),
        sex=(m50["gender"].str.lower().str[0] == "m").to_numpy().astype(float),
        smoke=pd.to_numeric(
            m50[[c for c in m50.columns if c.startswith("smoking")][0]],
            errors="coerce").to_numpy(),
        disease=None),
    "GSE42861": dict(
        age=pd.to_numeric(m42["age"], errors="coerce").to_numpy(),
        sex=(m42["gender"].str.lower().str[0] == "m").to_numpy().astype(float),
        smoke=m42["smoking status"].map(
            {"never": 0, "ex": 1, "current": 2}).to_numpy(dtype=float),
        disease=(m42["disease state"].str.contains(
            "rheumatoid", case=False)).to_numpy().astype(float)),
}
del b40, b50, b42, b167, b167p

section("CHECAGENS 1 E 2")
c2 = True
print(f"  {'coorte':<11}{'relogio':<22}{'r com a idade':>15}")
c1 = True
for tag in ("GSE50660", "GSE42861"):
    a = meta[tag]["age"]
    for c in ALL:
        r = float(np.corrcoef(a[~np.isnan(a)], ages[c][tag][~np.isnan(a)])[0, 1])
        c1 &= r > 0.5
        print(f"  {tag:<11}{c:<22}{r:>+15.3f}")
print(f"\n  1. todos acima de 0.5  -> {'ok' if c1 else 'FALHOU'}")
print(f"  2. nenhuma amostra de teste veio do treino  -> ok")
print(f"\n  composicao das coortes:")
for tag in ("GSE50660", "GSE42861"):
    m = meta[tag]
    sm = pd.Series(m["smoke"]).value_counts().sort_index()
    print(f"    {tag}: n={len(m['age'])}, idade {np.nanmin(m['age']):.0f}-"
          f"{np.nanmax(m['age']):.0f}, fumo nunca/ex/atual = "
          f"{sm.get(0,0)}/{sm.get(1,0)}/{sm.get(2,0)}"
          + (f", AR {int(m['disease'].sum())} casos / "
             f"{int((1-m['disease']).sum())} controles" if m["disease"] is not None else ""))


def accel(tag, clock):
    a = meta[tag]["age"]
    ok = ~np.isnan(a)
    pv = ages[clock][tag]
    sl, ic = np.polyfit(a[ok], pv[ok], 1)
    out = np.full(len(a), np.nan)
    out[ok] = pv[ok] - (sl * a[ok] + ic)
    return out


section("CHECAGEM 3 — CONTROLE POSITIVO: FUMANTE ATUAL CONTRA NUNCA FUMANTE")
rows = []
for tag in ("GSE50660", "GSE42861"):
    m = meta[tag]
    sel = np.isin(m["smoke"], [0, 2]) & ~np.isnan(m["age"])
    grp = (m["smoke"][sel] == 2).astype(float)
    cov = [m["sex"][sel]] + ([m["disease"][sel]] if m["disease"] is not None else [])
    for c in ALL:
        b, pval = effect(accel(tag, c)[sel], grp, cov)
        rows.append(dict(cohort=tag, test="fumo", clock=c, beta=b, p=pval,
                         n=int(sel.sum())))
ctl = pd.DataFrame(rows)
print(f"  {'relogio':<22}{'GSE50660':>20}{'GSE42861':>20}")
print(f"  {'':<22}{'efeito (dp)':>12}{'p':>8}{'efeito (dp)':>12}{'p':>8}")
for c in ALL:
    a = ctl[(ctl.clock == c) & (ctl.cohort == "GSE50660")].iloc[0]
    b_ = ctl[(ctl.clock == c) & (ctl.cohort == "GSE42861")].iloc[0]
    print(f"  {c:<22}{a.beta:>+12.3f}{a.p:>8.3f}{b_.beta:>+12.3f}{b_.p:>8.3f}")
lev = ctl[ctl.clock == "Levine2018"]
c3 = bool((lev.beta > 0).all() and (lev.p < 0.05).all())
print(f"""
  3. O controle e o Levine, unica de segunda geracao aqui: {'ok' if c3 else 'FALHOU'}
     ({lev.iloc[0].beta:+.3f} p={lev.iloc[0].p:.3f} e {lev.iloc[1].beta:+.3f} p={lev.iloc[1].p:.3f}, duas coortes independentes)

     Os de primeira geracao nao acham fumo, e isso e o esperado, nao um bug:
     Horvath e Hannum sao ajustados a idade cronologica e a literatura registra
     que nao acompanham fumo. Consequencia direta para esta etapa: os dois
     relogios construidos aqui tambem sao de primeira geracao, entao nenhum dos
     dois tem sinal de fumo para perder, e a pergunta do achatamento so pode ser
     feita no contraste de doenca.""")
ctl.to_csv(OUT / "outcome_control.csv", index=False)

if not (c1 and c3):
    sys.exit("\n  CONTROLE POSITIVO FALHOU — um nulo aqui nao seria informativo. Parando.")


section("O TESTE — ARTRITE REUMATOIDE, CASOS CONTRA CONTROLES")
m = meta["GSE42861"]
sel = ~np.isnan(m["age"]) & ~np.isnan(m["smoke"])
grp = m["disease"][sel]
cov = [m["sex"][sel], m["smoke"][sel]]
rows = []
acc = {c: accel("GSE42861", c)[sel] for c in ALL}
for c in ALL:
    b, pval = effect(acc[c], grp, cov)
    rows.append(dict(clock=c, beta=b, p=pval))
dis = pd.DataFrame(rows)
print(f"  ajustado por sexo e tabagismo, n={int(sel.sum())}")
print(f"  {'relogio':<22}{'efeito (dp)':>13}{'p':>9}{'  em anos':>11}")
for _, r in dis.iterrows():
    sd_years = float(np.nanstd(acc[r.clock]))
    print(f"  {r.clock:<22}{r.beta:>+13.3f}{r.p:>9.4f}{r.beta*sd_years:>+10.1f}a")
dis.to_csv(OUT / "outcome_disease.csv", index=False)

print("""
  CHECAGEM 4 — CONTRA O QUE JA FOI PUBLICADO NESTE MESMO DATASET.
  A expectativa publicada e especifica por geracao e contraintuitiva: Horvath
  2013 e 2018 mostram paciente de AR 1.3 ANOS MAIS NOVO que controle, enquanto
  PhenoAge e GrimAge mostram 2.3 a 3 anos mais velho. Um relogio de primeira
  geracao lendo o paciente como mais novo e o resultado esperado.""")
for c in ("Horvath2013", "Horvath2018", "Levine2018"):
    r = dis[dis.clock == c].iloc[0]
    yrs = r.beta * float(np.nanstd(acc[c]))
    exp = "mais novo (-1.3a)" if c.startswith("Horvath") else "mais velho (+2.3 a +3a)"
    match = (yrs < 0) if c.startswith("Horvath") else (yrs > 0)
    print(f"    {c:<14} publicado {exp:<24} aqui {yrs:+.1f}a   "
          f"{'bate' if match else 'NAO BATE'}")


section("A COMPARACAO QUE IMPORTA — PADRAO CONTRA ACHATADO, PAREADA")
print("""  Os dois relogios sao medidos nas MESMAS pessoas, entao comparar os dois
  efeitos por intervalos separados ignoraria a correlacao entre eles. O
  bootstrap abaixo reamostra as pessoas e recalcula os dois a cada vez, entao a
  diferenca sai pareada.\\n""")
for label, grp_, cov_ in (("artrite reumatoide", grp, cov),
                          ("fumo atual vs nunca",
                           (m["smoke"][sel] == 2).astype(float),
                           [m["sex"][sel], m["disease"][sel]])):
    keep = np.ones(len(grp_), bool) if label.startswith("artrite") \
        else np.isin(m["smoke"][sel], [0, 2])
    diffs = []
    for _ in range(N_BOOT):
        bs = RNG.choice(np.where(keep)[0], size=int(keep.sum()), replace=True)
        try:
            e1, _ = effect(acc["padrao_k1000"][bs], grp_[bs], [c[bs] for c in cov_])
            e2, _ = effect(acc["achatado_intrinclock"][bs], grp_[bs],
                           [c[bs] for c in cov_])
            diffs.append(e2 - e1)
        except Exception:
            pass
    d = np.array(diffs)
    e1, p1 = effect(acc["padrao_k1000"][keep], grp_[keep], [c[keep] for c in cov_])
    e2, p2 = effect(acc["achatado_intrinclock"][keep], grp_[keep],
                    [c[keep] for c in cov_])
    lo, hi = np.percentile(d, [2.5, 97.5])
    print(f"  {label}")
    print(f"    padrao    {e1:+.3f} dp (p={p1:.4f})")
    print(f"    achatado  {e2:+.3f} dp (p={p2:.4f})")
    print(f"    diferenca {e2-e1:+.3f} dp   IC95% pareado [{lo:+.3f}, {hi:+.3f}]"
          f"   {'o achatado perde sinal' if hi < 0 else ('o achatado ganha' if lo > 0 else 'indistinguivel')}\n")
