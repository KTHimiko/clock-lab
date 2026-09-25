#!/usr/bin/env python3
"""
Stage 46b — the saliva result, measured with an independent panel.

Stages 38-43 measured the saliva leftover with the first step of the EpiDISH
reference the correction was fitted on. That biases the measurement toward the
correction, so every saliva harm figure so far is a floor. 46a built a panel from
purified saliva fractions (GSE147318), sharing 3.3% of its probes with EpiDISH,
restricted to probes both arrays carry, and separating its own sorted fractions
without overlap. It is used as a relative immune score, not a proportion: the
absolute scale does not transfer between studies (46a), but the ordering does.

This re-runs the two saliva results that matter, unchanged except for the
measurement: the fit is still the nine-type EpiDISH hierarchy, the transports are
the same four directed pairs with GSE78874, the clocks are Levine 2018 and
Horvath 2018. The measured quantity is now the immune fraction from the
independent panel, and the composition term is scored on it.

PRE-REGISTERED (written before any independent-panel score was computed)
  1. the panel loads, every saliva cohort carries >= 90% of its probes, and the
     score correlates above 0.7 with the EpiDISH immune fraction — the two must
     be measuring the same axis for the comparison to mean anything. Hard stop.
  2. THE HARM SURVIVES AN INDEPENDENT MEASUREMENT: at matched n, at least 6 of 8
     (pair x clock) cells are harmful unpenalised. Stages 38 found 7 of 8 with
     the shared measurement.
  3. THE PENALTY STILL FAILS HERE: at alpha = 3, at least 4 of 8 cells remain
     harmful. Stage 38 found 7 of 8.
  Reported without a bar: the before-correction composition signal under the two
  measurements, which says how much the shared panel was flattering the
  correction, and the n = 40 draws.

Usage:  .venv/bin/python analysis/46b_saliva_independent_score.py
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
RES = ROOT / "results"; CACHE = RES / "cache"


from model.deconvolution import deconvolve

SAL = ["GSE232891", "GSE232332", "GSE78874"]
FIT_S = ["Epi", "Fib", "B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino"]
AGE = ["Levine2018", "Horvath2018"]
DATA = ROOT / "reference/data"


def section(t):
    print(f"\n{'='*74}\n{t}\n{'='*74}", flush=True)


REF = pd.read_csv(CACHE / "panel_saliva_indep.csv", index_col=0)
NEED = set(REF.index)
section("CHECAGEM 1 — COBERTURA DO PAINEL INDEPENDENTE")
OUT = CACHE / "saliva_ic_indep.csv"
if OUT.exists():
    IC = pd.read_csv(OUT, index_col=0).iloc[:, 0]
    print("  usando o cache existente")
else:
    def rows(path, **kw):
        parts = []
        for ch in pd.read_csv(path, index_col=0, chunksize=40_000, **kw):
            parts.append(ch.loc[ch.index.intersection(NEED)])
        return pd.concat(parts)

    pieces = []
    for tag in SAL + ["GSE149747"]:
        if tag == "GSE232891":
            hdr = pd.read_csv(DATA / f"{tag}_Processed_Beta_Values.csv.gz", nrows=0).columns
            b = rows(DATA / f"{tag}_Processed_Beta_Values.csv.gz",
                     usecols=[hdr[0]] + [c for c in hdr[1:] if not c.endswith("_detP")])
        elif tag == "GSE232332":
            import ast, gzip
            with gzip.open(DATA / f"{tag}_Matrix_Processed.txt.gz", "rt") as fh:
                head = ast.literal_eval(fh.readline().strip().strip('"')); first = fh.readline()
            nf = len(first.rstrip("\n").split(","))
            names = ["probe"] + (head[1:] if nf == len(head) else head[2:])
            b = rows(DATA / f"{tag}_Matrix_Processed.txt.gz", header=None, names=names, skiprows=1)
            b = b[[c for c in b.columns if not str(c).endswith("_detPValue") and "Replicate" not in str(c) and str(c).strip()]]
        elif tag == "GSE78874":
            raw = rows(DATA / f"{tag}_datSignal.csv.gz")
            ids = sorted({c[:-len("Methylated")] for c in raw.columns
                          if c.endswith("Methylated") and not c.endswith("Unmethylated")})
            M = raw[[i + "Methylated" for i in ids]].to_numpy(float)
            U = raw[[i + "Unmethylated" for i in ids]].to_numpy(float)
            b = pd.DataFrame(M / (M + U + 100), index=raw.index, columns=ids)
        else:
            b = rows(DATA / "GSE149747_MDL_Matrix_AverageBetas.csv.gz", encoding="utf-8-sig", low_memory=False)
            b = b[[c for c in b.columns if not c.startswith("Detection")]]
        b = b.apply(pd.to_numeric, errors="coerce")
        p = [x for x in REF.index if x in b.index]
        cov = len(p) / len(REF)
        sub = b.loc[p].dropna(axis=1, how="all")
        R = REF.loc[p]
        w = (R.IC - R.Epi).to_numpy(); w = w / np.linalg.norm(w)
        Z = sub.to_numpy(); Z = Z - np.nanmean(Z, axis=1, keepdims=True)
        sc = pd.Series(np.nan_to_num(Z).T @ w, index=sub.columns, name="IC_indep")
        pieces.append(sc)
        print(f"  {tag:<11}cobertura {cov:.1%}  escore dp {sc.std():.2f}"
              f" -> {'ok' if cov >= 0.90 else 'FALHOU'}", flush=True)
        if cov < 0.90:
            sys.exit("  parando.")
        del b
    IC = pd.concat(pieces)
    IC.to_frame().to_csv(OUT)

D = {}
for t in SAL:
    a = pd.read_csv(CACHE / f"{t}_ages.csv", index_col=0)
    ic = IC.reindex(a.index)
    D[t] = dict(chrono=a.chrono.to_numpy(float), y={c: a[c].to_numpy(float) for c in AGE},
                C=pd.read_csv(CACHE / f"{t}_compfit.csv", index_col=0).loc[a.index, FIT_S].to_numpy(),
                M_sh=pd.read_csv(CACHE / f"{t}_compmeas.csv", index_col=0).loc[a.index, ["Epi", "Fib", "IC"]].to_numpy(),
                M_in=ic.to_numpy().reshape(-1, 1))
    N = len(a)
    r = float(np.corrcoef(ic.to_numpy(), pd.read_csv(CACHE / f"{t}_compmeas.csv", index_col=0).loc[a.index, "IC"])[0, 1])
    print(f"  {t}: escore independente x fracao do EpiDISH, r = {r:.3f}"
          f" -> {'ok' if abs(r) > 0.7 else 'FALHOU'}")
    if abs(r) <= 0.7:
        sys.exit("  parando: os dois paineis nao medem o mesmo eixo.")


def r2(X, y):
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return 1 - (y - X @ beta).var() / y.var()


def inc_and_base(y, a, C):
    Xa = np.column_stack([np.ones(len(a)), a])
    base = r2(Xa, y)
    return r2(np.column_stack([Xa, C]), y) - base, base


def null_mean(y, a, C, n_perm=500, seed=46):
    rng = np.random.default_rng(seed)
    Xa = np.column_stack([np.ones(len(a)), a]); base = r2(Xa, y)
    return float(np.mean([r2(np.column_stack([Xa, C[rng.permutation(len(C))]]), y) - base
                          for _ in range(n_perm)]))


def partial_out(age, Y):
    Xa = np.column_stack([np.ones(len(age)), age])
    return Y - Xa @ np.linalg.pinv(Xa) @ Y


def ridge_coefs(age, C, y, alphas):
    Ct = partial_out(age, C[:, :-1]); yt = partial_out(age, y.reshape(-1, 1)).ravel()
    mu, sd = Ct.mean(0), Ct.std(0) + 1e-12
    Z = (Ct - mu) / sd
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    scale = float((s ** 2).mean()); Uty = U.T @ yt
    keep = s > np.finfo(float).eps * max(len(age), Z.shape[1]) * s.max()
    return {a: (Vt.T @ np.where(keep, s * Uty / (s ** 2 + a * scale), 0.0)) / sd for a in alphas}


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


NUL, BEF = {}, {}
section("SINAL ANTES DA CORRECAO, NAS DUAS MEDICOES")
print(f"  {'coorte':<11}{'relogio':<13}{'compartilhada':>15}{'independente':>15}")
for t in SAL:
    d = D[t]
    for c in AGE:
        for key, M in (("sh", d["M_sh"]), ("in", d["M_in"])):
            nul = null_mean(d["y"][c], d["chrono"], M)
            ib, base = inc_and_base(d["y"][c], d["chrono"], M)
            NUL[(t, c, key)] = nul; BEF[(t, c, key)] = ((ib - nul) / (1 - base), base)
        print(f"  {t:<11}{c:<13}{BEF[(t, c, 'sh')][0]:>+15.1%}{BEF[(t, c, 'in')][0]:>+15.1%}")


def after(dst, c, b, cbar, key):
    d = D[dst]; M = d["M_sh"] if key == "sh" else d["M_in"]
    yc = d["y"][c] - (d["C"][:, :-1] - cbar[:-1]) @ b
    ia, _ = inc_and_base(yc, d["chrono"], M)
    return (ia - NUL[(dst, c, key)]) / (1 - BEF[(dst, c, key)][1])


PAIRS = [("GSE232891", "GSE78874"), ("GSE78874", "GSE232891"),
         ("GSE232332", "GSE78874"), ("GSE78874", "GSE232332")]
section("CHECAGEM 2 E 3 — TRANSPORTE EM N PAREADO, MEDIDO DE FORMA INDEPENDENTE")
cells = []
for src, dst in PAIRS:
    s = D[src]; NS = len(s["chrono"]); n = min(NS, len(D[dst]["chrono"]))
    for c in AGE:
        res = {k: [] for k in ("ols_sh", "a3_sh", "ols_in", "a3_in")}
        for r in range(20):
            idx = (np.arange(NS) if n >= NS
                   else stratified_draw(s["chrono"], n, np.random.default_rng(4600 + 7 * r)))
            co = ridge_coefs(s["chrono"][idx], s["C"][idx], s["y"][c][idx], [0.0, 3.0])
            cb = s["C"][idx].mean(0)
            for key in ("sh", "in"):
                for lab, al in (("ols", 0.0), ("a3", 3.0)):
                    res[f"{lab}_{key}"].append(after(dst, c, co[al], cb, key) - BEF[(dst, c, key)][0])
            if n >= NS:
                break
        cells.append(dict(src=src, dst=dst, n=n, clock=c,
                          **{k: float(np.median(v)) for k, v in res.items()}))
P = pd.DataFrame(cells); P.to_csv(RES / "saliva_independent.csv", index=False)
print(f"  {'par':<20}{'relogio':<13}{'OLS compart.':>14}{'a=3 compart.':>14}{'OLS indep.':>12}{'a=3 indep.':>12}")
for r in P.itertuples():
    print(f"  {r.src[3:]+' -> '+r.dst[3:]:<20}{r.clock:<13}{r.ols_sh:>+14.1%}{r.a3_sh:>+14.1%}"
          f"{r.ols_in:>+12.1%}{r.a3_in:>+12.1%}")
k_o, k_a = int((P.ols_in > 0).sum()), int((P.a3_in > 0).sum())
c2, c3 = k_o >= 6, k_a >= 4
print(f"\n  nocivas medidas de forma independente: sem penalidade {k_o} de 8, alpha=3 {k_a} de 8")
print(f"  (medicao compartilhada: {int((P.ols_sh > 0).sum())} e {int((P.a3_sh > 0).sum())})")
section("FECHAMENTO")
print(f"  2. o dano sobrevive a medicao independente (barra >= 6): {k_o} -> {'ok' if c2 else 'FALHOU'}")
print(f"  3. a penalidade ainda falha (barra >= 4): {k_a} -> {'ok' if c3 else 'FALHOU'}")
