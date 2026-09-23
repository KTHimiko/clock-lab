#!/usr/bin/env python3
"""
The four figures of the paper.

Reads what the stages already saved (results/ncurve.csv from stage 17,
results/all_pairs.csv from stage 21) and recomputes the two things that were
printed but never written to disk: stage 17's permuted-coefficient reference
line, and stage 21's matched-n configurations. Both come off the stage 18 cache,
so nothing here reloads a series matrix.

COLOUR, validated rather than asserted. Every hex is taken unchanged from the
reference palette of the data-viz skill, and each combination was run through
that skill's validator (Node v24.21.0, installed locally for this):

  categorical, figures 1 and 3   #2a78d6 + #eb6834
      all checks pass — CVD ΔE 24.7 (protan), normal-vision ΔE 33.6, both far
      above the 8 and 15 floors; contrast >= 3:1
  ordinal ramp, figure 4         #5598e7 -> #184f95
      all checks pass — monotone lightness, ΔL gap clear, light end 2.91:1
      against the surface, hue spread 3°
  figure 2 uses one hue and carries identity on marker shape, so no categorical
  check applies to it

The light end of figure 4's ramp started at the palette's documented ordinal
floor, #86b6ef (2.06:1). That passes, but the validator's categorical run flags
it at chroma 0.097 and below 3:1, and this figure is going to be printed and
photocopied. Stepping it to #5598e7 buys 2.91:1 for nothing.

A fourth run — all four hues together as one categorical set, all pairs — FAILS,
and it is reported here so nobody re-derives it as a problem. No figure puts
those four on screen as identity: figure 4's two blues are an ordinal ramp,
which is the run above, and figures 1 and 3 use the categorical pair. The
failing combination does not exist in this paper.

TWO LANGUAGES. The data is computed once and rendered twice, into figures/en and
figures/pt, so the English manuscript and its Portuguese version never drift
apart in their numbers — only their labels differ.

GREYSCALE. Journals print in grey and reviewers photocopy. Identity is never
carried by hue alone in any of the four: figure 1 separates by line weight and
dash, figure 2 by marker shape, figures 3 and 4 by marker fill and by lightness
steps that survive desaturation.

Usage:  .venv/bin/python paper/figures.py
"""
import sys
from itertools import permutations, combinations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_extended import TYPES12
from model.clocks import CLOCKS
from model.deconvolution import TYPES as TYPES6

RES = ROOT / "results"; CACHE = RES / "cache"
OUT = ROOT / "paper" / "figures"
RNG = np.random.default_rng(20260927)

BLUE, ORANGE = "#2a78d6", "#eb6834"
BLUE_LIGHT, BLUE_DARK = "#5598e7", "#184f95"
INK, INK_2, MUTED = "#0b0b0b", "#52514e", "#8a8985"
GRID = "#e6e5e2"

CLOCKS3 = ["Horvath2013", "Levine2018", "Horvath2018"]
COHORTS = ["GSE40279", "GSE61151", "GSE50660", "GSE42861"]
SHAPES = {"GSE40279": "o", "GSE61151": "s", "GSE50660": "^", "GSE42861": "D"}
N_REPS = 20
ALPHA_FIX = 3.0

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.8,
    "xtick.color": INK_2, "ytick.color": INK_2,
    "axes.labelcolor": INK, "text.color": INK,
    "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "legend.frameon": False,
})


def style(ax, ylab=None, xlab=None, title=None):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, color=GRID, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    if ylab: ax.set_ylabel(ylab)
    if xlab: ax.set_xlabel(xlab)
    if title: ax.set_title(title, loc="left", pad=8, color=INK)


PCT_LANG = "en"          # set per render pass; the axis formatter has no other
                         # way to know which language it is drawing


def pct(x, _=None):
    return dec(f"{x:+.0%}", PCT_LANG) if x else "0"


def save(fig, name, lang):
    d = OUT / lang; d.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(d / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------- machinery off the cache --
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
    Sf, St = sigma(af, Cf), sigma(at, Ct)
    return float(np.trace(np.linalg.inv(Sf) @ St) / len(af))


def ridge_coefs(age, C, y, alphas):
    Ct = partial_out(age, C[:, :-1])
    yt = partial_out(age, y.reshape(-1, 1)).ravel()
    mu, sd = Ct.mean(axis=0), Ct.std(axis=0) + 1e-12
    Z = (Ct - mu) / sd
    U, s, Vt = np.linalg.svd(Z, full_matrices=False)
    scale = float((s ** 2).mean()); Uty = U.T @ yt
    return {a: (Vt.T @ (s * Uty / (s ** 2 + a * scale))) / sd for a in alphas}


print("carregando o cache ...", flush=True)
need = [CACHE / f"{t}_{w}.csv" for t in COHORTS for w in ("comp12", "comp6", "ages")]
if not all(f.exists() for f in need):
    sys.exit("cache ausente. Rode analysis/18_conditioning.py antes.")
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
NULL, BEFORE = {}, {}
for t in COHORTS:
    for c in CLOCKS3:
        NULL[(t, c)] = null_mean(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        ib, base = inc_and_base(D[t]["y"][c], D[t]["chrono"], D[t]["C6"])
        BEFORE[(t, c)] = ((ib - NULL[(t, c)]) / (1 - base), base)


def transport(src, dst, n, reps, alphas=(0.0,), permute=False):
    s, d = D[src], D[dst]
    rows = []
    for r in range(reps):
        idx = (np.arange(N[src]) if n >= N[src]
               else stratified_draw(s["chrono"], n, np.random.default_rng(700 + 7 * r + n)))
        af, Cf = s["chrono"][idx], s["C12"][idx]
        if permute:                      # stage 17's reference line
            Cf = Cf[RNG.permutation(len(Cf))]
        cbar = Cf.mean(axis=0)
        ti = transport_index(af, Cf, d["chrono"], d["C12"])
        for c in CLOCKS3:
            coefs = ridge_coefs(af, Cf, s["y"][c][idx], list(alphas))
            before, base = BEFORE[(dst, c)]
            for alpha, b in coefs.items():
                yc = d["y"][c] - (d["C12"][:, :-1] - cbar[:-1]) @ b
                ia, _ = inc_and_base(yc, d["chrono"], d["C6"])
                rows.append(dict(src=src, dst=dst, n=len(idx), clock=c, alpha=alpha,
                                 tindex=ti,
                                 delta=(ia - NULL[(dst, c)]) / (1 - base) - before))
    return rows


# ------------------------------------------------------------ the strings ---
T = {
 "en": dict(
  f1_title="Above zero, correcting is worse than not correcting",
  f1_y="composition left in the age residual\n(points of age-acceleration variance)",
  f1_x="samples in the fitting cohort",
  f1_med="median (3 clocks × 3 cohorts), IQR shaded",
  f1_coh="each test cohort, separately",
  f1_ref="reference: coefficients with no information",
  f1_refann="coefficients with no information: {v:+.1%}",
  f1_cross="the median crosses zero\nbetween {a} and {b}",
  f2_title="A quantity computed in advance orders the damage",
  f2_y="composition left, in points",
  f2_floor="below 0.05, none of\nthe 63 was harmful",
  f2_fit="fitted on {c}",
  f2_stats="Spearman ρ = {r:.3f}\n{n} configurations\np = {p:.0e}",
  f3_title="The index is asymmetric, and so is the damage",
  f3_sub="each pair in both directions, fitted at n = min(nA, nB) — called {k} of {m}",
  f3_x="composition left, in points",
  f3_lo="lower-index direction", f3_hi="higher-index direction (predicted worse)",
  f4_title="Penalising the coefficients zeroes all twelve transports",
  f4_sub="{k} of 12 pairs are harmful without a penalty; none with one",
  f4_x="composition left, in points",
  f4_ols="least squares (no penalty)", f4_ridge="ridge, α = 3"),
 "pt": dict(
  f1_title="Acima de zero, corrigir é pior do que não corrigir",
  f1_y="composição restante no resíduo de idade\n(pontos da variância da aceleração)",
  f1_x="amostras na coorte de ajuste",
  f1_med="mediana (3 relógios × 3 coortes), IQR sombreado",
  f1_coh="cada coorte de teste, separada",
  f1_ref="referência: coeficientes sem informação",
  f1_refann="coeficientes sem informação: {v:+.1%}",
  f1_cross="a mediana cruza zero\nentre {a} e {b}",
  f2_title="Uma conta feita de antemão ordena o estrago",
  f2_y="composição restante, em pontos",
  f2_floor="abaixo de 0,05 nenhuma\ndas 63 foi nociva",
  f2_fit="ajustado em {c}",
  f2_stats="Spearman ρ = {r:.3f}\n{n} configurações\np = {p:.0e}",
  f3_title="O índice é assimétrico, e o dano também",
  f3_sub="cada par nas duas direções, ajustado em n = min(nA, nB) — acertou {k} de {m}",
  f3_x="composição restante, em pontos",
  f3_lo="direção de índice menor", f3_hi="direção de índice maior (prevista pior)",
  f4_title="Penalizar os coeficientes zera os doze transportes",
  f4_sub="{k} dos 12 pares são nocivos sem penalidade; nenhum com ela",
  f4_x="composição restante, em pontos",
  f4_ols="mínimos quadrados (sem penalidade)", f4_ridge="ridge, α = 3"),
}
XLAB = {
 "en": (r"transport index    $(\sigma^2/n)\,\mathrm{tr}"
        r"(\Sigma_{\mathrm{fit}}^{-1}\Sigma_{\mathrm{test}})$"),
 "pt": (r"índice de transporte    $(\sigma^2/n)\,\mathrm{tr}"
        r"(\Sigma_{\mathrm{ajuste}}^{-1}\Sigma_{\mathrm{teste}})$"),
}


def dec(txt, lang):
    """Decimal separator follows the language, in the figures as in the text."""
    return txt.replace(".", ",") if lang == "pt" else txt


# ------------------------------------------------------------- the data -----
print("\ncalculando (uma vez para os dois idiomas) ...", flush=True)
cur = pd.read_csv(RES / "ncurve.csv"); V = cur[cur.verdict]
# one draw per cohort is a single realisation of a random quantity: stage 17
# reported +1.8% from four permutations. Twenty draws per cohort instead.
perm = pd.DataFrame(sum((transport("GSE40279", d, N["GSE40279"], 20, permute=True)
                         for d in ["GSE61151", "GSE50660", "GSE42861"]), []))
perm_line = float(perm.delta.median())
agg = V.groupby("n").delta.agg(["median", lambda s: s.quantile(.25),
                                lambda s: s.quantile(.75)])
agg.columns = ["med", "q25", "q75"]
cross = agg[agg.med < 0].index.min(); below = agg[agg.index < cross].index.max()

A = pd.read_csv(RES / "all_pairs.csv")
cfg = (A[A.alpha == 0.0].groupby(["src", "dst", "n"])
         .agg(indice=("tindex", "median"), dano=("delta", "median")).reset_index())
rho, pv = spearmanr(cfg.indice, cfg.dano)

rows3 = []
for a, b in combinations(COHORTS, 2):
    n = min(N[a], N[b])
    ga = pd.DataFrame(transport(a, b, n, N_REPS))
    gb = pd.DataFrame(transport(b, a, n, N_REPS))
    ia, ib = ga.tindex.median(), gb.tindex.median()
    da, db = ga.delta.median(), gb.delta.median()
    hi = (ia, da) if ia > ib else (ib, db)
    lo = (ib, db) if ia > ib else (ia, da)
    rows3.append(dict(par=f"{a[3:]} / {b[3:]}", n=n, hi=hi[1], lo=lo[1],
                      acerto=hi[1] > lo[1]))
R3 = pd.DataFrame(rows3).sort_values("hi")

rows4 = []
for src, dst in permutations(COHORTS, 2):
    n = min(N[src], N[dst])
    g = pd.DataFrame(transport(src, dst, n, N_REPS, alphas=(0.0, ALPHA_FIX)))
    rows4.append(dict(par=f"{src[3:]}→{dst[3:]}", n=n,
                      ols=g[g.alpha == 0.0].delta.median(),
                      ridge=g[g.alpha == ALPHA_FIX].delta.median()))
R4 = pd.DataFrame(rows4).sort_values("ols")
n_bad = int((R4.ols > 0).sum())


# ----------------------------------------------------------- the rendering --
def fig1(t, lang):
    fig, ax = plt.subplots(figsize=(6.6, 4.1))
    style(ax, ylab=t["f1_y"], xlab=t["f1_x"])
    ax.axhline(0, color=INK, linewidth=1.0, zorder=2)
    for dst, dash in (("GSE61151", (1, 1.6)), ("GSE50660", (3, 1.6)),
                      ("GSE42861", (5, 1.6))):
        g = V[V.cohort == dst].groupby("n").delta.median()
        ax.plot(g.index, g.values, color=MUTED, linewidth=1.0, dashes=dash, zorder=3)
    ax.fill_between(agg.index, agg.q25, agg.q75, color=BLUE, alpha=0.16, zorder=4, lw=0)
    ax.plot(agg.index, agg.med, color=BLUE, linewidth=2.0, zorder=6, marker="o",
            markersize=4.5, markerfacecolor="white", markeredgecolor=BLUE,
            markeredgewidth=1.4)
    ax.axhline(perm_line, color=ORANGE, linewidth=1.6, dashes=(4, 2), zorder=5)
    ax.annotate(dec(t["f1_refann"].format(v=perm_line), lang), xy=(300, perm_line),
                xytext=(0, 9), textcoords="offset points", ha="center",
                va="bottom", color=ORANGE, fontsize=8)
    ax.axvspan(below, cross, color=GRID, alpha=0.9, zorder=1)
    ax.annotate(t["f1_cross"].format(a=below, b=cross),
                xy=(np.sqrt(below * cross), 0.10), ha="center", va="bottom",
                fontsize=8, color=INK_2)
    ax.annotate(f"{agg.med.iloc[0]:+.0%}", xy=(agg.index[0], agg.med.iloc[0]),
                xytext=(6, 4), textcoords="offset points", fontsize=9,
                color=BLUE, fontweight="bold")
    ax.set_xscale("log")
    ticks = [40, 80, 160, 320, 656]     # 160 and 184 collide on a log axis
    ax.set_xticks(ticks); ax.set_xticklabels([str(n) for n in ticks])
    ax.minorticks_off()
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(pct))
    ax.set_ylim(-0.08, 0.30)
    ax.legend(handles=[
        Line2D([], [], color=BLUE, lw=2.0, marker="o", markersize=4.5,
               markerfacecolor="white", markeredgecolor=BLUE, label=t["f1_med"]),
        Line2D([], [], color=MUTED, lw=1.0, dashes=(3, 1.6), label=t["f1_coh"]),
        Line2D([], [], color=ORANGE, lw=1.6, dashes=(4, 2), label=t["f1_ref"]),
    ], loc="upper right", handlelength=2.6)
    ax.text(0, 1.06, t["f1_title"], transform=ax.transAxes, fontsize=10, color=INK)
    save(fig, "fig1_curva_n", lang)


def fig2(t, lang):
    fig, ax = plt.subplots(figsize=(6.6, 4.1))
    style(ax, ylab=t["f2_y"], xlab=XLAB[lang])
    ax.axhline(0, color=INK, linewidth=1.0, zorder=3)
    ax.axvline(0.05, color=MUTED, linewidth=1.0, dashes=(4, 2), zorder=2)
    ax.annotate(t["f2_floor"], xy=(0.048, 0.125), ha="right", va="top",
                fontsize=8, color=INK_2)
    size = {"o": 42, "s": 40, "^": 58, "D": 36}   # equal visual area, not equal s
    for src in COHORTS:
        g = cfg[cfg.src == src]
        ax.scatter(g.indice, g.dano, s=size[SHAPES[src]], marker=SHAPES[src],
                   facecolor=BLUE, edgecolor="white", linewidth=0.9, alpha=0.9,
                   zorder=5, label=t["f2_fit"].format(c=src))
    ax.set_xscale("log")
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(pct))
    ax.legend(loc="upper left", handletextpad=0.4)
    ax.text(0.03, 0.56, dec(t["f2_stats"].format(r=rho, n=len(cfg), p=pv), lang),
            transform=ax.transAxes, ha="left", va="top", fontsize=9, color=INK)
    ax.text(0, 1.06, t["f2_title"], transform=ax.transAxes, fontsize=10, color=INK)
    save(fig, "fig2_indice", lang)


def fig3(t, lang):
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    style(ax, xlab=t["f3_x"])
    ax.axvline(0, color=INK, linewidth=1.0, zorder=2)
    ypos = np.arange(len(R3))
    for y, r in zip(ypos, R3.itertuples()):
        ax.plot([r.lo, r.hi], [y, y], color=MUTED, linewidth=1.2, zorder=3)
    ax.scatter(R3.lo, ypos, s=52, marker="o", facecolor="white", edgecolor=BLUE,
               linewidth=1.6, zorder=5, label=t["f3_lo"])
    ax.scatter(R3.hi, ypos, s=52, marker="o", facecolor=ORANGE, edgecolor="white",
               linewidth=1.0, zorder=5, label=t["f3_hi"])
    # the one miss is explained in the caption: every in-figure placement
    # collided with either the data or the legend
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{r.par}   n={r.n}" for r in R3.itertuples()], fontsize=8)
    ax.set_ylim(-0.7, len(R3) - 0.3)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.05))
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(pct))
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2,
              handletextpad=0.4)
    ax.text(0, 1.08, t["f3_title"], transform=ax.transAxes, fontsize=10, color=INK)
    ax.text(0, 1.015, t["f3_sub"].format(k=int(R3.acerto.sum()), m=len(R3)),
            transform=ax.transAxes, fontsize=8, color=INK_2)
    save(fig, "fig3_assimetria", lang)


def fig4(t, lang):
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    style(ax, xlab=t["f4_x"])
    ax.axvline(0, color=INK, linewidth=1.0, zorder=2)
    ypos = np.arange(len(R4))
    for y, r in zip(ypos, R4.itertuples()):
        ax.annotate("", xy=(r.ridge, y), xytext=(r.ols, y),
                    arrowprops=dict(arrowstyle="-|>", color=MUTED, linewidth=1.1,
                                    shrinkA=5, shrinkB=5), zorder=3)
    ax.scatter(R4.ols, ypos, s=52, marker="o", facecolor=BLUE_LIGHT,
               edgecolor="white", linewidth=1.0, zorder=5, label=t["f4_ols"])
    ax.scatter(R4.ridge, ypos, s=52, marker="s", facecolor=BLUE_DARK,
               edgecolor="white", linewidth=1.0, zorder=5, label=t["f4_ridge"])
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{r.par}   n={r.n}" for r in R4.itertuples()], fontsize=8)
    ax.set_ylim(-0.7, len(R4) - 0.3)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.05))
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(pct))
    ax.legend(loc="lower right", handletextpad=0.4)
    ax.text(0, 1.06, t["f4_title"], transform=ax.transAxes, fontsize=10, color=INK)
    ax.text(0, 1.012, t["f4_sub"].format(k=n_bad), transform=ax.transAxes,
            fontsize=8, color=INK_2)
    save(fig, "fig4_conserto", lang)


for lang, t in T.items():
    globals()["PCT_LANG"] = lang
    for f in (fig1, fig2, fig3, fig4):
        f(t, lang)
    print(f"  figures/{lang}/ — quatro figuras", flush=True)

print(f"\n  fig 1: cruzamento entre n={below} e n={cross}")
print(f"  referencia embaralhada: mediana {perm_line:+.2%}, "
      f"IQR {perm.delta.quantile(.25):+.2%} a {perm.delta.quantile(.75):+.2%}, "
      f"{len(perm)} valores")
print(f"  fig 2: rho = {rho:.3f} sobre {len(cfg)} configuracoes")
print(f"  fig 3: {int(R3.acerto.sum())} de {len(R3)} direcoes acertadas")
print(f"  fig 4: {n_bad} de 12 nocivos com OLS, {int((R4.ridge > 0).sum())} com ridge")
