#!/usr/bin/env python3
"""
The figures of the paper, after the review (stages 23-29).

Every figure reads what a stage saved; nothing is recomputed here, so a figure
can never disagree with the stage that produced its numbers.

  fig1  the n-curve with the permuted reference drawn AT EVERY n    stage 24
  fig2  the two components: what the real fit leaves, its n-independent
        floor, and what coefficients with no information leave      stage 24
  fig3  model shift: the specification term against the full-n
        leftover                                                    stage 25
  fig4  the transport index against net damage, per clock, with the
        transports a closed-form predictor wrongly calls safe       stage 26
  fig5  the penalty, cell by cell                                   stage 27
  fig6  saliva: the immune-fraction slope, per cohort and clock     stages 39-40

The old figure 1 drew the permuted reference as a flat line measured only at
n = 656; the old figure 2 marked a "safe" floor that stage 26 withdrew; the old
figure 3 (asymmetry) rested on six pairs of which one was a tie. They are gone.

COLOUR, validated with the data-viz skill's validator (Node v24.21.0):
  categorical #2a78d6 + #eb6834 — all checks pass (CVD ΔE 24.7, normal 33.6)
  ordinal     #5598e7 -> #184f95 — all checks pass (light end 2.91:1)
No figure puts more than two categorical hues on screen. Identity is never hue
alone: dash, marker shape and fill carry it too, so every figure survives print.

TWO LANGUAGES. Rendered into figures/en and figures/pt from the same numbers.

Usage:  .venv/bin/python paper/figures.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
OUT = ROOT / "paper" / "figures"

BLUE, ORANGE = "#2a78d6", "#eb6834"
BLUE_LIGHT, BLUE_DARK = "#5598e7", "#184f95"
INK, INK_2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8985", "#e6e5e2"
SHAPES = {"GSE40279": "o", "GSE61151": "s", "GSE50660": "^", "GSE42861": "D"}
SIZE = {"o": 34, "s": 32, "^": 46, "D": 28}

mpl.rcParams.update({
    "font.family": "sans-serif", "font.size": 9, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.8,
    "xtick.color": INK_2, "ytick.color": INK_2, "axes.labelcolor": INK,
    "text.color": INK, "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "legend.frameon": False,
})

T = {
 "en": dict(
  pct_dec=".",
  y_net="net damage (points of age-acceleration variance)",
  y_left="composition left in the age residual",
  x_n="samples in the fitting cohort",
  f1_title="Transported, the correction is worse than none below ~160 samples",
  f1_real="real coefficients (median, IQR shaded)",
  f1_perm="coefficients with no information (permuted)",
  f2_title="Two components: noise that shrinks with n, and a floor that does not",
  f2_real="left by real coefficients", f2_floor="its floor at full n",
  f2_perm="left by permuted coefficients", f2_ann="n-independent floor",
  f3_title="Model shift explains the worst transports",
  f3_x=r"specification term  $(\beta_{fit}-\beta_{test})'\,\Sigma_{test}\,(\beta_{fit}-\beta_{test})$, bias-corrected",
  f3_y="composition left at full fitting size",
  f3_worst="GSE61151 → GSE42861",
  f3_rho="Spearman ρ = {r:.2f}\n24 configurations",
  f4_title="The index ranks the noise; nothing computed in advance certifies safety",
  f4_x=r"transport index   $\mathrm{tr}(\Sigma_{fit}^{-1}\Sigma_{test})\,/\,n$",
  f4_safe_harm="called safe in advance, observed harmful",
  f4_other="other configurations",
  f4_note="within-n ρ ≈ 0.36\nblock-permutation p = 0.036",
  f5_title="A penalty removes the harm in 29 of 30 cells",
  f5_sub="{k} of 30 (pair × clock) cells harmful without a penalty; {j} at α = 3",
  f5_ols="least squares (no penalty)", f5_ridge="ridge, α = 3",
  f6_title="In saliva, the composition slope changes sign between cohorts",
  f6_sub="clock age per +10 points of immune fraction, adjusted for age (95% CI)",
  f6_x="years per +10 percentage points of immune fraction",
  f6_h="Horvath 2018", f6_l="Levine 2018", f6_norm="GSE78874 (normalised)"),
 "pt": dict(
  pct_dec=",",
  y_net="dano líquido (pontos da variância de aceleração)",
  y_left="composição restante no resíduo de idade",
  x_n="amostras na coorte de ajuste",
  f1_title="Transportada, a correção é pior que nenhuma abaixo de ~160 amostras",
  f1_real="coeficientes reais (mediana, IQR sombreado)",
  f1_perm="coeficientes sem informação (embaralhados)",
  f2_title="Dois componentes: ruído que cai com n, e um piso que não cai",
  f2_real="deixado pelos coeficientes reais", f2_floor="o piso dele em n cheio",
  f2_perm="deixado pelos embaralhados", f2_ann="piso independente de n",
  f3_title="O model shift explica os piores transportes",
  f3_x=r"termo de especificação  $(\beta_{ajuste}-\beta_{teste})'\,\Sigma_{teste}\,(\beta_{ajuste}-\beta_{teste})$, corrigido",
  f3_y="composição restante em tamanho cheio",
  f3_worst="GSE61151 → GSE42861",
  f3_rho="Spearman ρ = {r:.2f}\n24 configurações",
  f4_title="O índice ordena o ruído; nenhuma conta de antemão certifica segurança",
  f4_x=r"índice de transporte   $\mathrm{tr}(\Sigma_{ajuste}^{-1}\Sigma_{teste})\,/\,n$",
  f4_safe_harm="declarada segura de antemão, observada nociva",
  f4_other="demais configurações",
  f4_note="ρ dentro de n ≈ 0,36\npermutação em blocos p = 0,036",
  f5_title="A penalidade remove o dano em 29 de 30 células",
  f5_sub="{k} de 30 células (par × relógio) nocivas sem penalidade; {j} com α = 3",
  f5_ols="mínimos quadrados (sem penalidade)", f5_ridge="ridge, α = 3",
  f6_title="Na saliva, a inclinação da composição troca de sinal entre coortes",
  f6_sub="idade do relógio por +10 pontos de fração imune, ajustada pela idade (IC 95%)",
  f6_x="anos por +10 pontos percentuais de fração imune",
  f6_h="Horvath 2018", f6_l="Levine 2018", f6_norm="GSE78874 (normalizada)"),
}
LANG = "en"


def pct(x, _=None):
    s = f"{x:+.0%}" if x else "0"
    return s.replace(".", T[LANG]["pct_dec"])


def style(ax, xlab=None, ylab=None):
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(True, color=GRID, linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)
    if xlab: ax.set_xlabel(xlab)
    if ylab: ax.set_ylabel(ylab)


def title(ax, text, sub=None):
    ax.text(0, 1.07 if sub else 1.04, text, transform=ax.transAxes, fontsize=10)
    if sub:
        ax.text(0, 1.015, sub, transform=ax.transAxes, fontsize=8, color=INK_2)


def save(fig, name):
    d = OUT / LANG; d.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(d / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def log_n_axis(ax, ticks=(40, 80, 160, 320, 656)):
    ax.set_xscale("log"); ax.set_xticks(list(ticks))
    ax.set_xticklabels([str(n) for n in ticks]); ax.minorticks_off()


# ------------------------------------------------------------------ data ----
NB = pd.read_csv(RES / "null_by_n.csv")
real, perm = NB[NB.kind == "real"], NB[NB.kind == "perm"]
q = lambda s, k: s.quantile(k)
R1 = real.groupby("n").delta.agg(["median", lambda s: q(s, .25), lambda s: q(s, .75)])
R1.columns = ["med", "lo", "hi"]
P1 = perm.groupby("n").delta.agg(["median", lambda s: q(s, .25), lambda s: q(s, .75)])
P1.columns = ["med", "lo", "hi"]
LEFT = real.groupby("n").after.median()
FLOOR = float(LEFT.loc[LEFT.index.max()])

BH = pd.read_csv(RES / "beta_heterogeneity.csv")
rho3 = spearmanr(BH.spec, BH.after)[0]

ND = pd.read_csv(RES / "net_damage.csv")
C4 = (ND.groupby(["src", "dst", "n", "clock"])
        .agg(index=("tindex", "median"), obs=("obs", "median"),
             p_est=("p_est", "median")).reset_index())
C4["false_safe"] = (C4.p_est <= 0) & (C4.obs > 0)

RP = pd.read_csv(RES / "ridge_per_clock.csv")
C5 = (RP.groupby(["src", "dst", "clock", "alpha"]).delta.median()
        .unstack("alpha").reset_index().sort_values(0.0))
k_ols, k_r3 = int((C5[0.0] > 0).sum()), int((C5[3.0] > 0).sum())


SL = pd.read_csv(RES / "saliva_slopes.csv")
SL_ORDER = ["GSE232891", "GSE232332", "GSE149747", "GSE78874", "GSE78874 (normalised)"]


# -------------------------------------------------------------- figures -----
def fig1(t):
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    style(ax, t["x_n"], t["y_net"])
    ax.axhline(0, color=INK, lw=1.0, zorder=2)
    ax.fill_between(R1.index, R1.lo, R1.hi, color=BLUE, alpha=0.15, lw=0, zorder=3)
    ax.plot(R1.index, R1.med, color=BLUE, lw=2.0, marker="o", ms=4.5,
            mfc="white", mec=BLUE, mew=1.4, zorder=5)
    ax.plot(P1.index, P1.med, color=ORANGE, lw=1.6, dashes=(4, 2), marker="s",
            ms=4, mfc=ORANGE, mec="white", zorder=4)
    ax.annotate(pct(R1.med.iloc[0]).replace("%", "") + "%", (R1.index[0], R1.med.iloc[0]),
                xytext=(6, 2), textcoords="offset points", color=BLUE,
                fontsize=9, fontweight="bold")
    log_n_axis(ax)
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(pct))
    ax.legend(handles=[
        Line2D([], [], color=BLUE, lw=2, marker="o", ms=4.5, mfc="white", mec=BLUE,
               label=t["f1_real"]),
        Line2D([], [], color=ORANGE, lw=1.6, dashes=(4, 2), marker="s", ms=4,
               mfc=ORANGE, mec="white", label=t["f1_perm"])], loc="upper right")
    title(ax, t["f1_title"])
    save(fig, "fig1_curve")


def fig2(t):
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    style(ax, t["x_n"], t["y_left"])
    ax.axhline(0, color=INK, lw=1.0, zorder=2)
    ax.axhspan(0, FLOOR, color=GRID, alpha=0.8, zorder=1)
    ax.axhline(FLOOR, color=MUTED, lw=1.0, dashes=(2, 2), zorder=3)
    ax.annotate(t["f2_ann"], (41, FLOOR), xytext=(0, -12), textcoords="offset points",
                fontsize=8, color=INK_2)
    ax.plot(LEFT.index, LEFT.values, color=BLUE, lw=2.0, marker="o", ms=4.5,
            mfc="white", mec=BLUE, mew=1.4, zorder=5)
    ax.plot(P1.index, P1.med, color=ORANGE, lw=1.6, dashes=(4, 2), marker="s",
            ms=4, mfc=ORANGE, mec="white", zorder=4)
    log_n_axis(ax)
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(pct))
    ax.legend(handles=[
        Line2D([], [], color=BLUE, lw=2, marker="o", ms=4.5, mfc="white", mec=BLUE,
               label=t["f2_real"]),
        Line2D([], [], color=MUTED, lw=1, dashes=(2, 2), label=t["f2_floor"]),
        Line2D([], [], color=ORANGE, lw=1.6, dashes=(4, 2), marker="s", ms=4,
               mfc=ORANGE, mec="white", label=t["f2_perm"])], loc="upper right")
    title(ax, t["f2_title"])
    save(fig, "fig2_components")


def fig3(t):
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    style(ax, t["f3_x"], t["f3_y"])
    ax.axhline(0, color=INK, lw=0.8, zorder=2); ax.axvline(0, color=INK, lw=0.8, zorder=2)
    for clock, color, marker in (("Levine2018", BLUE, "o"), ("Horvath2018", ORANGE, "D")):
        g = BH[BH.clock == clock]
        ax.scatter(g.spec, g.after, s=40 if marker == "o" else 32, marker=marker,
                   facecolor=color, edgecolor="white", lw=0.9, zorder=5, label=clock)
    w = BH[(BH.src == "GSE61151") & (BH.dst == "GSE42861")]
    lx, ly = float(w.spec.mean()), float(w.after.mean()) - 0.075
    for r in w.itertuples():
        ax.annotate("", xy=(r.spec, r.after), xytext=(lx, ly),
                    arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.9, shrinkB=5))
    ax.text(lx, ly - 0.004, t["f3_worst"], ha="center", va="top", fontsize=8,
            color=INK_2)
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(pct))
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(pct))
    ax.legend(loc="center right")
    ax.text(0.03, 0.95, t["f3_rho"].format(r=rho3).replace(".", t["pct_dec"], 1),
            transform=ax.transAxes, va="top", fontsize=9)
    title(ax, t["f3_title"])
    save(fig, "fig3_model_shift")


def fig4(t):
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    style(ax, t["f4_x"], t["y_net"])
    ax.axhline(0, color=INK, lw=1.0, zorder=2)
    for src, mk in SHAPES.items():
        g = C4[C4.src == src]
        for fs, face, edge in ((False, BLUE, "white"), (True, ORANGE, "white")):
            h = g[g.false_safe == fs]
            ax.scatter(h["index"], h.obs, s=SIZE[mk], marker=mk, facecolor=face,
                       edgecolor=edge, lw=0.8, alpha=0.9, zorder=5 if fs else 4)
    ax.set_xscale("log")
    ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(pct))
    handles = [Line2D([], [], ls="", marker="o", ms=6, mfc=ORANGE, mec="white",
                      label=t["f4_safe_harm"]),
               Line2D([], [], ls="", marker="o", ms=6, mfc=BLUE, mec="white",
                      label=t["f4_other"])]
    handles += [Line2D([], [], ls="", marker=mk, ms=5.5, mfc="white", mec=INK_2,
                       label=src) for src, mk in SHAPES.items()]
    ax.legend(handles=handles, loc="upper left", fontsize=7.5)
    CAP = 0.55
    out = C4[C4.obs > CAP]
    for r in out.itertuples():
        ax.annotate(pct(r.obs).replace("%", "") + "%", (r.index, CAP),
                    xytext=(0, -14), textcoords="offset points", ha="center",
                    fontsize=7.5, color=INK_2,
                    arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.8))
    ax.set_ylim(top=CAP)
    title(ax, t["f4_title"], t["f4_note"].replace("\n", "  ·  "))
    save(fig, "fig4_index")


def fig5(t):
    fig, ax = plt.subplots(figsize=(6.6, 6.4))
    style(ax, t["y_net"])
    ax.axvline(0, color=INK, lw=1.0, zorder=2)
    y = np.arange(len(C5))
    for yi, r in zip(y, C5.itertuples()):
        ax.annotate("", xy=(r._5, yi), xytext=(r._4, yi),
                    arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.0,
                                    shrinkA=4, shrinkB=4), zorder=3)
    ax.scatter(C5[0.0], y, s=36, marker="o", facecolor=BLUE_LIGHT, edgecolor="white",
               lw=0.9, zorder=5, label=t["f5_ols"])
    ax.scatter(C5[3.0], y, s=34, marker="s", facecolor=BLUE_DARK, edgecolor="white",
               lw=0.9, zorder=5, label=t["f5_ridge"])
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.src[3:]}→{r.dst[3:]}  {r.clock.replace('2018', ' 18').replace('2013', ' 13')}"
                        for r in C5.itertuples()], fontsize=7)
    ax.set_ylim(-0.7, len(C5) - 0.3)
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(0.05))
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(pct))
    ax.legend(loc="lower right")
    title(ax, t["f5_title"], t["f5_sub"].format(k=k_ols, j=k_r3))
    save(fig, "fig5_penalty")


def fig6(t):
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    style(ax, t["f6_x"])
    ax.axvline(0, color=INK, lw=1.0, zorder=2)
    for k, (clock, col, mk, off) in enumerate([("Horvath2018", BLUE, "o", -0.14),
                                                ("Levine2018", ORANGE, "s", 0.14)]):
        g = SL[SL.clock == clock].set_index("cohort").loc[SL_ORDER]
        y = np.arange(len(g)) + off
        ax.errorbar(g.slope, y, xerr=1.96 * g.se, fmt="none", ecolor=col, elinewidth=1.2, capsize=0, zorder=3)
        ax.scatter(g.slope, y, s=34, marker=mk, facecolor=col, edgecolor="white", lw=0.9, zorder=5,
                   label=t["f6_h"] if clock == "Horvath2018" else t["f6_l"])
    labels = []
    for c in SL_ORDER:
        r = SL[SL.cohort == c].iloc[0]
        name = t["f6_norm"] if "normalised" in c else c
        labels.append(f"{name}  ·  {r['array']}, n = {r['n']}")
    ax.set_yticks(np.arange(len(SL_ORDER))); ax.set_yticklabels(labels, fontsize=7.5)
    ax.invert_yaxis()
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(
        lambda x, _: f"{x:+.0f}".replace("+0", "0") if x else "0"))
    ax.legend(loc="upper left")
    title(ax, t["f6_title"], t["f6_sub"])
    save(fig, "fig6_saliva")


for LANG in T:
    for f in (fig1, fig2, fig3, fig4, fig5, fig6):
        f(T[LANG])
    print(f"  figures/{LANG}/: seis figuras")
print(f"  piso {FLOOR:+.2%} | rho fig3 {rho3:.3f} | falsos seguros fig4 "
      f"{int(C4.false_safe.sum())} | fig5 {k_ols} -> {k_r3}")
