#set page(paper: "a4", margin: (x: 2.3cm, y: 2.4cm), numbering: "1")
#set text(font: "Libertinus Serif", size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.62em, first-line-indent: 0pt, spacing: 0.9em)
#set heading(numbering: "1.1")
#show heading: it => block(above: 1.4em, below: 0.7em)[
  #set text(size: if it.level == 1 { 12pt } else { 10.5pt }, weight: "bold")
  #if it.numbering != none [#counter(heading).display(it.numbering)#h(0.6em)]
  #it.body
]
#show figure.caption: it => [
  #set text(size: 9pt)
  #set par(justify: true)
  *#it.supplement #context it.counter.display(it.numbering).* #it.body
]
#set figure(gap: 0.9em)
#show raw: set text(font: "DejaVu Sans Mono", size: 9pt)

#align(center)[
  #block(text(size: 15pt, weight: "bold")[
    Cell-composition adjustment of epigenetic age can inject\
    the confounding it removes
  ])
  #v(0.2em)
  #block(text(size: 11.5pt)[
    A transport-risk index, and a one-line fix
  ])
  #v(1.1em)
  #text(size: 10.5pt)[Luan Ivepe]
  #v(0.2em)
  #text(size: 9.5pt, style: "italic")[Independent researcher]
  #v(0.2em)
  #text(size: 9pt)[#link("mailto:luanivepe@gmail.com")[luanivepe\@gmail.com]]
  #v(0.2em)
  #text(size: 9pt)[Preprint draft — #datetime.today().display("[day] [month repr:long] [year]")]
]

#v(1.2em)

#block(inset: (x: 1.2em), [
  #text(weight: "bold")[Abstract] #h(0.6em)
  Epigenetic age acceleration in blood is routinely adjusted for immune cell
  composition, and the adjustment is almost always estimated and applied inside
  the same cohort. We show that when the coefficients of such an adjustment are
  estimated on one cohort and applied to another — the situation whenever a
  published correction is reused — the procedure can *add* composition signal
  rather than remove it. Fitting on forty samples and transporting, the residual
  composition term rises by 12.2 points of age-acceleration variance and the
  correction is worse than doing nothing in 89% of draws, while composition
  coefficients carrying no information at all sit at +0.5% with an interquartile
  range straddling zero: the damage comes from real coefficients estimated
  badly, not from noise. The failure is invisible within the fitting cohort,
  where the adjustment removes the composition term exactly by construction.
  Neither sample size nor the conditioning of the fitting cohort's composition
  matrix accounts for it; the cohort that fails worst is the best conditioned in
  our set. What does account for it is the standard excess-risk term for least
  squares under covariate shift, $(sigma^2 \/ n) dot tr(Sigma_"fit"^(-1)
  Sigma_"test")$, which we call the transport index. It ranks 63 directed
  configurations across four fitting cohorts at Spearman $rho = 0.907$, is
  asymmetric where symmetric distances are not, and predicts our worst observed
  failure to within half a point. Ridge-penalising the composition coefficients
  bounds the amplification directly and holds all twelve directed cohort pairs
  at or below zero, including a transport that costs +15.9 points unpenalised.
  Ordinary cross-validation on the fitting cohort does *not* size the penalty
  correctly on the cohort that fails. We report the index as a pre-computable
  diagnostic and the penalty as a default.
])

= Introduction

Blood is a mixture, and the mixture changes with age. Lymphocytes fall,
myeloid cells rise, and naive compartments give way to memory ones. Because
every leukocyte subtype carries its own methylome, an epigenetic clock applied
to whole blood reads partly the age of the cells and partly the census of the
tube @jaffe2014.

The size of that contamination is now well characterised. At twelve-type
resolution, immune composition explains 13% to 34% of age-acceleration variance
depending on the clock @zhang2024, and naive CD8 T cells read 15–20 years
younger than effector memory CD8 from the same donor @tomusiak2024. Two
responses exist. One designs the confounding out at the level of CpG selection
@tomusiak2024. The other — far more common, and the subject of this paper —
subtracts it after the fact: regress clock age on chronological age and
estimated cell proportions, and keep the residual. This is intrinsic epigenetic
age acceleration, and it is computed inside whatever dataset is at hand.

That last detail is the one that matters. *An adjustment estimated and evaluated
in the same cohort cannot be validated there.* The residual of a least-squares
fit is orthogonal to its predictors by construction, so the composition term is
driven to exactly zero whatever the coefficients are worth. Any assessment of
whether the adjustment works must therefore be made out of cohort — and the
field's own benchmark of eight cell-type correction methods does not test
cross-dataset transfer for any of them #footnote[Genome Biology 17:84 (2016),
"An evaluation of methods correcting for cell-type heterogeneity in DNA
methylation studies"; PMC4855979. Author list to be verified before
submission.].

Transfer is not a hypothetical. Whenever a published correction is reused, a
clinical cohort of a few dozen is adjusted with coefficients from a larger
study, or two groups compare adjusted estimates, coefficients cross a cohort
boundary. We ask what happens when they do.

= Methods

== Cohorts, reference panels and clocks

Four public whole-blood 450k series with chronological age: GSE40279
($n = 656$, ages 19–101), GSE61151 ($n = 184$), GSE50660 ($n = 464$) and
GSE42861 ($n = 689$). Cell proportions were estimated by constrained
non-negative least squares against two reference panels built here: a six-type
panel from Reinius et al. (GSE35069) and a twelve-type panel from the
FlowSorted.BloodExtended.EPIC reference (GSE167998) @salas2022. Our twelve-type
panel uses a simpler t-statistic probe selection than the published IDOL
library and recovers the known proportions of the twelve reconstructed mixtures
at $r = 0.79$, mean absolute error $0.027$; every composition number here should
be read as a lower bound on what the published library achieves. The panel's
monocyte channel runs high, and no per-type monocyte coefficient is interpreted.

Three clocks with published open coefficients carry the analysis: Horvath 2013,
Levine 2018 (PhenoAge) and Horvath 2018. Hannum 2013 was trained on GSE40279 and
is excluded wherever that cohort appears.

== The correction, and how it is scored

The correction is the standard one. On the fitting cohort we fit
$y = beta_0 + beta_1 a + C beta_c$, where $y$ is clock age, $a$ chronological
age and $C$ the estimated proportions (one column dropped, since they sum to
one), keep $beta_c$, and apply $y - (C_"test" - macron(C)_"fit") beta_c$ in the
test cohort. Chronological age stays in the fitting model — composition drifts
with age, and omitting it would strip real ageing along with the blood count —
but only the composition coefficients travel, so applying the correction needs
no age.

The outcome is the composition term left in the age residual of the test
cohort: the increment in $R^2$ from adding composition to a model of clock age
on chronological age, minus a permutation null, as a share of the uncorrected
age-acceleration variance. *Fit and measurement use different panels* — fitted
with twelve types, measured with six — because scoring a correction with its own
representation scores it against its own account of what it removed. We report
$Delta$ = after $-$ before; positive means the correction left the clock worse
than not correcting at all.

== The transport index

Write the estimation error $e = hat(beta)_c - beta_c$. The transported
correction leaves $-C_"test" e$, whose variance in the test cohort is
$e' Sigma_"test" e$. For least squares, $"Cov"(e) = sigma^2 \/ n dot
Sigma_"fit"^(-1)$, so

$ EE["damage"] prop (sigma^2 / n) dot tr(Sigma_"fit"^(-1) Sigma_"test") . $

We call $tr(Sigma_"fit"^(-1) Sigma_"test") \/ n$ the *transport index*. Both
second-moment matrices are taken after partialling out intercept and
chronological age, which is the space the coefficients live in.

*We claim no novelty for this quantity.* It is the standard excess-risk term for
least squares under covariate shift, and the failure mode it describes is named
in that literature: spectral inflation, directions carrying little variation in
training that carry more at evaluation @spectral2023. Theory for choosing ridge
regularisation under covariate shift exists as well @patil2024. What is new here
is the observation that a widely used epidemiological adjustment is exposed to
it, and the measurement of when it bites.

For ridge with penalty $lambda$ the corresponding term is
$tr(M Sigma_"fit" M Sigma_"test") \/ n$ with $M = (Sigma_"fit" + lambda I)^(-1)$.
Penalties are quoted as $alpha$, a scale-free multiplier of the mean eigenvalue
of the standardised cross-product, so the same $alpha$ means the same thing at
every $n$ and in every cohort. Coefficients are fitted by Frisch–Waugh–Lovell on
standardised columns and returned to the original scale; at $alpha = 0$ this
reproduces the unpenalised fit exactly.

== Subsampling

Where a fitting cohort is subsampled, draws are *stratified by age decile*. An
unstratified draw of forty from a cohort spanning 19 to 101 narrows and shifts
the age range, which would confound sample size with training range. Across the
grid the median age of a draw stays within 0.5 years of the full cohort's and
the 10–90 percentile spread within 2% of it.

= Results

== A transported correction can be worse than no correction

#figure(
  image("figures/fig1_curva_n.png", width: 95%),
  caption: [*The correction is harmful below roughly 160 fitting samples.*
  Composition left in the age residual after transporting, as a share of
  age-acceleration variance; above zero the correction left the clock worse than
  leaving it alone. Corrections fitted on age-stratified subsamples of GSE40279
  and carried to three external cohorts, 30 draws per point, three clocks. Grey
  lines are the three test cohorts separately. The orange rule is the same
  transport with composition coefficients fitted on permuted composition, which
  carry no information by construction.],
) <fig1>

Fitting the correction on subsamples of GSE40279 and transporting it to three
external cohorts (@fig1), the median $Delta$ crosses zero between $n = 160$ and
$n = 184$. Below that the standard adjustment is worse than doing nothing, and
at the bottom of the range it is not close: at forty samples $Delta = +12.2%$
and the correction is harmful in 89% of draws. Expressed against the
composition signal that was present, the correction at $n = 40$ does not remove
$0%$ of it — it *adds 2.7 times* what was there. At the full 656 it removes 61%
out of cohort.

The reference line matters for interpretation. Composition coefficients fitted
on *permuted* composition carry no information, and transporting them yields a
median $Delta$ of $+0.5%$ with an interquartile range of $-0.8%$ to $+2.1%$ over
180 draws — a distribution straddling zero. Meaningless coefficients are close
to harmless. It is real coefficients, estimated on too few samples, that do the
damage, and they do roughly twenty times more of it than noise does.

The median is not the whole story. Taking the spread across draws *within* each
clock-by-cohort cell, the median cell's upper quartile does not fall below zero
until $n approx 320$: the typical draw starts helping near 160–184, but an
unlucky one keeps hurting well past that.

== Neither sample size nor conditioning explains it

Sample size is a strong cause and not a sufficient one. Fitting on 184 samples
of GSE40279 and transporting gives $Delta = -0.9%$ — mildly beneficial — while
fitting on GSE61151's own 184 and transporting to GSE40279 gives $+5.0%$, with
all three clocks harmed. Something about the fitting cohort beyond its size is
doing the work.

The natural candidate is conditioning: collinear cell proportions give unstable
coefficients, and variance inflation factors above 100 are documented to flip
the sign of most coefficients in this setting @meredith2019. In our twelve-type
panel the largest variance inflation factor runs 215 to 302, worse than the
113.7 reported at six types. But conditioning barely moves with $n$ — the
condition number of the residualised composition block is 65.9 at $n = 40$ and
54.5 at $n = 656$ — and within each $n$-by-clock-by-cohort cell it carries no
information about the damage at all (median Spearman $rho = -0.025$, positive in
48% of cells, a coin flip).

Decisively: *GSE61151, the cohort whose transported correction fails worst, is
the best conditioned cohort in the set* (condition number 40.8, largest VIF
117.8, against 54.4–65.9 and 212–302 for every GSE40279 subsample including the
full cohort). Conditioning of the fitting cohort alone cannot be the mechanism.

== The transport index ranks the damage

#figure(
  image("figures/fig2_indice.png", width: 95%),
  caption: [*The transport index orders 63 directed configurations.* Each point
  is one (fitting cohort, test cohort, fitting size) configuration; marker shape
  gives the fitting cohort. Spearman $rho = 0.907$. No configuration below an
  index of 0.05 was harmful; 81% above it were.],
) <fig2>

The index is a joint property of both cohorts, which is precisely why measuring
$Sigma_"fit"$ alone returned nothing: a well-conditioned $Sigma_"fit"$ can still
have $Sigma_"fit"^(-1)$ amplify the directions in which $Sigma_"test"$ happens
to carry variance.

Across 63 directed configurations built from four fitting cohorts, the index
ranks the median damage at $rho = 0.907$ ($p = 2 times 10^(-24)$; @fig2), and it
does so within each fitting cohort separately — $rho$ = 0.858, 0.951, 0.907 and
0.897 — so the ranking is not a property of one cohort's subsample ladder. Damage
crossed zero near an index of 0.05: *none* of the 21 configurations below that
value was harmful, and 81% above it were. Between 0.05 and 0.16 both outcomes
occur, so the threshold is a one-sided floor rather than a crossing point.

One negative result belongs here. The index does *not* predict the damage of an
individual draw: within a fixed $n$ and test cohort, median $rho = 0.021$ at 56%
of cells, worse than the quantity it was meant to beat. The theory predicts this.
With $n$ and the test cohort held fixed, $Sigma_"fit"$ is nearly identical across
draws so the index barely varies, while the realised $e' Sigma_"test" e$ is
chi-square-like on eleven degrees of freedom. A near-constant predictor cannot
track an outcome dominated by realisation noise. The index predicts expected
damage across configurations, not the draw in front of you.

#figure(
  image("figures/fig3_assimetria.png", width: 95%),
  caption: [*The index is asymmetric, and so is the damage.* Each reversible pair
  fitted in both directions at $n = min(n_A, n_B)$, so only the fitting
  covariance differs. Orange marks the direction the index calls worse; it is to
  the right of the open marker in five of six. In the sixth
  (GSE40279 / GSE42861) the two indices are 0.0275 and 0.0271 — the index
  predicted indifference and was scored as though it had predicted a
  direction.],
) <fig3>

Symmetric distances between composition covariance structures cannot do this
job, and the asymmetry is the reason. Fitted at matched $n$ so that only the
fitting covariance differs, the index names the worse direction in five of six
reversible pairs (@fig3; binomial $p = 0.109$, which is not significance and is
not presented as such — with six pairs no threshold both clears 0.05 and
tolerates one miss). The miss is a pair whose two indices are 0.0275 and 0.0271:
the index predicted indifference and was scored as though it had predicted a
direction.

For the configuration that defeated our earlier attempts, the index for
GSE61151 $arrow$ GSE40279 is 0.174 against 0.038 for its own reverse, a ratio of
4.6. An index near 0.18 corresponds to about $+4.6%$ of damage in @fig2; the
measured value is $+5.0%$.

A simulation closes the mechanism. Drawing two synthetic cohorts with the *same*
true coefficient vector, the same noise and the same $n$, so that the only thing
that can go wrong is estimation carried across a covariance mismatch, swapping
which covariance does the fitting moves the residual from $+21.2%$ to $+84.6%$ at
$n = 40$. Observed residuals regress on predicted with slope 0.92 through the
origin, $r = 0.815$.

== Penalising the coefficients removes the failure

#figure(
  image("figures/fig4_conserto.png", width: 95%),
  caption: [*Ridge holds all twelve directed transports at or below zero.* Each
  row is one directed cohort pair at $n = min(n_A, n_B)$; the arrow runs from the
  unpenalised fit to $alpha = 3$.],
) <fig4>

Ridge replaces $Sigma_"fit"^(-1)$ with $(Sigma_"fit" + lambda I)^(-1)$ and bounds
the amplification directly. A penalty only counts as a fix if it removes the harm
*and* keeps the benefit: shrinking a correction to zero removes harm and is not a
method. Requiring, before the run, that a penalty be no worse than doing nothing
at $n = 40$ *and* retain at least 70% of the unpenalised benefit at $n = 656$,
exactly one value on our grid qualifies, $alpha = 3$. Its curve is flat — between
$-1.4%$ and $-2.0%$ across a sixteen-fold range of fitting size — which is the
point: the correction stops depending on how large the fitting cohort was.
$alpha = 10$ is the degenerate corner the criterion was written to catch:
coefficient norm 13, damage pinned near $-0.9%$ everywhere, harmless and useless.

The same $alpha$, arrived at independently from the real cohort that failed,
turns that transport from $+5.0%$ to $-0.5%$. Across all twelve directed pairs at
matched $n$ it holds the median at or below zero (@fig4), including
GSE61151 $arrow$ GSE42861, the worst transport we found, which goes from $+15.9%$
to $-2.5%$. It is not free: where the unpenalised correction already transports
well, the penalty costs some of the benefit — GSE40279 $arrow$ GSE42861 gives up
most of it, $-2.8%$ to $-0.7%$. Penalising is insurance and insurance has a
premium in the cases that did not need it.

*Cross-validation is not enough to size it.* Leave-one-out on the fitting cohort
selects sensibly on GSE40279 subsamples, and on the cohort that actually breaks
the transport it selects $alpha = 0.3$ and leaves $+2.9%$ of the damage standing.
Cross-validation optimises prediction inside the fitting cohort, and that
objective does not know the coefficients are about to be shipped elsewhere.

== The diagnostic survives two teams and two reference panels

Everything above uses one deconvolution for both ends of each transport, which
is not the situation the result is about. Estimating the fitting cohort's
composition from one reference and the test cohort's from another — the same six
labels, different donors, probes and array generation — the two panels agree
well (median per-type $r = 0.936$), and the index still ranks configurations at
$rho = 0.776$ against $rho = 0.793$ for matched panels. Two references disagreeing
about a cell type *is* a covariance mismatch, so the index sees it without being
told. Ridge at $alpha = 3$ holds all twelve pairs at or below zero in the
mismatched arm as well.

= Discussion

Three things follow for practice.

*First, a within-cohort adjustment cannot be validated within that cohort.* This
is arithmetic, not an empirical finding, and it applies to every paper that
reports an intrinsic age acceleration without an out-of-cohort check. The
apparent perfection of the adjustment at home is guaranteed and carries no
information.

*Second, small fitting cohorts are not merely noisy.* A coefficient vector with
no information does essentially nothing when transported. Real coefficients
estimated on forty samples add 2.7 times the composition signal that was
present. The failure is not a diluted version of success.

*Third, the risk is computable in advance.* The transport index needs only the
two cohorts' estimated proportions and the fitting size — no outcome, no second
validation cohort, no held-out labels. Below 0.05 we observed no harm in 63
configurations. We would not read that as a guarantee, but as a threshold worth
computing before reusing a published correction.

Our recommendation is narrower than "use ridge and cross-validate". Penalise —
and do not let cross-validation on the fitting cohort alone decide how much.

= Limitations

Four cohorts, all whole blood on 450k, all adult, and twelve directed pairs. The
asymmetry count is six. Both reference panels are built here with a cruder probe
selection than the published libraries, so the disagreement between them is a
lower bound on what two real laboratories would produce. $alpha = 3$ is this
panel with these clocks; what we would expect to transport is "penalise", not the
number. The index is derived for least squares with a correctly specified linear
composition term, and the simulation confirming its constant assumes exactly
that. And the index predicts expected damage across configurations, not the
individual draw — a negative result we report rather than bury.

Nothing here bears on whether epigenetic clocks measure anything real. This is a
statement about a correction applied to them.

= Data and code availability

All five series are public: GSE40279, GSE61151, GSE50660, GSE42861 (whole blood),
GSE35069 and GSE167998 (reference panels). No third-party data is redistributed.
Analysis code, the full stage-by-stage record including every conclusion later
overturned, and the figure scripts are at
#link("https://github.com/KTHimiko/clock-lab")[github.com/KTHimiko/clock-lab].

#v(0.8em)
#block(fill: luma(245), inset: 9pt, radius: 3pt, width: 100%, [
  #text(size: 9pt, weight: "bold")[Before submission — to verify]
  #text(size: 9pt)[
    #list(
      spacing: 0.45em,
      [The author list of the 2016 Genome Biology eight-method benchmark
       (PMC4855979) and of the covariate-shift preprint arXiv:2312.17463 were not
       confirmed and are marked as such rather than guessed.],
      [Repository is currently private; it must be public before the preprint
       cites it.],
      [A reproducibility statement and the exact package versions should be
       added.],
    )
  ]
])

#bibliography("refs.bib", title: "References", style: "nature")
