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

#align(center)[
  #block(text(size: 15pt, weight: "bold")[
    Transported cell-composition corrections of epigenetic age\
    can add the confounding they are meant to remove
  ])
  #v(0.2em)
  #block(text(size: 11.5pt)[Estimation noise, model shift, and a penalty as the safeguard])
  #v(1.1em)
  #text(size: 10.5pt)[Luan Ivepe]
  #v(0.2em)
  #text(size: 9.5pt, style: "italic")[Independent researcher]
  #v(0.2em)
  #text(size: 9pt)[#link("mailto:luanivepe@gmail.com")[luanivepe\@gmail.com]]
  #v(0.2em)
  #text(size: 9pt)[Preprint draft, revised after review — #datetime.today().display("[day] [month repr:long] [year]")]
]

#v(1.2em)

#block(inset: (x: 1.2em), [
  #text(weight: "bold")[Abstract] #h(0.6em)
  Epigenetic age acceleration in blood is routinely adjusted for immune cell
  composition by regressing clock age on estimated cell proportions. Within the
  cohort where it is fitted, the adjustment does what it claims. We study the case
  where its coefficients are applied to a different cohort. Across four public
  whole-blood cohorts, scoring only clocks never trained on the cohorts involved,
  a correction fitted on forty samples and transported increased the composition
  signal in 93% of draws for the age-estimating clocks, and on the naive/memory
  lymphocyte axis the increase was three to four times what a six-type panel
  shows; for a pace-of-ageing clock, DunedinPACE, the same small-sample transport
  was neutral. The error has two parts.
  Estimation noise, amplified by differences between the cohorts' composition
  covariance, falls as $1\/n$ and is reproduced almost exactly by coefficients
  that carry no information. Model shift — a composition effect that differs
  between cohorts — does not fall with fitting size and accounts for the worst
  transports. Because the second part depends on the target cohort's own
  coefficients, nothing we could compute beforehand certified a transport as
  safe: a closed-form predictor labelled 30 of 73 harmful transports as safe, and
  the transport index added only modest information beyond sample size
  (within-$n$ Spearman $rho approx 0.36$; block-permutation $p = 0.036$). A ridge
  penalty on the composition coefficients cut harmful (pair × clock) transports
  from half to 3%, across measurement and reference panels. Cross-validation on
  the fitting cohort chose penalties too weak to do the same.
])

= Introduction

Blood composition shifts with age: lymphocytes decline, myeloid cells increase,
and memory cells replace naive ones. Each leukocyte subtype has its own
methylation profile, so a clock applied to whole blood partly measures the
composition of the sample @jaffe2014. With twelve cell types resolved,
composition explains 13–34% of age-acceleration variance depending on the clock
@zhang2024, and naive CD8 T cells read 15–20 years younger than effector memory
CD8 cells from the same donor @tomusiak2024.

The usual remedy is to regress clock age on chronological age and estimated
proportions, and keep the residual. This is normally done inside the dataset at
hand, where it removes the linear composition term exactly. For the same reason
it cannot be checked there: a least-squares residual is orthogonal to its
predictors regardless of how good the coefficients are. Published evaluations of
cell-type adjustment are simulation-based and within-dataset @mcgregor2016.

We examine the case where the coefficients leave the cohort that produced them —
reuse of published coefficients, small cohorts borrowing from larger ones, or a
fixed correction applied to new samples.

= Methods

== Cohorts, reference panels and clocks

Four public whole-blood 450k series with chronological age: GSE40279
($n = 656$, ages 19–101), GSE61151 ($n = 184$), GSE50660 ($n = 464$, smoking
cohort) and GSE42861 ($n = 689$, rheumatoid arthritis case-control). Proportions
were estimated by constrained non-negative least squares against a six-type
panel (GSE35069) and a twelve-type panel resolving naive and memory lymphocytes
(GSE167998) @salas2022, both built here; the twelve-type panel recovers known
mixture proportions at $r = 0.79$, mean absolute error 0.027.

Clocks: Horvath 2013 @horvath2013, Levine 2018 @levine2018 and Horvath 2018
@horvath2018. A clock is excluded from any pair involving a cohort it was trained
on. Hannum 2013 and Horvath 2013 were both trained on GSE40279 (Horvath 2013 lists
it as training set 3; GSE42861 was a test set only). Levine 2018 (InCHIANTI) and
Horvath 2018 were trained on none of the four. Pairs involving GSE40279 are
therefore scored with Levine 2018 and Horvath 2018 only. As a clock of a different
kind we add DunedinPACE @belsky2022, which estimates the pace of ageing and was
trained on a cohort not in GEO; we reimplemented it from its package's published
model data and validated it (cohort means 0.93–1.05; current smokers +0.14 faster
than never smokers, $p = 3 times 10^(-11)$).

== Correction and scoring

On the fitting cohort we fit $y = beta_0 + beta_1 a + C beta_c$ (clock age $y$,
chronological age $a$, proportions $C$ with one column dropped) and apply
$y - (C_"test" - macron(C)_"fit") beta_c$ in the test cohort. The outcome is the
composition term remaining in the test cohort's age residual — the $R^2$
increment from composition over chronological age, minus a permutation null, as a
share of uncorrected age-acceleration variance. Net damage $Delta$ is this share
after correction minus before; positive $Delta$ means the transported correction
did worse than no correction. The correction is fitted with twelve types and
scored with six to avoid scoring it on its own representation; sensitivity
analyses score it on twelve types and on the naive/memory columns alone.

All counts are per clock: a pair where one clock is harmed and another helped
counts as two cells. Pooling clocks within a pair concealed a substantial part of
the harm in an earlier version of this analysis.

== Two components of transport error

Let $b$ be fitted on cohort A and $S_B$ be the composition covariance of cohort B,
both after removing intercept and age. The composition left in B is
$(beta_B - b)' S_B (beta_B - b)$. With $b = beta_A + e$ and
$"Cov"(e) = (sigma^2 \/ n) S_A^(-1)$, its expectation is a specification term
$(beta_A - beta_B)' S_B (beta_A - beta_B)$ plus an estimation term
$sigma^2 dot tr(S_A^(-1) S_B) \/ n$. We call $tr(S_A^(-1) S_B)\/n$ the transport
index. The estimation term is the standard excess risk of least squares under
covariate shift @eyre2024; the specification term is model shift @lei2021. We
claim neither as new.

A reference with no information is obtained by shuffling whole rows of the
fitting cohort's composition matrix before fitting, which preserves the
collinearity between cell types. It is measured at every fitting size.

== Penalty, subsampling and inference

Ridge coefficients are fitted on standardized, age-residualized composition with
penalty $alpha$ times the mean eigenvalue, so that $alpha$ is comparable across
cohorts and sizes; $alpha = 0$ reproduces least squares. Under distribution shift
the optimal penalty need not even be positive @patil2024; a fixed positive
penalty is used here as a conservative default. Subsamples are stratified by age
decile. Where configurations share cohorts, significance is assessed by block
permutation @winkler2015, swapping whole index profiles between cohort pairs.

= Results

== A transported correction can be worse than none

#figure(
  image("figures/en/fig1_curve.png", width: 95%),
  caption: [*Net damage of a transported correction by fitting size.* Fitted on
  age-stratified subsamples of GSE40279 and applied to three external cohorts;
  Levine 2018 and Horvath 2018; 30 draws per size. Above zero the correction did
  worse than none. Orange: the same procedure with composition rows shuffled
  before fitting.],
) <fig1>

Fitted on 40 samples and transported, the correction had a median $Delta$ of
+18.9% (IQR +8.2% to +34.2%) and was harmful in 93% of draws (@fig1). Across three
independent sets of draws the median at $n = 40$ ranged from +11.8% to +18.9%.
The median was near zero between 160 and 240 samples and −1.6% at the full 656,
where 33% of (clock × cohort) values were still harmful. Scored on twelve cell
types instead of six, the damage at $n = 40$ was +48.1%; scored on the
naive/memory columns alone, +30.6%. That measurement shares the fitting panel and
is biased toward the correction, so the six-type curve understates the harm.

== Two components

#figure(
  image("figures/en/fig2_components.png", width: 95%),
  caption: [*Composition left by real and by shuffled coefficients.* The real
  fit keeps a floor of about 4.5 points at full fitting size; above that floor,
  its excess tracks the shuffled reference.],
) <fig2>

Shuffled coefficients carry no information, yet at $n = 40$ they left +15.6% net
damage, falling with $n$ at a log-log slope of −0.93 (@fig2). Most of the
small-sample harm is therefore estimation noise. The real fit, however, still
left 4.5 points of composition at $n = 656$, where the shuffled reference left
1.0. Subtracting that floor, the remainder matched the shuffled reference closely
(9.6 vs 9.6 points at $n = 60$; 3.1 vs 3.1 at 120; 1.2 vs 1.3 at 240).

The floor is not an artefact of scoring with a different panel: applied within
GSE40279, the same correction removed 95% (Levine) and 85% (Horvath 2018) of the
six-type composition signal. Transported at full size, it removed between 92% and
−109% depending on the pair — in the worst case doubling the signal.

== Model shift

#figure(
  image("figures/en/fig3_model_shift.png", width: 95%),
  caption: [*Specification term against composition left at full fitting size.*
  24 (directed pair × clock) configurations. The term is bias-corrected for the
  estimation noise of both fits.],
) <fig3>

The bias-corrected specification term ranked the full-size leftover at Spearman
$rho = 0.633$ ($p = 0.0009$; @fig3). It accounts for the worst transport found,
GSE61151 into the arthritis cohort: predicted +26.6% and +51.8% for the two
clocks, observed +30.8% and +30.5%. Pairwise Wald tests of equal coefficients
rejected in 4 of 12 comparisons after Bonferroni correction, short of our
pre-set threshold of 6 (9 of 12 at nominal $p < 0.05$). The two measures weight
coefficient differences differently: Wald by estimation precision, the
specification term by variance in the target cohort. Within cohorts, Levine 2018
coefficients differed between arthritis cases and controls and between ever and
never smokers ($p = 0.013$ each); for Horvath 2018 this test did not reject
($p = 0.61$, $0.48$). Adjusting for disease and smoking left between-cohort
differences essentially unchanged.

Model shift is not a batch effect. Within GSE42861, where cases and controls share
study, array and laboratory, we fitted the correction on a random half of the
controls and applied it both to the other half and to the arthritis cases (30
paired splits). Applied to cases it left 7.4, 7.4 and 2.6 more points of
composition than applied to controls (Horvath 2013, Levine 2018, Horvath 2018),
although the transport index was lower for cases — so none of the excess is
estimation noise, and the estimated model-shift component was 8.7, 8.4 and 6.3
points. Horvath 2018 is affected even though the Wald test above did not detect a
coefficient difference. The choice of reference population also changes the
estimated disease effect: for Horvath 2018, the age-adjusted arthritis effect was
−0.74 years with a whole-cohort correction and −1.39 years with a controls-only
correction.

== What can be known before transporting

#figure(
  image("figures/en/fig4_index.png", width: 95%),
  caption: [*Transport index against net damage, per clock.* Orange: transports
  that the closed-form predictor labelled safe and that were harmful. Marker
  shape gives the fitting cohort.],
) <fig4>

Because the index scales as $1\/n$, part of its correlation with damage comes from
sample size alone. At the three sizes available for every pair, a block
permutation that keeps each pair's fitting sizes gave a null median of
$rho = 0.42$; the observed value was 0.579 ($p = 0.036$). Within a single fitting
size the index ranked damage at $rho$ = 0.41, 0.36 and 0.33. It carries real but
modest information about which pairs are risky (@fig4). Counted per clock, 12 of
40 transports with an index below 0.05 were harmful, so no index threshold marks
a safe region.

A closed-form estimate of net damage that assumes shared coefficients,
$2 hat(sigma)^2 dot "index" - b' S_B b$, needs only the fitting cohort and the
target's proportions. It predicted the sign in 71% of configurations
($rho = 0.633$), against 81% ($rho = 0.884$) for an oracle that knows the target's
own coefficients. It errs toward reassurance: 30 of the 73 transports it
labelled safe were harmful, concentrated where the specification term is large.

== A penalty removes most of the harm

#figure(
  image("figures/en/fig5_penalty.png", width: 95%),
  caption: [*Unpenalised against ridge-penalised transport, per (pair × clock)
  cell at matched fitting size* $n = min(n_A, n_B)$.],
) <fig5>

Without a penalty, 15 of 30 (pair × clock) cells were harmful. At $alpha = 3$, one
remained, at +0.2% (@fig5); at $alpha = 10$, none. The penalty has a cost: where
the unpenalised correction already helped, the median benefit fell from −4.0% to
−3.3% at $alpha = 3$ and to −1.5% at $alpha = 10$. The fixed penalty also held
when the two cohorts were deconvolved from different reference panels (14
harmful cells unpenalised, 1 at $alpha = 3$) and when scored on twelve types or
the naive/memory axis (no harmful cell at $alpha = 3$).

Choosing the penalty by leave-one-out cross-validation on the fitting cohort did
not work as well: it selected $alpha = 0.3$ in the median cell and left 30% of
cells harmful, against 3% at fixed $alpha = 3$. Cross-validation optimises fit
within the fitting cohort and cannot account for where the coefficients will be
used.

== A pace-of-ageing clock

For DunedinPACE, the transported correction fitted on 40 samples of GSE40279 was
neutral (median +0.3%, harmful in 51% of draws), although shuffled coefficients
still did +4.9% of damage: the noise component was present, but real
coefficients removed enough genuine signal to offset it, because DunedinPACE's
composition effect transports well from that cohort. The small-sample harm is
therefore a property of the clock and the cohort pair, not a constant. What
replicated was the rest: at matched sizes, 6 of 12 directed pairs were harmful
unpenalised and none at $alpha = 3$, and a correction fitted on controls left
10.5 more points of composition in arthritis cases than in other controls.

= Discussion

Within-cohort composition adjustment cannot be validated inside the cohort; this
follows from least squares, not from data. Transported, the adjustment's small-
sample damage is mostly estimation noise — shuffled coefficients do nearly as
much of it — while its large-sample damage comes from differences in the
composition effect between cohorts. That second component is invisible without
the target's own coefficients, and a target large enough to estimate them could
simply be adjusted within itself.

For practice this suggests three things. Where the target cohort is large enough,
fit the adjustment within it. Where coefficients must be transported, penalise
them with a fixed, substantial penalty rather than one chosen by cross-validation
on the source cohort. And treat any pre-transport diagnostic, including the
transport index, as a ranking of risk rather than a guarantee. The same caution applies
inside a single study: a correction fitted on controls and applied to patients is
a transport, and here it changed the estimated disease effect by up to a factor
of two.

= Limitations

Four adult whole-blood cohorts on one array, two of them defined by disease or
exposure; results for other tissues, ages or platforms are untested. Reference
panels were built here with simpler probe selection than published libraries.
Medians at small fitting sizes vary between independent sets of draws (+11.8% to
+18.9% at $n = 40$), and the small-sample harm itself depends on the clock and the
fitting cohort (neutral for DunedinPACE from GSE40279). The penalty value is specific to this panel and these
clocks. The decomposition assumes a linear composition effect. Sensitivity
scoring on twelve types shares the fitting panel. None of this bears on whether
epigenetic clocks measure biological ageing; it concerns one correction applied
to them.

= Data and code availability

All series are public (GSE40279, GSE61151, GSE50660, GSE42861, GSE35069,
GSE167998). Analysis code, the stage-by-stage record including every overturned
conclusion, and figure scripts are at
#link("https://github.com/KTHimiko/clock-lab")[github.com/KTHimiko/clock-lab]
(to be made public before submission).

#bibliography("refs.bib", title: "References", style: "nature")
