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
  #block(text(size: 11.5pt)[Estimation noise, model shift, and the reach of a penalty])
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
  composition by regressing clock age on estimated cell proportions. Within the
  cohort where it is fitted the adjustment does what it claims; we study what
  happens when its coefficients are applied to another cohort. In six public
  blood cohorts, scoring only clocks never trained on the cohorts involved, a
  correction fitted on forty samples and transported increased the composition
  signal in 83–95% of draws for the age-estimating clocks; for the
  pace-of-ageing clock DunedinPACE it was neutral. In the units a study reports,
  it moved the estimated effect of smoking or disease 0.91 years away from the
  within-cohort estimate, against 0.55 years for applying no correction at all,
  and reversed its sign in 3 of 24 configurations. The error has two parts:
  estimation noise, which falls as $1\/n$ and is reproduced by coefficients that
  carry no information, and model shift — a composition effect differing between
  cohorts — which did not fall with fitting size up to 2,639 samples. Because
  model shift depends on the target's own coefficients, nothing computable
  beforehand certified a transport as safe: of 73 transports a closed-form
  predictor called safe, 30 were harmful. A fixed ridge penalty cut harmful
  (pair × clock) transports from 23 of 72 to 1, where a penalty fading with
  sample size, or one chosen by cross-validation, did not. In saliva, where
  composition accounts for up to half of age acceleration, transports were
  harmful in 7 of 8 cells with or without the penalty, and pooling studies did
  not help. What separates the tissues is how far cohorts disagree relative to
  the average effect: the between-cohort standard deviation of the composition
  slope was 0.09–0.55 of the mean effect in blood and 0.91–5.34 in saliva. A
  penalty shrinks a coefficient toward zero, which is close to every cohort's own
  only in the first case.
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
fixed correction applied to new samples. Adjustment coefficients fitted in one
dataset and applied to others are already published practice: a saliva adaptation
of blood clocks fits them on eight pooled studies and applies them to held-out
ones @galkin2021, and its training set includes GSE78874, one of the cohorts
where we find transport fails. We test it in blood, and in saliva,
where composition is larger and a published adaptation already transports it.

= Methods

== Cohorts, reference panels and clocks

Four public whole-blood 450k series with chronological age: GSE40279
($n = 656$, ages 19–101), GSE61151 ($n = 184$), GSE50660 ($n = 464$, smoking
cohort) and GSE42861 ($n = 689$, rheumatoid arthritis case-control). Proportions
were estimated by constrained non-negative least squares against a six-type
panel (GSE35069) and a twelve-type panel resolving naive and memory lymphocytes
(GSE167998) @salas2022, both built here; the twelve-type panel recovers known mixture proportions at
$r = 0.79$, mean absolute error 0.027. Because both were built here and share
their construction, every blood result was repeated on published libraries: the
Salas et al. 2022 twelve-type reference on its IDOL-optimised 450k probe set for
the fit @salas2022 and the published seven-type blood reference for the
measurement @teschendorff2017, which share 20 probes of 600 and 333.

Clocks: Horvath 2013 @horvath2013, Levine 2018 @levine2018 and Horvath 2018
@horvath2018. A clock is excluded from any pair involving a cohort it was trained
on. Hannum 2013 and Horvath 2013 were both trained on GSE40279 (Horvath 2013 lists
it as training set 3; GSE42861 was a test set only). Levine 2018 (InCHIANTI) and
Horvath 2018 were trained on none of the four. Pairs involving GSE40279 are
therefore scored with Levine 2018 and Horvath 2018 only. A fifth cohort, GSE132203
($n = 795$, EPIC array, mostly African American), tests replication across array
generation and ancestry; clocks are scored on it only above 95% probe coverage
(Levine 2018, Horvath 2018). As a clock of a different
kind we add DunedinPACE @belsky2022, which estimates the pace of ageing and was
trained on a cohort not in GEO; we reimplemented it from its package's published
model data and validated it (cohort means 0.93–1.05; current smokers +0.14 faster
than never smokers, $p = 3 times 10^(-11)$). A sixth cohort, GSE55763
@lehne2015 (450k, London), serves as a large fitting cohort: after dropping all
72 technical-replicate arrays it has 2,639 unrelated adults (ages 24–75); it
postdates Horvath 2013 and is in no clock's training set, and every clock
cleared the coverage and age checks on it.

Saliva, a mixture of buccal epithelium and leukocytes, was tested in three adult
cohorts: GSE232891 (EPIC, $n = 552$; inflammatory bowel disease and controls),
GSE232332 (EPIC, $n = 265$ after removing technical replicates; oesophageal
cancer and controls) and GSE78874 (450k, $n = 259$; betas from raw signal).
Proportions came from the EpiDISH references @teschendorff2017 @zheng2018: a
nine-type hierarchical fit (epithelium, fibroblast, seven immune subtypes) and a
three-type measurement (epithelium, fibroblast, immune), which share their first
step and so bias measurement toward the correction. The first two cohorts come
from one group and their files carry no genotyping probes, so shared individuals
could not be excluded and they were never paired. Levine 2018 and Horvath 2018
cleared coverage and age checks in all three. Because that measurement is the
fit's own first step, saliva was also scored with an independent panel built from
GSE147318 @middleton2022, children's saliva sorted into immune and epithelial
fractions; its 300 probes share 3.3% with EpiDISH and it is used as a relative
immune score, since its absolute scale does not transfer between studies. A
fourth saliva cohort, GSE149747
(EPIC, another group), contributes its 44 baseline samples to the slope
comparison and to the pooled fits, but is too small for a transport curve.

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
penalty is used here as a conservative default. The value $alpha = 3$ was fixed
on the first four cohorts, before GSE132203 and GSE55763 entered the project, so
its behaviour on those two and on saliva is out-of-sample for the choice. Subsamples are stratified by age
decile. Where configurations share cohorts, significance is assessed by block
permutation @winkler2015, swapping whole index profiles between cohort pairs.

= Results

== A transported correction can be worse than none

#figure(
  image("figures/en/fig1_curve.png", width: 95%),
  caption: [*Net damage of a transported correction by fitting size.* Fitted on
  age-stratified subsamples of GSE40279 and applied to three external cohorts;
  Levine 2018 and Horvath 2018; 100 draws per size. Above zero the correction did
  worse than none. Orange: the same procedure with composition rows shuffled
  before fitting.],
) <fig1>

Fitted on 40 samples and transported, the correction had a median $Delta$ of
+16.1 p.p. (IQR +5.5 to +29.7) and was harmful in 88% of draws over 100 draws per
size (@fig1); three earlier sets of 30 draws gave medians from +11.8 to +18.9, so
single small sets are unstable at this size. The median was near zero between 160
and 240 samples and −1.6 p.p. at the full 656, where 33% of (clock × cohort)
values were still harmful. Scored on twelve cell types instead of six, the damage
at $n = 40$ was +48.1 p.p.; scored on the naive/memory columns alone, +30.6. That measurement shares the fitting panel and
is biased toward the correction, so the six-type curve understates the harm.

== Two components

#figure(
  image("figures/en/fig2_components.png", width: 95%),
  caption: [*Composition left by the real fit, against a floor plus noise.*
  Blue: what the transported correction leaves. Grey: the noise-free floor, the
  3.8 p.p. left at full fitting size after subtracting what shuffled coefficients
  still do there. Orange: that floor plus the damage shuffled coefficients do at
  each size. Both series are composition left, in the same units.],
) <fig2>

Shuffled coefficients carry no information, yet at $n = 40$ they did +17.4 p.p.
of net damage, falling with $n$ at a log-log slope of −1.16 (@fig2). Most of the
small-sample harm is therefore estimation noise. The real fit, however, still
left 4.5 p.p. of composition at $n = 656$, where shuffled coefficients still did
0.8; the difference, 3.8 p.p., is a floor that fitting size does not remove.
Adding that floor to the shuffled damage reproduces the real fit's curve within
0.6 p.p. at every size from 60 upward, and underestimates it by 1.2 at $n = 40$
(@fig2).

The floor is not an artefact of scoring with a different panel: applied within
GSE40279, the same correction removed 95% (Levine) and 85% (Horvath 2018) of the
six-type composition signal. Transported at full size, it removed between 92% and
−109% of the signal it found, depending on the pair — in the worst case doubling
it. (Percentages of the signal found are relative; $Delta$ and composition left
are in p.p. of age-acceleration variance throughout.)

== Model shift

#figure(
  image("figures/en/fig3_model_shift.png", width: 95%),
  caption: [*Specification term against composition left at full fitting size.*
  24 (directed pair × clock) configurations. The term is bias-corrected for the
  estimation noise of both fits.],
) <fig3>

The bias-corrected specification term ranked the full-size leftover at Spearman
$rho = 0.633$ ($p = 0.0009$; @fig3). It accounts for the worst transport at full fitting size,
GSE61151 into the arthritis cohort: predicted +26.6 and +51.8 p.p. for the two clocks, observed +30.8 and +30.5. Pairwise Wald tests of equal coefficients
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
−0.74 years (95% CI −1.05 to −0.37) with a whole-cohort correction and −1.39
(−2.08 to −0.80) with a controls-only one, a difference of −0.66 (−1.13 to −0.31)
by a bootstrap that refits both corrections in each resample.

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
target's proportions. It predicted the sign in 71% of configurations ($rho = 0.6332$ over 126 configurations — coincidentally near the specification term's $rho$ above, which is a different correlation over 24), against 81% ($rho = 0.884$) for an oracle that knows the target's own coefficients. It errs toward reassurance: 30 of the 73 transports it
labelled safe were harmful, concentrated where the specification term is large.

== A penalty removes most of the harm

#figure(
  image("figures/en/fig5_penalty.png", width: 95%),
  caption: [*Unpenalised against ridge-penalised transport, per (pair × clock)
  cell at matched fitting size* $n = min(n_A, n_B)$.],
) <fig5>

Without a penalty, 15 of 30 (pair × clock) cells were harmful. At $alpha = 3$, one remained, at +0.2 p.p. (@fig5); at $alpha = 10$, none. The penalty has a cost: where the unpenalised correction already helped, the median benefit fell from −4.0 p.p. to −3.3 at $alpha = 3$ and to −1.5 at $alpha = 10$. The fixed penalty also held
when the two cohorts were deconvolved from different reference panels (14
harmful cells unpenalised, 1 at $alpha = 3$) and when scored on twelve types or
the naive/memory axis (no harmful cell at $alpha = 3$).

Choosing the penalty by leave-one-out cross-validation on the fitting cohort did
not work as well: it selected $alpha = 0.3$ in the median cell and left 30% of
cells harmful, against 3% at fixed $alpha = 3$. Cross-validation optimises fit
within the fitting cohort and cannot account for where the coefficients will be
used.

The penalty's size is defined relative to the fitting cohort's composition
variance, so its proportional shrinkage does not fade with $n$. A conventional
fixed-$lambda$ penalty does fade; set equal to $alpha = 3$ at $n = 40$, it is
0.045 at $n = 2639$. Over all 30 directed pairs among the six cohorts at matched
size (72 age-clock cells), it left 11 cells harmful against 1 for $alpha = 3$ (23
unpenalised); the cells it missed were the model-shift transports. This matches
the distinction between covariate and regression shift in optimal ridge
regularisation @patil2024: an error that does not shrink with $n$ is not bounded
by a penalty that does.

== A pace-of-ageing clock

For DunedinPACE, the transported correction fitted on 40 samples of GSE40279 was neutral (median +0.3 p.p., harmful in 51% of draws), although shuffled coefficients still did +4.9 p.p. of damage: the noise component was present, but real
coefficients removed enough genuine signal to offset it. Fitted instead on the
fifth cohort, it was neutral again (−1.5 p.p., 42%), so this looks like a property of
the clock rather than of the fitting cohort. What
replicated was the rest: at matched sizes, 6 of 12 directed pairs were harmful
unpenalised and none at $alpha = 3$, and a correction fitted on controls left
10.5 more points of composition in arthritis cases than in other controls.

== Replication on another array and ancestry

Fitted on 40 samples of GSE132203 (EPIC, mostly African American) and transported
to the four 450k cohorts, the correction was harmful for the age clocks in 95% of draws (median +24.3 p.p.; Levine 2018 91%, Horvath 2018 99%), with shuffled coefficients at +16.4. Across the 8 directed pairs involving this cohort, 7 of 24 (pair × clock) cells
were harmful unpenalised and none at $alpha = 3$; the 24 count DunedinPACE
alongside the two age clocks (6 of its 16 age-clock cells were harmful).

== A fitting cohort four times larger

Fitted on GSE55763 and transported to the other five cohorts (13 age-clock
cells), the correction was harmful at $n = 40$ in 83% of draws (median +7.6 p.p.).
Between $n = 656$ and the full 2,639 the composition it left fell only from 1.9
to 1.3 points, a ratio of 0.72 against the 0.25 that $1\/n$ noise predicts; the
shuffled reference fell from 0.5 to 0.1 (log-log slope −1.03). The floor
therefore persists on median. It is not uniform: at full size it was +4.6 points
into the arthritis cohort and +3.0 into GSE40279, but zero into GSE50660 and
GSE61151, and it exceeded the shuffled reference in 8 of 13 cells, below our
pre-set 9. Model shift is a property of the pair of cohorts. With this many
fitting samples the unpenalised correction helped in 10 of 13 cells, and the
fixed penalty cost benefit (median −3.0 p.p. against −4.0) while removing the three harmful cells. Restricting GSE40279 to the fitting cohort's age range (24–75
years) left its floor unchanged (+3.5 to +3.3 points for Horvath 2018), so age
extrapolation does not explain it. The floor is confined to the cohort's
European-ancestry half (+4.3 and +3.9 points, against −0.4 and −0.2 in its
Hispanic half), which in this cohort was also processed on separate plates;
centring by plate did not remove it.

== Saliva

#figure(
  image("figures/en/fig6_saliva.png", width: 95%),
  caption: [*Immune-fraction slope of each clock in four saliva cohorts.* Clock
  age regressed on chronological age and immune fraction; slope per 10
  percentage points, 95% CI. GSE149747 at baseline.],
) <fig6>

Before correction, composition accounted for 10.8–51.6% of age-acceleration
variance for the age clocks in saliva, and 42–56% for DunedinPACE. Fitted on 40
samples and transported between saliva cohorts, the correction was harmful in 67% of draws (median +20.9 p.p.; shuffled +6.4). At matched size ($n = 259$), 7 of 8
(pair × clock) cells were harmful unpenalised, 7 of 8 at $alpha = 3$ and 6 of 8
with the fading penalty; for Horvath 2018 into GSE78874 the unpenalised
correction left 132 and 168 p.p. of composition where it had found 16 — eight to ten times the signal it was meant to remove. Within each
cohort, fitting on one random half and applying to the other was beneficial in 5 of 6 age-clock cells (−15 to −48 p.p.); the exception had 132 fitting samples and a
small initial signal. The failure is therefore in the transport. Every saliva pair
also crosses array and preprocessing, so technical and biological differences
cannot be separated, but quantile-normalising GSE78874 to the EPIC distribution
left 7 of 8 cells harmful, and so did a three-type fit whose composition matrix
is well conditioned (condition number ≈ 1, against 872–1,317 for nine types).
The cause is on the dominant axis (@fig6): adjusted for age, Horvath 2018 changed
by +1.0 and +1.9 years per 10 points of immune fraction in the EPIC cohorts and by
−0.9 in GSE78874, and Levine 2018, which keeps its sign in all four cohorts,
ranged from −0.2 to −4.7. Shrinkage moves a transported coefficient toward zero,
so it can only limit a correction that disagrees with the target this much: at
$alpha = 3$ the worst cell fell from +116 to +28 p.p., and Levine 2018's cells, sign intact, to a median of +3.2 rather than to zero. In a fourth EPIC saliva
cohort from another group (GSE149747, 44 adults at baseline) the slope was −0.82
(95% CI −1.71 to 0.07), on the side of the 450k cohort rather than the other EPIC
cohorts, which argues against the array as the explanation; by our pre-set rule
this comparison was inconclusive.

Pooling studies did
not rescue transport: fitted on the other saliva cohorts with a fixed effect per study and
applied to the one left out, the correction was harmful in 6 of 8 cells (5 of 8
at $alpha = 3$), because the pooled slope takes the sign of the pool's majority.

Blood differs in how far the cohorts disagree relative to the effect itself. In
a random-effects fit of the per-cohort slopes, the between-cohort standard
deviation over the absolute mean effect was 0.09–0.55 in the six blood cohorts
(naive CD8 and neutrophil axes, both clocks) and 0.91 for Levine 2018 and 5.34
for Horvath 2018 in saliva; $I^2$ was 43–84% against 97%. Sign reversal is the
extreme of that spread, and it is not the whole of it: with GSE78874 normalised,
Horvath 2018's saliva slope is −0.20 (95% CI −0.46 to +0.06), no longer clearly
opposite to the EPIC cohorts, and 7 of 8 cells remain harmful; Levine 2018 never
reverses and is still harmful at $alpha = 3$ in 3 of 4 cells.

== What it does to a reported association

Composition left is not what a study reports. For four exposures — smoking
(GSE50660), rheumatoid arthritis (GSE42861), inflammatory bowel disease and
oesophageal cancer (the two EPIC saliva cohorts) — we re-estimated the exposure
coefficient in (clock age) ~ age + exposure under three corrections: none, fitted
within the target, and transported. Taking the within-cohort correction as the
comparator, a correction fitted on 40 samples elsewhere moved the reported effect
a median of 0.91 years, against 0.55 years for applying no correction at all; at
$alpha = 3$, 0.48; at full fitting size, 0.25. Of 24 (target × source × clock)
cells, 13 moved by more than a year and 3 reversed sign; the largest was 4.19
years, for Levine 2018 carried from GSE132203 into the arthritis cohort, whose
within-cohort effect is +0.04 years. At $alpha = 3$, 6 of 24 still moved by more
than a year.

This metric also shows a tail that composition left, bounded above by
construction, cannot: at $n = 40$, 1.2% of 720 draws landed more than 10 years
from the within-cohort estimate, the worst at 37.6 years. Solved as exact least
squares, without dropping near-zero singular values, 10.1% exceeded 10 years, so
the tail's size is solver-dependent while its existence is not; it is
concentrated in saliva, where the composition matrix is nearly singular. At
$alpha = 3$ no draw exceeded 10 years.

== Robustness to the reference panels

Repeated on the published libraries, with nothing else changed, the blood results
hold: the correction fitted on 40 samples of GSE40279 was harmful in 77% of 600
draws (median +9.1 p.p., against 88% and +16.1 with the panels built here); over
the 60 (pair × clock) cells at matched size, 17 were harmful unpenalised and none
at $alpha = 3$; and fitted on all of GSE55763 the median composition left was
+1.9 p.p. The two twelve-type panels agree on what they both estimate
(neutrophils $r = 0.985$–0.997, naive CD8 $r = 0.852$–0.919 across the six
cohorts).

In saliva, the independent panel correlates 0.951–0.997 with the EpiDISH immune
fraction and sees almost the same composition signal before correction as
EpiDISH's own immune column (for example +1.5 against −0.1 p.p. in GSE232891 for
Levine 2018, +50.6 against +50.6 in GSE78874). What differs is the number of
axes, not their source: the three-type measurement sees +18.5 and +31.8 p.p. in
GSE232891 where the immune axis alone sees about zero, so in the EPIC cohorts
most of the composition signal in age acceleration lies on the
epithelial/fibroblast axes. Scored on the immune axis alone with the independent
panel, 6 of 8 cells were harmful unpenalised and 5 of 8 at $alpha = 3$, against 7
and 7 with the shared measurement.

= Discussion

Within-cohort composition adjustment cannot be validated inside the cohort; this
follows from least squares, not from data. Transported, the adjustment's small-
sample damage is mostly estimation noise — shuffled coefficients do nearly as
much of it — while its large-sample damage comes from differences in the
composition effect between cohorts. That second component is invisible without
the target's own coefficients, and a target large enough to estimate them could
simply be adjusted within itself.

In the units a study reports, the cost of borrowing a correction is about a year on a disease or exposure effect, and more than four years in the worst configuration we found — enough to change what a paper concludes.

For practice this suggests four things. Where the target cohort is large enough,
fit the adjustment within it. In blood, where coefficients must be transported,
penalise them with a fixed, substantial penalty rather than one chosen by
cross-validation on the source cohort or one that fades with $n$; with thousands
of fitting samples this costs about a point of benefit. In saliva, and plausibly
in any tissue where composition dominates, do not transport: no penalty and no
pooling made it safe. And treat any pre-transport diagnostic, including the
transport index, as a ranking of risk rather than a guarantee. The same caution
applies inside a single study: a correction fitted on controls and applied to
patients is a transport, and here it changed the estimated disease effect by up
to a factor of two.

The penalty's reach follows from how far cohorts disagree relative to the
effect being corrected. Shrinking toward zero replaces a transported coefficient
with a smaller one, which is close to every cohort's own coefficient when the
between-cohort spread is a fraction of the average effect, as in blood, and close
to none of them when the spread is as large as the effect, as in saliva. A
reversal of sign is the extreme case, not a separate mechanism, and we would
expect the penalty to fail wherever that ratio approaches one.

Transport is already published practice outside blood: a saliva adaptation of a
blood clock fits composition terms on about 960 pooled samples and applies them
to held-out studies @galkin2021. It was judged by accuracy against chronological
age, which cannot show composition left in acceleration; in our saliva cohorts,
pooling did not prevent harm. Composition-dependent heterogeneity within a
saliva cohort has also been reported @chan2026.

= Limitations

Six adult whole-blood cohorts, five on the 450k array and one on EPIC, two
defined by disease or exposure, and four adult saliva cohorts whose transport
pairs all cross array and preprocessing; other tissues and children are untested. The
panels built here were checked against published libraries, and the saliva
measurement against an independent panel, without changing any conclusion.
Medians at small fitting sizes are unstable across independent sets of 30 draws
(+11.8 to +18.9 p.p. at $n = 40$; the pooled figure from 100 draws is +16.1), and
the small-sample harm itself depends on the clock and the
fitting cohort (neutral for DunedinPACE from GSE40279). The penalty value is specific to this panel and these
clocks. The decomposition assumes a linear composition effect; squared terms for the four largest components add a median of 0.002 to within-cohort $R^2$ (significant in 4 of 12 cohort × clock cells), and a quadratic correction transports no better than the linear one (22 against 23 harmful cells of 60, same median). Sensitivity
scoring on twelve types shares the fitting panel. None of this bears on whether
epigenetic clocks measure biological ageing; it concerns one correction applied
to them.

= Data and code availability

All series are public (GSE40279, GSE61151, GSE50660, GSE42861, GSE132203,
GSE55763, GSE232891, GSE232332, GSE78874, GSE149747, GSE35069, GSE167998). Analysis code, the stage-by-stage record including every overturned
conclusion, and figure scripts are at
#link("https://github.com/KTHimiko/clock-lab")[github.com/KTHimiko/clock-lab],
private at the time of writing and available from the author on request.

#bibliography("refs.bib", title: "References", style: "nature")
