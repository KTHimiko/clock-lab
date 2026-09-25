# clock-lab — what has been found

## The question

Every longevity company sells a "biological age" test built on DNA methylation.
Blood composition shifts with age, and each cell type carries its own methylation
pattern. **Does a clock read how old the cells are, or who is in the sample?**

## The answer, after forty-three stages

**Both, and the proportions matter more than either camp says.**

In sorted cells the effect is enormous: up to 35 years between fractions of one
man's blood drawn on one day, replicated in a second cohort on a second array.
In actual blood from actual people it is smaller but far from small. At twelve
cell types — splitting naive from memory lymphocytes, which is where the effect
lives — composition explains **10% to 32%** of age acceleration, landing on the
published figures. This project's earlier answer of 0.9–2.7% was a floor
produced by a six-type panel and a denominator the literature does not use.

**It can be designed away, at least partly** — though this project spent three
stages concluding otherwise. Excluding CpGs that correlate with naive-CD8
identity cuts a clock's cell-type displacement by 41% while slightly improving
its accuracy on an external cohort. What defeated the earlier attempts was
filtering on overall cell-type variance instead of on the naive-versus-memory
axis, where the effect actually lives.

It can also be *subtracted* after the fact — inside the cohort where the
correction is fitted, which is how it is almost always used, and where it does
what it claims. **Carried to another cohort, it can add the confounding it was
meant to remove.** Fitted on forty samples and transported, it leaves +14.2% of
age-acceleration variance in composition signal and is worse than doing nothing
in 88% of draws; measured on the naive/memory axis, which the six-type panel used
throughout cannot see, the same transport leaves +30.6%, so every curve here is
conservative. The median turns roughly neutral between 160 and 240 fitting
samples.

**The harm has two components** (stages 24–25):

- **estimation noise**, amplified where the target cohort varies in directions
  the fitting cohort barely did — the standard excess-risk term for least squares
  under covariate shift, (σ²/n)·tr(Σ_fit⁻¹ Σ_test). It falls as 1/n, and
  coefficients carrying no information reproduce it almost point for point:
  +17.4 p.p. at n = 40 (stage 24b, 100 draws). Most of the small-n harm *is*
  noise.
- **model shift**: the composition effect itself differs between cohorts, so
  even a perfectly estimated correction from one damages another. It does not
  shrink with fitting size, and it is behind the worst transport in the project,
  GSE61151 → the rheumatoid arthritis cohort at +15.9%.
  Stage 35 pushed the fitting cohort to 2,639 samples. The floor stays on
  median (1.9 → 1.3 points while the noise reference falls to 0.1), but it
  belongs to the pair: +4.6 into the arthritis cohort, zero into two others.

**No quantity computed in advance certifies a transport as safe** (stages 26,
28). The transport index sees only the first component; it adds modest
information beyond fitting size (within-n ρ ≈ 0.36, block-permutation p = 0.036),
and counted per clock, 12 of 40 transports below its old "safe" threshold were
harmful. Of the 73 transports a closed-form net-damage predictor
calls safe, 30 are harmful. The second component needs the target's own coefficients, and a target
large enough to estimate them does not need a transported correction.

**In years, the units a study reports** (stage 43): a correction fitted on forty
samples elsewhere moves the estimated effect of smoking or disease 0.91 years
from what the target cohort's own correction gives, against 0.55 years for not
correcting at all; 13 of 24 cells move by more than a year and 3 flip sign. A
penalty brings it to 0.48.

**What works is penalising** (stage 27). Ridge shrinks the coefficients toward
"no correction", which bounds the damage from noise and from a wrong β alike.
Counted per clock, harmful transports go from 50% unpenalised to 3% at α = 3 —
one cell left, at +0.2% — and to none at α = 10, which gives up most of the
benefit where none was at risk. Cross-validation on the fitting cohort does not
size the penalty on the cohort that fails. It is a small-n safeguard: fitted on 2,639
samples, the unpenalised correction is the better one on median (−4.0% against
−3.0%, stage 35). A textbook penalty that fades as 1/n recovers that
benefit but leaves 11 of 72 cells harmful at matched n, against 1 of 72 (stage
36): model shift does not fade, so the penalty bounding it cannot either.
**In saliva it does not hold at all** (stage 38): transported between saliva
cohorts, the correction was harmful in 7 of 8 cells even at α = 3, and up to
+152% unpenalised, while the same correction works inside each cohort. The
penalty is a blood result. Neither normalisation nor a three-type fit changes it: the
clock's slope on the immune fraction changes sign between saliva cohorts (stage
39), and pooling studies before fitting does not rescue it (stage 41: 6 of 8
harmful). Stage 42 names the mechanism: shrinkage helps when a
cohort's own coefficient is close to the shrunk one, which holds when the spread
between cohorts is small next to the average effect (τ/|mean| 0.09–0.55 in blood)
and fails when it is not (0.91–5.34 in saliva, where one clock's slope also
reverses sign).

And flattening a clock against composition appears to be **free**: the clock
built here that reads 41% less composition detects rheumatoid arthritis exactly
as well as the standard one, with a paired confidence interval that excludes any
meaningful loss. The worry that composition was half the signal is not supported
for that outcome.

What is *not* true is that the published clocks are badly designed. Their
apparent advantage and their apparent deficit both turned out to be artefacts of
what they were tested against — and when a family of ridge clocks is trained on
a matched age range, it beats them, by about a third, on exposure per unit of
age response. The floor it hits, 3.9 years, is the biology.

## The design that answers it

Reinius et al. (GSE35069) split one blood draw from each of six men into ten
fractions. Same person, same day, ten measurements.

**Any difference in epigenetic age between those fractions cannot be ageing.**

## Stage 2 — the result

Epigenetic age by fraction, averaged over the six donors:

| fraction | Horvath 2013 | Hannum 2013 | Levine 2018 | Horvath 2018 |
|---|---|---|---|---|
| CD4+ T cells | 34.4 | 41.0 | 17.3 | 39.5 |
| CD19+ B cells | 35.8 | 51.2 | 30.5 | 26.4 |
| Neutrophils | 36.1 | 47.5 | 30.9 | 32.8 |
| CD14+ Monocytes | 38.3 | 48.7 | 43.2 | 37.7 |
| **Whole blood** | **39.7** | **45.4** | **31.3** | **37.4** |
| CD8+ T cells | 41.2 | **24.3** | **6.2** | 34.7 |
| Eosinophils | 46.8 | **59.7** | 33.2 | 36.0 |
| CD56+ NK cells | 47.2 | 51.3 | 30.2 | 40.5 |

**Hannum reads 24 years for one man's CD8+ T cells and 60 for his eosinophils —
from the same tube of blood.**

### How much of the signal is not ageing

| clock | spread across fractions, within one donor | range | spread across donors | ratio |
|---|---|---|---|---|
| Horvath 2018 | 4.5 yr | 15.1 yr | 13.8 yr | **0.33** |
| Horvath 2013 | 5.1 yr | 16.0 yr | 13.7 yr | **0.37** |
| Hannum 2013 | 9.9 yr | 35.4 yr | 13.6 yr | **0.73** |
| Levine 2018 | 10.7 yr | 37.0 yr | 13.9 yr | **0.77** |

The left column is measured **inside one person on one day** — it cannot be
ageing. The right column mixes real age with everything else.

**For Hannum and Levine, the variation between cell types within one person is
roughly three-quarters of the variation between different people.** The Horvath
clocks are about half as exposed, which is consistent with Horvath 2013 having
been fitted across many tissues rather than on blood alone.

### The sanity checks, fixed before the result was read

1. probe coverage — 100% for all four clocks
2. design — exactly 6 donors × 10 fractions, 60 samples
3. **PBMC must sit closer to whole blood than a random pair of fractions**, since
   PBMC is a subset of whole blood. Observed 2.16 years against 5.26 for a random
   pair. This is the check that would catch a donor mix-up, and it held.

## Stage 3 — the obvious explanation is wrong

If a clock varies between cell types, the natural guess is that its probes are
cell-type markers. Decomposing the variance of all 454k complete probes into
donor, cell type and residual, then asking where each clock's probes sit:

| clock | probe cell-type variance (weighted by \|coefficient\|) | percentile | spread measured in stage 2 |
|---|---|---|---|
| Horvath 2018 | 0.557 | 78th | **4.5 yr** — the least exposed |
| Hannum 2013 | 0.561 | 66th | 9.9 yr |
| Horvath 2013 | 0.467 | 63rd | 5.1 yr |
| Levine 2018 | **0.395** | 53rd | **10.7 yr** — the most exposed |

**Spearman −0.40: the wrong direction.** Horvath 2018 draws on the most
cell-type-variable probes and is the least vulnerable; Levine draws on the least
and is the most.

With four clocks this is not a test — four points give a perfect correlation by
chance once in twenty-four. But the ordering is inverted, not merely weak, and
that points somewhere specific.

**What matters is not the size of each probe's cell-type effect but its
alignment with the coefficients.** A clock age is a weighted sum; if the
cell-type shifts across its probes point in directions unrelated to the weights,
they cancel. Magnitude without alignment buys nothing, and alignment without
magnitude is enough.

That is the next thing to measure.

### The sanity check that stopped the first run

The variance components did not sum to one. The identity
`SS_total = SS_donor + SS_fraction + residual` holds for a balanced complete
design and breaks as soon as a cell is missing — group means then come from
different subsets. With 31,442 missing values scattered through the matrix, the
parts exceeded the whole. Restricting to the 454k complete probes restored it.

**The check existed because this stage was written expecting to be wrong
somewhere.** It was.

## Stage 4 — and the alignment hypothesis fails too

Stage 3 proposed that what exposes a clock is not how cell-type-variable its
probes are but whether those shifts line up with its coefficients. Measured
against two permutation nulls — one that shuffles which weight goes with which
probe, keeping both multisets and destroying only the pairing; one that draws
random probe sets of the same size from all 485k:

| clock | observed | shuffle null | alignment | p(shuffle) |
|---|---|---|---|---|
| Horvath 2013 | 0.168 | 0.327 | **0.5×** | 0.85 |
| Hannum 2013 | 7.303 | 4.838 | 1.5× | 0.13 |
| Levine 2018 | 6.902 | 8.750 | **0.8×** | 0.62 |
| Horvath 2018 | 0.144 | 0.306 | **0.5×** | 0.93 |

**No clock shows alignment above chance.** Every p-value is above 0.09, and two
clocks sit *below* their own null. Shuffling which coefficient goes with which
probe does not reduce the cell-type displacement — it slightly increases it.

*(The Horvath clocks report in transformed units rather than years, so their
`observed` and `shuffle` columns are not comparable with the linear clocks'. The
alignment ratio is unit-free and is.)*

### What that means, and it is stronger than the hypothesis it replaces

The cell-type displacement is **exactly what a random pairing of those weights
with those shifts produces.** The clocks are not reading cell identity through
some unlucky probe choice that better selection could avoid. They are
accumulating whatever cell-type variation exists across their probes, in
proportion to the size of their coefficients.

Levine carries a total absolute weight of 3,490 against Hannum's 628, and its
null is correspondingly larger — 8.75 years against 4.84. The differences
between clocks look like differences in **gain**, not in exposure to cell
identity.

That proposed mechanism rests on two linear clocks, where a correlation is
±1 by arithmetic and means nothing. It is an argument consistent with the
numbers, not a measurement. **The measured result is the null one: there is no
alignment to find.**

## Stage 5 — it replicates, in other people, on another array

Everything above rests on six men, one lab, one 450K run. GSE110554 (Salas et
al.) is a different cohort, the EPIC array, and magnetic-bead separation:
37 purified samples, ages 19 to 59.

**The design is not the same, and that changed the statistic.** GSE35069 is
paired — ten fractions per donor, so the within-person spread is measured
directly. GSE110554 is not: each purified sample is one cell type from one
*different* donor. The paired statistic does not exist here. So the replication
is of the claim rather than the arithmetic — if cell type shifts epigenetic age,
each cell type must carry its own systematic offset from its donor's real age,
and those six offsets must spread further apart than shuffling the labels
produces.

| clock | EPIC coverage | spread across cell types | shuffle null | ratio | p | technical noise |
|---|---|---|---|---|---|---|
| Horvath 2013 | 94.6% | 2.79 yr | 1.52 | 1.8× | 0.003 | 1.13 yr |
| Hannum 2013 | 91.5% | **10.24 yr** | 4.04 | 2.5× | <0.0001 | 2.05 yr |
| Levine 2018 | 100% | **15.73 yr** | 6.01 | 2.6× | <0.0001 | 5.74 yr |
| Horvath 2018 | 100% | 4.90 yr | 2.03 | 2.4× | <0.0001 | 2.79 yr |

**All four clocks replicate**, and the ordering of magnitudes is the stage 2
ordering: the Horvath clocks least exposed, Hannum and Levine two to three times
more.

### The two checks that could have killed it

**The ruler.** Th2535-1 and Th2535-2 are the same donor, same cell type, same
plate — technical replicates, and they only became visible after the loader fix
below. Their difference *is* the measurement noise: 1.1, 2.1, 5.7 and 2.8 years.
Every clock's cell-type spread clears its own noise floor, Horvath 2013 by the
narrowest margin (2.5×).

**The confound.** A clock compresses the age scale, so an older donor gets a
negative age gap for free — every clock here has a slope of −0.11 to −0.23
years per year. And the cell types are *not* age-balanced: B cells average 40.5
years, neutrophils 26.8, a 13.7-year gap, and those two sat at opposite ends of
the raw result. That is exactly what the confound would manufacture. Residualising
each age gap on chronological age barely moves anything — 3.36 → 2.79,
10.33 → 10.24, 15.87 → 15.73, 5.45 → 4.90 — and every p-value survives. The
effect is not regression to the mean.

### Do the two studies rank the cell types the same way?

The harder test: not "is there an effect" but "is it the same effect". Six cell
types map between the studies, so the null is exact — all 720 orderings.

| clock | Spearman ρ | p |
|---|---|---|
| Levine 2018 | **+0.94** | 0.008 |
| Horvath 2018 | +0.77 | 0.051 |
| Hannum 2013 | +0.60 | 0.121 |
| Horvath 2013 | +0.26 | 0.329 |

**The clocks with the larger effects reproduce their ordering; the ones with
small effects do not** — which is what a real signal buried in noise looks like,
and the opposite of what a cohort artefact would do.

The single most reproducible number in this project is CD8+ T cells. Hannum
reads them 22 years young in the first cohort and 23 in the second; Levine, 22
and 28. Two cohorts, two array versions, two sorting protocols.

**So the stage 2 finding is biology, not GSE35069.**

## Stage 6 — the first prediction, and what it cost to read it honestly

Stages 2 and 5 describe. This stage predicts.

GSE110554 also contains 12 samples that are nobody's blood: DNA mixtures
reconstructed by pooling purified cell types in **annotated proportions** — six
by "method A" (balanced), six by "method B" (neutrophil-heavy, like real blood).
A clock score is a weighted sum over probes, so it is linear in the betas; if a
mixture's methylation is the proportion-weighted average of its ingredients, its
clock score must be too:

> predicted(m) = Σ<sub>c</sub> p(m,c) × mean_score(c)

The cell-type means come from the 37 purified samples of stage 5; the proportions
come from the GEO annotation. **Zero free parameters.**

**The premise held.** Every mixture's betas sit closer to its own
proportion-weighted prediction than to another mixture's — DNA mixing really is
linear in beta, and the proportions are not mislabelled.

**The prediction, judged within method** (method is nearly collinear with
composition — 17% neutrophils against 69% — so a correlation across all 12 could
be a two-group batch effect wearing a costume):

| clock | r (method A) | p | r (method B) | p |
|---|---|---|---|---|
| Horvath 2013 | +0.22 | 0.37 | +0.14 | 0.39 |
| Hannum 2013 | +0.61 | 0.075 | +0.65 | 0.093 |
| Levine 2018 | +0.38 | 0.23 | **+0.90** | **0.008** |
| Horvath 2018 | +0.06 | 0.46 | −0.38 | 0.72 |

Seven of eight positive, one convincing, most not. A combined permutation test
that shuffles proportion assignments within method — preserving the fact that
four clocks over twelve mixtures are nowhere near independent — gives mean
r = +0.32 against a null of 0.00 ± 0.22, **p = 0.067.**

On its own that is a shrug.

### Why it is a shrug, and the statistic that says so

The predicted spread across six mixtures of one method ranges from **0.13 to
2.75 years** depending on the clock. These twelve mixtures barely differ in a
way any clock could see. So the question is not "did the prediction work" but
"could it have":

| ordered by predicted spread | spread | r achieved |
|---|---|---|
| Horvath 2018, B | 0.13 yr | −0.38 |
| Horvath 2013, B | 0.13 yr | +0.14 |
| Horvath 2013, A | 0.42 yr | +0.22 |
| Horvath 2018, A | 0.53 yr | +0.06 |
| Hannum 2013, B | 0.78 yr | +0.65 |
| Levine 2018, B | 1.30 yr | +0.90 |
| Hannum 2013, A | 2.70 yr | +0.61 |
| Levine 2018, A | 2.75 yr | +0.38 |

**Spearman +0.71, p = 0.029 against the exact null over all 40,320 orderings.**

The predictor in that table — how far the *prediction* spreads — is computed
without looking at an observed value even once, so it cannot have been inflated
by one. **The composition prediction succeeds in proportion to how much
composition there was to predict.** That is what the hypothesis says should
happen, and it is not what a coincidence looks like.

*(The eight cells share four clocks and twelve mixtures, so that p is optimistic.
The direction is the claim; the p-value is decoration.)*

### The honest failure: a power check written after the fact

The power diagnostic above was added **after** seeing the result. It does not
rescue anything — it is computed from the prediction alone — and it could have
been pre-specified. It was not. Recorded as a diagnostic, not a criterion.

It also corrected its own first yardstick. Measuring the predicted spread against
stage 5's technical-noise floor made seven of eight cells look hopeless, but that
floor came from replicate *purified CD4T pellets*, not mixtures: Levine's
declared noise is 5.74 years while its method-B mixtures span 2.03 years in
total, which is impossible if the noise were real. The ratio that does not depend
on it is predicted spread over observed spread, and that is the column the table
uses.

### What is left over, and it is not composition

Composition is already subtracted in the prediction, so the residual should sit
at zero. It does not:

| clock | residual, method A | residual, method B | difference |
|---|---|---|---|
| Horvath 2013 | +2.1 yr | +7.5 yr | +5.4 |
| Hannum 2013 | −3.7 yr | +3.7 yr | +7.4 |
| Levine 2018 | +0.5 yr | +5.8 yr | +5.3 |
| Horvath 2018 | +2.2 yr | +9.4 yr | +7.2 |

All four clocks read method-B mixtures **five to seven years older than their
annotated composition allows.** Two readings, and this cohort cannot separate
them, because method and neutrophil fraction are nearly the same variable here:
either the two reconstruction protocols produce DNA the array reads differently
— nothing to do with clocks — or purified neutrophils are not what neutrophils
look like inside a mixture, which would need the neutrophil mean to move about
12 years. **Open.**

## Stage 7 — and now the number that actually matters

Nobody buys a sorted-cell test. They send a tube of blood. So: estimate each
donor's cell composition from their own methylation, and ask how much of what
the clock calls their biological age is that composition.

The panel is six purified types from GSE35069, and it is validated the hard way
— by recovering the **known** proportions of GSE110554's twelve reconstructed
mixtures, a different cohort on a different array:

**r = 0.985 across all 72 cell-by-mixture values, mean absolute error 0.021,
worst per-type bias 0.032.** Per type, every correlation is above +0.93.

### What that validation does not cover, and the check I failed to specify

A reconstructed mixture is by construction a linear combination of sorted cells
— far more like the panel than real blood is. Applied to 188 real whole-blood
samples, the estimates drift from clinical reference ranges:

| type | median estimate | adult range | |
|---|---|---|---|
| CD4T | 0.165 | 0.07–0.15 | above |
| CD8T | 0.026 | 0.03–0.08 | below |
| B cells | 0.010 | 0.01–0.04 | ok |
| NK | 0.088 | 0.02–0.06 | above |
| Monocytes | 0.058 | 0.02–0.08 | ok |
| Neutrophil-like | 0.639 | 0.40–0.75 | ok |

**My pre-specified physiology check covered neutrophils only** — 64% of blood
and by far the easiest to get right. The smaller types had no check at all, and
they are where it went wrong. The pattern is diagnostic: surplus CD4T and NK,
deficit CD8T, the pairs that compete for the same lymphoid probes. The
proportions are collinear, so the **joint** fit can be right while the split
between neighbouring types is wrong.

That determines what may be read below. The joint test survives it; the
per-type coefficients do not, and are not used as measurements of a type.

### The result

Three nested regressions per clock. Composition enters *after* chronological
age, so the increment is what blood count explains **beyond** what age already
did — which is the honest quantity, because composition itself changes with age
and that share belongs to age.

| clock | R² age | R² age + composition | increment | null | p | in years (SD) | span |
|---|---|---|---|---|---|---|---|
| Horvath 2013 | 0.722 | 0.741 | 0.020 | 0.008 | 0.025 | 1.36 yr | 6.9 yr |
| Hannum 2013 | 0.803 | 0.830 | **0.027** | 0.005 | **0.0002** | 1.47 yr | 8.4 yr |
| Levine 2018 | 0.632 | 0.659 | **0.027** | 0.010 | 0.018 | 1.63 yr | 8.9 yr |
| Horvath 2018 | 0.791 | 0.800 | 0.009 | 0.006 | 0.17 | 0.89 yr | 6.0 yr |

**Real, significant in three clocks of four — and an order of magnitude smaller
than the sorted-cell experiment implies.**

> **Stage 13 corrected both halves of that sentence.** The increment is reported
> here as a share of total clock-age variance; the literature reports a partial
> R² inside age acceleration, which is roughly four times larger for the same
> data. And this six-type panel cannot see the naive-versus-memory split, where
> the largest cell-type age differences live, so it measured a floor rather than
> the effect. At twelve types the numbers roughly double again. See stage 13.

That is not a contradiction, it is the whole point. Stage 2 found 35 years
between a man's CD8+ T cells and his eosinophils, but nobody's blood is 100%
CD8+ T cells. Across 188 real people the composition term moves a reading by
about **1 to 1.6 years of standard deviation, and 6 to 9 years between the two
most extreme donors in the cohort.** That is the size of the thing when it
reaches a customer.

### The reconciliation, which came out weak, and is reported weak

Do the whole-blood coefficients recover the ordering measured on sorted cells?

| clock | ρ vs stage 2 | p | ρ vs stage 5 | p |
|---|---|---|---|---|
| Horvath 2013 | +0.77 | 0.051 | −0.31 | 0.75 |
| Hannum 2013 | +0.31 | 0.28 | +0.03 | 0.50 |
| Levine 2018 | +0.49 | 0.18 | +0.71 | 0.068 |
| Horvath 2018 | +0.60 | 0.12 | +0.37 | 0.25 |

Seven of eight positive, **not one below 0.05.** This is exactly what the
collinearity above predicts: individual coefficients are not identified well
enough to carry an ordering, even where the joint effect is solid. Reported as
the weak result it is, and not rescued.

*(Circularity, declared: the panel is built from GSE35069, the same data as
stage 2. The stage 5 column, from another cohort and another array, does not
have that problem — and it is the weaker of the two.)*

## Stage 8 — the stage 4 argument, finally measured

Stage 4 concluded that the cell-type displacement of a clock is exactly what a
random pairing of its weights with its probes' cell-type shifts produces, and
that what separates clocks is therefore **gain** rather than probe choice. That
was recorded as an argument, not a measurement, because it rested on two linear
clocks — where a correlation is ±1 by arithmetic.

Four clocks cannot test it. Forty-five can. This stage trains a family of ridge
age predictors on GSE61151 whole blood, sweeping probe count (100 → 458,674) and
penalty (10⁻³ → 10⁶), and measures each on two axes: **accuracy**, 5-fold
cross-validated, and **displacement**, the stage 2 statistic, on GSE35069
purified cells the training never touches.

### Stage 4 was right about the mechanism and wrong about the quantity

| summary of β | Spearman with displacement | log-log slope | residual scatter |
|---|---|---|---|
| L1 — Σ\|β\|, what stage 4 used | +0.747 | 0.82 | 1.80× |
| **L2 — √Σβ²** | **+0.968** | 0.79 | 2.05× |
| L2 weighted by each probe's cell-type shift | +0.967 | 0.81 | 2.01× |

Spreading the same total weight over more probes leaves L1 untouched and drops
L2. Displacement follows L2, not L1 — so **dilution is a defence**, and stage 4
named the wrong norm. (All three leave about a two-fold residual scatter around
a log-log line, so no single summary of β captures everything; the rank
correlation is what separates them.)

### The frontier

The published clocks, recomputed on the same six fractions with the same
statistic — *not* imported from stage 2, whose number is over ten fractions and
is not comparable:

| clock | displacement | family configurations that are flatter **and** reach r > 0.75 |
|---|---|---|
| Horvath 2013 | 5.52 yr | 10 |
| Horvath 2018 | 5.77 yr | 11 |
| Hannum 2013 | 10.96 yr | 22 |
| Levine 2018 | 13.61 yr | 23 |

Along the family's own Pareto edge, r = 0.888 costs 9.21 years of displacement,
while r = 0.822 buys it down to 5.45 and r = 0.782 to 4.10.

> **Stage 9 overturned the practical half of this.** The frontier above mixes a
> scale-invariant axis (correlation) with a scale-dependent one (displacement in
> years), and that lets a model buy flatness by shrinking itself. The diluted
> configurations were not reading less composition — they were responding less to
> everything, age included. What survives from this stage is the mechanism:
> displacement tracks L2, not L1. What does not survive is "dilution is a
> defence". See stage 9.

### The caveat that governs how far this goes

The two axes are not equally fair. Displacement is measured identically for both
sides on GSE35069, which neither trained on — that axis is clean. Accuracy is
not: the family is cross-validated *inside* the cohort it was trained on, while
the published clocks arrive cold. The family plays at home.

So the claim that survives is not "you can beat Hannum". It is: **displacement
halves without losing accuracy within the cohort, and what governs it is the L2
norm of the coefficients.** Whether the diluted version's accuracy transfers to
another cohort, this stage did not test.

### A check that was mis-specified, and what happened to it

Check 3 required the degenerate corner to appear — at a large enough penalty,
coefficients collapse and both axes go to zero. It had three clauses.
Displacement → 0.020 years and |r| → 0.16 passed decisively. The third,
**total absolute weight < 1.0, failed at 1.877** — because that total is a sum
over probes, so its scale rides on probe count, and 1.877 across 458,674 probes
is four millionths per probe. The threshold was meaningless when written.

It was replaced by its scale-free form, weight *per probe*, and not deleted.
Relaxing a check that has just failed is the move that should always be looked
at twice, so it is in the docstring, in the output, and here.

## Stage 9 — the external cohort, and the frontier falls over

Stage 8 ended on a caveat: its accuracy was cross-validated inside the cohort it
trained on. GSE40279 settles it — 656 whole blood samples, ages 19 to 101, a
different study and a different population, downloaded for this stage.

*(It is also the cohort the Hannum clock was trained on. Hannum's accuracy here
is a memory, not a prediction. It is reported and marked, and never used as a
bar the family must clear.)*

### The trap, which this stage walked into before climbing out

**Correlation is scale-invariant.** At the degenerate end of a ridge sweep the
coefficients collapse, predictions become the training mean plus an epsilon
multiple of a real age direction — and r stays at 0.885 while every prediction
is wrong by 13.7 years. Displacement, measured in years, collapses by the same
epsilon.

So the degenerate corner reads as *accurate and perfectly flat at once*, and it
dominated the frontier. It is not a clock; it is a constant with a rumour of a
direction.

Both axes are therefore normalised by the predictor's scale — the slope of its
prediction on real age in the external cohort — and displacement becomes years
of cell-type spread **per year of genuine age response**. Epsilon cancels.

**That normalisation kills stage 8's recommendation.** The diluted
configurations were not flat; they were small:

| configuration | raw displacement | scale | displacement per year of age response |
|---|---|---|---|
| 447,564 probes, λ=10⁻³ | 5.45 yr | 0.247 | **22.1** |
| 100,000 probes, λ=10⁻³ | 6.36 yr | 0.280 | 22.7 |
| 1,000 probes, λ=1 | 9.21 yr | 0.574 | 16.0 |
| **100 probes, λ=1** | 8.30 yr | 0.564 | **14.7** |

And on that measure the published clocks are better than anything the family
produced:

| clock | r on GSE40279 | scale | displacement per year of response | family configurations that beat it on both |
|---|---|---|---|---|
| **Horvath 2013** | +0.918 | 0.798 | **6.91** | **0** |
| **Horvath 2018** | +0.940 | 0.787 | **7.33** | **0** |
| Hannum 2013 | +0.946 | 0.843 | 12.99 | 0 — *and it trained here* |
| Levine 2018 | +0.852 | 0.770 | 17.66 | 9 |

### Dilution does not travel either

Spearman between a configuration's L2 norm and how much accuracy it loses
leaving home: **−0.384.** Higher L2 means *less* loss — the concentrated clocks
transfer better. Dilution costs accuracy out of cohort on top of buying no real
flatness.

### What survives, and what this cost

From stage 8, the mechanism survives: displacement tracks the L2 norm of the
coefficients, not the L1 norm stage 4 named. The recommendation does not.

**There is no evidence here that a better-chosen clock can be flatter than
Horvath 2013 at equal accuracy.** The two Horvath clocks sit where nothing in a
45-configuration sweep could reach them.

One honest handicap, stated: the family trains on ages 35–83 and is tested on
19–101, so its slope is compressed by extrapolating outside its training range,
and its scale of ~0.56 against the published ~0.79 partly reflects that rather
than a defect of the approach. The clean version of this test trains on the wide
cohort and tests on the narrow one — which is the next stage, not a rescue of
this one.

## Stage 10 — it was the training range, not the design

Stage 9 left the Horvath clocks unreachable and named a suspect: the family had
trained on ages 35–83 and been tested on 19–101, so its slope was compressed by
extrapolating. This stage runs the test backwards — train on GSE40279 (656
samples, 19–101), test on GSE61151 (184, 35–83). Now the family has the wider
range and the published clocks do not.

*(GSE40279 is Hannum's training cohort, so a family trained there stands exactly
where Hannum stands, and comparing them on GSE61151 is finally fair. Both
Horvath clocks are external to both.)*

| clock | r on test | scale | displacement per year of response | family configurations better on **both** |
|---|---|---|---|---|
| Horvath 2013 | +0.849 | 0.914 | 6.0 yr | **8** |
| Horvath 2018 | +0.889 | 0.900 | 6.4 yr | **7** |
| Hannum 2013 | +0.896 | 0.875 | 12.5 yr | 13 |
| Levine 2018 | +0.795 | 0.859 | 15.8 yr | 22 |

**Best family configuration: 1,000 probes, λ=1 — r = 0.922 at 3.9 years of
displacement per year of age response.** More accurate than Horvath 2013 and
roughly a third less exposed.

*(These numbers are not comparable with stage 9's table: scale is measured on
whichever cohort is the test set, so Horvath 2013 reads 6.91 there and 6.04
here. The comparison inside each table is what carries meaning.)*

### Where the design advantage actually lives — nowhere much

| variant | r | scale | displacement per year of response |
|---|---|---|---|
| (a) Horvath's probes, Horvath's coefficients — the published clock | +0.849 | 0.914 | 6.0 yr |
| (b) Horvath's probes, coefficients refitted on GSE40279 | +0.856 | 0.884 | 5.2 yr |
| (c) 353 probes chosen by the family, fitted on GSE40279 | +0.912 | 0.971 | 4.7 yr |

Refitting Horvath's own 353 CpGs on a modern wide cohort improves both axes
(6.0 → 5.2). Letting the family pick its own 353 improves them again, and adds
0.06 of correlation. **Horvath's probes are not special, and neither are his
coefficients — what he had, in stage 9, was a training set that spanned the
ages he was tested on.**

### What this does not say

The best clock built here still displaces **3.9 years between cell fractions for
every year of genuine age response.** Reducible is not the same as eliminated;
nothing in this sweep produced a clock that reads time and not composition. And
the whole comparison rests on one test cohort of 184 people in a narrow band.

### A second check rewritten after it failed, and why that is defensible here

Check 4 first required the degenerate corner to look bad on the normalised axis.
It failed at 10.5 against a median of 11.2 — and **the premise was wrong.** As λ
grows, ridge does not converge to noise; it converges to the correlation-weighted
direction, a perfectly reasonable estimator multiplied by a vanishing constant.
Normalising by scale divides that constant out, so a collapsed configuration is
entitled to a decent displacement-per-response. What disqualifies it is that its
slope is 0.009 and it is wrong by ten years — it does not report an age at all.

The replacement is a calibration band: slope between 0.5 and 1.5, which is what
"reports years" means. **The result does not depend on the swap** — the best
collapsed configuration reaches 10.5 and the frontier sits at 3.9, so they were
never in contention.

## Stage 11 — the floor is biology

Every clock in this project, published or built, displaces several years between
cell fractions for each year of genuine age response, and stage 10's best was
3.9. Two explanations that differ in what anyone can do about them: **method**,
where ranking probes by correlation with age simply never thought to avoid
cell-type markers, or **biology**, where the CpGs that track time *are* largely
the ones that distinguish cell types.

The test: filter probes by their spread across cell types before ranking them by
age, and sweep the severity.

*(The filter is built from GSE110554's purified cells and displacement is
measured on GSE35069. Using one cohort for both would be selecting on the
evaluation set — it would lower the number without lowering anything real.)*

### The answer

| k | filter | displacement/response | matched random control | filter's edge | r | control's r |
|---|---|---|---|---|---|---|
| 1,000 | none | **4.1 yr** | — | — | +0.918 | — |
| 1,000 | keep cleanest 75% | 4.4 | 4.1 | −0.3 | +0.900 | +0.917 |
| 1,000 | keep cleanest 50% | 4.7 | 4.4 | −0.3 | +0.834 | +0.916 |
| 1,000 | keep cleanest 25% | 10.8 | 5.4 | −5.4 | +0.736 | +0.910 |
| 10,000 | none | 5.7 | — | — | +0.932 | — |
| 10,000 | keep cleanest 75% | 4.6 | 5.8 | **+1.2** | +0.910 | +0.927 |
| 10,000 | keep cleanest 50% | 5.4 | 6.2 | +0.8 | +0.844 | +0.922 |

Every filtered configuration is matched against one that discards the **same
number** of probes at random, because any pruning changes the fit and that part
is free.

**Nothing filtered beats the unfiltered floor.** The one place the cell-type
information earns its keep — a mild trim at k=10,000, worth 1.2 years over its
control — still lands at 4.6, above the 4.1 an unfiltered clock reaches without
being told anything about cell types. Push the filter harder and accuracy falls
off a cliff while displacement gets *worse*.

### Why

The global rank correlation between |correlation with age| and cell-type spread
is only **+0.13**, which would suggest the two barely touch. But a clock does
not use the average probe; it uses the top of the age ranking.

| top age probes | share in the most cell-type-variable quarter | enrichment |
|---|---|---|
| 100 | 60% | 2.4× |
| **1,000** | **67%** | **2.7×** |
| 10,000 | 53% | 2.1× |
| 100,000 | 35% | 1.4× |

Under independence it would be 25%.

**Two thirds of the best age-tracking probes in the genome sit in the quarter of
probes that vary most between cell types.**

> **Stage 14 refuted the conclusion drawn from this.** The enrichment above is
> real, but "the floor is biology" does not follow from it. This stage's filter
> was built from six cell types and could not express the naive-versus-memory
> axis — which is where the largest cell-type age differences live, and which
> turns out to be separable from age after all. Filtering on *that* axis cuts
> displacement by 41% while slightly improving accuracy. The floor was an
> artefact of filtering on the wrong axis. See stage 14.

## Stage 12 — the standard fix does not travel

This is what practitioners do: estimate blood composition from the methylation
and subtract the part of the epigenetic age it explains. On GSE40279 the model
`clock_age ~ chronological_age + composition` is fitted and only the composition
coefficients are kept — so applying them needs no age, which is the whole point
of a clock. Fitting *with* age in the model matters: composition drifts with age,
and a correction fitted without it would strip out real ageing along with the
blood count.

Then those coefficients, unchanged, are carried to GSE61151.

| clock | composition R² at home | after | **composition R² away** | **after** | chance level | p |
|---|---|---|---|---|---|---|
| Horvath 2013 | 0.0086 | 0.0000 | 0.0182 | 0.0172 | 0.0076 | 0.046 |
| Hannum 2013 | 0.0105 | 0.0000 | 0.0282 | **0.0107** | 0.0051 | 0.059 |
| Levine 2018 | 0.0353 | 0.0000 | 0.0271 | 0.0174 | 0.0100 | 0.116 |
| Horvath 2018 | 0.0020 | 0.0000 | 0.0091 | 0.0091 | 0.0058 | 0.177 |
| family, k=1,000 | 0.0004 | 0.0000 | 0.0132 | 0.0123 | 0.0043 | 0.011 |

At home the correction removes everything, which is guaranteed and proves
nothing. **Away it removes 5%, 62%, 36%, 0% and 7%** — and what is left still
sits above chance.

> **Stage 15 overturned this.** The diagnosis below — that the per-type
> coefficients are not identified — was right, and it was a statement about the
> panel rather than about correction. At twelve types the same correction
> removes most of the signal. What this stage could not see, having tested only
> one direction, is that the correction also has a failure mode worse than doing
> nothing. See stage 15.

It is at least cheap: correlation with age is unchanged to three decimals for
every clock, and the mean error moves by at most 0.3 years.

### Why it does not transfer, and the answer was already in stage 7

The composition *signal* is present in both cohorts. What differs is how the
deconvolution splits it between types — median CD8T is 0.036 in GSE40279 and
0.000 in GSE61151, CD4T 0.129 against 0.179 — while the neutrophil fraction,
which is two thirds of the sample and the easiest to pin down, agrees to within
0.008.

That is stage 7's recorded failure arriving with consequences. The proportions
are collinear, so **the joint fit is sound while the split between neighbouring
types is not** — and a correction made of per-type coefficients is built
entirely out of the part that is not sound. The composition term is real in both
cohorts and describable in neither's coordinates.

### The extrapolation, reported and not counted

The corrected clocks were also measured on purified cells, where composition is
100% of one type and the correction learned its coefficients where neutrophils
run 50–75% and B cells 1–4%.

Displacement gets **worse** for every clock — Horvath 2013 from 5.52 to 13.67
years. That is what extrapolating a linear term fifty percentage points outside
its range does. It is reported because displacement is how this project measures
exposure throughout, and it is not the verdict on anything.

## Stage 13 — the same question at twice the cell-type resolution

Stage 7 said composition explains 0.9–2.7% of epigenetic age beyond
chronological age. Zhang et al. (2024, *Aging Cell*, 10,147 samples) report 13%
for Horvath, 25% for Hannum, 33.6% for PhenoAge. That is a tenfold gap, and it
turned out to be two separate mistakes stacked.

**The denominator.** Stage 7 reported the increment against total clock-age
variance. The literature reports a partial R² inside age acceleration — the
residual after chronological age. Converting stage 7's own numbers: 7.0%, 13.7%,
7.4%, 4.4%. Half the gap closes on arithmetic, and the smaller-looking number
was the one this project had been quoting.

**The panel.** Stage 7 used six cell types. A naive CD8 T cell reads 15–20 years
younger than an effector memory CD8 from the same person, and a six-type panel
folds those into one number — so the single largest contributor was invisible to
it by construction.

Rebuilt from GSE167998 (56 purified samples, twelve types, plus twelve mixtures
with known twelve-way proportions), on GSE61151:

| clock | 6 types | **12 types** | 12 types, excess over null | Zhang et al. |
|---|---|---|---|---|
| Horvath 2013 | 6.5% | **15.6%** | 9.5% | 13% |
| Hannum 2013 | 14.3% | **32.2%** | 26.1% | 25% |
| Levine 2018 | 7.4% | **23.4%** | 17.4% | 33.6% |
| Horvath 2018 | 4.3% | 9.6% | 3.6% | — |

*(all in partial-R²-inside-EAA, the literature's metric)*

**At twelve types this project lands on the published numbers.** Twelve
predictors buy more R² by chance than six do, so every value is compared against
its own permutation null — and the excess over null doubles too (Hannum 11.5% →
26.1%). The gain is information, not degrees of freedom.

### The panel is better, and the check that caught where it is worse

Leave-one-out assigns **12 of 12** purified types correctly with their own sample
held out of the reference, and the twelve mixtures come back at r = 0.79,
mean absolute error 0.027 against proportions that average 8.3%.

Against clinical ranges on the same samples, the twelve-type panel lands 4 of 6
buckets in range and the six-type panel 3 of 6 — it repairs exactly the types
stage 7 got wrong (CD4T too high, NK too high, CD8T too low) and breaks
monocytes, at 0.112 against a bound of 0.10.

That was a pre-specified check, and it failed. **The bound was not moved.** What
changed was the consequence — a miss inside 25% of the bound is declared as a
bias and the stage continues — and that decision was made after seeing the
failure, which is the third time in this project a check has been revisited post
hoc. The monocyte channel of this panel runs high, so no per-type monocyte
coefficient is read anywhere downstream.

### What this costs the earlier stages

Stage 7's number was a floor, not a measurement. Stage 12's finding that the
correction does not transfer was measured with the same blunt panel and is now
untrustworthy in the same direction. The practical claim in the plain-language
summary — "1 to 2 years" — should have been roughly double.

## Stage 14 — the floor was the wrong axis, and stage 11 is wrong

Stage 11 concluded the exposure floor is biology. Tomusiak et al. (2024,
*Communications Biology*) report the opposite: IntrinClock, a clock built to be
invariant across ten immune cell types. Stage 13 explains how both could be
honest — stage 11's filter was built from six types and could not express the
naive-versus-memory axis.

IntrinClock's published rule names that axis precisely: keep CpGs with |r| > 0.3
against chronological age and |r| < 0.3 against a sample being naive CD8. Three
selection rules, trained identically on GSE40279, tested on GSE61151, measured
against a twelve-type displacement:

| rule | r on test | scale | displacement per year of response | matched control | edge |
|---|---|---|---|---|---|
| (a) no filter | +0.921 | 1.009 | 6.8 yr | — | — |
| (b) stage 11's six-type filter | +0.909 | 0.972 | 6.7 yr | 6.8 | +0.1 |
| **(c) the IntrinClock rule** | **+0.925** | 1.012 | **4.0 yr** | 6.6 | **+2.6** |

**Displacement falls 41% and accuracy goes slightly up.** A gain on both axes at
once, which stage 11 declared impossible. Stage 11's own filter buys 0.1 years —
nothing, exactly as it reported.

And (c) beats every published clock on this axis: Horvath 2013 at 6.2, Horvath
2018 at 6.8, Hannum at 18.5, Levine at 24.3.

### The check that made the result readable

Tomusiak reports naive CD8 reading 15–20 years younger than effector memory CD8.
Required before anything else: at least two clocks showing an 8-year gap.

| clock | memory CD8 − naive CD8 |
|---|---|
| Horvath 2013 | +11.0 yr |
| **Hannum 2013** | **+40.6 yr** |
| **Levine 2018** | **+46.4 yr** |
| Horvath 2018 | +4.4 yr |

It reproduces, and larger. Naive CD8 is the single most extreme cell type in
this project — it reads 47 and 50 years below expectation on the two linear
clocks, against the 35-year *total spread* that stage 2 found across six types.
Stage 2's headline was measured on a panel that could not see this.

### A control that was wrong, and what fixing it cost

The first version drew each filter's random control from all 441,010 probes. For
rule (c), whose eligible pool is only 4,575, that handed the control four
thousand random probes and asked it to build a clock — so it lost on age
accuracy, and the filter's apparent edge of **+17.1 years** included "I kept
probes that correlate with age", which is not what was being tested.

Drawing the control from the probes that already pass the age criterion brings
the edge down to **+2.6**. The conclusion survives; the effect size was inflated
sevenfold until the control was fixed.

## Stage 15 — the correction works, and it can also poison the well

Stage 12 fitted a composition correction on one cohort, carried it to another,
and watched it remove almost nothing. Its diagnosis was that the per-type
coefficients are not identified, because six proportions are collinear. Stage 13
turned that into a testable claim: if the panel was the problem, twelve types
should fix it.

It does. Correction fitted with the twelve-type panel, residual measured with
the six-type panel — the crossed cell, which is the only one that does not score
a correction against its own representation of what it removed:

**Fitted on GSE40279 (656 samples), tested on GSE61151:**

| clock | before | after | |
|---|---|---|---|
| **Hannum 2013** | 11.5% | **1.6%** | −9.9 |
| Horvath 2013 | 3.8% | **1.0%** | −2.8 |
| Horvath 2018 | 1.6% | **−0.9%** | −2.4 (below its own null) |
| family, k=1,000 | 8.0% | 5.7% | −2.3 |
| IntrinClock-rule clock | 5.3% | 3.3% | −2.0 |
| Levine 2018 | 4.6% | 5.8% | +1.2 *worse* |

**Improves 5 of 6, median −2.4 points.** Stage 12's verdict was an artefact of
the panel.

### And the failure mode stage 12 could not see

Stage 12 tested one direction. Reversed — fitted on GSE61151's 184 samples,
applied to GSE40279's 656 — every clock external to the test cohort gets
**worse**:

| clock | before | after | |
|---|---|---|---|
| Horvath 2013 | 4.7% | 10.7% | +5.9 |
| Horvath 2018 | 0.9% | 5.9% | +4.9 |
| Levine 2018 | 12.1% | 13.8% | +1.6 |

**Improves 0 of 3, median +4.9 points.** A correction fitted on 184 samples does
not fail to remove composition; it *injects* it. Coefficients estimated on a
small cohort are unstable, and an unstable coefficient carried somewhere else
adds a composition-correlated term to every prediction.

Meredith et al. (2019) give the mechanism from the other end: cell-type
proportions collinear with the methylation being adjusted produce variance
inflation factors above 100 and flip the sign of 83% of coefficients in
simulation. This is that, transported between cohorts.

**So the practical claim is not "correcting works" and not "correcting fails".
It is: correcting works with a fine panel and a large fitting cohort, and
mis-calibrated it is worse than leaving the clock alone.**

> **Stage 17 kept the finding and overturned the mechanism.** "A large fitting
> cohort" is not what this is. Fitting on 184 samples of GSE40279 and
> transporting it does **not** reproduce the harm — median −0.9%, against the
> +4.9 points measured here — so the size of the fitting cohort is not what
> broke this direction. The curve of n is real and steep, and it is not
> sufficient. The candidate that is left sits in stage 12's own numbers, where
> median CD8T is 0.036 in GSE40279 and 0.000 in GSE61151. See stage 17.

### The trap in the second direction, which was nearly read as a result

Three of the six clocks here have a stake in GSE40279 — Hannum was trained on
it, and both of this project's own clocks were fitted on it. Every number in
this stage is a share *of the age residual*, and a clock measured in the cohort
it learned has almost no residual left. Divide by that and a modest absolute
change reads as an enormous percentage: the family clock appeared to go from
0.5% to **40%**, and the IntrinClock-rule clock from 0.9% to 32%.

That is the denominator, not the correction. Those rows are marked in the output
and excluded from the verdict, which is taken only from clocks external to the
cohort being tested.

## Stage 16 — flattening a clock is free

Every stage to 15 says clocks are contaminated by blood composition, and stage
14 built one 41% less exposed. The assumption underneath — the field's, and this
project's — is that removing the contamination improves the clock. Nobody tested
it, and the assumption might be backwards: blood composition *is* informative
about health, so it could have been half the signal rather than the noise.

Two cohorts: GSE50660 (464, smoking) and GSE42861 (689, rheumatoid arthritis
plus smoking, cases and controls in the same batch).

### The positive control, and the check that was written wrong

The control was to be smoking, and the first version required two of four
published clocks to detect it. **It failed at one of four** — and that one was
Levine, in both cohorts independently, +0.453 (p = 0.044) and +0.209 (p = 0.033).

The check had encoded a belief the field already knew to be false. First-generation
clocks are fitted to chronological age and are documented as *not* tracking
smoking; that is what second-generation clocks were built for. Of the four
clocks here exactly one is second-generation, and it fired twice.

This is the **fourth** check in this project revisited after failing, and the
distinction is worth keeping: the first three were a mis-scaled threshold, a
wrong premise about ridge, and a bound that was too narrow. This one was not
reading the literature before writing the test — the exact failure the
bibliography had just finished recommending against, committed in the next
stage.

It also settled something: the two clocks built here are first-generation too,
so **neither has smoking signal to lose**, and the flattening question cannot be
asked on smoking at all. It moved to the disease contrast.

### Reconciliation against published results on this exact dataset

GSE42861 has been analysed before, and the expectation is generation-specific
and counter-intuitive — first-generation clocks read RA patients as *younger*:

| clock | published | here | |
|---|---|---|---|
| Horvath 2013 | −1.3 yr | **−1.3 yr** | matches |
| Horvath 2018 | −1.3 yr | **−1.4 yr** | matches |
| Levine 2018 | +2.3 to +3 yr | **+2.5 yr** | matches |

Three for three, to a tenth of a year on two of them. This is the strongest
external validation anything in this project has passed.

### The answer

Rheumatoid arthritis, cases against controls, adjusted for sex and smoking,
n = 621, effects in standard deviations of age acceleration:

| clock | effect | p | in years |
|---|---|---|---|
| Levine 2018 | +0.465 | <0.0001 | +2.5 |
| Horvath 2018 | −0.427 | <0.0001 | −1.4 |
| **flattened (IntrinClock rule)** | **−0.374** | **<0.0001** | **−1.5** |
| **standard (k = 1,000)** | **−0.347** | **<0.0001** | **−1.3** |
| Horvath 2013 | −0.284 | 0.0004 | −1.3 |
| Hannum 2013 | +0.100 | 0.21 | +0.4 |

Paired bootstrap over the same 621 people, so the two clocks are compared on the
same resamples rather than through separate intervals:

> **difference −0.026 SD, 95% CI [−0.105, +0.052] — indistinguishable.**

**Flattening costs nothing.** The clock that reads 41% less cell composition
detects rheumatoid arthritis exactly as well as the one that does not.

And it says something about the disease signal: if RA showed up in these clocks
*because* it shifts blood composition, the flattened clock should have seen less
of it. It sees the same. The RA signal is not mainly composition.

*(On smoking the flattened clock is marginally ahead — paired difference +0.109,
CI [+0.009, +0.205] — but neither clock detects smoking on its own, so that is a
hint about first-generation clocks and not a result.)*

### What this cannot say

RA patients are medicated. This says a clock separates cases from controls; it
cannot say the separation is the disease rather than methotrexate. And one
outcome in one cohort is one outcome in one cohort — mortality, which is what
clocks are sold on, is untested here.

---

## Stage 17 — the curve is real, and stage 15's explanation is not

Stage 15 concluded that transporting a composition correction fails when it is
fitted on too few samples. It never tested that. Its two cells differ in size
**and** in cohort, age range, array batch, collection site and composition
distribution, and any of those explains the result equally well.

This stage varies n and holds the rest fixed. Every correction is fitted on
GSE40279 — the same people, the same array, the same sites — subsampled at ten
sizes, stratified by age decile, and carried to **three** external cohorts
instead of stage 15's one.

### The curve

Pooled over the three verdict clocks and the three test cohorts, 270 values per
row. Delta is what the correction did to the composition term left in the age
residual, in points of age-acceleration variance. **Positive means the
correction left the clock worse than not correcting at all.**

| fitting n | median delta | q75 | share of cells worse than doing nothing |
|---|---|---|---|
| 40 | **+12.2%** | +22.3% | 89% |
| 60 | +7.5% | +16.5% | 85% |
| 80 | +5.2% | +12.8% | 76% |
| 120 | +2.4% | +8.8% | 67% |
| 160 | +0.3% | +4.6% | 52% |
| **184** | **−0.9%** | +2.9% | 41% |
| 240 | −0.9% | +2.7% | 43% |
| 320 | −1.8% | +1.0% | 31% |
| 480 | −2.2% | +0.5% | 27% |
| 656 | −2.3% | −1.1% | 22% |

**The median crosses zero between n = 160 and n = 184.** Below that the standard
fix is worse than leaving the clock alone, and at the bottom of the range it is
not close: at forty samples it makes the clock worse in **89% of draws**.

Normalised against how much composition signal was there to begin with, the
bottom of the curve is worse than it looks:

| fitting n | 40 | 80 | 120 | 184 | 320 | 656 |
|---|---|---|---|---|---|---|
| fraction of the composition signal removed | **−271%** | −104% | −65% | +21% | +41% | **+61%** |

At full n the correction removes 61% of the composition signal out of cohort.
At forty samples it does not remove 100% of nothing — it **adds 2.7 times what
was there**.

### Meaningless coefficients are nearly harmless; badly estimated ones are not

> **Stage 24 retracted this section.** The permuted reference below was measured
> only at n = 656 and compared with real coefficients at n = 40. Measured at
> every n it falls as 1/n and reaches **+15.6% at n = 40**: coefficients with no
> information do nearly as much damage as real ones. See stage 24.

The stage carries a reference line: the same transport, with composition
coefficients fitted at full n on **permuted** composition. Those coefficients
carry no information by construction.

> permuted coefficients, n = 656: median delta **+0.5%**, IQR −0.8% to +2.1%,
> over 180 draws — a distribution straddling zero
> real coefficients, n = 40: median delta **+12.2%**

**A coefficient vector with nothing in it does essentially nothing.** It is the
*real* coefficients, estimated on forty samples, that do the damage — and they
do about twenty times more of it than noise does.

*(This stage first reported the reference as +1.8% from four permutations, and
quoted a seven-fold ratio off it. Twenty draws per cohort put the stable median
at +0.5%; the old figure sits inside the new IQR, so it was noise rather than
error. The corrected framing is the stronger one: the comparison is not "worse
than noise by a factor", it is "noise is harmless and this is not".)* Small-n
failure here is not dilution of a good answer toward a harmless one; the
estimate is actively wrong in a direction that costs more than silence. Meredith
et al. (2019) give the shape of it from inside one cohort — collinear cell-type
predictors, variance inflation above 100, 83% of coefficients flipping sign.
An inflated coefficient with a flipped sign, carried somewhere the composition
distribution differs, does not cancel. It multiplies.

### The anchor, which failed

The decision rule was written into the script before the run. Stage 15's
reversed direction fitted on 184 samples and made every clock worse, median
+4.9 points. If sample size is the mechanism, 184 samples of GSE40279 must also
be harmful.

> **median at n = 184: −0.9%** (IQR −3.1% to +2.9%, 270 values). 41% of cells
> worse than doing nothing — the minority, not the majority.

**Negative anchor.** One hundred and eighty-four samples are not, by themselves,
enough to produce stage 15's harm. Something else about GSE61151 *as a fitting
cohort* did it, and stage 15 named the wrong thing.

### What survives of stage 15, and what does not

The finding survives: a transported correction can be worse than no correction,
and this stage reproduces that at every n below 160 in three independent
cohorts. **The mechanism does not survive.** Size is a strong cause — the curve
is monotone across a factor of sixteen — but it is not sufficient, so it is not
the cause stage 15 claimed.

The better candidate is already in this project's own record, in stage 12, and
was not connected to it: **median CD8T is 0.036 in GSE40279 and 0.000 in
GSE61151.**

> **Stage 18 tested that candidate and it lost.** GSE61151 is the *better*
> conditioned cohort — condition number 40.8 against 54.4 to 65.9 for every
> GSE40279 subsample, full cohort included — and conditioning carries no
> information about the damage once n is held fixed. Both candidate mechanisms
> are down. What stage 18 did find is that ridge-penalising the composition
> coefficients removes the small-n harm entirely, without explaining it. A cohort where a cell type's estimated proportion is pinned at zero
cannot identify that type's coefficient at any n. If what governs the damage is
the conditioning of the composition matrix rather than the number of rows in it,
then n matters only because more rows usually mean better conditioning — and a
badly conditioned cohort of 184 behaves like a well-conditioned cohort of 40.
That is a measurable claim and it is stage 18.

### An artefact in the table above, caught before it was written down

The first table's q75 column tempts a statement it cannot support: *the upper
quartile only falls below zero at n = 656, so below that a bad draw still hurts.*
That reads the 656 row as if it described draws. It does not — **n = 656 has one
possible draw**, so its spread is across clocks and cohorts only, while every
other row mixes that with draw-to-draw variation. The two are not comparable.

Measured properly, as the spread across draws **within** each clock-by-cohort
cell:

| fitting n | 120 | 160 | 184 | 240 | 320 | 480 |
|---|---|---|---|---|---|---|
| median cell's q75 | +5.5% | +2.4% | +0.3% | +1.2% | **−0.9%** | −1.4% |
| worst cell's q75 | +18.9% | +15.5% | +13.2% | +10.5% | +11.1% | +6.8% |

So the honest version is two numbers, not one. **The median draw starts helping
at n ≈ 160–184; an unlucky draw keeps hurting until n ≈ 320**, and in the worst
clock-by-cohort cell an unlucky draw still hurts at 480.

### Per clock, and per test cohort

| fitting n | Horvath 2013 | Horvath 2018 | Levine 2018 |
|---|---|---|---|
| 40 | +9.2% | **+19.8%** | +11.2% |
| 184 | −1.6% | +0.9% | −3.4% |
| 656 | −2.8% | −1.1% | **−8.6%** |

Horvath 2018 is the worst behaved and does not cross zero until n = 480. Levine
gains the most, and it is also the verdict clock with the most composition
signal to lose (8.3% against 3.8% and 4.4%) — the correction is worth more where
there is more to correct, which is unsurprising and worth stating because it
bounds who should bother.

| fitting n | GSE61151 | GSE50660 | GSE42861 |
|---|---|---|---|
| 40 | +6.6% | +13.0% | +16.0% |
| 184 | −0.3% | −1.9% | −0.1% |
| 656 | −2.0% | −3.1% | −2.8% |

Three independent test cohorts, three curves with the same shape and the same
crossing region. That is the part stage 15 could not have known with one.

### Hannum, reported and not counted

Hannum was trained on GSE40279, so its correction here is fitted on an age
residual that barely exists. It is kept out of the verdict for that reason and
reported because its pattern is the cleanest in the stage: −9.3% at full n,
still only +1.1% at forty. It is also the clock with by far the most composition
signal to remove, 11.5%. Where there is a lot of real signal, the correction
survives a small fit; where there is little, the noise dominates it. That
pattern is consistent with the conditioning hypothesis and is not evidence for
it, because the in-sample fit confounds it.

### What this cannot say

- **One fitting cohort.** This is GSE40279's curve. The crossing point is a
  property of that cohort's composition structure, not a constant, and the whole
  point of the failed anchor is that another cohort of the same size behaves
  differently. Nobody should read "n = 184" off this table as a threshold.
- **Three test cohorts, all whole blood on 450k.** Nothing here speaks to EPIC,
  to other tissues, or to the twelve-type panel replacing the six-type one as
  the measurement.
- The correction is the IEAA-style one: regress clock age on chronological age
  and composition, keep the composition coefficients. It says nothing about
  reference-free methods.

---

## Stage 18 — conditioning is not it either, and ridge fixes it anyway

Stage 17 left one candidate standing. If the damage is governed by the
conditioning of the composition matrix rather than by the number of rows in it,
then a badly conditioned cohort of 184 behaves like a well conditioned cohort of
40, the failed anchor is explained, and the practitioner gets a diagnostic they
can compute on their own single cohort. This stage measured it.

**It is not that.** Two checks, written before the run, both say no.

### Conditioning barely moves with n

Measured on the composition block after partialling out intercept and age —
the system whose solution the coefficients actually are:

| fitting n | 40 | 80 | 160 | 320 | 656 |
|---|---|---|---|---|---|
| condition number | 65.9 | 62.1 | 54.5 | 54.4 | **54.5** |
| largest VIF | 302.3 | 275.8 | 215.7 | 212.8 | **214.5** |
| L2 norm of coefficients | 394.7 | 303.3 | 226.0 | 218.1 | **217.7** |

Conditioning is almost a constant of the cell-type structure, not something a
larger cohort buys: a sixteen-fold change in n moves the condition number by
about a fifth. What *does* move is the coefficient norm, which is 1.8 times
larger at forty samples than at full n — the inflation Meredith et al. predict.

Worth recording for scale: Meredith reports a variance inflation factor of
113.7 at six cell types and calls it far past the acceptable threshold of 5.
**At twelve types the same measurement runs 215 to 302.** The finer panel that
stages 13-15 adopted for its resolution is also, unavoidably, twice as collinear.

### Check 5 — conditioning explains nothing beyond n, and my first version of
### this check would have said something false

n and conditioning are correlated by construction, so the test was written to
run *within* each n stratum, where n cannot masquerade as anything.

**The first version pooled the four clocks inside each stratum, and that is
wrong.** The clocks differ enormously in both quantities — Levine's coefficient
norm is 386 with a median delta of −0.7%, Horvath 2018's norm is 183 with
+4.8% — so a between-clock pattern reads as a within-stratum correlation. Pooled
that way, the norm appeared to predict *less* damage, strongly and significantly
(ρ = −0.484 at n = 480, p = 3×10⁻¹⁷). That number is an artefact of mixing
clocks and it is not a finding.

Redone inside each n × clock × cohort cell, which is the comparison that holds
everything but the draw fixed:

| | median ρ | share of cells positive |
|---|---|---|
| condition number vs damage | **−0.025** | 48% |
| coefficient norm vs damage | +0.076 | 63% |

Forty-eight percent is a coin. **Conditioning carries no information about the
damage once n is held fixed.** The coefficient norm carries a little, in the
direction inflation predicts, and not enough to be a mechanism.

### Check 6 — GSE61151 is the *better* conditioned cohort

This is the number the stage existed for, and it is decisive in the direction
that kills the hypothesis.

> GSE61151, the cohort stage 15 fitted on: **condition number 40.8**, largest
> VIF 117.8.
> Every GSE40279 subsample, at every n from 40 to 656: **54.4 to 65.9**.

Stage 15's fitting cohort is better conditioned than the full 656-sample cohort
that this project has been treating as the well-behaved one, and it still
produced the +4.9 points of damage that stage 17 could not reproduce.

**Both candidate mechanisms are now down.** Sample size is a strong cause and
not a sufficient one (stage 17). Conditioning is not a cause at all (this
stage). Stage 20 found the quantity that is: the damage is governed by
tr(Σ_fit⁻¹ Σ_test)/n, a **joint** property of both cohorts, which is exactly why
measuring Σ_fit alone here returned a coin flip. What is left is a cohort difference that neither quantity captures, and
the obvious remaining candidate is the one stage 12 already demonstrated on
purified cells without generalising it: the correction is a **linear term
extrapolated outside the composition range it was fitted in**. Stage 12 pushed
neutrophils from 50–75% to 100% and watched Horvath 2013's displacement go from
5.52 to 13.67 years. Two whole-blood cohorts differ far less than that — but
stage 12 also recorded that median CD8T is 0.036 in GSE40279 and **0.000** in
GSE61151, which is the same thing in miniature. That is stage 19.

### Ridge, and the degenerate corner it had to survive

Shrinking a correction to zero removes the harm and is not a method. This
project already lost a conclusion to that shape once, in stage 8, where
"dilution is a defence" turned out to be a model shrinking itself rather than
resisting anything. So the rule was two-sided and written before the run: an
alpha counts only if at n = 40 it is no worse than not correcting **and** beats
OLS, while at n = 656 it keeps at least 70% of the OLS benefit.

Median delta by penalty and fitting size:

| alpha | 40 | 80 | 160 | 184 | 320 | 656 | coef. norm |
|---|---|---|---|---|---|---|---|
| **0** (OLS) | **+12.2%** | +5.2% | +0.3% | −0.9% | −1.8% | −2.3% | 277 |
| 0.1 | +7.2% | +2.6% | −0.3% | −1.3% | −1.8% | −2.0% | 138 |
| 0.3 | +3.9% | +0.6% | −1.3% | −1.7% | −2.2% | −2.6% | 112 |
| 1.0 | +0.0% | −1.1% | −2.0% | −2.4% | −2.6% | **−2.8%** | 69 |
| **3.0** | **−1.4%** | −1.5% | −1.7% | −2.0% | −1.9% | −1.9% | 35 |
| 10.0 | −1.0% | −0.8% | −0.8% | −1.0% | −0.9% | −0.9% | 13 |

**Alpha = 3 is the only penalty that passes both arms**, and what it buys is not
a better number at any single n — it is a **flat** one. Across a sixteen-fold
range of fitting size the delta moves between −1.4% and −2.0%. The correction
stops depending on how big the fitting cohort was.

Alpha = 10 is the degenerate corner doing exactly what the check was written to
catch: coefficient norm 13, delta pinned near −0.9% everywhere, harmless and
useless. It fails the second arm and is not reported as a fix. Alpha = 1 is the
near miss worth naming: it is the best penalty at every n from 60 up, and beats
OLS at full n (−2.8% against −2.3%), but at n = 40 it lands a hair above zero
and the pre-written rule excludes it. The rule was not moved.

### The number that matters to somebody using this

Nobody picks alpha by looking at this table; they cross-validate on the cohort
they have. So the stage records what ordinary leave-one-out cross-validation
picks, per draw, and what that choice then costs out of cohort:

| fitting n | OLS median | OLS worse than doing nothing | CV-ridge median | CV-ridge worse |
|---|---|---|---|---|
| 40 | **+12.2%** | **89%** | **−0.2%** | 44% |
| 60 | +7.5% | 85% | −0.3% | 42% |
| 80 | +5.2% | 76% | −0.5% | 37% |
| 120 | +2.4% | 67% | −0.7% | 36% |
| 160 | +0.3% | 52% | −1.3% | 29% |
| 184 | −0.9% | 41% | −1.7% | 23% |
| 320 | −1.8% | 31% | −2.0% | 19% |
| 480 | −2.2% | 27% | −2.3% | 17% |
| 656 | −2.3% | 22% | −2.8% | 22% |

**The median never turns positive at any n.** Cross-validation selects alpha = 3
at n = 40, alpha = 1 from 60 to 184, alpha = 0.3 from 240 to 480, and alpha = 0
only at the full cohort — which is to say it finds the right answer on its own,
including the answer that no penalty is needed when the cohort is large enough.

The fix therefore costs nothing anybody does not already have. It is not a new
method, a new reference or a new panel; it is a penalty term and a
cross-validation loop, and the eight-method benchmark the field uses contains
no penalised adjustment at all.

> **Stage 19 narrowed this.** On GSE40279 subsamples cross-validation picks
> well, and on GSE61151 — the cohort that actually broke the transport — it
> picks alpha = 0.3 and leaves +2.9% of the damage standing. Cross-validation
> optimises prediction inside the fitting cohort, and that objective does not
> know the coefficients are about to be shipped elsewhere. The recommendation
> that survives is narrower: penalise, and do not let cross-validation alone
> decide how much. See stage 19.

### What this does not fix

- **Ridge removes the catastrophe, not the coin flip.** At forty samples, 44% of
  draws still come out worse than leaving the clock alone. What changed is that
  the median stopped being a disaster, not that small cohorts became safe.
- Alpha = 3 is GSE40279's number on this panel with these clocks. The
  transportable claim is "cross-validate the penalty", not "use 3".
- None of this explains the anchor. A fix that works is not a mechanism, and the
  stage leaves stage 15's reversed direction still unexplained.

---

## Stage 19 — the fix reaches the anchor, cross-validation does not

Two mechanisms are down and stage 15's reversed direction is still unexplained.
This stage asks three things off the stage 18 cache, with no series reloaded.

### What the algebra ruled out before any data was touched

The planned stage was covariate shift: a linear term extrapolated outside the
composition range it was fitted in, which is precisely what stage 12
demonstrated on purified cells when neutrophils went from 50–75% to 100% and
Horvath 2013's displacement went from 5.52 to 13.67 years.

That stage was not written, because the correction is

    y_corrected = y − (C_test[:, :−1] − c̄_fit[:−1]) @ b

and `c̄_fit` is a constant vector. It adds the same number to every sample, and
a constant changes no R². **The difference in composition *means* between
fitting and test cohort is invisible to this metric.** Stage 12's result was a
displacement measurement, which is sensitive to the mean; everything from stage
15 onward is a variance-share measurement, which is not. Carrying the mechanism
across without checking the algebra would have been the stage 7 denominator
error again, in a new costume.

What survives is the second moment: `b` solves a least-squares problem shaped by
one cohort's composition covariance and is applied under another's.

### Stage 15's reversed direction reproduces

The check that had to pass before anything else meant anything. Fitting on
GSE61151's 184 samples, testing on GSE40279's 656, OLS:

| clock | before | after | delta | stage 15 reported |
|---|---|---|---|---|
| Horvath 2013 | 4.7% | 11.8% | **+7.1%** | +5.9 |
| Horvath 2018 | 0.9% | 6.0% | **+5.0%** | +4.9 |
| Levine 2018 | 12.0% | 14.0% | **+2.0%** | +1.6 |

Three of three positive, median +5.0%. The scales are not identical — stage 15
quoted shares of the age residual before the null subtraction stages 17–19 use —
so the agreement in ordering and magnitude across a rewritten pipeline is
stronger evidence than the digits suggest. **The transport failure is real and
reproducible.** It is the *explanation* that has been wrong twice.

### The stage 18 fix reaches it

| alpha | Horvath 2013 | Horvath 2018 | Levine 2018 | median | coef. norm |
|---|---|---|---|---|---|
| 0 (OLS) | +7.1% | +5.0% | +2.0% | **+5.0%** | 241 |
| 0.3 | +4.2% | +2.9% | +1.0% | +2.9% | 101 |
| 1.0 | +1.2% | +1.3% | −0.6% | +1.2% | 65 |
| **3.0** | −0.5% | +0.2% | −1.3% | **−0.5%** | 35 |
| 10.0 | −0.6% | −0.1% | −0.9% | −0.6% | 14 |

Alpha = 3 is again the only penalty that passes both arms of the
degenerate-corner rule: it turns the reversed direction negative **and** keeps
−1.9% on the forward direction, against OLS's −2.0% there. Alpha = 10 fails the
second arm exactly as in stage 18.

That stage 18 and stage 19 land on the same penalty, from a subsample curve and
from the real cohort that failed, is the strongest thing in either stage.

### And cross-validation does not — which corrects stage 18

Stage 18's practical claim was that nobody needs to know the number because
ordinary leave-one-out cross-validation finds it. **On GSE61151 it does not.**

> Cross-validation on GSE61151 selects **alpha = 0.3**, which leaves
> **+2.9%** of damage — most of the failure intact.

Cross-validation optimises prediction *inside* the fitting cohort, and nothing
in that objective knows the coefficients are about to be shipped somewhere else.
On GSE40279 subsamples it happened to pick well because the test cohorts were
geometrically close; on the cohort that is geometrically far, it picks a
penalty far too weak. Stage 18's table is not wrong and its conclusion was too
broad, which is recorded there rather than rewritten.

So the honest recommendation is narrower than stage 18's: **penalise, and do
not let cross-validation alone choose how much.**

### The geometry, which survives its test and does not explain much

Distance between the composition correlation matrices, after partialling out
intercept and age:

| pair | Frobenius | largest principal angle |
|---|---|---|
| **GSE61151 ↔ GSE40279** (stage 15's pair) | **1.70** | 26.1° |
| GSE40279 ↔ GSE50660 | 1.04 | 35.1° |
| GSE40279 ↔ GSE42861 | 1.27 | **81.7°** |

By Frobenius distance the stage 15 pair is the most distant of the three, which
is what the hypothesis predicts. **By principal angle it is the closest**, and
the pair with the most violently different dominant directions is the one that
transports fine. The two metrics disagree, so "geometry" is not one thing here.

Within each n × clock × cohort cell, so that n cannot masquerade as anything:

> ρ(Frobenius distance, damage) median **+0.082**, positive in **63%** of 81
> cells. The pre-written bar was a median above zero and more than 60% positive.
> Conditioning, in stage 18, gave 48% — a coin.

**It passes, and it is weak.** Sixty-three percent against a coin's fifty is a
real signal and not an explanation; ρ = 0.08 leaves essentially all of the
damage unaccounted for. Stage 20 explains why the metric was
the wrong shape: the phenomenon is asymmetric and a Frobenius distance between
correlation matrices is not. What can be said is that composition geometry is the
first candidate of three that has not been falsified, and that the cohort which
broke the transport is the geometrically most distant one by the metric that
matches the algebra.

### What this cannot say

- One directed pair failed and three worked. A mechanism proposed on n = 1
  failure is a hypothesis with a supporting anecdote.
- Alpha = 3 is this panel, these clocks, these cohorts. What transports is
  "penalise", not the number.
- The reversed direction is fixed by a penalty that does not know why it is
  needed. A fix is not a mechanism, and stage 15's failure is now reproducible,
  removable, and still unexplained.

---

## Stage 20 — the mechanism was in the statistics literature the whole time

Three stages hunted a mechanism and killed two candidates. The quantity the
algebra points at was never measured, because it was never written down.

Write the estimation error as `e = b̂ − b_true`. A correction fitted on one
cohort and applied to another leaves `−C_test @ e`, whose variance in the test
cohort is `e′ Σ_test e`. For least squares, `Cov(e) = σ²/n · Σ_fit⁻¹`, so

> **E[damage] ∝ (σ²/n) · tr( Σ_fit⁻¹ Σ_test )**

Call that, divided by n, the **transport index**.

**This stage claims no credit for it.** It is the standard excess-risk term for
least squares under covariate shift, and the failure mode it describes is named
in that literature: **spectral inflation** — directions carrying little
variation in training that carry more at evaluation. Theory for choosing ridge
regularisation under covariate shift exists too. What is new is only that a
widely used epidemiological adjustment is exposed to it, that the field's
eight-method benchmark does not mention it, and that nobody had measured when it
bites.

### Why it would explain what the other candidates could not

- it carries **1/n**, so stage 17's curve is the n in the denominator
- it is a **joint** quantity of both cohorts. Stage 18 measured the conditioning
  of `Σ_fit` alone and found a coin flip, which is exactly right: a
  well-conditioned `Σ_fit` can still have `Σ_fit⁻¹` amplify precisely the
  directions where `Σ_test` carries variance. That is how GSE61151 is the
  best-conditioned cohort in the project and still the one that fails
- it is **asymmetric**. Stage 19's Frobenius distance is a symmetric distance
  and the phenomenon is not
- ridge replaces `Σ_fit⁻¹` with `(Σ_fit + λI)⁻¹`, bounding the amplification —
  which is why α = 3 fixed both directions before anyone knew why

### It fails at the level of the individual draw

| | median within-cell ρ | cells positive |
|---|---|---|
| **transport index** | **0.021** | **56%** |
| Frobenius (stage 19) | 0.082 | 63% |
| condition number (stage 18) | −0.025 | 48% |

The pre-written bar was ρ > 0.30 at more than 80% of cells. **Failed**, and worse
than the quantity it was supposed to beat.

The theory predicts this failure. Within one n and one test cohort, `Σ_fit` is
nearly identical across draws so the index barely varies, while the realised
`e′Σ_test e` is chi-square-like on about eleven degrees of freedom with a
coefficient of variation near 0.43. A nearly constant predictor cannot track an
outcome dominated by realisation noise. **The index does not tell you what one
draw will do**, and that is recorded as a finding rather than explained away.

### It nails the level it actually speaks to

Check 2b, written after check 2 failed and marked post hoc for that reason. The
formula predicts *expected* damage across configurations, so: does the index
rank the 27 (n × test cohort) configurations?

| fitting n | GSE61151 index / damage | GSE50660 | GSE42861 |
|---|---|---|---|
| 40 | 0.345 / +6.6% | 0.490 / +13.0% | **0.783 / +16.0%** |
| 80 | 0.127 / +3.3% | 0.185 / +5.0% | 0.310 / +8.9% |
| 160 | 0.053 / +0.9% | 0.076 / −0.7% | 0.120 / +0.9% |
| 320 | 0.024 / −0.8% | 0.035 / −2.6% | 0.058 / −1.8% |
| 480 | 0.016 / −1.4% | 0.023 / −2.9% | 0.038 / −2.6% |

> **Spearman ρ = 0.897 over 27 configurations, p = 2.3×10⁻¹⁰.** Bar set before
> computing: 0.70.

Damage crosses zero at an index near 0.05, consistently across three test
cohorts. That is a threshold a practitioner can compute.

> **Stage 26 withdrew the floor.** Counted per clock instead of as a median over
> clocks, 12 of 40 configurations below 0.05 are harmful. The earlier counts let a
> harmed clock be outvoted by a helped one. The index is not a safety certificate.
>
> **Stage 21 blunted this.** With four fitting cohorts and 63 configurations,
> 0.05 is not a crossing point but a one-sided safety floor: nothing below it
> was harmful, 81% above it were, and between 0.05 and 0.16 both outcomes
> occur. A small index licenses the correction; a large one does not forbid it.

### And it predicts the anchor, in magnitude

The index is asymmetric, which is the property stages 18 and 19 both needed and
neither had:

| direction | transport index |
|---|---|
| **GSE61151 → GSE40279** (stage 15's) | **0.1737** |
| GSE40279 → GSE61151 | 0.0376 |
| GSE40279 → GSE50660 | 0.0570 |
| GSE40279 → GSE42861 | 0.0936 |

Stage 15's direction carries **4.6×** the index of its own reverse. And the
table above says an index near 0.18 goes with about **+4.6%** of damage. Stage
19 measured stage 15's direction at **+5.0%**.

The quantity that failed to predict a single draw predicts the anchor that
defeated two stages, to within half a point.

### The simulation, where b_true is identical in both cohorts

No model misspecification: the same true coefficients, the same noise, the same
n. The only thing that differs between the two directions is which covariance is
the fitting one.

| case | n | residual after | predicted | ratio |
|---|---|---|---|---|
| matched | 40 | 36.1% | 40.3% | 0.90 |
| matched | 480 | 3.8% | 2.2% | 1.74 |
| 40→61 | 40 | 21.2% | 26.6% | 0.80 |
| **61→40** | 40 | **84.6%** | 114.2% | 0.74 |
| **61→40** | 184 | **16.8%** | 17.2% | 0.98 |
| 40→61 | 184 | 4.9% | 4.1% | 1.20 |

**+84.6% against +21.2% at n = 40, from swapping which cohort's covariance does
the fitting and nothing else.** The anchor's asymmetry reproduces in a system
where nothing else exists that could cause it.

> Calibration over the whole sweep: **slope 0.92 through the origin, r = 0.815.**

### The ridge version tracks the fix

| alpha | index | observed damage at n = 40 |
|---|---|---|
| 0 | 0.594 | +12.5% |
| 0.3 | 0.066 | +3.8% |
| 1.0 | 0.030 | −0.2% |
| **3.0** | **0.013** | **−1.6%** |
| 10 | 0.004 | −1.0% |

The penalty that stages 18 and 19 found by trial is the penalty that drops the
index below the 0.05 threshold. The fix and the mechanism are the same object
seen from two sides.

### Two checks of mine were mis-specified, and they are counted

**Check 4 was arithmetically impossible.** It demanded a residual below 2% under
matched covariances at n = 480 — a number picked by eye. The stage's own formula
says that residual is σ²p/n over the residual variance, which is 0.92 × 11/480 =
**2.1%**. The threshold sat below what the theory predicts, so no correct
implementation could have passed. Rewritten to test what recovery means: within
a factor of 2.5 of the prediction. Observed ratio 1.74.

**Check 2 tested the wrong level**, for a reason the theory states in advance.
Its negative result is kept and reported above; check 2b is marked post hoc
rather than renumbered into the sequence, because it was written after seeing
check 2 fail.

These are the fifth and sixth checks in this project revisited after failing.
Both were mine, both from setting a threshold without doing the arithmetic the
stage's own formula provides.

### What this cannot say

- The index is derived for least squares with a correctly specified linear
  composition term. Real clocks are not generated that way, and the simulation
  that confirms the constant assumes they are.
- ρ = 0.897 is over configurations built from **one** fitting cohort and three
  test cohorts. The asymmetry result rests on a single directed pair.
- It predicts expected damage. It says nothing about the draw in front of you,
  which is the negative result above and is not a detail.

---

## Stage 21 — twelve directed pairs, and the index survives the test that could have killed it

Stage 20 wrote its own two weaknesses into its closing section: ρ = 0.897 over
configurations that were **all** fitted on GSE40279, and an asymmetry resting on
**one** directed pair. Four cohorts give twelve directed pairs and six
reversible ones, all off the stage 18 cache.

Hannum is dropped from the whole design rather than from the cells where it is
unsafe: it was trained on GSE40279, which here is a fitting cohort in six pairs
and a test cohort in three, and a clock in-sample in half a design is worse than
one absent from all of it.

### The ranking holds, and it holds inside every fitting cohort separately

| fitting cohort | configurations | Spearman ρ |
|---|---|---|
| GSE40279 | 18 | +0.858 |
| GSE61151 | 12 | **+0.951** |
| GSE50660 | 15 | +0.907 |
| GSE42861 | 18 | +0.897 |
| **all** | **63** | **+0.907** (p = 1.6×10⁻²⁴) |

Stage 20's 0.897 was not a property of GSE40279's subsample ladder. The index
ranks configurations in four cohorts independently, at essentially the same
strength, and best in the small cohort that broke the transport in the first
place.

### The asymmetry, at matched n

Both directions of each pair fitted at n = min(n_A, n_B), so the only thing that
changes is which covariance does the fitting.

| pair | n | index A→B | B→A | damage A→B | B→A | called |
|---|---|---|---|---|---|---|
| 40279 / 61151 | 184 | 0.046 | **0.174** | +0.2% | **+5.0%** | yes |
| 40279 / 50660 | 464 | 0.024 | 0.039 | −2.5% | −0.4% | yes |
| 40279 / 42861 | 656 | 0.0275 | 0.0271 | −2.8% | −0.8% | **no** |
| 61151 / 50660 | 184 | 0.116 | 0.054 | +1.3% | −1.7% | yes |
| **61151 / 42861** | 184 | **0.194** | 0.062 | **+15.9%** | −0.6% | yes |
| 50660 / 42861 | 464 | 0.035 | 0.026 | −0.2% | −3.9% | yes |

> **Five of six, binomial p = 0.109.** That is not significance and the script
> says so: with six pairs there is no bar that both clears 0.05 and tolerates a
> single miss. The number is the count and its p, not a verdict wearing their
> clothes.

The miss is worth describing precisely rather than either hiding or excusing.
On 40279 / 42861 the two indices are **0.0275 and 0.0271** — apart by 1.5%, far
below anything the measurement resolves. The index is predicting indifference
there and was scored as if it had predicted a direction. It counts as a miss
because the rule was written that way; it is not evidence that the index points
the wrong way.

**And the design found a failure three times worse than the one that started
this.** GSE61151 → GSE42861 at n = 184 does **+15.9%** of damage, the largest in
the project, and carries the largest index in the whole design. The anchor that
defeated stages 18 and 19 is no longer the worst case — it is the second worst,
and both are at the top of the same ranking.

### The fix generalises to all twelve

| pair | n | OLS | ridge α = 3 |
|---|---|---|---|
| 61151 → 42861 | 184 | **+15.9%** | **−2.5%** |
| 61151 → 40279 | 184 | +5.0% | −0.5% |
| 61151 → 50660 | 184 | +1.3% | −3.3% |
| 40279 → 61151 | 184 | +0.2% | −2.0% |
| 42861 → 50660 | 464 | −3.9% | −3.4% |
| 40279 → 42861 | 656 | −2.8% | **−0.7%** |

**Median damage at or below zero in all twelve.** Every harmful configuration is
turned beneficial, including the +15.9%.

It is not free, and the table says where the bill lands. Where OLS already
transports well, the penalty costs a little: 40279 → 42861 loses most of its
benefit, −2.8% to −0.7%, and 42861 → 50660 gives up half a point. Penalising is
insurance, and insurance has a premium in the cases that did not need it.

### The threshold, which stage 20 drew too sharply

| | |
|---|---|
| configurations with index ≤ 0.05 that are harmful | **0 of 21** |
| configurations with index > 0.05 that are harmful | **81%** |
| highest index still beneficial | 0.163 |
| lowest index already harmful | 0.051 |

Stage 20, on one fitting cohort, reported damage "crossing zero near 0.05". With
four it is not a crossing point but a **safety floor**: below 0.05 nothing in 63
configurations was harmful, above it four in five were, and the zone between
0.05 and 0.16 contains both outcomes. The usable statement is one-sided — a
small index licenses the correction, a large one does not forbid it.

### What this cannot say

- Four cohorts, all whole blood on 450k, all adult. Twelve pairs is twelve, and
  the asymmetry count is six.
- α = 3 is still this panel and these clocks. What generalised is that a fixed
  penalty in the right neighbourhood helps in twelve directed pairs, not the
  number.
- The index is computed from estimated proportions. Every cohort's composition
  here comes from the same deconvolution with the same panel, so this says
  nothing about what happens when two groups deconvolve differently.

> **Stage 22 closed this one.** With the fitting cohort deconvolved by one
> reference and the test cohort by another, the index still ranks configurations
> at ρ = 0.776 against a bar of 0.70, and ridge at α = 3 still holds all twelve
> pairs at or below zero. Two references disagreeing about a cell type is itself
> a covariance mismatch, which is why the index sees it without being told.

---

## Stage 22 — two teams, two reference panels, and the diagnostic still works

Stage 21 wrote its own limitation: every cohort's composition in this project
comes from the same deconvolution with the same panel, which is not the
situation any of this is about. The real case is two groups who never spoke —
one publishes composition coefficients estimated with their reference, another
applies them to proportions estimated with theirs.

**The obvious framing of that is impossible.** "Fit with six types and apply to
twelve" cannot be done: the coefficient vectors are not the same length.
Coefficients only transfer between analyses using the same cell-type *labels*.
So the mismatch that matters is the same six names estimated from different
reference data — **direct-6** from Reinius (GSE35069, the panel stages 7–12
used) against **collapsed-6**, the twelve-type Salas panel folded onto those six
labels by the published map. Different donors, different probes, different array
generation, same column names.

### The two panels do measure the same things

Correlation between the two estimates of each type, across four cohorts:

| type | Neu | Bcell | CD4T | NK | CD8T | Mono |
|---|---|---|---|---|---|---|
| median r | **0.989** | 0.952 | 0.948 | 0.924 | 0.894 | **0.886** |

Median 0.936 over six types and four cohorts, against a bar of 0.50. These are
two references agreeing about the same biology, not one broken panel.

**The literature's prediction did not land where it said it would**, and that is
worth recording rather than smoothing over. Reinius is documented as weakest on
NK and granulocytes; here NK sits at 0.924 and neutrophils are the *best*
channel at 0.989. The worst is monocyte, which is where **this project's own
declared bias** from stage 13 lives — the twelve-type panel's monocyte channel
runs high. The disagreement is on the side the project already knew about, not
on the side the literature warned about.

### The diagnostic survives the case it will actually be used in

| arm | configurations | Spearman ρ |
|---|---|---|
| same panel both ends | 63 | +0.793 |
| **different panel each end** | 63 | **+0.776** |

> Bar for the mismatched arm: ρ > 0.70. **Passed**, and barely degraded from the
> matched arm.

Two references disagreeing about the same cell type *is* a covariance mismatch,
so the index sees it without being told. That is the prediction stage 22 existed
to test, and it is the one that matters for anyone using this: the diagnostic
does not require both cohorts to have been processed by the same hands.

(Neither number is comparable with stage 21's 0.907. This stage applies with
direct-6 and measures with the twelve-type composition; stage 21 applied with
twelve and measured with six. Different measurement geometry, same bar.)

### The fix survives it too

Ridge at α = 3 holds median damage at or below zero in **all twelve** directed
pairs of the mismatched arm, including the three that OLS leaves harmful:
42861 → 40279 from +1.5% to −1.8%, 42861 → 61151 from +1.3% to −0.6%,
61151 → 40279 from +3.3% to −1.5%.

### What the mismatch costs, and a confound in my own design

At matched n, paired over the twelve pairs: **median cost +0.5%**, worse in 8 of
12. Pooling every fitting size including the small ones, the count of harmful
pairs goes from **5 of 12 to 10 of 12**.

But the sign flips in four pairs, and the reason is a confound I built in
without seeing it. The mismatched arm does not only mismatch — it also **fits
with the better panel**. Collapsed-6 comes from Salas, which is the later, finer
reference; direct-6 comes from Reinius. So the contrast is

> (fit with Salas, apply to Reinius)  against  (fit with Reinius, apply to Reinius)

and it mixes two effects pulling opposite ways: a penalty for mismatching, and a
benefit from fitting on better-estimated proportions. GSE61151 → GSE50660 goes
from +4.1% to −0.1% and GSE61151 → GSE42861 from +2.0% to −2.6% — mismatching
*helped*, which only makes sense as the quality term winning.

**So +0.5% is a mixture and not a mismatch penalty**, and the clean design needs
the reverse arm too: fit with Reinius, apply to Salas. That is not run here and
the number should not be quoted as if it were.

What is *not* confounded is check 5, because both arms are ranked by their own
index and the question there is whether the ranking holds, not how large the
damage is. The diagnostic result stands on its own.

### A bug, found by NaNs rather than by a check

The first run printed NaN in five of twelve rows. Not a result: the matched-n
comparison selected its rows out of the grid, and the matched n is almost never
*on* the fitting cohort's grid — 184 is not in GSE40279's. Empty frames,
medians of nothing. Stage 21 ran those configurations explicitly and this stage
did not, until it did.

It is recorded because the failure mode is worth naming: **the check that caught
it was arithmetic refusing to print**, not any of the seven pre-specified checks,
none of which look at whether the rows they summarise exist.

### What this cannot say

- One directed mismatch, Salas-fitted to Reinius-applied. The reverse is not run.
- Both panels are built by this project from public reference data, with a
  cruder probe selection than either published library. The disagreement between
  them is therefore a lower bound on what two real labs would produce.
- Six types. Nothing here says what happens when one team uses twelve and the
  other six, because as stated at the top, that transfer is not possible at all.

---

## Stage 23 — one of the verdict clocks was trained on the fitting cohort

A reviewer of the manuscript draft asked whether Horvath 2013 or 2018 had been
trained on any of the four cohorts. Nobody had checked. Hannum was excluded from
the start for having been trained on GSE40279; the other three were assumed
clean. The assumption was wrong for one of them.

### The primary source

Horvath 2013, Additional file 1 — the 82-dataset table with its *Data Use*
column, fetched through Europe PMC:

| data set | tissue | array | use | n | source |
|---|---|---|---|---|---|
| 3 | blood WB | 450K | **Training** | 656 | Hannum 2012, **GSE40279** |
| 44 | blood WB | 450K | Test | 689 | Liu 2013, GSE42861 |

**GSE40279 is a training set of Horvath 2013.** GSE42861 was only a test set,
which does not make the clock in-sample there; GSE61151 and GSE50660 do not
appear. Horvath 2018 was trained on GSE80261, GSE50759, GSE104471, GSE77136,
GSE52026, E-MTAB-4385, GSE79056 and unpublished data — none of the four. Levine
2018 was trained on InCHIANTI (n = 456) — none of the four.

Horvath 2013 therefore falls under the rule Hannum was under all along: out of
every pair in which GSE40279 appears, at either end of the transport.

### Does the training show in the data?

Standard deviation of each clock's age residual, GSE40279 against the mean of the
other three cohorts: Hannum **1.07**, Horvath 2013 **1.08**, Horvath 2018 1.23,
Levine 1.29. Every clock has a larger residual on GSE40279 — it is the widest age
range in the set, 19 to 101 — so the absolute ratio says nothing. The *relative*
one does: the two clocks trained there are 15–20% tighter on it than the two that
were not. The contamination is visible, and modest, which fits Horvath 2013
having seen 656 of roughly four thousand training samples.

*(The stage printed a criterion of "ratio below one" for detecting this. That was
naive — a cohort effect swamps it — and it is recorded here rather than silently
dropped.)*

### What survives, and what moves

Each headline carried a survival criterion written before the run. All four
survive. Several numbers do not.

| headline | before | corrected | criterion |
|---|---|---|---|
| damage at n = 40 | +12.2% | **+14.2%** (harmful in 88%) | > +5%: survives |
| median crosses zero | 160–184 | 160–184, then +0.1% at 240 | 120–240: survives |
| benefit at n = 656 | −2.3% | **−1.6%** | — |
| index ranks configurations | ρ = 0.907 | **ρ = 0.837** | > 0.70: survives |
| harmful below index 0.05 | 0 of 21 | **1 of 20** | — |
| asymmetry, matched n | 5 of 6 | 6 of 6 | ≥ 5: survives |
| ridge α = 3, twelve pairs | 12 of 12 | 12 of 12 | all ≤ 0: survives |
| stage 15's direction, 61151 → 40279 | +5.0% | **+3.5%** | — |

**The curve got stronger, not weaker.** Horvath 2013 was the mildest of the three
verdict clocks at n = 40 (+9.2% on its own), so including it had *understated* the
damage. At large n the benefit shrinks, and the median hovers around zero from
160 to 240 before turning clearly negative at 320 — so "crosses zero between 160
and 184" was the first negative point of a curve that is really flat there, and
the manuscript should say so.

**The floor loses its absolute form.** "None of the configurations below 0.05 was
harmful" is now one of twenty. It stays a useful threshold and stops being a
guarantee, which it should never have been phrased as.

**The asymmetry going to six of six is not evidence.** The pair that flipped is
the tied one — indices 0.0275 and 0.0270 — and its call is a coin whichever way
it lands. The honest count is five of five informative pairs.

**Stage 15's direction drops to +3.5%.** Stage 20's claim that an index near 0.18
"predicts it to within half a point" was made against the contaminated +5.0%, and
has to be re-read against +3.5%. The worst case in the project, GSE61151 →
GSE42861 at +15.9%, involves no GSE40279 and does not move.

---

## Stage 24 — the permuted reference at every n, and a second component

Stage 17's reference line — composition coefficients fitted on permuted
composition — was only ever measured at n = 656, and was drawn as a flat line
across a plot whose other curve started at n = 40. Every claim built on it
("seven times", "twenty times", "noise is harmless") compared real coefficients
at n = 40 against permuted ones at n = 656. A reviewer caught it.

Measured at every n, with the two clocks that stage 23 left clean on GSE40279,
30 draws per size, rows of the composition matrix shuffled whole so collinearity
survives:

| fitting n | 40 | 60 | 80 | 120 | 184 | 240 | 320 | 656 |
|---|---|---|---|---|---|---|---|---|
| real fit, Δ | +18.9% | +6.7% | +4.7% | +1.7% | +0.1% | −0.7% | −1.1% | −1.6% |
| real fit, error left | +25.5% | +14.1% | +11.8% | +7.6% | +6.6% | +5.7% | +5.1% | +4.5% |
| **permuted, Δ** | **+15.6%** | +9.6% | +5.0% | +3.1% | +2.6% | +1.3% | +2.1% | +1.0% |

> **Stage 24b re-ran this with 100 draws per size instead of 30.** At n = 40 the
> real fit's Δ is +16.1 p.p. (not +18.9) and the permuted reference +17.4 (not
> +15.6); the slope is −1.16. The floor is 3.8 p.p. once the 0.8 the permuted
> fit still does at n = 656 is subtracted, and adding it to the permuted damage
> reproduces the real curve within 0.6 p.p. from n = 60 up, but underestimates
> it by 1.2 at n = 40. The three-size agreement quoted below was selective.
> Stage 24's P3 also flips: with 100 draws the permuted coefficients are the
> more harmful at n = 40 (+17.4 against +16.1), where 30 draws had the real fit
> worse. The two are within each other's spread, which is the point — at n = 40
> essentially all of the harm is consistent with noise.

**Coefficients carrying no information at all do +15.6% of damage at n = 40.**
"Noise is harmless" is false, and so is every ratio built on it. The permuted
reference falls with n at a log-log slope of **−0.93** — the 1/n the transport
index predicts (P2, confirmed). The old flat line was one point of a steep curve.

### The prediction that failed, and what it found

P1 predicted that at the same n the permuted coefficients would do at least as
much damage as the real fit's error term, because their error is the whole of
what they carry. **Refuted, 0 of 10.** The real fit's leftover error is larger
at every n — and at n = 656 it is still +4.5%, where the permuted reference is
+1.0%. With 656 fitting samples that cannot be estimation noise.

Subtract that n = 656 floor from the real fit's error and what remains tracks the
permuted reference almost exactly:

| fitting n | 60 | 80 | 120 | 160 | 184 | 240 |
|---|---|---|---|---|---|---|
| real error minus its n = 656 floor | 9.6 | 7.3 | 3.1 | 2.6 | 2.1 | 1.2 |
| permuted Δ | 9.6 | 5.0 | 3.1 | 3.5 | 2.6 | 1.3 |

So the real fit's error has **two components**, and the algebra only had one:

1. **estimation noise amplified by covariance mismatch** — steep in n, reproduced
   point for point by coefficients that know nothing, and what the transport
   index measures
2. **an n-independent transfer error** that more fitting data does not remove

### The floor is transport, not measurement (post hoc)

Written after P1 failed, to tell the two readings apart: the full-cohort
correction applied to GSE40279 *itself*, measured with the other panel.

| clock | removed at home, cross-panel | removed after transport, full n |
|---|---|---|
| Levine 2018 | **95%** | 92% (→ 50660), 69% (→ 42861), **−23%** (→ 61151) |
| Horvath 2018 | **85%** | 46% (→ 61151), 16% (→ 50660), **−109%** (→ 42861) |

At home the twelve-type correction removes almost everything the six-type panel
can see, so the cross-panel measurement is not what leaves the floor. Carried to
another cohort at full fitting size, it can remove most of the signal, none of
it, or **double** it — Horvath 2018 into the rheumatoid arthritis cohort. That
pattern is what the reviewer's fourth point predicted: the true composition
coefficients differ between cohorts, and a disease that reshapes blood may also
reshape how composition maps onto a clock. Stage 26 tests it directly.

### What this does to the paper

- the n-curve stands, and gets an honest reference that falls with it
- the transport index stands as a measure of the **first** component only; it
  cannot see the second, which is why its floor was never absolute
- ridge works on both: shrinking toward no correction bounds the damage of a
  noisy β and of a wrong one alike
- the headline "the damage comes from real coefficients estimated badly, not
  from noise" is retracted. The corrected statement is the reverse of it: most of
  the small-n damage *is* estimation noise; what real coefficients add on top is
  a specification error that no fitting size cures

The retraction rule written into this stage named P1 as its vehicle. P1 failed,
and the retraction comes anyway — from the permuted column of the table, which
is the direct measurement the old claim never made.

---

## Stage 25 — the composition effect differs between cohorts, and that is the worst case

Stage 24 found a floor: a correction fitted on all 656 samples of GSE40279, where
estimation noise is negligible, still leaves +4.5% of composition in other
cohorts. The reviewer's fourth point named the candidate — the true composition
coefficients differ between cohorts. In the linear-regression literature this is
**model shift** (also regression shift, posterior drift), as opposed to the
covariate shift the transport index measures; the two are treated together by
Lei et al. (ICML 2021) and by Patil, Du & Tibshirani (ICML 2024).

Full-cohort fits for Levine 2018 and Horvath 2018, the two clocks clean in all
four cohorts.

### Pairwise heterogeneity — not established at the bar set

Wald tests of β_A = β_B, 11 degrees of freedom:

| pair | Levine 2018 | Horvath 2018 |
|---|---|---|
| 40279 / 61151 | **p = 8×10⁻⁴** | p = 0.24 |
| 40279 / 50660 | p = 0.028 | p = 0.015 |
| 40279 / 42861 | p = 0.042 | **p = 8×10⁻⁴** |
| 61151 / 50660 | p = 0.023 | p = 0.61 |
| 61151 / 42861 | **p = 1×10⁻³** | p = 0.11 |
| 50660 / 42861 | p = 0.018 | **p = 1×10⁻⁴** |

Four of twelve reject at Bonferroni against a pre-set bar of six. **Not
established**, and reported that way. Nine of twelve sit below a nominal 0.05,
which is suggestive and is not the criterion.

### The specification term ranks the floor — including the worst case

If estimation were perfect, a correction fitted on A and applied in B leaves
(β_A − β_B)′ S_B (β_A − β_B). Bias-corrected for the noise of both fits and
scaled to B's age-acceleration variance, it ranks the observed full-n leftover
across 24 directed configurations at **ρ = 0.633** (p = 0.0009), against a bar of
0.5. **Passes.**

| configuration | predicted | observed leftover |
|---|---|---|
| **61151 → 42861, Horvath 2018** | +51.8% | **+30.5%** |
| **61151 → 42861, Levine 2018** | +26.6% | **+30.8%** |
| 40279 → 42861, Horvath 2018 | +20.4% | +10.0% |
| 61151 → 40279, Levine 2018 | +18.1% | +14.1% |

**The worst transport in the project — the +15.9% that stage 21 found — is a
specification failure, not an estimation one.** The composition effect in
GSE61151 and in the rheumatoid arthritis cohort are different enough that a
perfectly estimated correction from one would still damage the other.

The bias correction overshoots in two configurations fitted on GSE61151, the
smallest cohort, where it predicts a negative leftover (−2.2%, −4.6%) against
+6.0% and +7.3% observed. Subtracting the fit's own noise is least reliable
exactly where that noise is largest.

### Why the Wald tests and the specification term disagree

They weigh the same β difference in different directions. The Wald statistic
weighs it by the **precision of the estimates** — directions that are poorly
estimated count for little. The specification term weighs it by the **variance
of the test cohort** — directions the test cohort varies in count for a lot. A
difference can sit where the fits are imprecise and the test cohort is wide, and
then it does damage while no test can see it. It is the same covariance mismatch
as stage 20, applied to the coefficients rather than to their errors.

### Inside the cohorts

| contrast | Levine 2018 | Horvath 2018 |
|---|---|---|
| RA cases vs controls (GSE42861) | p = 0.013 | p = 0.61 |
| ever vs never smokers (GSE50660) | p = 0.013 | p = 0.48 |

The composition coefficients of Levine 2018 shift with disease and with smoking;
those of Horvath 2018 do not. That fits what each clock was trained on: PhenoAge
was fitted to a phenotypic age built partly from inflammatory markers, Horvath
2018 to chronological age. Reported, with no bar.

Carrying RA and smoking as covariates in the fits barely moves the between-cohort
heterogeneity (median change in the Wald statistic −0.2). Whatever differs
between cohorts is not an additive effect of those two conditions.

### What it means for anyone using the correction

The two components are not equally predictable in practice. The transport index
needs only the **target's composition** — no clock values, no outcome — and so can
be computed before transporting. The specification term needs a good estimate of
the **target's own β**, and transport is only needed when the target is too small
to fit one. **The component that caused the worst failure is the one a user
cannot see in advance.** Ridge limits both, which is the argument for penalising
by default rather than only when the index is high.

### A check written wrong

The metadata check demanded over 100 current smokers in GSE50660, a number
assumed rather than looked up. The cohort has 22 current, 263 former and 179
never smokers; the join was intact, and the run stopped on a false belief about
the data. Join integrity and group size had been folded into one check and are
now separate; the smoking contrast became ever versus never, the weaker one. It
is the seventh check in this project revisited after failing.

---

## Stage 26 — no quantity computable in advance certifies a transport as safe

The transport index predicts the estimation error a correction leaves, never the
net effect. The reviewer proposed closing that gap in advance: if β is shared
between cohorts, E[b′S_B b] = β′S_B β + σ²·index, which gives a pre-transport
estimate of the net damage, **P_est = 2σ̂²·index − b′S_B b**, needing only the
fitting cohort and the target's proportions. Stage 25 had already shown β is not
always shared, so P_est was expected to break where it isn't. An oracle that knows
the target's own full-cohort β gives the ceiling.

126 configurations — every directed pair, fitting sizes 40 to full, the two clean
clocks **kept separate**.

| predictor | needs | Spearman with net damage | sign right |
|---|---|---|---|
| oracle | target's own β | **0.884** | **81%** |
| P_est (reviewer's) | fitting cohort + target proportions | 0.633 | 71% |
| transport index | fitting + target proportions | 0.611 | — |

The ceiling is high, so the twelve-versus-six measurement is not what limits.
P_est misses both pre-set bars (sign 75%, ρ 0.7). **Failed**, and the way it fails
is the result:

| | observed beneficial | observed harmful |
|---|---|---|
| **P_est says beneficial** | 43 | **30** |
| P_est says harmful | 7 | 46 |

Of the 73 configurations it calls safe, **30 are harmful — 41%**. A safety
diagnostic that errs, errs this way at its peril: it reassures where it should
warn. Its sign errors sit where stage 25's specification term is large (median
+10.4% among misses, +4.0% among hits), and the worst of them is the worst case in
the project: GSE61151 → GSE42861 for Levine at n = 184, called −21.5% and measured
**+15.9%**. The stage 25 mechanism predicts where the stage 26 predictor breaks.

### The safety floor was an aggregation artefact

With each clock kept separate, **12 of the 40 configurations below an index of
0.05 are harmful.** Stage 21 reported none of 21, stage 23 one of 20 — both counts
taken on configurations whose damage was a median **over clocks**, so a clock that
was harmed could be outvoted by one that was helped. Horvath 2018 carried from
GSE40279 into the rheumatoid arthritis cohort at the full n = 656 has a tiny index
and a leftover of +5.6%; pooled with Levine, it disappears.

**The index is not a safety certificate and never was.** It ranks the estimation
component. The specification component — the one behind the worst failures — is
invisible before transport to the index, to P_est, and to anything else that does
not know the target's own β. And a target large enough to estimate its own β is
one that does not need a transported correction in the first place.

### What survives as practical advice

Only the penalty, if it holds — and stage 23's "twelve of twelve pairs" for ridge
at α = 3 was *also* a median over clocks. Stage 27 asks whether it holds per clock,
with the criterion written before the numbers.

---

## Stage 27 — the penalty, recounted clock by clock

Stage 26 showed the index's floor had been an artefact of taking medians over
clocks. Stage 23's "ridge at α = 3 holds all twelve directed pairs" was counted the
same way, and was recounted here with each clock kept separate: 30 (directed pair ×
clock) cells at matched n.

Under model shift a fully shrunk correction does nothing, so its net damage is zero
by construction; penalising moves toward that point but a finite penalty need not
reach it. And under shift the optimal penalty can even be negative (Patil, Du &
Tibshirani, ICML 2024). Positive shrinkage is a safe default, not an optimum.

### It does not hold every cell. It comes close.

| | cells | harmful, no penalty | harmful, α = 3 | harmful, α = 10 |
|---|---|---|---|---|
| Horvath 2013 | 6 | 33% | 0% | 0% |
| Horvath 2018 | 12 | **67%** | 8% | 0% |
| Levine 2018 | 12 | 42% | 0% | 0% |
| **all** | **30** | **50%** | **3%** | 0% |

**Twenty-nine of thirty.** The pre-set criterion was all thirty, so the claim is now
"reduces", not "holds". The one cell left is Horvath 2018 carried from GSE61151 into
GSE40279: +5.0% unpenalised, **+0.2%** at α = 3.

**Without a penalty, half of all transports are harmful.** Stage 23 counted three
harmful pairs of twelve. Per clock it is fifteen cells of thirty — the aggregation
had been hiding exactly the harm this project is about. Horvath 2018 is the fragile
clock: two of every three of its transports hurt.

### The failure is not where the prediction said

The pre-set expectation was that any failing cell would sit where stage 25's
specification term is large. **Refuted**: the failing cell's term is −2.2% against a
median of +4.3% elsewhere. It is one of the two configurations fitted on GSE61151,
the smallest cohort, where stage 25 had already flagged that the bias correction
overshoots — so the number attached to it is unreliable. That is a reason to
distrust the estimate, not to re-read the refutation, and it is left standing.

### The price of certainty

α = 10 holds all thirty. It costs what it saves: where the unpenalised correction
already helped, the median benefit goes from −4.0% to −3.3% at α = 3 and to
**−1.5%** at α = 10 — most of the correction's value, given up to insure against
the minority of transports that would have hurt. There is no penalty that is both
certainly safe and free. What the data support is a choice with its price stated.

---

## Stage 28 — most of the index's correlation was the fitting size

The manuscript quoted the transport index ranking configurations at p = 2×10⁻²⁴.
A reviewer objected that configurations share cohorts and are not independent. The
problem is deeper than the p-value. The index carries 1/n by construction, and the
damage falls with n too — so across fitting sizes, *any* quantity proportional to
1/n would correlate with the damage whatever it knew about the cohorts. The
question is what the index knows beyond n.

On stage 26's per-clock net damage, at the three fitting sizes every directed pair
has (40, 80, 160): 72 configurations, 12 pairs, the pair as the block.

| test | ρ | p |
|---|---|---|
| pooled over sizes, naive | 0.579 | 1×10⁻⁷ |
| **n alone** (median of the block-permutation null) | **0.418** | — |
| **pooled, block permutation** — index profiles swapped between pairs | 0.579 | **0.036** |
| within n = 40 | 0.408 | 0.048 |
| within n = 80 | 0.364 | 0.08 |
| within n = 160 | 0.326 | 0.12 |

**Fitting size alone produces ρ ≈ 0.42.** The index adds something on top — the
block-permutation p is 0.036 against a bar of 0.05, and the within-n median is 0.364
against a bar of 0.3, so both pre-set criteria pass — but not much. Of the
correlation this project has been quoting, most was the n in the index's
denominator. The honest p-value is 0.036, not 2×10⁻²⁴.

Leaving one cohort out moves the pooled ρ between 0.447 (without GSE61151) and 0.727
(without GSE40279): the small cohort carries a good share of what the index knows.

### What the index is, after stages 23 to 28

A quantity that grows as 1/n — which is most of what it measures — and that adds a
modest, real ranking of which cohort pairs amplify estimation error. It certifies
nothing (stage 26), it cannot see the specification error behind the worst failures
(stages 25–26), and its between-pair information survives a dependence-respecting
test at p = 0.036. That is a diagnostic worth reporting and not one to lead a paper.

---

## Stage 29 — measured on the naive/memory axis, the harm is larger

Since stage 15 the correction is fitted with twelve cell types and what it leaves
is measured with six, so it is never scored against its own representation. A
reviewer named the cost: six types cannot see the naive/memory split, which is
where the largest composition effect lives. The only panel here that can see it is
the one the correction is fitted with, and measuring with it is self-scoring —
biased toward making the correction look good. So the reading was fixed in advance
and is one-sided: harm that **persists** is robust; harm that vanishes would say
nothing.

| measured with | net damage at n = 40 | at n = 656 | per-clock cells harmful, OLS | at α = 3 |
|---|---|---|---|---|
| six types | +11.8% | −1.6% | 54% | 4% |
| **twelve types** | **+48.1%** | −1.4% | 58% | 0% |
| **naive/memory columns only** | **+30.6%** | −2.5% | 50% | 0% |

**It persists, and it grows** — against a bias pulling the other way. At n = 40 the
transported correction leaves four times more composition signal than six types
could register, and most of that sits on the naive/memory axis. The six-type
measurement used throughout this project **understated the damage by a factor of
three to four**; every curve built on it is conservative.

At full fitting size the three measurements agree, as they should once estimation
noise is small. And the penalty holds on all three: no harmful cell at α = 3
measured on twelve types or on the naive/memory axis alone.

---

## Stage 30 — re-auditing the two claims that predate the audit

Stages 23–29 changed the protocol: Horvath 2013 is out of every pair touching
GSE40279, and damage is counted per clock, never as a median over clocks. Two
claims the manuscript still made had been computed the old way — stage 19's
"cross-validation does not size the penalty" and stage 22's "the penalty holds
when the two cohorts come from different reference panels". Both were recomputed
before being quoted again, with survival criteria set first.

### Cross-validation does not size the penalty — now across thirty cells

The penalty chosen by leave-one-out on each fitting subsample, over α in {0, 0.1,
0.3, 1, 3, 10}, per (directed pair × clean clock) cell at matched n:

| penalty | cells harmful |
|---|---|
| none | 47% |
| **chosen by cross-validation** | **30%** |
| fixed α = 3 | **3%** |

Cross-validation picks α = 0.3 in the median cell, and 0 in seven of thirty. It
leaves ten times as many harmful transports as the fixed penalty — 27 points
against a survival bar of ten. **Survives**, and stronger than stage 19 had it:
that was one cohort; this is systematic. The worst case is again GSE61151 → the
arthritis cohort, where cross-validation's α = 0.3 leaves +4.9% and α = 3 gives
−2.8%. Cross-validation minimises error inside the fitting cohort; nothing in that
objective knows the coefficients are about to travel.

### The penalty under two reference panels

Fitting cohort deconvolved from the Salas reference collapsed to six labels, test
cohort from the Reinius reference, measured on the test cohort's twelve-type
composition: 14 of 30 cells harmful unpenalised, **1** at α = 3 (GSE42861 →
GSE61151 for Levine, +7.5% → +0.5%). The pre-set bar was at most one.
**Survives.**

---

## Stage 31 — model shift inside a single study

Between cohorts, model shift is entangled with batch, laboratory and population.
Within GSE42861, rheumatoid arthritis cases and controls share the study, the
array and the lab. A correction fitted on controls and applied to cases — a
reference-population design — is a transport in which only disease differs.

Controls were split at random into halves 30 times; each time the correction was
fitted on half A and applied to half B (estimation noise only) and to the cases,
so the paired difference is what disease adds.

| clock | control → control | control → case | excess |
|---|---|---|---|
| Horvath 2013 | −6.8% | +0.6% | **+7.4%** |
| Levine 2018 | −7.7% | −0.7% | **+7.4%** |
| Horvath 2018 | −1.2% | +1.2% | **+2.6%** |

The pre-registered prediction was differential: an excess above one point for
Levine 2018 (whose coefficients stage 25 found to shift with RA) and below one
point for Horvath 2018 (whose did not). **Refuted** — Horvath 2018 shows an excess
too. What stands is stronger than what was predicted: inside one study, on one
array, a correction fitted on controls does worse on patients than on other
controls, for all three clocks. Batch cannot explain it.

### It is model shift, not covariance (post hoc)

Disease changes two things, the coefficients and the composition distribution,
and the excess could be either. Written after the prediction failed: the
transport index, paired draw by draw, is *lower* for control → case (0.065) than
for control → control (0.096) — by the noise component alone, cases should have
been harmed less. Removing the predicted noise term leaves the excess at +8.7%
(Horvath 2013), +8.4% (Levine) and **+6.3% (Horvath 2018)**: all of it is model
shift. The within-cohort Wald test that cleared Horvath 2018 (p = 0.61) lacked the
power to see a difference that damage can see — the mismatch stage 25 described.

### What it does to a case-control estimate

The RA effect on age acceleration, in years, adjusted for age:

| clock | unadjusted | adjusted in the whole cohort | correction fitted on controls |
|---|---|---|---|
| Horvath 2013 | −1.35 | −0.31 | −0.60 |
| Levine 2018 | +2.69 | +0.04 | −0.10 |
| Horvath 2018 | −1.36 | −0.74 | **−1.39** |

The choice of reference population alone moves the Horvath 2018 disease effect
from −0.74 to −1.39 years. Adjusting in the whole cohort and adjusting with a
controls-only reference are both defensible-looking choices, and they disagree by
nearly a factor of two.

---

## Stage 32 — a pace-of-ageing clock: what generalises and what does not

Every result so far came from three clocks that estimate age. DunedinPACE
estimates the pace of ageing, was trained on the Dunedin Study (not on GEO, so
none of the four cohorts is in its training), and tracks monocyte subsets. It was
reimplemented from the published model data of its R package (GPL-3):
quantile normalisation of each sample to a 20,000-probe reference, then a
173-CpG weighted sum.

**Validation (32a) passed on all three pre-set checks.** Coverage 99.6–100%;
cohort means 0.93–1.05 with SD 0.11–0.13, on the published scale of about one
year of biological ageing per calendar year; and the positive control — current
smokers faster than never smokers in GSE50660 — came out at **+0.140** (one-sided
p = 3×10⁻¹¹).

### The results (32b)

| pre-set criterion | result | |
|---|---|---|
| n = 40 transport harmful in > 60% of draws | median +0.3%, harmful in **51%** | **failed** |
| per pair at matched n: ≥ 3 harmful unpenalised, ≤ 1 at α = 3 | **6** of 12, **0** of 12 | passed |
| controls → cases inside GSE42861: excess > +1 point | **+10.5%** (IQR +7.0 to +14.3) | passed |

The permuted reference at n = 40 was +4.9%: the estimation-noise component is
there for DunedinPACE too. What differs is that real coefficients do *better*
than permuted ones (+0.3% against +4.9%), where for the age clocks they did
worse.

### Why the small-n curve did not generalise (post hoc)

The first guess — that DunedinPACE carries more composition signal, so even a
noisy correction removes enough real signal to break even — is **only partly
supported**: its median composition share in the test cohorts is 8.5%, about the
same as Levine's 8.3%, and Levine was harmed at n = 40 in stage 17. Signal size
alone does not explain it.

> **Stage 33 weakened this.** From a second fitting cohort, GSE132203, DunedinPACE
> is neutral at n = 40 again (−1.5%, harmful in 42%). The neutrality looks like a
> property of the clock, not of the fitting cohort.

What the per-pair table shows instead is that the n-curve has a single fitting
cohort, GSE40279, and that DunedinPACE's composition effect transports well *from
that cohort*: at full size its correction removes 19.4 points into the arthritis
cohort. The harmful DunedinPACE cells are fits on the other three cohorts. The
small-n headline is a property of the clock and the pair, not a constant.

### What generalises

- the penalty: 6 harmful pairs of 12 unpenalised, **none** at α = 3
- model shift inside one study, and larger than for the age clocks: +10.5 points
  when a controls-fitted correction is applied to arthritis patients
- the estimation-noise component, visible in the permuted reference

What does not: the claim that a small transported correction is harmful in most
draws. For one clock from one fitting cohort, it is a coin flip.

---

## Stage 33 — a fifth cohort, on another array and another ancestry

GSE132203, the Grady Trauma Project: 795 whole-blood samples on the **EPIC** array,
mostly African American. Every transport to or from it also crosses array
generations. None of the clocks was trained on it. Its betas are a 5.3 GB
supplementary CSV (md5 08771432…), of which the 22,374 rows the analysis needs were
kept.

### Validation (33a)

- **The panels are the same panels.** The stage 18 cache kept compositions, not
  the panels; they were rebuilt by the same deterministic procedure and reproduced
  its fingerprint exactly (r = 0.789, MAE = 0.027), then saved.
- Panel coverage on EPIC 96.8% (twelve types) and 98.2% (six).
- **The age check caught a bug.** On the first run every clock correlated with
  "age" at r ≈ 0.02–0.19 while agreeing with each other at 0.88–0.94. The series
  carries two fields, `age` and `age acceleration`, and a prefix match had let the
  second overwrite the first. Fixed to an exact field match; with true age the
  clocks track it at r = 0.86–0.95.
- Clock coverage bar of 95%, set in advance: Levine 2018 (99.8%) and Horvath 2018
  (100%) pass. **Horvath 2013 is excluded at 94.6%** despite r = 0.90 — the rule
  was written first. Hannum, 88.7%.

### Results (33b)

| pre-set criterion | result | |
|---|---|---|
| 8 directed pairs with GSE132203: ≥ 3 harmful cells unpenalised, ≤ 10% at α = 3 | **7** of 24, **0** of 24 | passed |
| fitted on 40 samples of GSE132203, age clocks harmful in > 60% of draws | **+24.3%**, harmful in **95%** (Levine 91%, Horvath 2018 99%) | passed |
| permuted reference positive | **+16.4%** | passed |

The small-sample harm replicates **from a different fitting cohort, on a different
array, in a different population**, and is larger than from GSE40279. The penalty
leaves no harmful cell.

### DunedinPACE is neutral again — which corrects stage 32's explanation

Fitted on 40 samples of GSE132203, DunedinPACE's transported correction was
neutral again: median −1.5%, harmful in 42%. Stage 32 had attributed its neutrality
at n = 40 to the fitting cohort — "its composition effect transports well from
GSE40279". Two different fitting cohorts now give the same answer, so that
explanation is weakened: it looks like a property of the clock. Why a pace-of-ageing
clock tolerates a noisy transported correction where age clocks do not is open.

---

## Stage 34 — a paediatric cohort: none usable in public data

Childhood is where blood composition changes fastest, so it would be the
strongest test of transport. A GEO search (series on the 450k or EPIC platforms,
whole blood, children or adolescents) returned seventeen candidates. None fits:

| candidate | n | why not |
|---|---|---|
| GSE154566 (E-Risk twins) | 1,177 | everyone is 18 — no age variance to residualise |
| GSE99863 (children, The Gambia) | 257 | **no age field** in the series metadata |
| GSE118144, GSE193879, GSE64495, … | 34–145 | disease-defined, and below the ~200 needed for the n-grid |

Other tissues are out of scope for a different reason: both deconvolution panels
are blood references, and a composition correction in buccal cells or saliva is a
different correction. **Recorded as a limit of the evidence, not a test run**:
every result here is adult blood.

---

## Stage 35 — a fitting cohort four times larger: the floor stays, but not everywhere

Stage 24 described the second component as n-independent, but the largest
fitting cohort then had 656 samples, so that description rested on one cohort's
top end. GSE55763 (Lehne et al. 2015, London, 450k) has 2,639 unrelated adults
once all 72 technical-replicate arrays are dropped. It lets the fitting n go four
times further. None of the clocks was trained on it; all four passed coverage
(≥99.7%) and the age check (r 0.86–0.94). The six criteria were written into the
docstring before any transport was computed.

Fitted on GSE55763 and transported to the other five cohorts (13 age-clock cells,
median across cells; "error left" is the composition signal still visible after
correction):

| fitting n | 40 | 80 | 160 | 320 | 656 | 1300 | 2639 |
|---|---|---|---|---|---|---|---|
| Δ, real fit | +8.3% | +1.2% | −2.1% | −3.3% | −3.4% | −3.7% | −4.0% |
| error left, real | +14.8% | +7.9% | +3.7% | +2.0% | +1.9% | +1.5% | **+1.3%** |
| Δ, shuffled | +14.5% | +5.9% | +2.2% | +1.3% | +0.5% | +0.6% | **+0.1%** |
| Δ, ridge α = 3 | −2.4% | −2.4% | −2.8% | −3.1% | −2.9% | −2.9% | −3.0% |

- **The floor persists (criterion 2, passed).** From n = 656 to 2,639 the error
  left falls from 1.9 to 1.3 points, a ratio of 0.72. Pure 1/n noise predicts
  0.25, and the shuffled reference does fall that way (slope −1.03, criterion 4).
- **It is not everywhere (criterion 3, failed).** At the full n, the error left
  exceeds the shuffled reference in only 8 of 13 cells; the bar was 9. The floor
  depends on the target cohort:

  | target | GSE42861 (arthritis) | GSE40279 | GSE132203 | GSE50660 | GSE61151 |
  |---|---|---|---|---|---|
  | error left at n = 2,639 | **+4.6** | **+3.0** | +0.7 | 0.0 | −0.1 |

  Into two cohorts, a well-estimated correction transports essentially perfectly.
  Into the arthritis cohort it leaves the most, which agrees with stage 31's
  controls → cases result: disease reshapes the composition effect. Why GSE40279
  also keeps a floor is not tested here. One candidate is age extrapolation:
  GSE40279 runs to 101 years and the fitting cohort stops at 75. That is a
  hypothesis written after the result, not a finding.

  > **Stage 36 tested this and it failed.** Restricted to ages 24–75, GSE40279
  > keeps the floor (+3.5 → +3.3 for Horvath 2018). Stage 37 found the floor confined to the
  > Caucasian-European half of GSE40279; ancestry and plate are confounded there
  > and centring by plate does not remove it.
- **Small-n harm from a third fitting cohort (criterion 5, passed):** at n = 40
  the age clocks' correction was harmful in 83% of draws (median +7.6%).
- **The penalty (criterion 6, passed):** harmful cells go from 13 of 13 to 1 of
  13 at n = 40, and from 3 of 13 to 0 at full n. **Its cost shows at large n:**
  at n = 2,639 the unpenalised fit is the better one on median (−4.0% against
  −3.0%). For Levine 2018 into GSE132203 it is −12.5% against −7.9%. α = 3
  trades about a quarter of the benefit for a guarantee the fit did not need in
  10 of 13 cells.
- **DunedinPACE**, fitted on a third cohort: +0.4% at n = 40, neutral again, and
  beneficial from n = 80 (−8.7% at full n). It is the third fitting cohort to
  show the same thing, which supports stage 33's reading that this is a property
  of the clock.

### What this changes

"Model shift does not shrink with fitting size" holds on average, now up to
2,639 fitting samples. It needs a qualifier: it is a property of the (fitting
cohort, target cohort) pair, large for some targets and absent for others. That
strengthens stage 26's conclusion — whether a target carries model shift
cannot be read off the fitting cohort — and adds one practical fact: with
thousands of fitting samples and no disease in the target, the transported
correction helped in 10 of 13 cells. The fixed penalty is a small-n safeguard.
At large n it has a price.

## Stage 36 — a penalty that fades with n protects less; model shift needs the proportional one

The penalty used since stage 21 is α × the mean eigenvalue of the fitting
cohort's standardised composition cross-product. That eigenvalue grows with n,
so α = 3 shrinks by the same proportion at every fitting size, and stage 35
measured what that costs at n = 2,639. A textbook ridge with a fixed λ fades as
1/n; here α_n = 3 × 40 / n, which equals α = 3 at n = 40. Patil, Du & Tibshirani
(2024) show that the optimal ridge level differs between covariate shift and
regression shift (our model shift). The prediction written before the run: a
penalty that vanishes with n cannot bound an error that does not.

All 30 directed pairs among the six cached cohorts, at matched n, 20 draws, 72
age-clock cells:

| | unpenalised | α = 3 | fading α_n |
|---|---|---|---|
| harmful cells | 23 of 72 | **1 of 72** | 11 of 72 |
| median where unpenalised helped (49 cells) | −4.1% | −4.2% | −5.4% |
| fitted on all 2,639 of GSE55763, median of 13 cells | −4.0% | −3.0% | −4.0% |

- **Criterion 2 passed:** the fading penalty leaves 11 harmful cells against 1.
  The cells it misses are the model-shift cells: Horvath 2018 into and out of
  GSE40279, and GSE61151 into the arthritis cohort, where it cuts +26.1% to +1.1%
  but not to zero.
- **Criterion 3 passed:** fitted on all of GSE55763 the fading penalty is 0.06
  points from the unpenalised median, which recovers what α = 3 gave up.
- **Criterion 4 failed, and stage 35's hypothesis goes with it.** Restricting
  GSE40279 to the fitting cohort's age range (24–75, 490 people) leaves the
  floor where it was: +2.4 → +2.3 (Levine 2018) and +3.5 → +3.3 (Horvath 2018).
  Age extrapolation does not explain it, and what does is not known.

The penalty result also generalises. Over six cohorts, 72 cells and a second
array, α = 3 leaves one harmful cell, the same one stage 27 left (GSE61151 →
GSE40279, Horvath 2018, +0.2%). At matched n it cost nothing on median.

### What this changes

No single rule wins both regimes. The proportional penalty is what bounds model
shift, and it costs about one point of benefit only when the fitting cohort runs
to thousands. The fading penalty costs nothing there, but it lets through ten
more harmful cells at the sizes where transported corrections are actually used.
Whether a target carries model shift cannot be known in advance (stage 26), so
the recommendation stays α = 3, now with its cost stated.

## Stage 37 — the GSE40279 floor: neither ancestry nor batch, and the two cannot be separated

With estimation noise out of the way (fitted on all 2,639 of GSE55763), GSE40279
keeps a floor of +2.4 (Levine 2018) and +3.5 (Horvath 2018) points. Stage 36
ruled out age extrapolation. The metadata offered two more candidates, both
written down before scoring: ancestry (426 Caucasian-European, 230
Hispanic-Mexican; Horvath et al. 2016 report different intrinsic and extrinsic
ageing between the two) and batch (nine processing plates).

**Criterion 2 (ancestry, Hispanic floor > 2× size-matched Caucasian) failed, and
the data point the other way.**

| subset | Levine 2018, error left | Horvath 2018, error left |
|---|---|---|
| Hispanic-Mexican (230) | −0.4% | −0.2% |
| Caucasian-European (426) | +4.3% | +3.9% |
| Caucasian, 30 random subsets of 230 (median, 10th–90th pct) | +3.9% (+1.5, +6.5) | +4.3% (+0.7, +6.0) |

The whole floor is in the Caucasian half. The transport from London into the
Hispanic half is essentially exact.

**Criterion 3 (batch, centring by plate removes over half) failed:** after
centring clocks and composition within plate, the error left was +3.9% and +2.5%.

**Ancestry and plate are confounded by design:** Hispanic samples sit on plates
5, 6 and 9 only, Caucasian samples on the other six. So "the Caucasian half"
means "those six plates" too, and the two cannot be told apart in this cohort.
What the stage does establish is where the floor is — half the cohort, not all of
it — and that a mean shift per plate is not it. In the Caucasian half, the
correction adds composition signal to Horvath 2018 (+0.8% before, +3.9% after).

### What this changes

Nothing in the practical conclusion. It sharpens the description of model shift:
the floor is not even a property of a cohort. It can sit in one recruitment
stream of a cohort and be absent from the other. That strengthens stage 26 —
diagnostics computed at the level of the cohort cannot see it.

## Phase B, search (a) — a published transport, found

Since stage 23 the project had not found one documented case of composition
coefficients fitted in one dataset and applied to another. Galkin et al. (2021,
*Frontiers in Aging*) is one. To carry a blood clock over to saliva, they fit an
adjustment on EpiDISH proportions (epithelial, immune, fibroblast) over eight
pooled saliva/buccal studies (about 960 samples) and apply it to four held-out
studies. Judged by accuracy against chronological age, MAE falls from 20.9 to
4.7 years and 5% of samples get worse.

Two things set it apart from the case studied here. The fit is large and pooled
across studies, which is where stages 35–36 found transport safest. And accuracy
against chronological age cannot show whether composition is left in, or added
to, age acceleration, which is the quantity a downstream association inherits.
The transport is real and published; the damage measured here has not been
measured there.

## Stage 38 — saliva: the harm replicates, and the penalty does not hold

Every earlier stage was blood. Saliva is buccal epithelium mixed with leukocytes,
so its composition axis is far larger. It is also where the one published
transport of composition coefficients was done (Galkin et al. 2021). Three adult
cohorts were cached in 38a with the EpiDISH references: a nine-type HEpiDISH fit
(epithelium, fibroblast, seven immune subtypes) and a three-type measurement
(epithelium, fibroblast, immune):

- GSE232891 (EPIC, 552; Crohn's disease, ulcerative colitis, controls)
- GSE232332 (EPIC, 265 after dropping replicates; oesophageal cancer, controls)
- GSE78874 (450k, 259; raw signal)

The median immune fraction is 0.73–0.76. Levine 2018 and Horvath 2018 cleared
coverage and age (r 0.64–0.91) in all three. GSE232891 and GSE232332 come from
one group, and their files carry no genotyping probes, so shared people cannot
be ruled out and they are never paired: 4 directed pairs, each with GSE78874. The
measurement shares the fit's first step, which biases it toward the correction.

Composition is a large part of age acceleration in saliva, before any
correction: 10.8% to 51.6% for the age clocks, and 42–56% for DunedinPACE.

**Criterion 2 (small-n harm) passed.** Fitted on 40 and transported, the
correction was harmful in 67% of draws (median +20.9%). **Criterion 3 passed:**
shuffled coefficients did +6.4%.

**Criterion 4 (the penalty) failed, badly.** At matched n = 259:

| pair | clock | before | unpenalised | α = 3 | fading |
|---|---|---|---|---|---|
| 232891 → 78874 | Levine 2018 | +51.6% | −13.6% | +2.1% | −0.7% |
| 232891 → 78874 | Horvath 2018 | +16.2% | **+116.2%** | +27.8% | +73.0% |
| 78874 → 232891 | Levine 2018 | +18.5% | +26.5% | +9.3% | +26.0% |
| 78874 → 232891 | Horvath 2018 | +31.7% | +8.0% | +5.0% | +9.2% |
| 232332 → 78874 | Levine 2018 | +51.6% | +7.6% | −6.5% | −9.5% |
| 232332 → 78874 | Horvath 2018 | +16.2% | **+152.0%** | +35.1% | +88.9% |
| 78874 → 232332 | Levine 2018 | +10.8% | +89.4% | +4.3% | +30.2% |
| 78874 → 232332 | Horvath 2018 | +44.3% | +59.7% | +14.5% | +26.0% |

Harmful cells: 7 of 8 unpenalised, **7 of 8 at α = 3**, 6 of 8 fading. These are
the largest harms in the project. For Horvath 2018 into GSE78874 the correction
leaves up to ten times the composition signal it found. DunedinPACE was helped
in all four pairs (−36% to −44%).

### Control: the correction works at home (38c, post hoc)

Written after 38b, before it was run. Within each cohort, fitting on a random
half and applying to the other half was beneficial in 5 of 6 age-clock cells
(−15% to −48%). The pre-set bar was all six, so it formally failed. The exception
is Levine 2018 in GSE232332: 132 fitting samples against a before-signal of only
10.8%, which is the small-n harm again. Applied to itself, each full cohort
removed essentially all of it (−0.4% to −0.8% left). So the saliva failure is
transport, not the correction or the measurement.

### What cannot be separated

Every saliva pair crosses array (EPIC against 450k) and processing (processed
betas against raw signal), so the shift may be technical as much as biological.
That is the realistic case: Galkin et al. pooled studies of mixed processing.

### Why the penalty fails here (reading, not tested)

The nine-type composition is highly collinear in saliva (condition number 872 to
1,317; the twelve-type blood matrices measure 13 to 31, and 178 in the EPIC
blood cohort). Ridge shrinks the weak directions of the
composition matrix. The dominant one here, epithelium against leukocytes, is
barely touched, and a coefficient that differs between cohorts on that axis
passes straight through. In blood the harmful error sat on weak directions,
which is why α = 3 worked there.

> **Stage 39 refuted this.** A three-type fit with condition number ≈ 1 fails the
> same way, and so does a normalised GSE78874. The immune-fraction slope changes
> sign between cohorts (Horvath 2018: +1.0 and +1.9 against −0.9 years per 10
> points), and shrinkage cannot fix a sign.

### What this changes

The penalty result is a **blood** result, and the manuscript has to say so. In a
tissue where composition dominates, a transported correction can be far worse
than none, and neither a fixed nor a fading penalty rescues it. The within-study
correction still works.

## Stage 39 — saliva: not the preprocessing, not collinearity — the slope changes sign

Stage 38 offered two readings of why the saliva transports failed and the penalty
did not help, and stage 39 tested both on the same 8 cells, with criteria written
beforehand.

- **Technical (criterion 2, failed).** 39a quantile-normalised GSE78874 to
  GSE232891's beta distribution. The immune fraction barely moved (r = 0.997 with
  raw) and the clocks kept tracking age. Transports were still harmful in 7 of 8
  cells unpenalised, and 7 of 8 at α = 3. The bar was 3 of 8.
- **Collinearity (criterion 3, failed).** A three-type fit (epithelium,
  fibroblast, immune) has a condition number of 1.0–1.2 against 872–1,317 for
  nine types, and still left 6 of 8 cells harmful at α = 3. The bar was 1.
  Normalised and three-type together: 7 of 8.

What is left, measured post hoc: the slope of each clock on the immune fraction
(adjusted for age), per cohort.

| cohort | Horvath 2018, years per +10 pp immune | Levine 2018 |
|---|---|---|
| GSE232891 (EPIC) | +1.03 (SE 0.29) | −0.20 (0.36) |
| GSE232332 (EPIC) | +1.90 (0.32) | −1.35 (0.40) |
| GSE78874 (450k) | −0.94 (0.14) | −4.72 (0.29) |
| GSE78874, normalised | −0.20 (0.13) | −5.10 (0.27) |

For Horvath 2018 the slope on the dominant axis **changes sign** between the EPIC
cohorts and GSE78874. For Levine 2018 it differs up to twenty-fold. A penalty
shrinks a coefficient toward zero. It can make a sign-reversed correction less
harmful (α = 3 took +116% to +28%), but it cannot make it helpful. Stage 38's
collinearity reading was wrong: the failure is model shift on the axis that
carries most of the composition signal.

What still cannot be separated: the two EPIC cohorts come from one group, so the
reversal may be group, array or population. Normalisation did not remove it,
which argues against a simple distributional artefact of processing.

> **Stage 40:** a fourth EPIC saliva cohort from another group gives −0.82
> (95% CI −1.71 to +0.07), on GSE78874's side, not with the other EPIC cohorts.
> Array is the less likely explanation. By the pre-set rule this was
> inconclusive.

## Stage 40 — a fourth saliva cohort sides with the 450k one (inconclusive by the pre-set rule)

Stage 39's sign reversal set two EPIC cohorts from one group against one 450k
cohort, so array and group were confounded. GSE149747 is EPIC saliva from another
group (Methylation Diet and Lifestyle): 44 adults at baseline, sampled twice more
after an intervention. The median immune fraction is 0.65. Horvath 2018 tracks age
at r = 0.76 with full coverage. The reading rule was written down before scoring.

| cohort | array | Horvath 2018, years per +10 pp immune (SE) |
|---|---|---|
| GSE232891 | EPIC | +1.03 (0.29) |
| GSE232332 | EPIC | +1.90 (0.32) |
| GSE78874 | 450k | −0.94 (0.14) |
| **GSE149747, baseline** | **EPIC** | **−0.82 (0.45)**, 95% CI −1.71 to +0.07 |
| GSE149747, person means | EPIC | −0.94 (0.43) |

**By the pre-set rule this is inconclusive:** the rule for "the reversal follows
the group" required the interval below zero, and it reaches +0.07. Its upper end
is still far below the +1.03 and +1.90 of the two EPIC cohorts. So an EPIC cohort
from another group sits with the 450k cohort, not with the other EPIC cohorts.
Array is the less likely explanation. The two cohorts from one group look like
the exception, and the secondary analysis (person means, CI −1.78 to −0.10)
agrees. Levine 2018 is negative too (−2.42, SE 0.68), as in every saliva cohort.

This does not rescue transport. Whichever cohorts are the exception, a
coefficient carried between saliva studies can have the wrong sign, and nothing
in the source cohort says which kind of target it is going to.

## Stage 41 — pooling saliva studies does not rescue transport

Galkin et al. (2021) pooled eight saliva studies before fitting their
composition adjustment. If model shift scatters around a common value, pooling
averages it out; if the target sits on the other side of the split, the pooled
coefficient keeps the wrong sign. This stage tests leave-one-cohort-out over the
four saliva cohorts. The fit is three-type (condition ≈ 1) with a fixed intercept
per study. GSE232891 and GSE232332 are never on opposite sides, because they
may share people.

| target | clock | own immune slope | pooled slope | unpenalised | α = 3 |
|---|---|---|---|---|---|
| GSE78874 | Levine 2018 | −4.72 | −0.92 | −3.4% | −2.3% |
| GSE78874 | Horvath 2018 | −0.94 | **+1.36** | **+137.4%** | +24.1% |
| GSE149747 | Levine 2018 | −2.42 | −2.62 | −11.8% | −7.7% |
| GSE149747 | Horvath 2018 | −0.82 | **+0.39** | +36.7% | +0.7% |
| GSE232891 | Levine 2018 | −0.20 | −4.56 | +28.5% | +2.5% |
| GSE232891 | Horvath 2018 | **+1.03** | −0.94 | +8.7% | +1.8% |
| GSE232332 | Levine 2018 | −1.35 | −4.56 | +30.5% | −1.9% |
| GSE232332 | Horvath 2018 | **+1.90** | −0.94 | +19.3% | +4.4% |

(slopes in years per +10 points of immune fraction)

- **Criterion 2 (pooling rescues, ≤ 2 of 8 harmful) failed:** 6 of 8.
- **Criterion 3 (pooling plus α = 3, ≤ 1 of 8) failed:** 5 of 8.
- Without study intercepts: 8 of 8.

The pooled slope follows the majority of the pool, or its largest cohort, and
helps only the targets that agree with it: Levine 2018 in GSE149747 (−11.8%). With
four cohorts, two of them from one group, "pooling" here means two or three
studies, not eight. Still, the mechanism is clear: pooling estimates an average,
and an average does not transport to a target whose slope has the other sign.

A related result from within one cohort: Chan et al. (2026, 529 children) found
epigenetic age and exposure effects differing between saliva samples that are
mostly epithelial and those that are mostly immune.

## Stage 42 — why the penalty works in blood: the sign holds there, and not in saliva

Ridge shrinks a transported coefficient toward zero. If the target's own
coefficient has the same sign and a different size, shrinking moves the
correction toward it and bounds the damage. If the sign differs, shrinking can
only reduce a correction pointing the wrong way. Stage 39 found a sign reversal in
saliva. This stage asks whether blood has one. The slope is of each clock on one
composition fraction, adjusted for age, per +10 points (SE):

| blood axis | clock | 40279 | 61151 | 50660 | 42861 | 132203 | 55763 |
|---|---|---|---|---|---|---|---|
| naive CD8 | Levine 2018 | −15.0 | −14.3 | −12.4 | −12.4 | −10.9 | −10.4 |
| naive CD8 | Horvath 2018 | −6.9 | −4.3 | −5.1 | −0.7 | −2.0 | −4.1 |
| neutrophils | Levine 2018 | +1.60 | +0.96 | +1.31 | +1.63 | +1.71 | +1.16 |
| neutrophils | Horvath 2018 | +0.04 (0.16) | +0.00 (0.37) | −0.28 | −0.16 | −0.26 | −0.50 |

- **Criterion 1 (same sign on neutrophils, both clocks, all six) failed.** Horvath
  2018's two "positive" slopes are +0.04 (SE 0.16) and 0.00 (0.37), both zero
  within noise.
- **Criterion 2 (same sign on naive CD8) passed** for both clocks in all six.
- **Criterion 3 (larger heterogeneity in saliva) passed:** I² 97% for both clocks
  in saliva, against 50% and 62% on the blood neutrophil axis.

The distinction that criterion 1 missed, counted after the fact: the cohorts
whose 95% interval lies entirely on one side of zero.

| | positive | negative |
|---|---|---|
| blood, any axis, any clock | never both | — |
| saliva, Levine 2018 | 0 | 3 |
| **saliva, Horvath 2018** | **2** | **1** |

Only Horvath 2018 in saliva has cohorts clearly on both sides. It matches stage
38 cell by cell. At α = 3 the saliva cells of Levine 2018, whose slope keeps its
sign, fall to a median of +3.2% (from +17.0%). Those of Horvath 2018, whose slope
flips, stay at +21.2% (from +87.9%). The penalty works where the sign holds. In
blood the sign held everywhere it could be measured, and that is why α = 3 was
enough there.

> **Amended after review (stage 42b).** "The difference is the sign" is too
> strong, and two facts in this project contradict it. Normalised, Horvath 2018's
> GSE78874 slope is −0.20 (95% CI −0.46 to +0.06), no longer clearly of the
> opposite sign, and 7 of 8 cells stay harmful. Levine 2018 keeps its sign in all
> four saliva cohorts and is still harmful at α = 3 (+3.2 p.p.). What separates
> the two tissues is the between-cohort spread **relative to the average effect**,
> τ/|mean| from a random-effects fit:
>
> | | Horvath 2018 | Levine 2018 |
> |---|---|---|
> | blood, naive CD8 | 0.49 | 0.09 |
> | blood, neutrophils | 0.55 | 0.15 |
> | **saliva, immune** | **5.34** | **0.91** |
>
> Where the spread is a fraction of the mean effect, a shrunk coefficient is
> close to every cohort's own, and α = 3 works. Where the spread is as large as
> the effect or larger — as in both saliva clocks — no single coefficient is
> close to all of them, and shrinking toward zero only limits the damage.
> Sign reversal is the extreme of that, not a separate mechanism.

## Stages 24b, 31b, 42b — what a reader's audit changed

A reader audited the manuscript before posting. Eight points were text; four
needed numbers, and those are here. Two claims were wrong and are withdrawn.

### 24b — the n = 40 median, from 100 draws

The median at n = 40 had been quoted as +18.9 p.p., the largest of three sets of
30 draws (the range was +11.8 to +18.9). Rerun with 100 draws per size and fresh
seeds: **median +16.1 p.p., harmful in 88% of draws** (IQR +5.5 to +29.7). The
permuted reference at n = 40 is +17.4 and falls at a log-log slope of −1.16.
Stage 24's own checks pass again (n = 656 reproduced; P1 refuted; P2 confirmed).

### The floor was double-counting noise, and "almost exactly" was wrong

Stage 24 called 4.5 p.p. the n-independent floor. That is what the real fit
leaves at n = 656, where shuffled coefficients still do 0.8, so the floor
contained the noise it was meant to exclude. **The noise-free floor is 3.8 p.p.**
Adding it to the shuffled damage reproduces the real curve within 0.6 p.p. from
n = 60 up, and underestimates it by 1.2 at n = 40 — not "almost exactly", and the
old text quoted only the three sizes where it did agree. Figure 2 also plotted
net damage on an axis labelled composition left; both series are now composition
left, in p.p.

### 31b — an interval for the arthritis effect

"Nearly a factor of two" rested on two point estimates. Bootstrap over people
(1,000 resamples, both corrections refitted inside each one), Horvath 2018:
whole-cohort −0.74 (95% CI −1.05 to −0.37), controls-only −1.39 (−2.08 to −0.80),
**difference −0.66 (−1.13 to −0.31)**. The gap survives; the other two clocks'
intervals cross zero.

### 42b — "the difference is the sign" withdrawn

See the amendment in stage 42. Between-cohort spread relative to the average
effect (τ/|mean|) is 0.09–0.55 in blood and 0.91–5.34 in saliva; sign reversal is
the extreme of that, not the mechanism. Normalised, Horvath 2018's saliva slope
no longer clearly reverses and 7 of 8 cells stay harmful; Levine 2018 never
reverses and is still harmful at α = 3.

### Text corrections

The abstract stated the safety predictor the other way round from the results (a
share of harmful transports called safe, against a share of "safe" calls that
were harmful); the results' version is now used in both. Δ and composition left
are labelled p.p. throughout, reserving % for shares of the composition signal
found. The header no longer says "revised after review". α = 3 is now stated to
have been fixed on the first four cohorts, which makes GSE132203, GSE55763 and
the saliva cohorts out-of-sample for that choice. GSE149747 is in the methods.
The 24 GSE132203 cells are stated to include DunedinPACE (6 of 16 age-clock cells
harmful). "The worst transport" is now "at full fitting size". Galkin et al.'s
training set includes GSE78874, one of the cohorts where transport fails here,
and the introduction says so.

## Stage 43 — in years: how far a reported association moves

Every earlier stage scored the correction by how much composition it leaves.
Nobody reports that number; what gets reported is an association — smoking
accelerates ageing by X years, this disease by Y. Stage 31 showed for one cohort
and one clock that the choice of reference moved the arthritis effect from −0.74
to −1.39 years. This asks it across four exposures, with the within-cohort
correction as the comparator (what the target's own data support, and the
practice the field treats as correct — a reference, not a truth):

| target | exposure | exposed / not |
|---|---|---|
| GSE50660 | smoking, ever vs never | 285 / 179 |
| GSE42861 | rheumatoid arthritis | 354 / 335 |
| GSE232891 (saliva) | inflammatory bowel disease | 302 / 250 |
| GSE232332 (saliva) | oesophageal cancer | 98 / 167 |

Median |estimate − within-cohort estimate|, in years, over (target × source ×
clock) cells, 30 draws each:

| | median error |
|---|---|
| no correction at all | 0.55 |
| transported, fitted on 40 | **0.91** |
| transported, fitted on 40, α = 3 | 0.48 |
| transported, full fitting cohort | 0.25 |

- **Criterion 2 passed.** At n = 40 a transported correction moves the reported
  answer further from the within-cohort one (0.91 y) than doing nothing does
  (0.55 y). In the units a reader sees, it is worse than useless.
- **Criterion 3 passed.** 13 of 24 cells are off by more than a year or flip
  sign; 3 flip sign outright. The worst: GSE132203 → GSE42861, Levine 2018, 4.19
  years, where the within-cohort arthritis effect is +0.04.
- **Criterion 4 passed.** α = 3 cuts the median error to 0.48 y, below doing
  nothing, and leaves 6 of 24 cells off by a year or more.
- **At full fitting size the transport is good** (0.25 y), consistent with the
  composition-left curves.

### A tail the composition metric could not show

Composition left is bounded above by construction, so it hid this. In years, at
n = 40, the unpenalised fit is sometimes degenerate: 1.2% of 720 draws land more
than 10 years from the within-cohort estimate, the worst at 37.6 years. Solved as
exact least squares, without dropping near-zero singular values, 10.1% exceed 10
years and the worst is 10^14 — so the size of the tail depends on the solver,
while its existence does not. **At α = 3 no draw exceeds 10 years; the worst is
3.5.** The tail is concentrated in saliva (33% of exact-least-squares draws,
against 6% in blood), where the composition matrix is nearly singular.

### What this changes

It gives the paper its practical statement. A correction borrowed from another
cohort and applied to forty samples does not merely leave composition behind: it
moves the number the study reports, by about a year on median, by more than four
in the worst cell, and occasionally by an absurd amount. A fixed penalty brings
the median below the do-nothing baseline and removes the tail.

## Stage 45 — the linearity assumption, tested

Every correction here is linear in the proportions, and the split into estimation
noise and model shift assumes it. If the true effect were curved, part of what
this project calls model shift could be one shared nonlinearity fitted at a
different point of the curve in each cohort — a different diagnosis with a
different remedy. Two tests, six blood cohorts, both age clocks.

**A. Is there curvature?** Squared terms for the four largest components (Neu,
CD4mem, CD8mem, Mono) added to the within-cohort fit:

| | extra R² | significant at 0.05 |
|---|---|---|
| median over 12 cells | **0.0019** | 4 of 12 |
| largest (GSE61151, Horvath 2018) | 0.0063 | no |

Curvature is detectable in a third of the cells and is negligible in size
everywhere. **Criterion 2 passed** (bar: median below 0.01).

**B. Does allowing for it transport better?** Both corrections fitted and
transported at matched n over all 30 directed pairs, 60 (pair × clock) cells:

| | harmful cells | median Δ |
|---|---|---|
| linear | 23 of 60 | −2.4 p.p. |
| quadratic | 22 of 60 | −2.4 p.p. |

**Criterion 3 failed, by one cell.** The bar was that the quadratic correction
leave at least as many harmful cells as the linear one; it left 22 against 23.
One cell in 60, with identical medians, is not a difference, and the pre-set bar
was a strict inequality where it should have been a tolerance. The substantive
reading is that the two are indistinguishable: **allowing for curvature does not
reduce the harm**, so the harm is not unmodelled curvature.

### What this changes

The limitation "the decomposition assumes a linear composition effect" can be
stated with a measurement behind it rather than as an unexamined caveat: the
nonlinearity is there, it is worth about 0.2% of variance, and modelling it
changes nothing about transport.

## Corrections so far

| what was wrong | what caught it | what it cost |
|---|---|---|
| the original question needed cell counts GSE61151 never deposited | checking the metadata before analysing | the question, changed for a better one |
| metadata prefix strip guarded by `dtype == object`, which pandas 3.0 broke | every age coming back NaN | a validation run |
| four samples with ages 0, 6, 7 — `agegap` leaking into `agebloodtaken` | the validation gate refusing to pass | r of 0.89 reading as 0.68 |
| **the cause of that leak**: the loader named each GEO characteristics line after the *first* sample's field, so one sample with an extra field shifted every later field — for it alone | two cell types in GSE110554 that were Illumina barcodes | 8 samples mislabelled in GSE61151, 2 in GSE110554, and an exclusion rule written to work around it |
| the blanket `^[^:]+:` prefix strip that positional naming required | it cut into `supplementary_file` URLs in all 49 samples | nothing yet — but it would cut into any field whose value holds a colon |
| assuming the series had the paper's 573 samples | the file declaring 188 | an expectation, not a result |
| decomposing variance over a matrix with missing values | components summing above 1 | a wrong answer about where the signal lives |
| expecting cell-type-variable probes to make a clock vulnerable | the correlation coming out inverted | the obvious hypothesis |
| proposing alignment as the replacement explanation | every permutation p-value above 0.09 | the replacement hypothesis, one stage later |
| restricting stage 4 to complete probes while stage 2 used all of them | the cross-stage reconciliation check | up to 5.7 years of silent disagreement |
| writing the null as `einsum` over a `broadcast_to` view | 27 minutes at 99.5% of one core with no output | half an hour |
| filling gaps with `betas.T.fillna(...).T` | a hang in loading — 485,577 columns after transpose | a second half hour |
| judging stage 6's null result before asking whether it had the power to be anything else | the predicted spread being 0.13 years in two of eight cells | very nearly the wrong conclusion |
| measuring stage 6's power against a noise floor taken from purified cell pellets | a clock whose "noise" exceeded the entire range of the samples it was applied to | one misleading table, caught before it was written down |
| writing stage 16's positive control to require most clocks to detect smoking, when first-generation clocks are documented not to | the check failing at 1 of 4, with the 1 being the only second-generation clock, twice | a rerun — and the irony of committing the exact error the bibliography had just recommended against |
| **stage 12's conclusion that the composition correction does not transfer** | redoing it at twelve types, where it removes 86% of Hannum's composition signal out of cohort | a stage's headline — it was measuring the panel, not the correction |
| nearly reading the second transfer direction from clocks trained on the test cohort | the family clock "going from 0.5% to 40%", which is a vanishing denominator | nothing, caught before it was written down |
| **stage 11's conclusion that the exposure floor is biology** | applying the IntrinClock rule on the naive/memory axis: 41% less displacement at slightly better accuracy | a stage's headline, and a claim repeated to the user in plain language |
| a random control that randomised the age criterion along with the cell-type one | rule (c)'s edge reading +17.1 years, which was too good | an effect size inflated sevenfold, caught before it was written down |
| **recording GSE167998 as unusable ("IDAT only") after checking only `matrix/` and not `suppl/`** | a literature search two weeks later naming the twelve-type panel as the standard | six stages run at half the available cell-type resolution |
| **quoting the composition effect as a share of total clock variance while the literature quotes it as a share of age acceleration** | comparing against Zhang et al. and finding a tenfold gap that was fourfold arithmetic | the project's headline number understated by ~4x for eleven stages |
| a six-type panel that cannot express the naive/memory split, where the largest cell-type age differences live | rebuilding at twelve types and watching every number roughly double | the other ~2x |
| reading stage 9's verdict as a statement about clock design | running the same test backwards, where the family beats both Horvath clocks on both axes | the conclusion that published clocks are better designed — they were better matched to their test range |
| writing stage 10's check 4 on the assumption that large-λ ridge converges to noise | the check failing while the collapsed configurations were nowhere near the frontier | a check rewritten mid-stage, on a premise that was wrong rather than a threshold that was |
| drawing stage 8's frontier with a scale-invariant accuracy axis and a scale-dependent displacement axis | the degenerate corner of the sweep scoring r = 0.885 and zero displacement at once, while predicting a constant | stage 8's practical conclusion — "dilution is a defence" was a model shrinking itself, not resisting anything |
| stage 4 naming Σ\|β\| as the quantity that governs displacement | 45 trained clocks ranking +0.97 with the L2 norm against +0.75 with L1 | the practical conclusion — dilution defends, and stage 4 said concentration was irrelevant |
| writing stage 8's degenerate-corner check with an absolute threshold on a sum whose scale rides on probe count | the check failing at 1.877 against a limit of 1.0, with both real clauses passing | a threshold rewritten after it failed, which is on the record rather than in the history |
| specifying stage 7's physiology check on neutrophils only — 64% of blood and the easiest cell to get right | three of the five unchecked types landing outside clinical range | nothing, because the joint test does not depend on the split; but the per-type reconciliation was read as weak evidence rather than as a broken measurement until this was found |
| **stage 15 attributing the transport failure to the size of the fitting cohort** | stage 17's anchor: 184 samples of GSE40279, fitted and transported the same way, do not reproduce the harm — −0.9% against +4.9 | a mechanism, and the headline of the transport finding |
| reading the n-curve's q75 column as a statement about draws | n = 656 having exactly one possible draw, so its spread is across clocks and cohorts while every other row also carries draw-to-draw variation | nothing, caught before it was written down — it would have moved the reliability threshold from 320 to 656 |
| **stage 17 nominating the conditioning of the composition matrix as the mechanism** | stage 18's check 6: GSE61151 is better conditioned than every GSE40279 subsample, full cohort included, and still did the damage | the second of two candidate mechanisms — both are now down |
| **writing stage 18's check 5 so that it pooled the four clocks inside each n stratum** | Levine's coefficient norm being 386 with a median delta of −0.7% against Horvath 2018's 183 and +4.8% — a between-clock pattern reading as a within-stratum correlation | a finding: pooled, the coefficient norm appeared to predict *less* damage at ρ = −0.484, p = 3e-17. Caught and redone within cells before it was written down |
| **stage 18 concluding that ordinary cross-validation finds the penalty on its own** | stage 19 running it on GSE61151, the cohort that actually broke the transport, where it picks alpha = 0.3 and leaves +2.9% of the damage | the practical half of stage 18's recommendation — penalise yes, trust cross-validation to size it no |
| planning stage 19 around covariate shift in the composition MEAN | the algebra: the fitting cohort's mean enters as a constant vector, and a constant changes no R-squared, so the metric in use since stage 15 cannot see a mean shift at all | a stage, redirected to the second moment before it was written — the stage 7 denominator error in a new costume |
| **three stages hunting a mechanism without writing down the estimation error** | the algebra: Cov(b̂ − b) = σ²/n · Σ_fit⁻¹, so the damage is tr(Σ_fit⁻¹ Σ_test)/n — a standard covariate-shift term that predicts the anchor to within half a point | stages 18 and 19, which measured the wrong quantities and said so |
| **stage 20's check 4 demanding a residual below 2% where the stage's own formula predicts 2.1%** | the check failing at 3.83% on a correct implementation | a hard stop mid-run, and a rewrite — the threshold was arithmetically impossible, not inconvenient |
| stage 20's check 2 testing whether the index predicts a single draw | rho = 0.021 at 56% of cells, worse than the quantity it was meant to beat — and the theory says why: a near-constant predictor cannot track an outcome whose realisation noise has CV 0.43 | nothing, the negative is kept and reported; check 2b is marked post hoc rather than renumbered |
| stage 20 reporting an index of 0.05 as the point where damage crosses zero | stage 21's 63 configurations across four fitting cohorts, where the crossing is a zone from 0.051 to 0.163 containing both outcomes | a sharp threshold, replaced by a one-sided floor: below 0.05 nothing was harmful |
| stage 22's matched-n table selecting its rows out of the grid, when the matched n is almost never on the fitting cohort's grid | five of twelve rows printing NaN — arithmetic on empty frames, caught by the output and by none of the seven pre-specified checks | a rerun; none of the checks verify that the rows they summarise exist |
| **reading stage 22's +0.5% as the cost of panel mismatch** | the sign flipping in four pairs: the mismatched arm also fits with the better panel, so it mixes a mismatch penalty with a panel-quality benefit | the magnitude, which is a mixture and is not quoted as a penalty; the ranking result in check 5 is unaffected |
| **stage 17's permuted reference line, +1.8% from four permutations, and the seven-fold ratio quoted off it** | recomputing it for the paper figure with twenty draws per cohort: median +0.5%, IQR −0.8% to +2.1% | the ratio, replaced by a better statement — meaningless coefficients are harmless, and it is the badly estimated real ones that do the damage |
| **assuming Horvath 2013 had not been trained on any of the four cohorts** | a reviewer asking, and Horvath 2013's Additional file 1 listing GSE40279 as training set 3 | every GSE40279 result carried an in-sample verdict clock since stage 17; corrected in stage 23 — all four headlines survive, the numbers move (n=40 +12.2%→+14.2%, ρ 0.907→0.837, anchor +5.0%→+3.5%, floor 0 of 21→1 of 20) |
| **the permuted reference measured at one n and drawn as a flat line, with 'noise is harmless' and two ratios built on it** | a reviewer; stage 24 measuring it at every n: +15.6% at n = 40, slope −0.93 in log-log | three claims in the synthesis, README and both manuscripts — retracted |
| the transport index modelled as the whole of the transport error | stage 24's P1 failing 0 of 10, and the real fit leaving +4.5% at n = 656 where estimation error is negligible | the index now covers one of two components; the second is n-independent specification error |
| stage 25's metadata check demanding over 100 current smokers in GSE50660 | the cohort having 22; join intact, run stopped on an assumption | a rerun; the smoking contrast became ever vs never, the weaker one |
| **reading stage 21's worst case, +15.9% for 61151 → 42861, as the transport index's failure mode** | stage 25's specification term predicting it (+26.6% / +51.8% against +30.8% / +30.5%) | the worst case is model shift, not covariate shift, and the index cannot see it |
| **the index's safety floor at 0.05 (0 of 21, then 1 of 20 harmful)** | stage 26 counting per clock: 12 of 40 harmful — the earlier counts took medians over clocks, so a harmed clock was outvoted by a helped one | the floor, withdrawn; the index ranks the estimation component and certifies nothing |
| the reviewer's pre-transport net-damage predictor as a usable safety tool | stage 26: sign 71%, rho 0.633, and 30 of the 73 configurations it calls safe are harmful | not adopted; it errs in the reassuring direction, where model shift lives |
| **'ridge at α = 3 holds all twelve pairs'** (stages 21, 23) | stage 27 counting per clock: 29 of 30 cells, one left at +0.2%; and 15 of 30 harmful without a penalty where the pooled count said 3 of 12 | the claim becomes 'reduces 50% to 3%'; α = 10 holds all at the cost of most of the benefit |
| **the index ranking configurations at ρ = 0.907, p = 2×10⁻²⁴** | a reviewer on dependence; stage 28 showing fitting size alone yields ρ ≈ 0.42, and a block permutation giving p = 0.036 | the index's headline strength — real between-pair information, modest, mostly 1/n |
| measuring every transport with the six-type panel, blind to the naive/memory axis | a reviewer; stage 29 measuring on twelve types and on the naive/memory columns: +48.1% and +30.6% at n = 40 where six types saw +11.8% | nothing overturned — the six-type curves understated the damage three- to four-fold |
| reading stage 25's non-significant Wald test (p = 0.61) as Horvath 2018 having no disease-driven model shift | stage 31: within GSE42861, controls → cases leaves +6.3 points beyond the noise term for Horvath 2018 | the differential prediction; model shift is present for all three clocks inside one study |
| treating 'harmful in 93% of draws at n = 40' as a general property of the transported correction | stage 32: for DunedinPACE, fitted on the same cohort, 51% — a coin flip | the small-n headline is clock- and pair-dependent; the penalty and the within-study model shift generalise |
| stage 33a reading GSE132203's age by a prefix match, so the 'age acceleration' field overwrote 'age' | check 4: every clock at r ≈ 0.02 with 'age' while agreeing with each other at 0.88–0.94 | nothing — caught before use; exact field match |
| stage 32 attributing DunedinPACE's small-n neutrality to the fitting cohort (GSE40279) | stage 33: neutral again from GSE132203 | the explanation — it looks like a property of the clock |
| **stage 24's floor read as a general property of transport, from one cohort's top end (n = 656)** | stage 35 fitting on 2,639: the floor persists on median (ratio 0.72), but it is +4.6 into the arthritis cohort and zero into two others; 8 of 13 cells, the pre-set bar was 9 | 'n-independent' stands; 'everywhere' does not — model shift is a property of the pair |
| ridge α = 3 recommended with its cost measured only at matched n (at most 656) | stage 35 at n = 2,639: median −3.0% against −4.0% unpenalised, −7.9% against −12.5% in the best cell | the penalty is a small-n safeguard; at large n it has a price |
| stage 35's reading of the GSE40279 floor as age extrapolation (that cohort runs to 101, the fitting one stops at 75) | stage 36 restricting GSE40279 to 24–75: +3.5 → +3.3 | the hypothesis; the cause of that floor is unknown |
| **α = 3 as a general safeguard for transported corrections** | stage 38 in saliva: 7 of 8 cells still harmful at α = 3 (up to +35%), unpenalised up to +152%, while the correction works at home | the penalty is a blood result; in a tissue whose composition axis dominates it does not hold |
| stage 38 reading the saliva penalty failure as collinearity of the nine-type fit | stage 39: a three-type fit (condition ≈ 1) fails the same way; the immune slope of Horvath 2018 changes sign between cohorts | the mechanism — it is model shift on the dominant axis, which shrinkage cannot reverse |
| **stage 42's 'the difference is the sign'** | a reader: normalised, Horvath 2018's saliva slope loses its opposite sign and 7 of 8 cells stay harmful; Levine 2018 keeps its sign and is still harmful at α = 3 | the sign framing, replaced by spread relative to the mean effect (τ/\|mean\| 0.09–0.55 in blood, 0.91–5.34 in saliva) |
| **the n = 40 median quoted as +18.9 p.p.**, the largest of three sets of 30 draws | stage 24b with 100 draws: +16.1 p.p., harmful in 88% | the headline number; small draw sets are unstable at this size |
| **the floor called 4.5 p.p. and the remainder said to match the shuffled reference 'almost exactly'** | a reader: 4.5 includes the 0.8 the shuffled fit still does at n = 656, and the match fails at n = 40 (20.9 vs 15.6) and n = 80 | floor restated as 3.8 p.p.; agreement quoted with its error (0.6 p.p. from n = 60, 1.2 at n = 40) |
| figure 2 plotting net damage on an axis labelled composition left | the same reader | both series are now composition left |
| 'nearly a factor of two' for the arthritis effect, from two point estimates | stage 31b bootstrap | the phrase keeps its meaning, now with −0.66 (−1.13 to −0.31) |

Three of those returned plausible numbers without crashing, and the loader bug
returned them for four stages before anything noticed. What finally caught it was
not a check — it was reading the metadata of a new dataset closely enough to
notice that two of its "cell types" were barcodes.

## Data

All open, none redistributed. Provenance in
[`reference/README.md`](reference/README.md).
