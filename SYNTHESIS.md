# clock-lab — what has been found

## The question

Every longevity company sells a "biological age" test built on DNA methylation.
Blood composition shifts with age, and each cell type carries its own methylation
pattern. **Does a clock read how old the cells are, or who is in the sample?**

## The answer, after twelve stages

**Both, and the proportions matter more than either camp says.**

In sorted cells the effect is enormous: up to 35 years between fractions of one
man's blood drawn on one day, replicated in a second cohort on a second array.
In actual blood from actual people it is small — composition explains **0.9% to
2.7%** of epigenetic age beyond chronological age, worth about **1 to 1.6 years
of standard deviation** and 6 to 9 years between the extremes of a cohort.

It cannot be designed away: two thirds of the best age-tracking probes in the
genome sit in the quarter that varies most between cell types, so removing one
removes the other, and no filtered clock ever beat an unfiltered one. It cannot
reliably be subtracted away either: a correction fitted on one cohort removes
most of nothing in another, because the deconvolution's split between
neighbouring cell types does not transport even when its total does.

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
probes that vary most between cell types.** That is the floor. The selection is
not being careless — the probes that mark time are, substantially, the probes
that mark which cell you are looking at. Removing one removes the other.

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

---

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
| **recording GSE167998 as unusable ("IDAT only") after checking only `matrix/` and not `suppl/`** | a literature search two weeks later naming the twelve-type panel as the standard | six stages run at half the available cell-type resolution |
| **quoting the composition effect as a share of total clock variance while the literature quotes it as a share of age acceleration** | comparing against Zhang et al. and finding a tenfold gap that was fourfold arithmetic | the project's headline number understated by ~4x for eleven stages |
| a six-type panel that cannot express the naive/memory split, where the largest cell-type age differences live | rebuilding at twelve types and watching every number roughly double | the other ~2x |
| reading stage 9's verdict as a statement about clock design | running the same test backwards, where the family beats both Horvath clocks on both axes | the conclusion that published clocks are better designed — they were better matched to their test range |
| writing stage 10's check 4 on the assumption that large-λ ridge converges to noise | the check failing while the collapsed configurations were nowhere near the frontier | a check rewritten mid-stage, on a premise that was wrong rather than a threshold that was |
| drawing stage 8's frontier with a scale-invariant accuracy axis and a scale-dependent displacement axis | the degenerate corner of the sweep scoring r = 0.885 and zero displacement at once, while predicting a constant | stage 8's practical conclusion — "dilution is a defence" was a model shrinking itself, not resisting anything |
| stage 4 naming Σ\|β\| as the quantity that governs displacement | 45 trained clocks ranking +0.97 with the L2 norm against +0.75 with L1 | the practical conclusion — dilution defends, and stage 4 said concentration was irrelevant |
| writing stage 8's degenerate-corner check with an absolute threshold on a sum whose scale rides on probe count | the check failing at 1.877 against a limit of 1.0, with both real clauses passing | a threshold rewritten after it failed, which is on the record rather than in the history |
| specifying stage 7's physiology check on neutrophils only — 64% of blood and the easiest cell to get right | three of the five unchecked types landing outside clinical range | nothing, because the joint test does not depend on the split; but the per-type reconciliation was read as weak evidence rather than as a broken measurement until this was found |

Three of those returned plausible numbers without crashing, and the loader bug
returned them for four stages before anything noticed. What finally caught it was
not a check — it was reading the metadata of a new dataset closely enough to
notice that two of its "cell types" were barcodes.

## Data

All open, none redistributed. Provenance in
[`reference/README.md`](reference/README.md).
