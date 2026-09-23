# clock-lab — what has been found

## The question

Every longevity company sells a "biological age" test built on DNA methylation.
Blood composition shifts with age, and each cell type carries its own methylation
pattern. **Does a clock read how old the cells are, or who is in the sample?**

## The answer, after nineteen stages

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

It can also be *subtracted* away after the fact, but only under conditions
narrow enough to matter. With a twelve-type panel and the full 656-sample
fitting cohort, the correction removes 61% of the composition signal out of
cohort. Fitted on forty samples it **adds 2.7 times what was there** — seven
times more damage than composition coefficients with no information in them —
and it is worse than doing nothing in 89% of draws. The median turns beneficial
somewhere around 160 to 184 samples and an unlucky draw keeps hurting to about
320. Sample size is not the whole story: 184 samples of one cohort help while
184 of another did the damage that stage 15 blamed on size — and conditioning,
the other candidate, turned out to carry no information about it at all.

**The failure can be removed without being explained.** Ridge-penalising the
composition coefficients, with the penalty chosen by ordinary cross-validation
on the fitting cohort, holds the median delta negative at every fitting size
from forty samples to six hundred: the +12.2% catastrophe at n = 40 becomes
−0.2%, and the benefit at full n is kept. No penalised adjustment appears among
the eight methods the field benchmarks.

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

### Seven times worse than no information at all

The stage carries a reference line: the same transport, with composition
coefficients fitted at full n on **permuted** composition. Those coefficients
carry no information by construction.

> permuted coefficients, n = 656: median delta **+1.8%** (range −3.8% to +4.2%)
> real coefficients, n = 40: median delta **+12.2%**

A correction fitted on forty samples is **seven times more damaging than a
coefficient vector with nothing in it**. That is the part worth keeping. Small-n
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
stage). What is left is a cohort difference that neither quantity captures, and
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
damage unaccounted for. What can be said is that composition geometry is the
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

Three of those returned plausible numbers without crashing, and the loader bug
returned them for four stages before anything noticed. What finally caught it was
not a check — it was reading the metadata of a new dataset closely enough to
notice that two of its "cell types" were barcodes.

## Data

All open, none redistributed. Provenance in
[`reference/README.md`](reference/README.md).
