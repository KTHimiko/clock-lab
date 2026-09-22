# clock-lab — what has been found

## The question

Every longevity company sells a "biological age" test built on DNA methylation.
Blood composition shifts with age, and each cell type carries its own methylation
pattern. **Does a clock read how old the cells are, or who is in the sample?**

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

Three of those returned plausible numbers without crashing, and the loader bug
returned them for four stages before anything noticed. What finally caught it was
not a check — it was reading the metadata of a new dataset closely enough to
notice that two of its "cell types" were barcodes.

## Data

All open, none redistributed. Provenance in
[`reference/README.md`](reference/README.md).
