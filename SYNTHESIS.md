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

---

## Corrections so far

| what was wrong | what caught it | what it cost |
|---|---|---|
| the original question needed cell counts GSE61151 never deposited | checking the metadata before analysing | the question, changed for a better one |
| metadata prefix strip guarded by `dtype == object`, which pandas 3.0 broke | every age coming back NaN | a validation run |
| four samples with ages 0, 6, 7 — `agegap` leaking into `agebloodtaken` | the validation gate refusing to pass | r of 0.89 reading as 0.68 |
| assuming the series had the paper's 573 samples | the file declaring 188 | an expectation, not a result |
| decomposing variance over a matrix with missing values | components summing above 1 | a wrong answer about where the signal lives |
| expecting cell-type-variable probes to make a clock vulnerable | the correlation coming out inverted | the obvious hypothesis |
| proposing alignment as the replacement explanation | every permutation p-value above 0.09 | the replacement hypothesis, one stage later |
| restricting stage 4 to complete probes while stage 2 used all of them | the cross-stage reconciliation check | up to 5.7 years of silent disagreement |
| writing the null as `einsum` over a `broadcast_to` view | 27 minutes at 99.5% of one core with no output | half an hour |
| filling gaps with `betas.T.fillna(...).T` | a hang in loading — 485,577 columns after transpose | a second half hour |

Two of those returned plausible numbers without crashing.

## Data

All open, none redistributed. Provenance in
[`reference/README.md`](reference/README.md).
