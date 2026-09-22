# clock-lab

Do epigenetic clocks measure ageing, or the cell composition of blood?

Every longevity company sells a "biological age" test built on DNA methylation.
The serious, unresolved objection: blood composition shifts with age —
lymphocytes fall, myeloid cells rise — and each cell type carries its own
methylation pattern. A clock fitted on whole blood may be reading **who is in the
sample** rather than **how old those cells are**.

## The design that makes it answerable

Reinius et al. (GSE35069) sequenced **six donors, ten fractions each**: whole
blood, PBMC, granulocytes, CD4+ T, CD8+ T, monocytes, B cells, NK cells,
neutrophils, eosinophils.

Same person, same day, ten measurements. **Any difference in epigenetic age
between those fractions cannot be ageing.** It is composition, and its size
measures the confound directly.

That is a cleaner test than adjusting for composition statistically: it removes
the mixture instead of correcting for it.

## Status

Stage 0 — feasibility — is done, and it already changed the question.

The original plan needed a cohort with methylation *and* measured cell counts.
GSE61151 was the candidate: 573 samples, and the paper reports flow cytometry on
95 of them. **Those counts were never deposited.** Fifteen metadata fields, none
of them cellular. That is the same wall the previous project hit four times — an
aggregated table answers questions about probes, and every question about people
needs the cohort it was built from.

Checking first cost twenty minutes. Finding out later would have cost a night.

## Clocks

Four, all with published open coefficients, none refitted here:

| clock | CpGs | transform |
|---|---|---|
| Horvath 2013 | 353 | inverse-log below adulthood |
| Hannum 2013 | 71 | linear |
| Levine 2018 (PhenoAge) | 513 | linear |
| Horvath 2018 (skin & blood) | 391 | inverse-log below adulthood |

Coefficients come from the `dnaMethyAge` package's data files, which carry the
values published with each paper.

## Layout

```
model/      the clocks
analysis/   numbered stages
scripts/    data loading
reference/  third-party data (not versioned)
```

## Data

All open. Provenance, licence and download for each set in
[`reference/README.md`](reference/README.md). Nothing is redistributed here.
