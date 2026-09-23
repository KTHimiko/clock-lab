# clock-lab

Do epigenetic clocks measure ageing, or the cell composition of blood?

Every longevity company sells a "biological age" test built on DNA methylation.
The serious, unresolved objection: blood composition shifts with age —
lymphocytes fall, myeloid cells rise — and each cell type carries its own
methylation pattern. A clock fitted on whole blood may be reading **who is in
the sample** rather than **how old those cells are**.

## The answer, after seventeen stages

**Both, and the proportions matter more than either camp says.**

In sorted cells the effect is enormous: up to 35 years between fractions of one
donor's blood drawn on one day, replicated on a second cohort and a second
array — and that was measured on a panel too coarse to see the largest component
of it. A naive CD8 T cell reads 47 years below expectation where a memory CD8
reads near zero.

In whole blood from real people, composition explains **10% to 32% of age
acceleration**, depending on the clock, which lands on the published figures.

It can be **partly designed away**: excluding CpGs that correlate with naive-CD8
identity cuts a clock's cell-type displacement by 41% while slightly improving
its accuracy on an external cohort. It can be **subtracted** away after the fact
only with a fine panel and a large fitting cohort: fitted on 656 samples the
correction removes 61% of the composition signal in a cohort it has never seen,
while the same correction fitted on forty *adds* nearly three times what was
there — and is worse than doing nothing in 89% of draws.

Full argument, every number and every correction:
[`SYNTHESIS.md`](SYNTHESIS.md).

## How it is built

Each stage is one script in `analysis/`, numbered, with its **sanity checks
written into the docstring before the result is read**. Where a check failed,
what happened to it is on the page rather than in the history — three checks in
this project were revisited after failing, and each of those is recorded with
the reasoning that justified it.

`SYNTHESIS.md` carries a corrections table. It is the most useful thing in the
repository: several stages reached conclusions that later stages overturned, and
the overturned sections are **amended in place** rather than rewritten, so the
path stays readable.

## Data

Four cohorts, all public, none redistributed here. Provenance, sizes, checksums
and download commands in [`reference/README.md`](reference/README.md).

| accession | what it is |
|---|---|
| GSE35069 | 6 donors × 10 purified fractions — the within-donor design |
| GSE61151 | 188 whole blood samples with age, 35–83 |
| GSE110554 | EPIC: 37 purified + 12 mixtures with known proportions |
| GSE40279 | 656 whole blood samples with age, 19–101 |
| GSE167998 | 12 leukocyte subtypes, naive and memory split |

## Clocks

Four, all with published open coefficients, none refitted except where a stage
says so explicitly:

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
model/      the clocks and the deconvolution
analysis/   numbered stages, one question each
scripts/    data loading
reference/  third-party data provenance (the data itself is not versioned)
```

Analysis outputs land in `results/` and are not versioned either — every one of
them is reproducible from the scripts and the download commands.
