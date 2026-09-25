# clock-lab

Do epigenetic clocks measure ageing, or the cell composition of blood?

Every longevity company sells a "biological age" test built on DNA methylation.
The serious, unresolved objection: blood composition shifts with age —
lymphocytes fall, myeloid cells rise — and each cell type carries its own
methylation pattern. A clock fitted on whole blood may be reading **who is in
the sample** rather than **how old those cells are**.

## The answer, after forty-seven stages

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
its accuracy on an external cohort.

It can be **subtracted** inside the cohort where the correction is fitted — the
usual practice, and one that does what it claims. **Carried to another cohort,
the correction can add the confounding it was meant to remove.** Fitted on forty
samples and transported, it is worse than doing nothing in 88% of draws, and the
damage measured on the naive/memory axis is three to four times what a six-type
panel registers.

The harm has two parts. Estimation noise, amplified by the difference between
the two cohorts' composition covariance — a textbook covariate-shift term, which
coefficients carrying no information reproduce almost exactly. And **model
shift**: the composition effect itself differs between cohorts, does not shrink
with more fitting data (tested up to 2,639 fitting samples), and is behind the
worst transports found. It belongs to the pair of cohorts: large into the
arthritis cohort, absent into two others.

**No quantity computed in advance certifies a transport as safe**, because the
second part needs the target cohort's own coefficients. **What works is a
penalty**: ridge-shrinking the coefficients takes harmful transports from half of
all (pair × clock) cells to 3%, at the cost of some benefit where none was at
risk. With thousands of fitting samples it costs about a point of benefit; a
penalty that fades with n avoids that cost but lets model shift through (11 of 72
cells harmful, against 1 of 72). **In saliva, where composition dominates, the penalty does
not hold:** 7 of 8 transports between saliva cohorts stayed harmful at α = 3.

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

DunedinPACE, a pace-of-ageing clock of a different kind, is reimplemented in
`model/dunedinpace.py` from the model data its package publishes, and validated
against published behaviour before use (stage 32).

## Layout

```
model/      the clocks and the deconvolution
analysis/   numbered stages, one question each
scripts/    data loading
reference/  third-party data provenance (the data itself is not versioned)
```

Analysis outputs land in `results/` and are not versioned either — every one of
them is reproducible from the scripts and the download commands.

## Paper

`paper/manuscript.typ` (English, the reference version) and
`paper/manuscript_pt.typ` (Portuguese), compiled with Typst; figures are
regenerated by `paper/figures.py` from the stage outputs.

## How this was made

This project was carried out with extensive use of an AI system, Claude
(Anthropic), through Claude Code. The AI wrote the code, designed and ran the
analyses — each with its checks written into the docstring before the result
was read — searched the literature and drafted the text. The author set the
aims, chose which questions to pursue and what to report, brought in critical
readings of the drafts, and is responsible for the content. The corrections
table in `SYNTHESIS.md` records every claim that was later overturned, whoever
caught it.

## Licence

Code: MIT (see `LICENSE`). Text, synthesis and figures: CC BY 4.0. The data are
third-party and are not redistributed; `reference/README.md` records where each
series comes from.
