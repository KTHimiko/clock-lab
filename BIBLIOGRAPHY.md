# Bibliography

What the literature says about each claim this project makes, including where it
corrected the project. Marked **[full]** where the paper was read in full,
**[abstract]** where only the abstract, summary or preprint was reachable.

---

## The objection itself

**Jaffe AE, Irizarry RA (2014).** Accounting for cellular heterogeneity is
critical in epigenome-wide association studies. *Genome Biology* 15:R31.
[abstract] — <https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r31>

The founding statement of the problem, ten years before the current wave. Blood
is a heterogeneous collection of cell types each with a very different
methylation profile, and an age association measured in bulk blood may be a
composition association.

---

## How large the effect is — where this project was wrong

**Zhang Z, Salas LA, Christensen BC, et al. (2024).** Deciphering the role of
immune cell composition in epigenetic age acceleration. *Aging Cell* 23:e14071.
[full] — <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10928575/>

10,147 blood samples, twelve-type deconvolution. Immune composition explains
**13% (Horvath), 25% (Hannum), 33.6% (PhenoAge)** of age-acceleration variance,
ranking second after chronological age and first for Hannum.

> **This corrected stage 7.** This project had been quoting 0.9–2.7%, which was
> a different denominator (share of total clock variance, not of age
> acceleration) measured on a six-type panel. Stage 13 redid it and landed on
> 15.6% / 32.2% / 23.4%.

**Tong H, Teschendorff AE, et al.** Immune cell composition is an important
contributor to epigenetic age variation. [full] —
<https://pmc.ncbi.nlm.nih.gov/articles/PMC11980466/>

Reports that roughly **39% of the age-prediction accuracy of blood** is driven
by lymphocyte composition, and that naive CD4 proportion is the strongest single
contributor for Hannum, PhenoAge and DunedinPACE. Also notes that clocks still
reflect cell-intrinsic ageing independent of composition — which is the position
this project ended at.

---

## The reference panel this project used

**Salas LA, Koestler DC, Christensen BC, et al. (2022).** Enhanced cell
deconvolution of peripheral blood using DNA methylation for high-resolution
immune profiling. *Nature Communications* 13:761. [full] —
<https://pmc.ncbi.nlm.nih.gov/articles/PMC8828780/>

The twelve-type panel: neutrophils, eosinophils, basophils, monocytes, naive and
memory B, naive and memory CD4 T, naive and memory CD8 T, NK, Treg. Reference
data deposited as **GSE167998**, which stages 13–15 build from directly.

Two points that bear on this project:

- They select probes with the **IDOL** iterative algorithm and reach an average
  R² near 1 on artificial mixtures with a 1,200-probe library. Stages 13–15 use
  a simpler t-statistic selection over 600 probes and reach r = 0.79, mean
  absolute error 0.027. **The panel here is cruder than the published one**, and
  every twelve-type number in this repository should be read as a lower bound on
  what the method can do.
- They state that "distinguishing hierarchically close cell types such as CD4mem
  and Tregs is challenging as the pool of potential specific markers is
  limited", and that residual heterogeneity remains in the CD8 memory
  compartment. **That is exactly where this project's panel has its largest
  errors** — CD4mem 0.057, CD4nv 0.050, Treg 0.043, against 0.002 for
  neutrophils. Independent confirmation of the specific weak points.

**Salas LA, et al. (2018).** An optimized library for reference-based
deconvolution of whole-blood biospecimens assayed using the Illumina EPIC
BeadArray. *Genome Biology* 19:64. [abstract] —
<https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1448-7>

**Koestler DC, et al. (2016).** Improving cell mixture deconvolution by
identifying optimal DNA methylation libraries (IDOL). *BMC Bioinformatics*
17:120. [abstract]

---

## Designing a clock around the confound — where this project was wrong again

**Tomusiak A, et al. (2024).** Development of an epigenetic clock resistant to
changes in immune cell composition. *Communications Biology* 7:934. [abstract;
full text behind a redirect] — <https://www.nature.com/articles/s42003-024-06609-4>
· preprint <https://www.biorxiv.org/content/10.1101/2023.03.01.530561v1>

IntrinClock, selected to be invariant across ten immune cell types. Published
selection rule: keep CpGs with **|r| > 0.3 against chronological age and
|r| < 0.3 against a sample being naive CD8**. Reports naive CD8 T cells reading
15–20 years younger than effector memory CD8 from the same individual.

> **This corrected stage 11**, which had concluded the exposure floor was
> biology. Stage 14 implemented the published rule and reproduced the result:
> displacement falls 41% while accuracy slightly improves. Stage 11's filter was
> built from six cell types and could not express the axis the effect lives on.
>
> The naive/memory gap also reproduced, and larger: **+11.0 (Horvath 2013),
> +40.6 (Hannum), +46.4 (Levine)** years between memory and naive CD8.

---

## Correcting for composition after the fact — where stage 15 lands

**Horvath S, et al.** Intrinsic epigenetic age acceleration (IEAA): the residual
of clock age regressed on chronological age *and* blood cell estimates.
[abstract] — the standard correction, and precisely what stage 15 tests by
transporting the composition coefficients between cohorts.

**Meredith JM, Melton PE, Huang RC, et al. (2019).** In epigenomic studies,
including cell-type adjustments in regression models can introduce
multicollinearity, resulting in apparent reversal of direction of association.
*Frontiers in Genetics* 10:816. [full] —
<https://pmc.ncbi.nlm.nih.gov/articles/PMC6746958/>

812 adolescents. Without cell-type adjustment, lower CDKN2A methylation
associated with higher BMI (β = −0.004, p = 0.003); adding six cell types
**reversed the sign** (β = +0.004). Granulocyte proportion correlated with
methylation at r = −0.783, giving a variance inflation factor of **113.7**
against an acceptable threshold of 5. In simulation, at that correlation 83% of
coefficients flipped sign.

> **This is the closest prior work to stage 15's negative direction.** Stage 15
> finds that a correction fitted on 188 samples and applied to 656 *injects*
> composition signal rather than removing it — Horvath 2013 from 4.7% to 10.7%.
> Meredith et al. give the mechanism: the composition predictors are collinear
> with the methylation being corrected, so the coefficients are unstable, and an
> unstable coefficient carried to a new cohort adds variance.

**Zheng SC, Breeze CE, Beck S, Teschendorff AE (2018).** Identification of
differentially methylated cell types in epigenome-wide association studies /
comparative analysis of cell-type adjustment methods. *Briefings in
Bioinformatics*. [abstract] —
<https://academic.oup.com/bib/advance-article/doi/10.1093/bib/bby068/5066710>

**Cell-type deconvolution in EWAS: a review and recommendations.**
*Epigenomics* (2017). [abstract] —
<https://www.tandfonline.com/doi/full/10.2217/epi-2016-0153>

Standing recommendation: adjust with measured counts where available; estimated
counts are acceptable for exposures other than age. Note the exception — **age
is the phenotype where cell counts matter most**, which is the phenotype this
whole project is about.

---

## Independent confirmation of stage 2

**Epigenetic age estimates in blood shift by decades with cell population: a
within-donor variance partition across two platforms and two study designs.**
*Research Square* preprint, 29 July 2026. [abstract] —
<https://www.researchsquare.com/article/rs-10503047/v1>

Reports within-donor medians of **15.6 years (Horvath), 36.7 (Hannum), 36.7
(PhenoAge)**, reaching 54.0 in one donor. Stage 2 of this project, run
independently, found 16.0 / 35.4 / 37.0 on the same design.

---

## Related, not load-bearing here

- **Epigenetic clocks and inflammaging: pitfalls caused by ignoring cell-type
  heterogeneity.** *GeroScience* (2025). —
  <https://link.springer.com/article/10.1007/s11357-025-01677-8>
- **Dissecting the impact of differentiation stage, replicative history, and
  cell type composition on epigenetic clocks.** *Stem Cell Reports* (2024). —
  <https://www.cell.com/stem-cell-reports/fulltext/S2213-6711(24)00218-2>
- **Guidelines on optimizing DNA methylation reference panels for cell-type
  deconvolution.** *Communications Biology* (2026). —
  <https://www.nature.com/articles/s42003-026-09745-1>

---

## A note on method

Both of this project's reversals — stages 7 and 11 — came from reading the
literature, not from its own pre-specified checks. Those checks caught six
implementation bugs and not one framing error, because an internal check cannot
know the whole ruler is wrong. **Search the literature at the start of a block,
not at the end of one.** Twenty minutes before stage 7 would have saved six
stages run at half the available cell-type resolution.
