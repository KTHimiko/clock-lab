# Third-party data — provenance

Nothing in this directory is versioned. Every set is public and downloadable
from the source below.

## GSE35069 — Reinius et al. 2012, purified blood cell types

The project's main dataset. Six donors, ten fractions each, Illumina 450K.

- GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE35069>
- Paper: Reinius LE et al., *Differential DNA methylation in purified human
  blood cells*, PLoS ONE 2012;7:e41361
- 485,577 probes × 60 samples, 118 MB compressed

```bash
curl -sLO "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE35nnn/GSE35069/matrix/GSE35069_series_matrix.txt.gz"
```

Fractions: whole blood, PBMC, granulocytes, CD4+ T, CD8+ T, CD14+ monocytes,
CD19+ B, CD56+ NK, neutrophils, eosinophils — six of each.

## GSE61151 — Breakthrough Generations Study, whole blood

Used for validating the clock implementations, because it carries chronological
age. The paper reports 573; the series matrix carries 188, and that
discrepancy is on the corrections table in `SYNTHESIS.md`.

- GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61151>
- 485,577 probes × 188 samples, 710 MB compressed

```bash
curl -sLO "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE61nnn/GSE61151/matrix/GSE61151_series_matrix.txt.gz"
```

**The white blood cell counts are not here.** The paper reports automated flow
cytometry on 95 samples; the values were not deposited with the series. The
metadata carries gender, age, plate, well, chip and batch — nothing cellular.
Recorded so the check is not repeated.

## GSE40279 — Hannum et al., 656 whole blood samples with age

Downloaded 2026-09-22 for stage 9, as an **external** cohort: a clock trained on
GSE61151 and tested here has never seen any of these people.

- GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40279>
- 473,034 probes × 656 samples, 1,237,679,628 bytes compressed
- md5 `5506575f4598dd4b88cbcdbdd90ac826`
- ages 19 to 101, whole blood, two collection sites

```bash
curl -sLO "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40279/matrix/GSE40279_series_matrix.txt.gz"
```

**This is the cohort the Hannum clock was trained on.** So it is a fair external
test for anything trained elsewhere, and an in-sample one for Hannum. Any table
that puts Hannum's accuracy here beside another clock's is comparing a memory
with a prediction, and says so.

## GSE167998 — FlowSorted.BloodExtended.EPIC, twelve leukocyte subtypes

Downloaded 2026-09-22 for stage 13. 56 purified samples across twelve types —
the six of the earlier panel split into naive and memory lymphocytes, plus
regulatory T cells, eosinophils and basophils — and twelve reconstructed
mixtures whose twelve-way proportions are known.

- GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167998>
- 865,859 probes × 68 samples, 819,745,924 bytes compressed
- md5 `160eef349c185198ff3c2f3a1ff2252c`

```bash
curl -sLO "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167998/suppl/GSE167998_matrix_processed.txt.gz"
curl -sLO "https://raw.githubusercontent.com/immunomethylomics/FlowSorted.BloodExtended.EPIC/master/inst/extdata/Pheno.csv"
```

**An earlier pass recorded this accession as unusable — "IDAT only, would need
R/minfi".** That was wrong, and it cost the project six stages of analysis at
half the available cell-type resolution. GEO carries the processed betas as a
*supplementary* file rather than as a series matrix, so a check that only looked
at `matrix/` found nothing. The lesson is in the corrections table: when an
accession looks empty, list `suppl/` before writing it off.

The phenotype table is the package's own `inst/extdata/Pheno.csv`, which joins
to the matrix on the sentrix barcode and carries each purified sample's cell
type and each mixture's known proportions.

## GSE42861 — rheumatoid arthritis, 689 peripheral blood leukocyte samples

Downloaded 2026-09-22 for stage 16. Carries **both** variables that stage needs
in the same people, same batch, same platform: a hard outcome and the positive
control.

- GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42861>
- 450K (GPL13534), 689 samples, 2,698,668,510 bytes compressed
- md5 `413f7958e443b89dea7cd3f273c47691`
- characteristics: `disease state` (rheumatoid arthritis / normal), `subject`
  (Patient / Normal), `smoking status` (never / ex / current), `age`, `gender`

```bash
curl -sLO "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42861/matrix/GSE42861_series_matrix.txt.gz"
```

**Two confounds to carry, not to forget.** Smoking is a risk factor for
rheumatoid arthritis, so the two variables are correlated inside this cohort and
each test must adjust for the other. And RA patients are medicated —
methotrexate, steroids — which is its own methylation exposure. This cohort can
say whether a clock separates cases from controls; it cannot say the separation
is the disease rather than its treatment.

## GSE50660 — smoking, 464 peripheral blood samples

Downloaded 2026-09-22 for stage 16, as an independent replication of the
positive control alone.

- GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50660>
- 450K (GPL13534), 464 samples, 604,598,383 bytes compressed
- md5 `b7a926ba9814a573dea5a6383c09001a`
- characteristics: `smoking` coded 0/1/2 for never/former/current, `age`,
  `gender`

```bash
curl -sLO "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50660/matrix/GSE50660_series_matrix.txt.gz"
```

The never/former/current coding in both cohorts allows a gradient rather than a
two-group split, which is harder to produce by accident.

## Clock coefficients — `clocks/`

From the `dnaMethyAge` R package, which stores the values published with each
original paper.

- Repository: <https://github.com/yiluyucheng/dnaMethyAge>
- Files: `HorvathS2013.rda`, `HannumG2013.rda`, `LevineM2018.rda`,
  `HorvathS2018.rda`

```bash
curl -sLO "https://raw.githubusercontent.com/yiluyucheng/dnaMethyAge/main/data/HorvathS2013.rda"
```

Read with `pyreadr`; no R needed. `DunedinPACE.rda` fails to parse and is not
used.
