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
age for 573 samples.

- GEO: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61151>
- 485,577 probes × 573 samples, 710 MB compressed

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
