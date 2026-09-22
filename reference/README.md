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
