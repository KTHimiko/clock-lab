#!/usr/bin/env python3
"""
Read the FlowSorted.BloodExtended.EPIC reference data (GSE167998).

GEO carries this one as a supplementary processed matrix rather than a series
matrix, which is why an earlier pass recorded it as IDAT-only and unusable. It
is not: the processed file has beta values for 68 samples, interleaved with a
detection p-value column per sample, and columns named by sentrix barcode
instead of GSM. The barcodes join to the phenotype table shipped in the
package's `inst/extdata/Pheno.csv`, which carries the cell type of each purified
sample and the known proportions of each reconstructed mixture.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

TYPES12 = ["Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem", "CD8nv",
           "Eos", "Mono", "Neu", "NK", "Treg"]

# how the twelve collapse onto the six of the earlier panel, so the two can be
# compared on the same footing rather than across two different partitions
COLLAPSE = {"CD4T": ["CD4nv", "CD4mem", "Treg"], "CD8T": ["CD8nv", "CD8mem"],
            "Bcell": ["Bnv", "Bmem"], "NK": ["NK"], "Mono": ["Mono"],
            "Neu": ["Neu", "Eos", "Bas"]}


def read_extended(matrix_path, pheno_path, max_probes=None):
    header = pd.read_csv(matrix_path, sep="\t", nrows=0)
    cols = list(header.columns)
    beta_cols = [c for c in cols[1:] if "Detection Pval" not in c]
    betas = pd.read_csv(matrix_path, sep="\t", index_col=0,
                        usecols=[cols[0]] + beta_cols, nrows=max_probes,
                        na_values=["NA", "null", ""], low_memory=False)
    betas = betas.astype(np.float32)

    pheno = pd.read_csv(pheno_path).set_index("Sample_ID")
    keep = [c for c in betas.columns if c in pheno.index]
    return betas[keep], pheno.loc[keep]


if __name__ == "__main__":
    b, p = read_extended(sys.argv[1], sys.argv[2],
                         int(sys.argv[3]) if len(sys.argv) > 3 else None)
    print(f"sondas {b.shape[0]:,}  amostras {b.shape[1]}")
    print(p.CellType.value_counts().to_string())
