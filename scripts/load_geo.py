#!/usr/bin/env python3
"""
Read a GEO series matrix into a beta matrix plus sample metadata.

Series matrix files are a header of !Sample_* lines followed by a probe-by-
sample table. They are large — GSE61151 is 710 MB compressed — so the table is
read in chunks and stored as float32, which halves the memory and is well below
the precision that matters for a beta value bounded in [0, 1].
"""
import gzip
import sys
from pathlib import Path
import numpy as np
import pandas as pd


def read_series_matrix(path, max_probes=None):
    path = Path(path)
    meta_rows, samples = {}, None
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            key = parts[0][8:]
            vals = [v.strip('"') for v in parts[1:]]
            if key == "geo_accession":
                samples = vals
            elif key == "characteristics_ch1":
                # each characteristics line is one field, named by its own prefix
                name = vals[0].split(":")[0].strip() if vals else "unknown"
                meta_rows.setdefault(name, vals)
            else:
                meta_rows.setdefault(key, vals)

    meta = pd.DataFrame({k: pd.Series(v) for k, v in meta_rows.items()})
    if samples is not None:
        meta.insert(0, "gsm", pd.Series(samples))
    for c in meta.columns:
        if meta[c].dtype == object:
            meta[c] = meta[c].astype(str).str.replace(r"^[^:]+:\s*", "", regex=True)

    betas = pd.read_csv(path, sep="\t", comment="!", index_col=0,
                        nrows=max_probes, na_values=["NA", "null", ""],
                        low_memory=False)
    betas = betas.astype(np.float32)
    return betas, meta


if __name__ == "__main__":
    b, m = read_series_matrix(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else None)
    print(f"sondas {b.shape[0]:,}  amostras {b.shape[1]}  "
          f"memoria {b.memory_usage(deep=True).sum()/1e6:.0f} MB")
    print(f"campos de metadados: {list(m.columns)}")
