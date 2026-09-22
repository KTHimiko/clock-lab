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


def _line_name(vals):
    """Name a characteristics line by the prefix most of its cells agree on."""
    names = [c.split(":")[0].strip() for c in vals if ":" in c]
    return max(set(names), key=names.count) if names else "unknown"


def read_series_matrix(path, max_probes=None):
    path = Path(path)
    meta_rows, char_rows, samples = {}, [], None
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
                # Name each cell by ITS OWN prefix rather than by the first
                # sample's. GEO writes one characteristics line per field only
                # when every sample carries the same fields in the same order;
                # a sample with one extra field shifts every later field by one
                # column, for that sample alone. Naming the whole line after
                # vals[0] then files those shifted values under the wrong names
                # and nothing complains. In GSE110554 two CD4T samples carry a
                # duplicated `sample_name`, which put their barcode under
                # `cell type` and their cell type under `cd4t`.
                char_rows.append(vals)
            else:
                meta_rows.setdefault(key, vals)

    # Split each characteristics cell on its own first colon, so the field name
    # travels with the value instead of with the column position. This also
    # replaces the blanket prefix strip that used to run over every column: that
    # strip was needed because the names came from position, and it would have
    # cut into any ordinary field whose value happened to contain a colon.
    n = len(samples) if samples is not None else max((len(v) for v in char_rows), default=0)
    chars = {}
    for vals in char_rows:
        for j, cell in enumerate(vals[:n]):
            name, sep, value = cell.partition(":")
            if not sep:
                # no colon: keep the line's own shape, filed under the line's
                # dominant name, so nothing is silently dropped
                name, value = _line_name(vals), cell
            name = name.strip()
            chars.setdefault(name, [None] * n)
            if chars[name][j] is None:
                chars[name][j] = value.strip()

    meta = pd.DataFrame({k: pd.Series(v) for k, v in meta_rows.items()})
    for k, v in chars.items():
        if k not in meta.columns:
            meta[k] = pd.Series(v, dtype=object)
    if samples is not None:
        meta.insert(0, "gsm", pd.Series(samples))

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
