#!/usr/bin/env python3
"""
Stage 32a — DunedinPACE for the four cohorts, validated before it is used.

DunedinPACE measures pace of ageing, not age, and was trained on the Dunedin
Study, which is not on GEO — so none of the four cohorts here is in its training
data. Its values are cached next to the stage 18 cache so later stages need not
reload the series.

PRE-REGISTERED VALIDATION (the implementation is ours, from the package's model
data, so it must reproduce the known behaviour before anything is built on it)
  1. background coverage: at least 80% of the 20,000 reference probes and of the
     173 model probes present in every cohort. Hard stop.
  2. scale: the cohort mean of DunedinPACE between 0.8 and 1.2 and its standard
     deviation between 0.04 and 0.2, in every cohort — the published scale is
     about one year of biological ageing per calendar year.
  3. POSITIVE CONTROL: current smokers faster than never smokers in GSE50660,
     one-sided p < 0.05, adjusted for age. Smoking is the most replicated
     correlate of DunedinPACE.

Usage:  .venv/bin/python analysis/32a_dunedinpace_cache.py
"""
import sys, gzip
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix
from model.dunedinpace import predict

CACHE = ROOT / "results" / "cache"
out, cov = {}, {}
for tag in ["GSE40279", "GSE61151", "GSE50660", "GSE42861"]:
    print(f"carregando {tag} ...", flush=True)
    b, _ = read_series_matrix(ROOT / "reference/data" / f"{tag}_series_matrix.txt.gz")
    s, c = predict(b)
    out[tag], cov[tag] = s, c
    s.to_frame().to_csv(CACHE / f"{tag}_dpace.csv")
    print(f"  cobertura de fundo {c:.1%} | media {s.mean():.3f} dp {s.std():.3f}", flush=True)
    del b

print("\nCHECAGEM 1 — cobertura >= 80% em todas:",
      "ok" if all(v >= 0.8 for v in cov.values()) else "FALHOU")
c2 = all(0.8 <= s.mean() <= 1.2 and 0.04 <= s.std() <= 0.2 for s in out.values())
print("CHECAGEM 2 — escala (media 0.8-1.2, dp 0.04-0.2):", "ok" if c2 else "FALHOU")

samples, smk = None, None
with gzip.open(ROOT / "reference/data/GSE50660_series_matrix.txt.gz", "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("!series_matrix_table_begin"):
            break
        p = [v.strip('"') for v in line.rstrip("\n").split("\t")]
        if p[0] == "!Sample_geo_accession":
            samples = p[1:]
        elif p[0] == "!Sample_characteristics_ch1" and p[1].startswith("smoking"):
            smk = [float(x.partition(":")[2]) for x in p[1:]]
        elif p[0] == "!Sample_characteristics_ch1" and p[1].startswith("age"):
            age = [float(x.partition(":")[2]) for x in p[1:]]
df = pd.DataFrame({"smk": smk, "age": age}, index=samples).join(out["GSE50660"])
g = df[df.smk.isin([0, 2])]
X = np.column_stack([np.ones(len(g)), g.age, (g.smk == 2).astype(float)])
beta, res, *_ = np.linalg.lstsq(X, g.DunedinPACE, rcond=None)
resid = g.DunedinPACE - X @ beta
se = np.sqrt(resid.var(ddof=3) * np.linalg.inv(X.T @ X)[2, 2])
t = beta[2] / se; p = stats.t.sf(t, len(g) - 3)
c3 = beta[2] > 0 and p < 0.05
print(f"CHECAGEM 3 — fumante atual - nunca: {beta[2]:+.3f} (p unilateral {p:.2g}) ->",
      "ok" if c3 else "FALHOU")
