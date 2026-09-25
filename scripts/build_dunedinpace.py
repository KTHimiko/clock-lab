#!/usr/bin/env python3
"""
Build reference/data/clocks/dunedinpace.npz, which model/dunedinpace.py reads.

The source is DunedinPACE.rda from the dnaMethyAge repository (see
reference/README.md), holding two objects:
  coefs                a table of probe and coefficient, whose first row is the
                       intercept
  gold_standard_means  the Dunedin Study's mean beta at 20,000 background probes,
                       the target of the quantile normalisation

It is read with the `rdata` package; `pyreadr` cannot parse it. The model
probes' means are the gold-standard means at those probes.

Usage:  .venv/bin/python scripts/build_dunedinpace.py
"""
from pathlib import Path
import numpy as np
import rdata

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "reference/data/clocks/DunedinPACE.rda"
OUT = ROOT / "reference/data/clocks/dunedinpace.npz"

r = rdata.read_rda(str(SRC))
coefs, gold = r["coefs"], r["gold_standard_means"]
is_int = coefs["Probe"].astype(str) == "Intercept"
assert is_int.sum() == 1, "esperava exatamente uma linha de intercepto"
intercept = float(coefs.loc[is_int, "Coefficient"].iloc[0])
body = coefs.loc[~is_int]
probes = np.array([str(x) for x in body["Probe"]])   # dtype <U, loadable without pickle
weights = body["Coefficient"].astype(float).to_numpy()
gold_probes = np.array([str(x) for x in gold.coords[gold.dims[0]].values])
gold_means = gold.values.astype(float)
pos = {p: i for i, p in enumerate(gold_probes)}
missing = [p for p in probes if p not in pos]
assert not missing, f"sondas do modelo fora do fundo: {missing[:5]}"
model_means = gold_means[[pos[p] for p in probes]]
np.savez(OUT, model_probes=probes, weights=weights, weight_probes=probes,
         intercept=np.array(intercept), model_means=model_means,
         gold_probes=gold_probes, gold_means=gold_means)
print(f"escrito {OUT.name}: {len(probes)} sondas do modelo, {len(gold_probes)} de fundo, "
      f"intercepto {intercept:+.6f}")
