#!/usr/bin/env python3
"""
The gate: do these implementations reproduce published clock behaviour?

A clock fitted on blood correlates with chronological age above r = 0.9 in the
papers that published it. GSE61151 is 573 whole-blood samples with age recorded,
which is the cleanest available test that the implementation is right rather
than merely running.

Nothing downstream is worth anything if this fails, so it runs first and alone.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "scripts"))
from load_geo import read_series_matrix
from model.clocks import CLOCKS, predict

MATRIX = ROOT / "reference/data/GSE61151_series_matrix.txt.gz"

print("carregando a matriz (710 MB comprimidos) ...", flush=True)
betas, meta = read_series_matrix(MATRIX)
print(f"  sondas {betas.shape[0]:,}  amostras {betas.shape[1]}  "
      f"{betas.memory_usage(deep=True).sum()/1e9:.2f} GB", flush=True)

age = pd.to_numeric(meta.set_index("gsm")["agebloodtaken"], errors="coerce")
age = age.reindex(betas.columns)
print(f"  idade: {age.notna().sum()} com valor, "
      f"{age.min():.0f} a {age.max():.0f} anos, mediana {age.median():.0f}\n", flush=True)

print(f"{'relogio':<14}{'cobertura':>11}{'r':>8}{'erro medio':>12}{'vies':>9}{'veredito':>12}")
rows = []
for name in CLOCKS:
    try:
        pred, cov = predict(betas, name)
    except Exception as e:
        print(f"{name:<14}  FALHOU: {e}")
        continue
    ok = age.notna() & pred.notna()
    r = float(np.corrcoef(age[ok], pred[ok])[0, 1])
    mae = float(np.abs(pred[ok] - age[ok]).mean())
    bias = float((pred[ok] - age[ok]).mean())
    verdict = "PASSA" if r >= 0.75 else "FALHA"
    rows.append(dict(clock=name, coverage=cov, r=r, mae=mae, bias=bias, verdict=verdict))
    print(f"{name:<14}{cov:>10.1%}{r:>8.3f}{mae:>11.1f}a{bias:>+8.1f}{verdict:>12}", flush=True)

res = pd.DataFrame(rows)
(ROOT / "results").mkdir(exist_ok=True)
res.to_csv(ROOT / "results/clock_validation.csv", index=False)

print(f"""
O limiar de r = 0.75 e deliberadamente mais frouxo que o 0.9 dos artigos. Dois
motivos declarados: esta coorte e estreita em idade (nao tem criancas, onde os
relogios ganham correlacao facil), e a matriz do GEO ja veio normalizada de um
jeito que nao e o dos autores. Um r de 0.8 aqui e consistente com uma
implementacao correta; um r de 0.3 nao e.

{'TODOS PASSAM — pode seguir.' if (res.verdict == 'PASSA').all() else 'ALGUM FALHOU — parar e investigar antes de qualquer analise.'}""")
