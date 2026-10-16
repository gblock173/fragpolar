import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import numpy as np, pandas as pd
D = CFG.DATA
OUT = CFG.work('tumour_hepatitis_age')
RES, FIG = OUT + '/results', OUT + '/figures'
SEED = 20261014
METRICS = ['opp', 'strictM', 'strictU', 'mixed', 'B', 'U', 'mu', 'BV']

def load():
    """wide: index (sid), columns MultiIndex (stratum, metric); samples table."""
    L = pd.read_parquet(RES + '/sample_strata_long.parquet')
    W = L.pivot_table(index='sid', columns=['stratum', 'metric'], values='value', dropna=False)
    nm = L.pivot_table(index='sid', columns=['stratum', 'metric'], values='n_markers')
    # B/mu(1-mu): 층 평균 B / 층 평균 V (V=U+B)
    strata = W.columns.get_level_values(0).unique()
    for s in strata:
        W[(s, 'BV')] = W[(s, 'B')] / W[(s, 'V')]
    W = W.sort_index(axis=1)
    S = pd.read_csv(RES + '/sample_table.tsv', sep='\t').set_index('AnalysisID')
    return W, nm, S

def boot_idx(n, B, rng):
    return rng.integers(0, n, size=(B, n))

def ols(X, y):
    Xd = np.column_stack([np.ones(len(y)), X]); ok = np.isfinite(y) & np.isfinite(Xd).all(1)
    beta = np.linalg.lstsq(Xd[ok], y[ok], rcond=None)[0]
    return beta, int(ok.sum())
