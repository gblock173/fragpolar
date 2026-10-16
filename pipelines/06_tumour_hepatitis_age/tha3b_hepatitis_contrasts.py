"""hepatitis step (included hepatitis samples only): three-group contrasts HCC-H, hep-H, HCC-hep for the 59 strata + 4 non-island strata.

RECONSTRUCTED SCRIPT (no original code was kept for this table). It rebuilds results/D1_contrasts_hep_included.tsv
(columns stratum, metric, contrast, d, lo, hi), which is read by figures/main/fig456_prep.py (Fig. 5c) and
figures/supp/s10.py (Supplementary Fig. 10), using the tha3_hepatitis.py conventions: group mean difference of the stratum values
(healthy = all healthy, HCC = all HCC, hepatitis = sample_lists.hepatitis_included), 95% CI from group-wise multinomial
bootstrap (2,000 replicates, seed 20261014). See docs/CODE_STRUCTURE_KO.md for the agreement with the stored table.
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import numpy as np, pandas as pd
from tha_common import *
W, nm, S = load(); NB = 2000
rng = np.random.default_rng(SEED)
MET = ['B', 'U', 'mixed', 'mu', 'opp', 'strictM']
HEP = set(CFG.ids('hepatitis_included'))
g = S.group
grp = {'H': [x for x in S.index[g == 'healthy'] if x in W.index], 'HCC': [x for x in S.index[g == 'liver'] if x in W.index],
       'hep': [x for x in S.index[g == 'hepatitis'] if x in HEP and x in W.index]}
print({k: len(v) for k, v in grp.items()})
cols = [c for c in W.columns if c[1] in MET]
X = W[cols]


def bmeans(Xg):
    v = Xg.to_numpy(float); val = np.isfinite(v); f = np.where(val, v, 0.0); n = len(v)
    cnt = rng.multinomial(n, np.full(n, 1 / n), size=NB).astype(float)
    return f.sum(0) / np.maximum(val.sum(0), 1), (cnt @ f) / np.maximum(cnt @ val.astype(float), 1)


est = {k: bmeans(X.loc[v]) for k, v in grp.items()}
rows = []
for cn, (a, b) in {'HCC-H': ('HCC', 'H'), 'hep-H': ('hep', 'H'), 'HCC-hep': ('HCC', 'hep')}.items():
    d = est[a][0] - est[b][0]; lo, hi = np.nanpercentile(est[a][1] - est[b][1], [2.5, 97.5], axis=0)
    for k, (s, m) in enumerate(cols):
        rows.append(dict(stratum=s, metric=m, contrast=cn, d=d[k], lo=lo[k], hi=hi[k]))
T = pd.DataFrame(rows)
T['_c'] = T.contrast.map({'HCC-H': 0, 'hep-H': 1, 'HCC-hep': 2})
T = T.sort_values(['stratum', 'metric', '_c']).drop(columns='_c')
T.to_csv(RES + '/hep_contrasts_included.tsv', sep='\t', index=False)
print(T[T.stratum == 'All|All'].to_string())
