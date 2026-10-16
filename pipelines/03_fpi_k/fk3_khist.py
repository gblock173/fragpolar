"""FPI_K step: engine(earlier engine run/<sid>/meta.json)에 저장된 표본 단위 K 히스토그램(K>=4 fragment, raw)으로
표본별 E[1/K]와 1/mean_K를 비교해 Jensen 보정 계수 c = E[1/K] * mean(K) (>=1)를 구한다. 읽기 전용."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import json, os
import numpy as np, pandas as pd
A = CFG.WORK
s = pd.read_csv(A + '/group_effects/results/agg_samples.tsv', sep='\t')
k = np.arange(64.0); rows = []
for r in s.itertuples():
    f = f'{CFG.inp("prior_engine_dir")}/{r.AnalysisID}/meta.json'
    if not os.path.exists(f): continue
    m = json.load(open(f)); h = np.array(m['K_hist']['raw']); h[:4] = 0; N = h.sum()
    e1 = (h[4:] / k[4:]).sum() / N; mk = (h * k).sum() / N
    rows.append(dict(AnalysisID=r.AnalysisID, group=r.group, split=r.split, n_frag=N, mean_K=mk, inv_meanK=1 / mk, E_invK=e1, c=e1 * mk, frac_K_ge63=h[63] / N,
                     p50=np.searchsorted(np.cumsum(h) / N, .5), var_K=((h * k * k).sum() / N - mk ** 2)))
d = pd.DataFrame(rows); d.to_csv(A + '/fpi_k/results/fpik_khist_correction_by_sample.tsv', sep='\t', index=False)
g = d.groupby(['group', 'split']).agg(n=('c', 'size'), c_median=('c', 'median'), c_min=('c', 'min'), c_max=('c', 'max'), mean_K=('mean_K', 'mean'), E_invK=('E_invK', 'mean'), inv_meanK=('inv_meanK', 'mean'), frac_K_ge63=('frac_K_ge63', 'mean')).reset_index()
g.to_csv(A + '/fpi_k/results/fpik_khist_correction_summary.tsv', sep='\t', index=False); print(g.to_string())
