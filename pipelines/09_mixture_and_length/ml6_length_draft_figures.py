import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
E = CFG.work('mixture_and_length')
man = pd.read_csv(E + '/data/length_manifest.tsv', sep='\t'); L = pd.read_csv(E + '/results/len_length_metrics_per_sample.tsv', sep='\t')
lab = {'healthy': 'healthy', 'liver': 'HCC', 'pbmc': 'PBMC', 'hepatitis': 'hepatitis'}; col = {'healthy': '#3b6ea5', 'liver': '#c0392b', 'pbmc': '#7f8c8d', 'hepatitis': '#e08e0b'}
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
for g in ['healthy', 'liver', 'pbmc', 'hepatitis']:
    hs = [np.load(f'{E}/results/length_hist/{s}.npz')['h_all'] for s in man[man.group == g].sid]
    h = np.mean([x / x.sum() for x in hs], 0); ax[0].semilogy(np.arange(2000)[:501], h[:501], c=col[g], label=f'{lab[g]} (n={len(hs)})')
ax[0].set_xlabel('fragment length (bp)'); ax[0].set_ylabel('mean fraction'); ax[0].legend(); ax[0].set_title('length distribution (cpg_calls, <=500 bp)')
M = pd.read_csv(E + '/results/len_sample_master.tsv', sep='\t'); 
Wl = pd.read_parquet(CFG.work('tumour_hepatitis_age', 'results', 'sample_strata_long.parquet'))
y = Wl[(Wl.stratum == 'All|All') & (Wl.metric == 'mixed')].set_index('sid').value
for g, c in [('healthy', col['healthy']), ('HCC', col['liver'])]:
    d = M[M.grp == g]; ax[1].scatter(d.frac300 * 100, y.reindex(d.AnalysisID).values * 100, s=6, c=c, alpha=.5, label=g)
ax[1].set_xlabel('% fragments > 300 bp'); ax[1].set_ylabel('All|All strict-mixed (%)'); ax[1].legend(); ax[1].set_title('homogenisation vs long-fragment fraction')
R = pd.read_csv(E + '/results/len_adjust_retention_summary.tsv', sep='\t'); R = R[(R.set == 'all') & (R.model.isin(['len1', 'len3', 'alpha', 'full']))]
comps = ['strict-mixed', 'strictM', 'WFV', 'sOUF', 'BFV', 'FPI']; w = .2
for k, m in enumerate(['len1', 'len3', 'alpha', 'full']):
    q = R[R.model == m].set_index('component').loc[comps]; ax[2].bar(np.arange(6) + k * w, q.AllAll_retention, w, yerr=[q.AllAll_retention - q.ret_lo, q.ret_hi - q.AllAll_retention], label=m)
ax[2].axhline(1, c='k', lw=.7); ax[2].set_xticks(np.arange(6) + .3); ax[2].set_xticklabels(comps, rotation=20); ax[2].set_ylabel('HCC-healthy effect retained'); ax[2].legend(fontsize=7); ax[2].set_title('All|All, adjusted / unadjusted')
plt.tight_layout(); plt.savefig(E + '/figures/len_length_and_retention.png', dpi=150)
