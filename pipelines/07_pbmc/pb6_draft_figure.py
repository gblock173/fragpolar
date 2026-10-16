"""보충 그림 초안: a 군별 All-marker 값, b 층별 효과 산점도(PBMC-healthy vs HCC-healthy), c 연결 쌍."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
PB_DIR = CFG.work('pbmc')
L = pd.read_parquet(f'{PB_DIR}/results/sample_strata_long.parquet'); ms = pd.read_csv(f'{PB_DIR}/results/plasma_sample_table.tsv', sep='\t'); tk = pd.read_csv(f'{PB_DIR}/data/pbmc_tasks.tsv', sep='\t')
E = pd.read_csv(f'{PB_DIR}/results/effects.tsv', sep='\t'); CR = pd.read_csv(f'{PB_DIR}/results/effect_correlation.tsv', sep='\t'); P = pd.read_csv(f'{PB_DIR}/results/paired_diffs.tsv', sep='\t')
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})
col = {'Healthy': '#4C78A8', 'HCC': '#E45756', 'PBMC': '#54A24B'}
grp = {'Healthy (Train)': ('plasma', ms[(ms.group == 'healthy') & (ms.split == 'Train')].AnalysisID, col['Healthy']),
       'HCC (Train)': ('plasma', ms[(ms.group == 'liver') & (ms.split == 'Train')].AnalysisID, col['HCC']),
       'HCC (Test)': ('plasma', ms[(ms.group == 'liver') & (ms.split == 'Test')].AnalysisID, '#F2A19F'),
       'PBMC': ('PBMC', tk.pbmc_lab_id, col['PBMC'])}
fig = plt.figure(figsize=(7.2, 6.4)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.05], hspace=.45, wspace=.55)
rng = np.random.default_rng(1)
for k, (m, lab) in enumerate([('opp', 'Opposite strict-uniform\nfraction'), ('mixed', 'Strict-mixed fraction'), ('B', 'Between-fragment B')]):
    ax = fig.add_subplot(gs[0, k]); d = L[(L.metric == m) & (L.stratum == 'All|All')].set_index('sid').value
    for i, (g, (_, ids, c)) in enumerate(grp.items()):
        v = d.reindex(list(ids)).dropna().to_numpy(); ax.scatter(i + rng.uniform(-.2, .2, len(v)), v, s=3 if len(v) > 100 else 9, color=c, alpha=.5 if len(v) > 100 else .9, lw=0)
        ax.hlines(v.mean(), i - .32, i + .32, color='k', lw=1.2)
    ax.set_xticks(range(len(grp))); ax.set_xticklabels(list(grp), rotation=40, ha='right'); ax.set_ylabel(lab); ax.set_title('abcdef'[k] if False else '', fontsize=8)
    if k == 0: ax.text(-.5, 1.12, 'a', transform=ax.transAxes, fontsize=11, fontweight='bold')
ax = fig.add_subplot(gs[1, :2]); ee = E[E.metric == 'opp'].pivot(index='stratum', columns='contrast', values='effect')
x = ee['HCCTrain_minus_HTrain'] * 100; y = ee['PBMC_minus_HTrain'] * 100; ax.scatter(x, y, s=12, color='#555', alpha=.8)
ax.scatter(x.loc['All|All'], y.loc['All|All'], s=40, facecolor='none', edgecolor='k')
ax.axhline(0, color='#bbb', lw=.6); ax.axvline(0, color='#bbb', lw=.6)
c = CR[(CR.metric == 'opp') & (CR.hcc_contrast == 'HCCTrain_minus_HTrain')].iloc[0]
ax.text(.03, .97, f"59 strata\nPearson r = {c.pearson:.2f} [{c.pearson_lo:.2f}, {c.pearson_hi:.2f}]\nSpearman = {c.spearman:.2f} [{c.spearman_lo:.2f}, {c.spearman_hi:.2f}]\ncircle: All markers", transform=ax.transAxes, va='top', fontsize=7)
ax.set_xlabel('HCC (Train) - healthy (Train), opposite fraction (pp)'); ax.set_ylabel('PBMC - healthy plasma (Train),\nopposite fraction (pp)'); ax.text(-.12, 1.02, 'b', transform=ax.transAxes, fontsize=11, fontweight='bold')
ax = fig.add_subplot(gs[1, 2]); q = P[(P.metric == 'opp') & (P.stratum == 'All|All') & P.plasma_value.notna()]
for r in q.itertuples():
    tr = r.split == 'Train'; ax.plot([0, 1], [r.plasma_value, r.pbmc_value], '-o', color='#4C78A8' if tr else '#aaa', ms=3, lw=.8, mfc='#4C78A8' if tr else 'w')
ax.set_xticks([0, 1]); ax.set_xticklabels(['Same-person\nhealthy plasma', 'PBMC']); ax.set_xlim(-.3, 1.3); ax.set_ylabel('Opposite strict-uniform fraction')
ax.text(.02, .02, 'filled: Train (n=%d)\nopen: Test (n=%d)' % (q[q.split == 'Train'].pbmc.nunique(), q[q.split == 'Test'].pbmc.nunique()), transform=ax.transAxes, fontsize=6.5)
ax.text(-.35, 1.02, 'c', transform=ax.transAxes, fontsize=11, fontweight='bold')
for ext in ['png', 'pdf']: fig.savefig(f'{PB_DIR}/figures/FigS_PBMC_draft.{ext}', dpi=300, bbox_inches='tight')
print('fig ok')
