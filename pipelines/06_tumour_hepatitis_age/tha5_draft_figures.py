import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from tha_common import *
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'font.family': 'NanumGothic', 'axes.unicode_minus': False})
W, nm, S = load()
# ---- tumour-fraction step ----
fs = pd.read_csv(RES + '/tf_f_samples.tsv', sep='\t').set_index('AnalysisID'); TF = pd.read_csv(RES + '/tf_dose_response.tsv', sep='\t')
hc = S[(S.group == 'healthy') & (S.split == 'Train')].index
panels = [('All|All', 'opp'), ('U|All', 'opp'), ('M|All', 'opp'), ('All|All', 'mixed'), ('All|All', 'B'), ('All|All', 'BV')]
fig, ax = plt.subplots(2, 3, figsize=(9, 5.6))
for a, (s, e) in zip(ax.ravel(), panels):
    med = W.loc[hc, (s, e)].median(); y = (W.loc[fs.index, (s, e)] - med).to_numpy(); f = fs.f.to_numpy()
    col = np.where(fs.split == 'Train', '#1f77b4', '#d95f02'); a.scatter(f, y, c=col, s=14)
    ok = np.isfinite(y); b = np.polyfit(f[ok], y[ok], 1); xx = np.linspace(0, f.max(), 50); a.plot(xx, b[0] * xx + b[1], 'k-', lw=1); a.axhline(0, color='gray', lw=.6, ls=':')
    r = TF[(TF.stratum == s) & (TF.endpoint == e) & (TF.subset == 'All46') & (TF.ref == 'median')].iloc[0]
    a.set_title(f'{s} {e}\nρ={r.rho:.2f} [{r.rho_lo:.2f},{r.rho_hi:.2f}] n={int(r.n)}', fontsize=8); a.set_xlabel('종양분율 f (CNV)'); a.set_ylabel('HCC - 건강인 Train 중앙값')
fig.suptitle('종양분율 용량-반응 (파랑 Train 28, 주황 Test 18)', fontsize=9); fig.tight_layout(); fig.savefig(FIG + '/tf_dose_response.png', dpi=200); plt.close(fig)
# ---- hepatitis step ----
M = pd.read_csv(RES + '/hep_contrasts.tsv', sep='\t'); E = M[(M.subset == 'full')]
key = ['All|All', 'U|All', 'mid|All', 'M|All', 'U|Island']; mets = ['opp', 'strictM', 'strictU', 'mixed', 'B', 'U']
fig, ax = plt.subplots(1, 3, figsize=(10, 3.6), gridspec_kw={'width_ratios': [2.2, 1, 1]})
a = ax[0]; w = .16; x0 = np.arange(len(mets)); s0 = 'All|All'
for k, (c, col) in enumerate([('HCC_minus_healthy', '#c0392b'), ('hepatitis_minus_healthy', '#2e86c1'), ('HCC_minus_hepatitis', '#7f8c8d')]):
    d = E[(E.stratum == s0) & (E.contrast == c)].set_index('metric').loc[mets]
    a.bar(x0 + (k - 1) * w * 1.1, d.diff_SD, w, color=col, label=c.replace('_', ' '))
    a.errorbar(x0 + (k - 1) * w * 1.1, d.diff_SD, yerr=[(d['diff'] - d.lo) / d.healthy_SD, (d.hi - d['diff']) / d.healthy_SD], fmt='none', ecolor='k', lw=.6)
a.set_xticks(x0); a.set_xticklabels(mets); a.axhline(0, color='k', lw=.6); a.set_ylabel('차이 / 건강인 SD'); a.set_title('All|All 세 대비 (full)'); a.legend(fontsize=6)
for a, m in zip(ax[1:], ['opp', 'B']):
    w_ = E[E.metric == m].pivot(index='stratum', columns='contrast', values='diff_SD').loc[[s for s in E.stratum.unique() if '|nonIsland' not in s]]
    a.scatter(w_.HCC_minus_healthy, w_.hepatitis_minus_healthy, s=9, c='k'); a.axhline(0, color='gray', lw=.5); a.axvline(0, color='gray', lw=.5)
    a.set_xlabel('HCC - 건강인 (SD)'); a.set_ylabel('간염 - 건강인 (SD)'); a.set_title(f'59층 {m}', fontsize=8)
fig.tight_layout(); fig.savefig(FIG + '/hep_three_groups.png', dpi=200); plt.close(fig)
# ---- ageing step ----
I = pd.read_csv(RES + '/age_island_vs_nonisland.tsv', sep='\t')
fig, ax = plt.subplots(1, 3, figsize=(11, 3.6))
d = I[(I.metric == 'mixed') & (I.comparator == 'U|nonIsland')]
subs = ['Train', 'Test']; ser = ['Island', 'comparator', 'Island_minus_comparator']; cols = ['#1b7837', '#762a83', 'k']
for j, sub in enumerate(['Train', 'Test']):
    for k, sr in enumerate(ser):
        r = d[(d.subset == sub) & (d.series == sr)].iloc[0]; ax[0].errorbar(k + j * 0.25 - .125, r.coef_per10y * 100, yerr=[[(r.coef_per10y - r.lo) * 100], [(r.hi - r.coef_per10y) * 100]], fmt='o' if j == 0 else 's', color=cols[k], capsize=2, ms=4)
ax[0].axhline(0, color='gray', lw=.6); ax[0].set_xticks(range(3)); ax[0].set_xticklabels(['U|Island', 'U|nonIsland', 'Island-nonIsland']); ax[0].set_ylabel('10년당 strict-mixed 기울기 (%p)'); ax[0].set_title('D2 건강인 (원 Train / 사각 Test)', fontsize=8)
for a, m in zip(ax[1:], ['mixed', 'opp']):
    J = pd.read_csv(RES + f'/age_direction_by_stratum_Test_{m}.tsv', sep='\t')
    band = J.stratum.str.split('|').str[0]; cmap = {'All': '#555', 'U': '#1f77b4', 'mid': '#999', 'M': '#d62728'}
    a.scatter(J.z_age, J.z_hcc, c=band.map(cmap), s=12); a.axhline(0, color='gray', lw=.5); a.axvline(0, color='gray', lw=.5)
    a.set_xlabel('노화 효과 (건강인 10년당, SD 단위)'); a.set_ylabel('HCC 효과 (SD 단위)'); a.set_title(f'59층 {m} (Test)', fontsize=8)
fig.tight_layout(); fig.savefig(FIG + '/age_age_island_and_direction.png', dpi=200)
