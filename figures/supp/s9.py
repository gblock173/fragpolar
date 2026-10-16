from sc import *
R = AN / 'pbmc/results/'
MN = {'opp': 'sOUF (%)', 'mixed': 'strict-mixed (%)', 'B': 'BFV (×10$^{-2}$)'}
g = pd.read_csv(R / 'group_means.tsv', sep='\t'); g = g[(g.stratum == 'All|All') & g.metric.isin(MN) & g.group.isin(['H_Train', 'HCC_Train', 'HCC_Test', 'PBMC'])]
wt(g[['metric', 'group', 'n_samples_valid', 'mean', 'lo', 'hi']], 'S9a.tsv')
e = pd.read_csv(R / 'effects.tsv', sep='\t'); e = e[e.metric.isin(MN) & e.contrast.isin(['PBMC_minus_HTrain', 'HCCTrain_minus_HTrain', 'HCCTest_minus_HTrain'])]
ep = e.pivot_table(index=['metric', 'stratum'], columns='contrast', values='effect').reset_index(); ep.columns.name = None; wt(ep, 'S9b.tsv')
cr = pd.read_csv(R / 'effect_correlation.tsv', sep='\t'); cr = cr[cr.metric.isin(MN) & cr.hcc_contrast.isin(['HCCTrain_minus_HTrain', 'HCCTest_minus_HTrain'])]
wt(cr[['metric', 'hcc_contrast', 'n_strata', 'pearson', 'pearson_lo', 'pearson_hi', 'spearman', 'n_same_sign']], 'S9b_corr.tsv')
p = pd.read_csv(R / 'paired_diffs.tsv', sep='\t'); p = p[(p.stratum == 'All|All') & p.metric.isin(MN) & p.split.isin(['Train', 'Test'])].dropna(subset=['diff'])
wt(p[['pbmc', 'split', 'metric', 'pbmc_value', 'plasma_value', 'diff']].assign(pbmc=lambda d: 'pair' + d.pbmc.astype('category').cat.codes.astype(str)), 'S9e.tsv')
A, B, C, P = rt('S9a.tsv'), rt('S9b.tsv'), rt('S9b_corr.tsv'), rt('S9e.tsv'); rng = np.random.default_rng(2)
sc_ = {'opp': 100, 'mixed': 100, 'B': 100}; GL = [('H_Train', 'Healthy\nTrain', HEALTHY, 'o'), ('HCC_Train', 'HCC\nTrain', HCC, '^'), ('HCC_Test', 'HCC\nTest', HCC, '^'), ('PBMC', 'PBMC', PBMC, 's')]
cv = SC(112, 9)
cv.panel('a', 1, 1); cv.title(8, 2.2, 'Group means, All markers (95% CI)')
for k, m in enumerate(MN):
    ax = cv.ax(12 + k * 29, 12, 24, 32)
    for i, (gk, lb, col, mk) in enumerate(GL):
        r = A[(A.metric == m) & (A.group == gk)].iloc[0]; f = sc_[m]; ax.plot([r.lo * f] * 2 if False else [i, i], [r.lo * f, r.hi * f], color=col, lw=1.0); ax.scatter([i], [r['mean'] * f], s=14, marker=mk, zorder=3, facecolor=col if gk != 'HCC_Test' else 'white', edgecolor=col, linewidths=1.0)
    ax.set_xticks(range(4)); ax.set_xticklabels([x[1] for x in GL], fontsize=5.5); ax.set_ylabel(MN[m], labelpad=1); ax.set_xlim(-.6, 3.6)
cv.panel('b', 98, 1); cv.title(105, 2.2, 'PBMC − plasma of the same person')
for k, m in enumerate(MN):
    ax = cv.ax(110 + k * 24, 12, 18, 32)
    for i, sp in enumerate(('Train', 'Test')):
        v = P[(P.metric == m) & (P.split == sp)]['diff'].values * 100; ax.scatter(i + rng.uniform(-.12, .12, len(v)), v, s=8, color=SPLIT_C[sp], lw=0, alpha=.85)
        ax.plot([i - .25, i + .25], [np.mean(v)] * 2, color=INK, lw=1.1)
    ax.axhline(0, color=PALE, lw=1.0, zorder=0); ax.set_xticks([0, 1]); ax.set_xticklabels([f'Train\n(n = {(P[(P.metric == "opp") & (P.split == "Train")]).shape[0]})', f'Test\n(n = {(P[(P.metric == "opp") & (P.split == "Test")]).shape[0]})'], fontsize=5.5); ax.set_xlim(-.6, 1.6)
    ax.set_ylabel('Δ ' + MN[m], labelpad=1, fontsize=5.5)
cv.panel('c', 1, 55); cv.title(8, 56.2, 'Stratum effects: HCC (Train) versus PBMC, both minus healthy Train')
for k, m in enumerate(MN):
    ax = cv.ax(14 + k * 58, 66, 36, 34); g_ = B[B.metric == m]; f = 100
    ax.scatter(g_.HCCTrain_minus_HTrain * f, g_.PBMC_minus_HTrain * f, s=6, color=INK, alpha=.6, lw=0); ax.axhline(0, color=GRID, lw=.8, zorder=0); ax.axvline(0, color=GRID, lw=.8, zorder=0)
    ax.set_xlabel('HCC − healthy: Δ ' + MN[m], labelpad=1, fontsize=5.5); ax.set_ylabel('PBMC − healthy' if k == 0 else '', labelpad=1, fontsize=5.5)
    r = C[(C.metric == m) & (C.hcc_contrast == 'HCCTrain_minus_HTrain')].iloc[0]
    ax.set_title(f'Pearson r = {r.pearson:.2f} ({r.pearson_lo:.2f} to {r.pearson_hi:.2f})', fontsize=5.5, pad=2, loc='left')
cv.save(); print(C)
