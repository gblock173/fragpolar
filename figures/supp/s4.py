from sc import *
R = AN / 'fpi_k/results/'
c = pd.read_csv(R / 'fpik_FPI_vs_FPIK_by_stratum.tsv', sep='\t')[['split', 'stratum', 'delta_FPI', 'q_BH_FPI', 'delta_FPIK', 'q_BH_FPIK', 'h_healthy', 'delta_h', 'ratio']]
wt(c, 'S4a.tsv')
s = pd.read_csv(R / 'fpik_sample_AllAll_values.tsv', sep='\t'); s = s[s.group.isin(['healthy', 'liver'])].copy(); s['group'] = s.group.map({'healthy': 'Healthy', 'liver': 'HCC'})
wt(s[['AnalysisID', 'group', 'split', 'h_floor', 'FPI', 'FPI_K']].rename(columns={'AnalysisID': 'sample', 'h_floor': 'h'}), 'S4b.tsv')
e = pd.read_csv(R / 'fpik_effect_summary.tsv', sep='\t'); lab = {'FPI': 'FPI', 'FPI_K': 'FPI_K', 'FPI_K_cglobal': 'FPI_K, global c', 'FPI_K_csample': 'FPI_K, sample c (upper bound)', 'FPI_K_markermean': 'FPI_K, marker mean', 'FPI_markermean': 'FPI, marker mean'}
e = e[e.metric.isin(lab)].copy(); e['variant'] = e.metric.map(lab)
wt(e[['split', 'variant', 'n_strata', 'n_delta_pos', 'n_q05_pos', 'allall_delta', 'allall_ci_lo', 'allall_ci_hi']], 'S4c.tsv')
j = pd.read_csv(R / 'fpik_raw_jensen_by_stratum_mean.tsv', sep='\t'); wt(j, 'S4d.tsv')
a, b, cc, dd = rt('S4a.tsv'), rt('S4b.tsv'), rt('S4c.tsv'), rt('S4d.tsv')
cv = SC(116, 4)
cv.panel('a', 1, 1); cv.title(8, 2.2, 'FPI and FPI_K effects per stratum')
ax = cv.ax(12, 10, 36, 36)
for sp, mk in (('Train', 'o'), ('Test', 's')):
    g = a[a.split == sp]; ax.scatter(g.delta_FPI, g.delta_FPIK, s=6, marker=mk, color=SPLIT_C[sp], alpha=.75, lw=0, label=sp)
ax.plot([-.01, .09], [-.01, .09], color=PALE, ls='--', lw=1.0, zorder=0); ax.plot([-.01, .09], [-.012, .109], color=HCC, lw=.8, zorder=0)
ax.axhline(0, color=GRID, lw=.8, zorder=0); ax.axvline(0, color=GRID, lw=.8, zorder=0); ax.set_xlabel('Δ FPI', labelpad=1); ax.set_ylabel('Δ FPI_K', labelpad=1)
ax.legend(frameon=False, fontsize=5.5, loc='upper left', handletextpad=.1, borderaxespad=0)
med = a.groupby('split').ratio.median(); cv.text(12, 57, f'Dashed, identity; line, slope 1/(1 − h) = 1.22.\nMedian Δ FPI_K / Δ FPI = {med["Train"]:.2f} (Train), {med["Test"]:.2f} (Test)', 5.5, color=MUTED)
cv.panel('b', 62, 1); cv.title(69, 2.2, 'K floor h does not explain FPI')
ax = cv.ax(72, 10, 36, 36)
for g, mk, col in (('Healthy', 'o', HEALTHY), ('HCC', '^', HCC)):
    x = b[b.group == g]; ax.scatter(x.h, x.FPI, s=3 if g == 'Healthy' else 5, marker=mk, color=col, alpha=.4, lw=0, label=g)
ax.set_xlabel('h = 1 / mean K (sample)', labelpad=1); ax.set_ylabel('FPI (All markers)', labelpad=1); ax.legend(frameon=False, fontsize=5.5, loc='upper left', handletextpad=.1, borderaxespad=0, markerscale=1.5)
cv.text(62, 57, 'Healthy FPI vs h: Spearman −0.26 (Train, n = 673).\nGroup difference in h ≈ 0.', 5.5, color=MUTED)
cv.panel('c', 122, 1); cv.title(129, 2.2, 'All-marker effect, variants')
ax = cv.ax(148, 10, 28, 36)
vs = ['FPI', 'FPI_K', 'FPI_K, global c', 'FPI_K, sample c (upper bound)']; ypos = {v: i for i, v in enumerate(vs)}
for sp, off, col in (('Train', -.14, SPLIT_C['Train']), ('Test', .14, SPLIT_C['Test'])):
    g = cc[(cc.split == sp) & cc.variant.isin(vs)]
    for _, r in g.iterrows():
        y = ypos[r.variant] + off; ax.plot([r.allall_ci_lo, r.allall_ci_hi], [y, y], color=col, lw=1.0); ax.scatter([r.allall_delta], [y], s=9, color=col, lw=0, zorder=3)
ax.set_yticks(range(4)); ax.set_yticklabels(['FPI', 'FPI_K', 'FPI_K\nglobal c', 'FPI_K\nsample c'], fontsize=5.5); ax.invert_yaxis(); ax.set_xlim(.02, .05); ax.set_xticks([.02, .035, .05]); ax.set_xlabel('Δ (HCC − healthy)', labelpad=1)
cv.text(122, 57, 'Bars, 95% CI; Train, Test.', 5.5, color=MUTED)
cv.panel("d", 1, 68); cv.title(8, 69.2, 'Approximation error of h = 1/mean K')
ax = cv.ax(12, 80, 80, 28); dd = dd.sort_values('diff_mean').reset_index(drop=True)
ax.errorbar(range(len(dd)), dd.diff_mean, yerr=dd.diff_sd_between_samples, fmt='o', ms=2, color=INK, ecolor=PALE, elinewidth=.8, lw=0)
ax.axhline(0, color=PALE, lw=1.0, ls='--'); ax.set_xlabel('59 strata (sorted)', labelpad=1); ax.set_ylabel('FPI_K exact − FPI_K (1/mean K)', labelpad=1); ax.set_xticks([])
aa = dd[dd.stratum == 'All|All'].iloc[0]
cv.text(100, 82, f'Exact E[1/K] lowers FPI_K by {-dd.diff_mean.median():.3f} (median over strata;\nAll markers {-aa.diff_mean:.3f}; 3 healthy samples).\nStratum rank Spearman 0.996.\nBars, SD across the 3 samples.', 5.5, color=MUTED)
cv.save()
