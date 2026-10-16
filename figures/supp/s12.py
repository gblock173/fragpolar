from sc import *
THA_DIR = AN / 'tumour_hepatitis_age/results/'
ML = {'opp': 'sOUF', 'mixed': 'strict-mixed', 'U': 'WFV', 'B': 'BFV', 'BV': 'FPI', 'strictM': 'fully methylated'}
d = pd.read_csv(THA_DIR / 'age_and_HCC_effects_all_strata.tsv', sep='\t'); d = d[d.metric.isin(ML) & d.stratum.isin(strata_order()) & d.subset.isin(['Train', 'Test'])].copy(); d['metric_name'] = d.metric.map(ML)
d['ci_excludes_0'] = ~((d.zlo <= 0) & (d.zhi >= 0)); d['q'] = np.where(d.ci_excludes_0, 0.0, 1.0)
wt(d[['kind', 'subset', 'stratum', 'metric_name', 'coef', 'z', 'zlo', 'zhi', 'ci_excludes_0', 'n']], 'S12abc.tsv')
i = pd.read_csv(THA_DIR / 'age_island_vs_nonisland.tsv', sep='\t'); i = i[(i.metric == 'mixed') & (i.series == 'Island_minus_comparator') & i.subset.isin(['Train', 'Train_strictAvail', 'Test'])].copy()
i['subset'] = i.subset.replace({'Train_strictAvail': 'Train'}); i = i.drop_duplicates(['subset', 'comparator'])  # M/mid comparators exist only as Train_strictAvail (= full Train after strict recovery)
wt(i.assign(comparator_label=i.comparator.map(lambda c: stratum_label(c, long=True)))[['subset', 'comparator', 'comparator_label', 'coef_per10y', 'lo', 'hi', 'n']], 'S12d.tsv')
D, I = rt('S12abc.tsv'), rt('S12d.tsv'); D['q'] = np.where(D.ci_excludes_0, 0.0, 1.0); D['metric'] = D.metric_name
print(D.groupby('kind').z.describe())
cv = SC(232, 12); mets = list(ML.values())
specs = [('a', 'age_per10y', 'Train', 'Train: age effect per 10 years in healthy (n = 673)', 14, .4), ('b', 'age_per10y', 'Test', 'Test: age effect per 10 years in healthy (n = 336)', 74, .4), ('c', 'HCC_vs_healthy_adj', 'Test', 'Test: HCC − healthy, adjusted for age, sex, depth (same scale)', 134, 3.0)]
for lt, kd, sp, ttl, y0, vm in specs:
    cv.panel(lt, 1, y0 - 12.5); cv.title(8, y0 - 11.3, ttl)
    g = D[(D.kind == kd) & (D.subset == sp)].copy(); cm, norm, xs = heat(cv, g, 20, y0, ch=3.0, vmax=vm, vcol='z', qcol='q', metrics=mets, catbars=False)
    cbar(cv, cm, norm, 22, y0 + 37, 24, 'Effect in healthy SD' + (' per 10 y' if kd.startswith('age') else ''), [-vm, 0, vm])
cv.ov.scatter([62], [14 + 38.5], s=1.4, color=INK, lw=0); cv.ov.text(64, 14 + 38.5, '95% CI excludes 0 (age panels: age slope; panel c: HCC coefficient)', fontsize=5.5, va='center')
cv.panel('d', 1, 184); cv.title(8, 185.2, 'Unmethylated-baseline CpG island minus comparator: strict-mixed slope per 10 years')
ax = cv.ax(40, 193, 40, 9); cmp_ = ['U|nonIsland', 'U|Shore', 'U|Enhancer', 'U|TSS', 'U|ReprPC', 'U|promoter', 'U|intergenic', 'M|nonIsland', 'mid|nonIsland']
ax.remove(); ax = cv.ax(40, 193, 44, 28)
for j, c in enumerate(cmp_):
    for sp, off in (('Train', -.17), ('Test', .17)):
        r = I[(I.subset == sp) & (I.comparator == c)]
        if len(r): r = r.iloc[0]; ax.plot([r.lo * 100, r.hi * 100], [j + off] * 2, color=SPLIT_C[sp], lw=1.0); ax.scatter([r.coef_per10y * 100], [j + off], s=8, color=SPLIT_C[sp], lw=0, zorder=3)
ax.axvline(0, color=PALE, ls='--', lw=1.0, zorder=0); ax.set_yticks(range(len(cmp_))); ax.set_yticklabels([stratum_label(c, long=True) for c in cmp_], fontsize=5.5); ax.invert_yaxis(); ax.set_xlabel('CpG island − comparator (%p per 10 y, 95% CI)', labelpad=1); ax.legend(handles=[Line2D([], [], color=SPLIT_C[k], marker='o', ms=3, label=k) for k in ('Train', 'Test')], frameon=False, fontsize=5.5, loc='lower right')
cv.save()
