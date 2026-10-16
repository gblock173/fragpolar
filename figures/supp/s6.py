from sc import *
R = AN / 'generators/results/'
cs = pd.read_csv(R / 'gen_cosine_by_stratum_rebuilt.tsv', sep='\t'); cs = cs[cs.seedset == 'TrainSeed20']
CN = {'Train_matched': 'Train matched', 'Test_paired': 'Test matched', 'Train_unmatched': 'Train unmatched', 'Test_unmatched': 'Test unmatched'}
rows = []
for _, r in cs.iterrows():
    for g in ('G1', 'G2', 'G3'):
        for f in ('0.01', '0.02', '0.05', '0.1', '0.2'): rows.append((CN[r.contrast], r.stratum, g, float(f), r[f'{g}_{f}']))
wt(pd.DataFrame(rows, columns=['comparison', 'stratum', 'generator', 'f', 'cosine']), 'S6ac.tsv')
p = pd.read_csv(R / 'gen_polarization_summary_rebuilt.tsv', sep='\t'); p = p[(p.seedset == 'TrainSeed20') & (p.f_rule != 'fixed')].copy(); p['comparison'] = p.contrast.map(CN)
wt(p[['comparison', 'f_rule', 'n_strata', 'polarization', 'within_fragment']], 'S6b.tsv')
c0 = pd.read_csv(R / 'gen_cosine_by_stratum_rebuilt.tsv', sep='\t'); c0 = c0[c0.contrast == 'Test_paired'][['seedset', 'stratum', 'G1_0.05', 'G2_0.05', 'G3_0.05']]
c0['seedset'] = c0.seedset.map({'TrainSeed20': 'Train-derived base', 'TestSeed20_orig': 'Test-derived base'}); wt(c0, 'S6d.tsv')
A, B, D = rt('S6ac.tsv'), rt('S6b.tsv'), rt('S6d.tsv'); rng = np.random.default_rng(3)
GC = {'G1': '#C76B43', 'G2': '#4a5568', 'G3': '#2C7DA0'}; comps = list(CN.values()); LS = {'Train matched': '-', 'Test matched': '--', 'Train unmatched': ':', 'Test unmatched': '-.'}
cv = SC(112, 6)
cv.panel('a', 1, 1); cv.title(8, 2.2, 'Generator–observed cosine by fraction f')
ax = cv.ax(12, 11, 50, 42)
for g in GC:
    for cmp_ in comps[:2]:
        x = A[(A.generator == g) & (A.comparison == cmp_)].groupby('f').cosine.median(); ax.plot(x.index * 100, x.values, ls=LS[cmp_], color=GC[g], marker='o', ms=2.5)
ax.set_xscale('log'); ax.set_xticks([1, 2, 5, 10, 20]); ax.set_xticklabels(['1', '2', '5', '10', '20']); ax.axhline(0, color=PALE, lw=1.0); ax.set_ylim(-.8, .9)
ax.set_xlabel('Mixing fraction f (%)', labelpad=1); ax.set_ylabel('Median cosine over 59 strata', labelpad=1)
ax.legend(handles=[Line2D([], [], color=GC[g], label=g) for g in GC] + [Line2D([], [], color=MUTED, ls='-', label='Train matched'), Line2D([], [], color=MUTED, ls='--', label='Test matched')], frameon=False, fontsize=5.5, loc='center right', bbox_to_anchor=(1.02, .66), handlelength=1.6, labelspacing=.25)
cv.panel('b', 72, 1); cv.title(79, 2.2, 'Strata judged polarized (of 59)')
ax = cv.ax(82, 11, 40, 42); fs = ['0.01', '0.02', '0.05', '0.1', '0.2', 'best']
for cmp_, col in zip(comps, [SPLIT_C['Train'], SPLIT_C['Test'], '#a0aec0', '#f2c46d']):
    g = B[B.comparison == cmp_].set_index(B[B.comparison == cmp_].f_rule.astype(str)).reindex(fs)
    ax.plot(range(6), g.polarization.values, color=col, marker='o', ms=3, label=cmp_, ls='-' if 'matched' in cmp_ and 'un' not in cmp_ else '--')
ax.set_xticks(range(6)); ax.set_xticklabels(['1', '2', '5', '10', '20', 'best f'], fontsize=5.5); ax.set_xlabel('f (%)', labelpad=1); ax.set_ylabel('Polarized strata', labelpad=1); ax.set_ylim(0, 27)
ax.legend(frameon=False, fontsize=5.5, loc='upper left', bbox_to_anchor=(0, 1.02), ncol=2, handlelength=1.5, columnspacing=.6)
cv.panel('c', 131, 1); cv.title(138, 2.2, 'Cosine at f = 5%')
ax = cv.ax(141, 11, 38, 42); pos = []; data = []; cols = []
for i, cmp_ in enumerate(comps):
    for k, g in enumerate(GC):
        v = A[(A.generator == g) & (A.comparison == cmp_) & (A.f == .05)].cosine.values; data.append(v); pos.append(i + (k - 1) * .26); cols.append(GC[g])
box(ax, data, pos, cols, w=.22, rng=rng, s=1.8, alpha=.5); ax.axhline(0, color=PALE, lw=1.0, zorder=0)
ax.set_xticks(range(4)); ax.set_xticklabels(['Train\nmatched', 'Test\nmatched', 'Train\nunmatched', 'Test\nunmatched'], fontsize=5.5, rotation=90, ha='center'); ax.set_ylabel('Cosine (59 strata)', labelpad=1)
cv.panel('d', 1, 66); cv.title(8, 67.2, 'Base-sample independence (Test matched, f = 5%)')
for k, g in enumerate(GC):
    ax = cv.ax(14 + k * 36, 78, 28, 28); x = D[D.seedset == 'Train-derived base'].set_index('stratum')[f'{g}_0.05']; y = D[D.seedset == 'Test-derived base'].set_index('stratum')[f'{g}_0.05'].reindex(x.index)
    ax.scatter(x, y, s=5, color=GC[g], alpha=.7, lw=0); lo, hi = min(x.min(), y.min()) - .05, max(x.max(), y.max()) + .05; ax.plot([lo, hi], [lo, hi], color=PALE, ls='--', lw=1.0, zorder=0)
    ax.set_xlabel(f'{g} cosine, Train-derived base', labelpad=1, fontsize=5.5); ax.set_ylabel('Test-derived base' if k == 0 else '', labelpad=1, fontsize=5.5)
cv.text(126, 80, 'Base healthy samples: 20 Train or 20 Test.\nMean response vectors agree\n(stratum cosine ≥ 0.999 at every f).\nOne point per stratum.', 5.5, color=MUTED)
cv.save()
