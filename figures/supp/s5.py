from sc import *
rows = []; per = []
for root, sp_in, sp in (('train_rebuild', 'Train(strict보유)', 'Train'), ('group_effects', 'Test', 'Test')):
    m = pd.read_csv(AN / root / 'results/ge_margin_summary.tsv', sep='\t'); m = m[m.split == sp_in].copy(); m['split'] = sp; rows.append(m)
    r = pd.read_csv(AN / root / 'results/ge_margin_reduction_by_stratum.tsv', sep='\t'); r = r[(r.split == sp_in) & (r.baseline == 'full')].copy(); r['split'] = sp; per.append(r)
M = pd.concat(rows)[['baseline', 'margin', 'split', 'n_strata_eval', 'n_dOpp_pos', 'n_dOpp_q05_pos', 'n_dB_pos', 'n_dB_q05_pos', 'median_reduction_pct', 'n_reduction_ge50', 'n_red_evaluable']]
wt(M, 'S5ab.tsv'); P = pd.concat(per)[['split', 'stratum', 'margin', 'opp', 'B', 'Bdel', 'Bsub', 'reduction']]; P['reduction_pct'] = P.reduction * 100; wt(P.drop(columns='reduction'), 'S5c.tsv')
M, P = rt('S5ab.tsv'), rt('S5c.tsv'); full = M[M.baseline == 'full']; rng = np.random.default_rng(1)
cv = SC(72, 5)
cv.panel('a', 1, 1); cv.title(8, 2.2, 'Reduction of ΔBFV on removing sOUF')
ax = cv.ax(12, 11, 36, 42)
for sp in ('Train', 'Test'):
    g = full[full.split == sp]; ax.plot(g.margin, g.median_reduction_pct, marker='o', ms=3, color=SPLIT_C[sp], label=sp)
    g2 = M[(M.baseline == 'trH') & (M.split == sp)]; ax.plot(g2.margin, g2.median_reduction_pct, ls=':', color=SPLIT_C[sp], lw=1.0)
ax.set_xlabel('Excluded margin around 0.5', labelpad=1); ax.set_ylabel('Median reduction (%)', labelpad=1); ax.set_ylim(40, 90); ax.set_xticks([0, .1, .2, .3]); ax.legend(frameon=False, fontsize=5.5, loc='upper right', handletextpad=.3)
cv.text(12, 59, 'Solid, baseline from all Train healthy;\ndotted, Train healthy only.', 5.5, color=MUTED)
cv.panel('b', 62, 1); cv.title(69, 2.2, 'Strata with increase (of 59)')
ax = cv.ax(72, 11, 36, 42)
for sp in ('Train', 'Test'):
    g = full[full.split == sp]
    ax.plot(g.margin, g.n_dOpp_q05_pos, marker='o', ms=3, color=SPLIT_C[sp]); ax.plot(g.margin, g.n_dB_q05_pos, marker='s', ms=3, color=SPLIT_C[sp], ls='--')
ax.set_ylim(30, 60); ax.set_xticks([0, .1, .2, .3]); ax.set_xlabel('Excluded margin around 0.5', labelpad=1); ax.set_ylabel('Strata, BH q < 0.05 and Δ > 0', labelpad=1)
ax.legend(handles=[Line2D([], [], color=MUTED, marker='o', ms=3, label='ΔsOUF'), Line2D([], [], color=MUTED, marker='s', ms=3, ls='--', label='ΔBFV')], frameon=False, fontsize=5.5, loc='lower left', handletextpad=.3)
cv.panel('c', 122, 1); cv.title(129, 2.2, 'Per-stratum reduction')
ax = cv.ax(132, 11, 44, 42); ms = sorted(P.margin.unique()); data = []; pos = []; cols = []
for i, mg in enumerate(ms):
    for k, sp in enumerate(('Train', 'Test')):
        v = P[(P.split == sp) & (P.margin == mg)].reduction_pct.clip(-50, 150).dropna().values; data.append(v); pos.append(i + (k - .5) * .38); cols.append(SPLIT_C[sp])
box(ax, data, pos, cols, w=.3, rng=rng, s=2.5, alpha=.45); ax.axhline(100, color=PALE, ls='--', lw=1.0, zorder=0)
ax.set_xticks(range(len(ms))); ax.set_xticklabels([f'{m:g}' for m in ms]); ax.set_xlabel('Excluded margin', labelpad=1); ax.set_ylabel('Reduction (%, clipped −50 to 150)', labelpad=1)
cv.save()
