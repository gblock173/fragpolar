"""Supplementary Figure 18: ageing in healthy plasma (y ~ age/10 + sex + log depth; model fixed in Train 673, refitted in Test 336).
a, age slope of sOUF and mean methylation in 59 strata grouped by healthy baseline; b, unmethylated-baseline CpG island vs non-CpG island
strict-mixed slopes; c, CpG island - comparator differences for unmethylated-baseline comparators."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
from sc import *
CDR = AN / 'tumour_hepatitis_age/results'
A = pd.read_csv(CDR / 'age_and_HCC_effects_all_strata.tsv', sep='\t')
A = A[(A.kind == 'age_per10y') & A.subset.isin(['Train', 'Test']) & A.metric.isin(['opp', 'mu', 'B']) & A.stratum.isin(strata_order())].copy()
A['metric_name'] = A.metric.map({'opp': 'sOUF', 'mu': 'mean methylation', 'B': 'BFV'}); A['stratum_label'] = A.stratum.map(stratum_label)
assert A.groupby(['subset', 'metric']).stratum.nunique().eq(59).all()
wt(A[['subset', 'metric_name', 'stratum', 'stratum_label', 'coef', 'lo', 'hi', 'z', 'zlo', 'zhi', 'n']].rename(columns={'coef': 'coef_per10y', 'lo': 'ci_lo', 'hi': 'ci_hi'}), 'S18a.tsv')
S = pd.read_csv(CDR / 'age_direction_summary.tsv', sep='\t'); S = S[S.subset.isin(['Train_strictAvail', 'Test']) & S.metric.isin(['opp', 'mu', 'B'])].copy()
S['subset'] = S.subset.replace({'Train_strictAvail': 'Train'}); S['metric_name'] = S.metric.map({'opp': 'sOUF', 'mu': 'mean methylation', 'B': 'BFV'})
wt(S[['subset', 'metric_name', 'strata_set', 'n_strata', 'age_up', 'age_down', 'age_up_CI', 'age_down_CI']], 'S18a_counts.tsv')
I = pd.read_csv(CDR / 'age_island_vs_nonisland.tsv', sep='\t'); I = I[(I.metric == 'mixed') & I.subset.isin(['Train', 'Test']) & I.comparator.str.startswith('U|')].copy()
SER = {'Island': 'CpG island', 'comparator': 'comparator', 'Island_minus_comparator': 'CpG island minus comparator'}
I = I.assign(series=I.series.map(SER), comparator=I.comparator.map(lambda s: stratum_label(s, long=True)), slope_pp_per10y=I.coef_per10y * 100, ci_lo_pp=I.lo * 100, ci_hi_pp=I.hi * 100)
wt(I[['subset', 'comparator', 'series', 'slope_pp_per10y', 'ci_lo_pp', 'ci_hi_pp', 'n']], 'S18bc.tsv')
A, S, I = rt('S18a.tsv'), rt('S18a_counts.tsv'), rt('S18bc.tsv')
for sp in ('Train', 'Test'):
    print(sp, S[(S.subset == sp) & (S.strata_set == 'baseline_M')][['metric_name', 'n_strata', 'age_up', 'age_down', 'age_up_CI', 'age_down_CI']].to_string())

cv = SC(160, 18); ov = cv.ov
strata, xs, span = strata_x(A.stratum.unique(), bases=('All', 'U', 'mid', 'M'), gap=1.6)
cv.panel('a', 1, 1, 'Age slope per 10 years in 59 strata (healthy; coefficient and 95% CI)')
X0, XW, top, hh = 22, 138, 15, 24
for i, (nm, yl, key, arrow) in enumerate((('sOUF', 'sOUF (%p per 10 y)', 'age_down', '↓'), ('mean methylation', 'Mean methylation\n(%p per 10 y)', 'age_up', '↑'))):
    ax = cv.ax(X0, top + i * (hh + 6), XW, hh)
    for sp, off in (('Train', -.2), ('Test', .2)):
        g = A[(A.subset == sp) & (A.metric_name == nm)].set_index('stratum').loc[strata]; xx = np.array([xs[s] + off for s in strata])
        ax.vlines(xx, g.ci_lo * 100, g.ci_hi * 100, color=SPLIT_C[sp], lw=.6, zorder=2); ax.scatter(xx, g.coef_per10y * 100, s=4.5, marker=SPLIT_MK[sp], color=SPLIT_C[sp], lw=0, zorder=3)
    zeroline(ax); ax.set_xlim(-1, max(xs.values()) + 1); ax.set_ylabel(yl, labelpad=1)
    for b, (a0, a1) in span.items():
        if a0 > 0: ax.axvline(a0 - .8, color=GRID, lw=.6, zorder=0)
        c = S[(S.metric_name == nm) & (S.strata_set == f'baseline_{b}')].set_index('subset')
        ax.text((a0 + a1) / 2, 1.02, f'{arrow} Train {c.loc["Train", key]}/{c.loc["Train", "n_strata"]}, Test {c.loc["Test", key]}/{c.loc["Test", "n_strata"]}',
                transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=5.5, color=MUTED)
        if i == 0: ax.text((a0 + a1) / 2, 1.13, BASE_HEAD[b], transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=5.8, fontweight='bold')
    ax.set_xticks([xs[s] for s in strata]); ax.set_xticklabels([ann_label(s.split('|')[1]) for s in strata] if i == 1 else [], rotation=90, fontsize=5.5)
    c = S[(S.metric_name == nm) & (S.strata_set == 'all59')].set_index('subset')
    ax.text(1.01, .5, f'{"Decrease" if arrow == "↓" else "Increase"},\n59 strata:\nTrain {c.loc["Train", key]}\nTest {c.loc["Test", key]}', transform=ax.transAxes, fontsize=5.5, va='center', color=INK2)
hs = [Line2D([], [], marker=SPLIT_MK[k], ls='-', color=SPLIT_C[k], ms=3, lw=.8, label=f'{k} (n = {n})') for k, n in (('Train', 673), ('Test', 336))]
cv.fig.legend(handles=hs, loc='upper right', bbox_to_anchor=(178 / W, 1 - 1 / cv.H), ncol=2, fontsize=5.5, frameon=False, handletextpad=.3)
# ---- b
B0 = 89
cv.panel('b', 1, B0, 'Unmethylated-baseline strict-mixed slope')
ax = cv.ax(36, B0 + 10, 26, 30); sers = [('CpG island', HEALTHY), ('comparator', PALE), ('CpG island minus comparator', INK)]
labs = ['CpG island', 'Non-CpG island', 'CpG island −\nnon-CpG island']
g0 = I[I.comparator == 'U|non-CpG island']
for j, (s_, col) in enumerate(sers):
    for sp, off in (('Train', -.16), ('Test', .16)):
        r = g0[(g0.subset == sp) & (g0.series == s_)].iloc[0]; ci_row(ax, r.slope_pp_per10y, r.ci_lo_pp, r.ci_hi_pp, j + off, col, SPLIT_MK[sp], ms=3)
        ax.text(1.03, j + off, f'{r.slope_pp_per10y:+.2f} [{r.ci_lo_pp:.2f}, {r.ci_hi_pp:.2f}]'.replace('-', '−'), transform=ax.get_yaxis_transform(), fontsize=5.5, va='center')
refline(ax, 0, 'v'); ax.set_yticks(range(3)); ax.set_yticklabels(labs); ax.set_ylim(2.6, -.6); ax.tick_params(axis='y', length=0)
ax.set_xlim(-.2, .3); ax.set_xlabel('%p per 10 years', labelpad=1)
cv.note(2, B0 + 47, 'Upper point, Train (circle); lower point, Test (square).\nSame samples; difference fitted within sample.')
# ---- c
cv.panel('c', 86, B0, 'CpG island versus unmethylated-baseline comparators (strict-mixed)')
comps = ['U|non-CpG island', 'U|CpG shore', 'U|Enhancer', 'U|TSS', 'U|Promoter', 'U|Intergenic', 'U|Polycomb-repressed']
for k, (s_, ttl, xl) in enumerate((('comparator', 'Comparator slope', '%p per 10 y'), ('CpG island minus comparator', 'CpG island − comparator', '%p per 10 y'))):
    ax = cv.ax(118 + k * 31, B0 + 12, 24, 40)
    for j, c in enumerate(comps):
        for sp, off in (('Train', -.17), ('Test', .17)):
            r = I[(I.comparator == c) & (I.series == s_) & (I.subset == sp)].iloc[0]
            ci_row(ax, r.slope_pp_per10y, r.ci_lo_pp, r.ci_hi_pp, j + off, SPLIT_C[sp], SPLIT_MK[sp], ms=3)
    if k == 0:
        for sp, ls in (('Train', '-'), ('Test', (0, (2, 1.5)))):
            r = I[(I.comparator == comps[0]) & (I.series == 'CpG island') & (I.subset == sp)].iloc[0]
            ax.axvline(r.slope_pp_per10y, color=HEALTHY, lw=.9, ls=ls, zorder=0)
    refline(ax, 0, 'v'); ax.set_yticks(range(len(comps))); ax.set_yticklabels([c.split('|')[1] for c in comps] if k == 0 else [])
    ax.set_ylim(len(comps) - .4, -.6); ax.tick_params(axis='y', length=0); ax_title(ax, ttl, pad=2); ax.set_xlabel(xl, labelpad=1)
    ax.axhspan(len(comps) - 1.5, len(comps) - .5, color=tint(TUMOUR, .12), zorder=-1, lw=0)
cv.note(86, B0 + 60, 'Blue lines, CpG island slope (solid Train, dashed Test). Shaded row: Polycomb-repressed\nregions, where the age-related increase in mixing is larger than in CpG islands.')
cv.save()
