"""Supplementary Figure 15: reduced within-fragment mixing in HCC (secondary observation).
a, 59-strata HCC - healthy differences of strict-mixed, healthy-direction uniform, fully methylated fractions and WFV (Train, Test, 95% CI);
b, all-marker composition change; c, retention in the subset with the same collection year and processing conditions."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
from sc import *
MET = [('mixed', 'strict-mixed', 100, 'Δ strict-mixed (%p)'), ('conc', 'healthy-direction uniform', 100, 'Δ healthy-direction\nuniform (%p)'),
       ('strictM', 'fully methylated', 100, 'Δ fully methylated (%p)'), ('U', 'WFV', 1000, 'Δ WFV (×10$^{-3}$)'), ('opp', 'sOUF', 100, 'Δ sOUF (%p)')]
STRICT = {'mixed', 'conc', 'strictM', 'opp'}


def eff(root, sp):
    e = pd.read_csv(AN / root / 'results/ge_strata_effects_all_metrics.tsv', sep='\t'); out = []
    for m, nm, f, _ in MET:
        lab = 'Train(strict보유)' if (sp == 'Train' and m in STRICT) else sp
        x = e[(e.split == lab) & (e.metric == m)].copy(); x['split'] = sp; x['metric_name'] = nm; out.append(x)
    x = pd.concat(out); assert x.stratum.nunique() == 59
    if sp == 'Train': assert (x.n_healthy == 673).all() and (x.n_HCC == 104).all()
    else: assert (x.n_healthy == 336).all() and (x.n_HCC == 73).all()
    return x


E = pd.concat([eff('train_rebuild', 'Train'), eff('group_effects', 'Test')])
E['stratum_label'] = E.stratum.map(stratum_label)
wt(E[['split', 'metric_name', 'stratum', 'stratum_label', 'n_HCC', 'n_healthy', 'mean_healthy', 'mean_HCC', 'delta', 'ci_lo', 'ci_hi', 'q_BH']], 'S15ab.tsv')
# c: subset with the same collection year and processing conditions (labels only; no site names)
S = pd.read_csv(AN / 'mixture_and_length/results/len_same_institution_AllAll.tsv', sep='\t')
S['comparison'] = np.select([S.comparison.str.startswith('FULL'), S.comparison.str.startswith('SAME')], ['full comparison', 'same collection year and processing'], 'other')
S = S[S.comparison != 'other'].copy(); S['component'] = S.component.replace({'strictM': 'fully methylated'})
R = S.pivot(index='component', columns='comparison', values='diff_mean'); R['retained_pct'] = R['same collection year and processing'] / R['full comparison'] * 100
Sc = S.merge(R[['retained_pct']].reset_index(), on='component')
wt(Sc[['comparison', 'component', 'n_A', 'n_B', 'diff_mean', 'ci_lo', 'ci_hi', 'retained_pct']].rename(columns={'n_A': 'n_HCC', 'n_B': 'n_healthy'}), 'S15c.tsv')
print(R.round(4))
E, Sc = rt('S15ab.tsv'), rt('S15c.tsv')

cv = SC(191, 15); ov = cv.ov
strata, xs, span = strata_x(E.stratum.unique(), bases=('All', 'U', 'mid', 'M'), gap=1.6)
cv.panel('a', 1, 1, 'HCC − healthy difference in 59 strata (mean, 95% CI)')
rows = MET[:4]; top = 13; hh = 21; X0, XW = 22, 138
for i, (m, nm, f, yl) in enumerate(rows):
    ax = cv.ax(X0, top + i * (hh + 3.2), XW, hh)
    for sp, off in (('Train', -.2), ('Test', .2)):
        g = E[(E.split == sp) & (E.metric_name == nm)].set_index('stratum')
        xx = np.array([xs[s] + off for s in strata]); v = g.loc[strata]
        ax.vlines(xx, v.ci_lo * f, v.ci_hi * f, color=SPLIT_C[sp], lw=.6, zorder=2)
        ax.scatter(xx, v.delta * f, s=4.5, marker=SPLIT_MK[sp], color=SPLIT_C[sp], lw=0, zorder=3)
        n_dn = int((v.ci_hi < 0).sum()); n_up = int((v.ci_lo > 0).sum())
        ax.text(1.005, .78 if sp == 'Train' else .55, f'{sp}: ↑{n_up} ↓{n_dn}', transform=ax.transAxes, fontsize=5.5, color=SPLIT_C[sp], va='center')
    zeroline(ax); ax.set_xlim(-1, max(xs.values()) + 1); ax.set_ylabel(yl, labelpad=1)
    for b, (a0, a1) in span.items():
        if a0 > 0: ax.axvline(a0 - .8, color=GRID, lw=.6, zorder=0)
    if i == 0:
        for b, (a0, a1) in span.items():
            ax.text((a0 + a1) / 2, 1.06, BASE_HEAD[b], transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=5.8, fontweight='bold')
    if i == len(rows) - 1:
        ax.set_xticks([xs[s] for s in strata]); ax.set_xticklabels([ann_label(s.split('|')[1]) for s in strata], rotation=90, fontsize=5.5)
    else:
        ax.set_xticks([xs[s] for s in strata]); ax.set_xticklabels([])
cv.text(X0 + XW + 1, top - 4.2, 'CI excl. 0', 5.5, color=MUTED)
# b: all-marker composition change
B0 = 129
cv.panel('b', 1, B0, 'All markers: composition change')
ax = cv.ax(14, B0 + 9, 60, 34); comps = [('strict-mixed', 100), ('healthy-direction uniform', 100), ('sOUF', 100)]
for j, (nm, f) in enumerate(comps):
    for sp, off in (('Train', -.19), ('Test', .19)):
        r = E[(E.split == sp) & (E.metric_name == nm) & (E.stratum == 'All|All')].iloc[0]
        ax.bar(j + off, r.delta * f, .34, color=tint(SPLIT_C[sp], .55), edgecolor=SPLIT_C[sp], lw=.8)
        ax.plot([j + off] * 2, [r.ci_lo * f, r.ci_hi * f], color=INK, lw=1.0)
        yv = r.delta * f; ax.text(j + off, (r.ci_hi * f + .25) if yv > 0 else (r.ci_lo * f - .25), f'{yv:+.2f}'.replace('-', '−'), ha='center',
                                    va='bottom' if yv > 0 else 'top', fontsize=5.5)
zeroline(ax); ax.set_xticks(range(3)); ax.set_xticklabels(['Strict-mixed', 'Healthy-direction\nuniform', 'sOUF']); ax.set_ylim(-4.2, 3.6)
ax.set_ylabel('HCC − healthy (%p)', labelpad=1); ax.set_xlim(-.6, 2.6)
ax.legend(handles=[Patch(facecolor=tint(SPLIT_C[k], .55), edgecolor=SPLIT_C[k], label=f'{k} ({n})') for k, n in (('Train', '104 vs 673'), ('Test', '73 vs 336'))],
          loc='upper left', fontsize=5.5, frameon=False, handlelength=1.0)
# c: retention in the subset with the same collection year and processing conditions
cv.panel('c', 86, B0, 'Same collection year and processing conditions')
ax = cv.ax(118, B0 + 9, 46, 34)
order = [('sOUF', 'sOUF', STATE['sOUF']), ('BFV', 'BFV', STATE['sOUF']), ('FPI', 'FPI', MUTED), ('strict-mixed', 'Strict-mixed', STATE['strict-mixed']),
         ('WFV', 'WFV', STATE['strict-mixed']), ('fully methylated', 'Fully methylated', STATE['strict-mixed'])]
Rr = Sc[Sc.comparison == 'same collection year and processing'].set_index('component')
for j, (k, lab, col) in enumerate(order):
    v = Rr.loc[k, 'retained_pct']; ax.barh(j, v, .62, color=tint(col, .6), edgecolor=col, lw=.8)
    ax.text(103, j, f'{v:.0f}%', va='center', fontsize=5.5)
ax.axvline(100, color=INK, lw=1.0, ls='--'); ax.set_yticks(range(len(order))); ax.set_yticklabels([o[1] for o in order]); ax.invert_yaxis()
ax.set_xlim(0, 118); ax.set_xticks([0, 50, 100]); ax.set_xlabel('Subset effect / full effect (%)', labelpad=1)
nH, nh = int(Rr.n_HCC.iloc[0]), int(Rr.n_healthy.iloc[0]); nHf = int(Sc[Sc.comparison == 'full comparison'].n_HCC.iloc[0]); nhf = int(Sc[Sc.comparison == 'full comparison'].n_healthy.iloc[0])
cv.note(86, B0 + 51, f'Mean HCC − healthy difference, all markers: subset {nH} vs {nh}\ndivided by the full comparison {nHf:,} vs {nhf:,} (Train + Test).')
cv.save()
