"""Supplementary Figure 14: 12 summaries x 59 strata, full (unmatched) group comparison, Test 73 vs 336 and Train 104 vs 673.
Same row/column order and colour scale as Fig. 2a (matched pairs); agreement with the matched comparison in panel c."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
from sc import *
from scipy.stats import pearsonr
d = CFG.unmap_features(pd.read_csv(AN / 'generators/results/gen_observed_cells_deltaZ_rebuilt.tsv', sep='\t'))
d['metric'] = d.feature.map(DISP)
full = d[d.contrast.isin(['Test_unmatched', 'Train_unmatched'])].copy(); full['split'] = full.contrast.str.split('_').str[0]
NH = {'Test': 336, 'Train': 673}
full['n_HCC'] = full.n.astype(int); full['n_healthy'] = full.split.map(NH)
wt(full[['split', 'stratum', 'metric', 'delta_z', 'p', 'q', 'n_HCC', 'n_healthy']], 'S14ab.tsv')
# agreement with the matched comparison (Test 73 pairs = Fig. 2a; Train 104 pairs = Supplementary Fig. 3a)
rows = []; sc_rows = []
for sp, um, mt in (('Test', 'Test_unmatched', 'Test_paired'), ('Train', 'Train_unmatched', 'Train_matched')):
    a = d[d.contrast == um].set_index(['stratum', 'feature']).delta_z; b = d[d.contrast == mt].set_index(['stratum', 'feature']).delta_z.loc[a.index]
    rows.append(dict(split=sp, n_cells=len(a), pearson_r=pearsonr(a, b)[0], same_sign_frac=float((np.sign(a) == np.sign(b)).mean()),
                     n_q05_full=int((d[d.contrast == um].q < .05).sum()), n_q05_matched=int((d[d.contrast == mt].q < .05).sum())))
    sc_rows.append(pd.DataFrame({'split': sp, 'stratum': a.index.get_level_values(0), 'metric': a.index.get_level_values(1).map(DISP),
                                 'delta_z_full': a.values, 'delta_z_matched': b.values}))
C = wt(pd.DataFrame(rows), 'S14c_agreement.tsv'); wt(pd.concat(sc_rows), 'S14c.tsv'); print(C)
D = rt('S14ab.tsv'); SC_ = rt('S14c.tsv'); C = C.set_index('split')
print('cells beyond |1|:', (D.delta_z.abs() > 1).sum())
cv = SC(178, 14); ov = cv.ov
for lab, sp, tp in (('a', 'Test', 2), ('b', 'Train', 64)):
    g = D[D.split == sp]; nh, nc = int(g.n_healthy.iloc[0]), int(g.n_HCC.max())
    cv.panel(lab, 1, tp - .5); cv.text(8, tp + .7, f'{sp}: HCC − healthy effect on 12 summaries across 59 strata (full comparison, {nc} vs {nh})', 6.3, 'bold')
    cm, norm, xs = heat(cv, g, 20, tp + 12, vmax=1.0)
cbar(cv, cm, norm, 22, 124.5, 26, 'Signed Δz (HCC − healthy)', [-1, 0, 1]); ov.scatter([62], [125.4], s=1.4, color=INK, lw=0)
ov.text(64, 125.4, 'BH q < 0.05', fontsize=5.5, va='center')
ov.text(82, 125.4, 'Same rows, columns and colour scale as Fig. 2a; values beyond ±1 shown at the end colour', fontsize=5.5, va='center', color=MUTED)
# c: full vs matched
cv.panel('c', 1, 132); cv.text(8, 133.2, 'Full comparison agrees with the matched comparison', 6.3, 'bold')
cat = {m[1]: m[2] for m in METRICS}
for k, (sp, ref) in enumerate((('Test', 'Fig. 2a, 73 pairs'), ('Train', 'Suppl. Fig. 3a, 104 pairs'))):
    ax = cv.ax(14 + k * 52, 141, 30, 28); g = SC_[SC_.split == sp]
    ax.scatter(g.delta_z_matched, g.delta_z_full, s=2.6, c=[CAT[cat[m]] for m in g.metric], alpha=.6, lw=0)
    ax.plot([-1.2, 1.2], [-1.2, 1.2], color=PALE, ls='--', lw=1.0, zorder=0); zeroline(ax); zeroline(ax, 'v')
    ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2); ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1])
    ax.set_xlabel(f'Matched Δz ({ref})', labelpad=1); ax.set_ylabel('Full-comparison Δz', labelpad=1)
    r = C.loc[sp]; ax_title(ax, f'{sp}: r = {r.pearson_r:.2f}, same sign {r.same_sign_frac * 100:.0f}%', pad=2)
hs = [Line2D([], [], ls='', marker='s', ms=3.4, mew=0, color=CAT[c], label=c) for c in dict.fromkeys(m[2] for m in METRICS)]
cv.fig.legend(handles=hs, loc='upper left', bbox_to_anchor=(118 / W, 1 - 141 / cv.H), fontsize=5.5, frameon=False, handletextpad=.3, labelspacing=.5)
cv.text(118, 167, 'One point per summary × stratum cell (708 cells).\nPearson r and sign agreement across cells.', 5.5, color=MUTED)
cv.save()
