import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
from sc import *
d = CFG.unmap_features(pd.read_csv(AN / 'generators/results/gen_observed_cells_deltaZ_rebuilt.tsv', sep='\t'))
d = d[d.contrast.isin(['Train_matched', 'Test_paired'])].copy(); d['metric'] = d.feature.map(DISP); d['split'] = d.contrast.map({'Train_matched': 'Train', 'Test_paired': 'Test'})
wt(d[['split', 'stratum', 'metric', 'delta_z', 'p', 'q', 'n']].rename(columns={'n': 'n_pairs'}), 'S3a.tsv')
D = rt('S3a.tsv'); print(D.groupby('split').n_pairs.max(), D.delta_z.abs().quantile(.99))
cv = SC(172, 3); ov = cv.ov
for lab, sp, y0, tp in (('a', 'Train', 14, 2), ('b', 'Test', 76, 64)):
    cv.panel(lab, 1, tp - .5); n = int(D[D.split == sp].n_pairs.max())
    cv.text(8, tp + .7, f'{sp}: HCC − healthy effect on 12 summaries across 59 strata ({n} matched pairs)', 6.3, 'bold')
    g = D[D.split == sp]; cm, norm, xs = heat(cv, g.rename(columns={'delta_z': 'delta_z'}), 20, y0 + (0 if lab == 'a' else 0), vmax=1.0)
cbar(cv, cm, norm, 22, 123, 26, 'Signed Δz (HCC − healthy)', [-1, 0, 1]); ov.scatter([62], [124], s=1.4, color=INK, lw=0); ov.text(64, 124, 'BH q < 0.05', fontsize=5.5, va='center')
# c: Train vs Test agreement
cv.panel('c', 1, 131); cv.text(8, 132.2, 'Train and Test agree for every summary', 6.3, 'bold')
ax = cv.ax(12, 138, 46, 28); T = D[D.split == 'Train'].set_index(['stratum', 'metric']); E_ = D[D.split == 'Test'].set_index(['stratum', 'metric'])
cat = {m[1]: m[2] for m in METRICS}; rows = []
for m in [x[1] for x in METRICS]:
    st = [s for s in T.index.get_level_values(0).unique() if (s, m) in T.index and (s, m) in E_.index]
    a = T.loc[[(s, m) for s in st], 'delta_z'].values; b = E_.loc[[(s, m) for s in st], 'delta_z'].values
    ax.scatter(a, b, s=3, color=CAT[cat[m]], alpha=.6, lw=0)
    from scipy.stats import spearmanr
    qa = T.loc[[(s, m) for s in st], 'q'].values; qb = E_.loc[[(s, m) for s in st], 'q'].values
    rows.append((m, cat[m], len(st), spearmanr(a, b)[0], (np.sign(a) == np.sign(b)).sum(), (qa < .05).sum(), (qb < .05).sum()))
S = pd.DataFrame(rows, columns=['metric', 'category', 'n_strata', 'spearman_train_test', 'n_same_sign', 'n_q05_train', 'n_q05_test']); wt(S, 'S3c.tsv')
ax.plot([-1, 1], [-1, 1], color=PALE, ls='--', lw=1.0, zorder=0); ax.axhline(0, color=GRID, lw=.8, zorder=0); ax.axvline(0, color=GRID, lw=.8, zorder=0)
ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_xlabel('Train Δz', labelpad=1); ax.set_ylabel('Test Δz', labelpad=1)
ax2 = cv.ax(78, 138, 42, 28); ax2.barh(range(12), S.spearman_train_test, color=[CAT[c] for c in S.category], height=.7); ax2.set_yticks(range(12)); ax2.set_yticklabels(S.metric); ax2.invert_yaxis()
ax2.set_xlabel('Spearman ρ, Train vs Test (59 strata)', labelpad=1); ax2.set_xlim(-.2, 1); ax2.axvline(0, color=PALE, lw=1.0)
ax3 = cv.ax(138, 138, 38, 28); yy = np.arange(12)
ax3.barh(yy - .19, S.n_q05_train, height=.36, color=SPLIT_C['Train']); ax3.barh(yy + .19, S.n_q05_test, height=.36, color=SPLIT_C['Test'])
ax3.set_yticks(yy); ax3.set_yticklabels(S.metric); ax3.invert_yaxis(); ax3.set_xlim(0, 59); ax3.set_xlabel('Strata with BH q < 0.05 (of 59)', labelpad=1)
ax3.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=SPLIT_C[k]) for k in ('Train', 'Test')], labels=['Train', 'Test'], frameon=False, fontsize=5.5, loc="upper left", bbox_to_anchor=(0, 1.22), ncol=2, handlelength=.8)
cv.save(); print(S)
