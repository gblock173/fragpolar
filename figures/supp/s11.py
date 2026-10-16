import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
from sc import *
e = pd.read_csv(CFG.inp('prior_cnv_effect_summary'), sep='\t')
SCN = {'full_primary': ('All markers (default)', 'CNV'), 'panel_variable': ('10,000-marker panel', 'CNV'), 'fixed_depth': ('Strict: 4 CpG, 10 fragments', 'CNV'),
       'fixed_balanced': ('Strict, balanced strata', 'CNV'), 'panel_balanced': ('Panel, balanced strata', 'CNV'), 'legacy_panel': ('Panel, earlier matched controls', 'CNV'),
       'negative_panel': ('Panel', 'No CNV'), 'negative_fixed': ('4 CpG, 10 fragments', 'No CNV'), 'negative_strict_fixed': ('Strict, balanced strata', 'No CNV')}
x = e[(e.level == 'segment') & (e.stratum == 'All') & e.scenario.isin(SCN)].copy()
x['condition'] = x.scenario.map(lambda s: SCN[s][0]); x['patients'] = x.scenario.map(lambda s: SCN[s][1])
MM_ = {'strict_opposite': 'sOUF', 'strict_mixed': 'strict-mixed', 'B_predelete': 'BFV before removal', 'B_deleted': 'BFV after sOUF removal', 'B_random': 'BFV after random removal'}
x = x[x.metric.isin(MM_)]; x['state'] = x.metric.map(MM_)
out = x[['patients', 'condition', 'state', 'contrast', 'n_patients', 'n_healthy', 'mean', 'lo', 'hi']].rename(columns={'n_healthy': 'n_healthy_matched'})
wt(out, 'S11.tsv'); X = rt('S11.tsv')
cv = SC(122, 11); order = [v[0] + '|' + v[1] for v in SCN.values()]
def labels(df):
    d = df.drop_duplicates(['patients', 'condition']); return [(r.patients, r.condition, int(r.n_patients)) for r in d.itertuples()]
def pan(ax, state, contrast, factor=100):
    g = X[(X.state == state) & (X.contrast == contrast)]; labs = []
    for i, (p, c) in enumerate([(v[1], v[0]) for v in SCN.values() if v[1] and (X[(X.state == state) & (X.contrast == contrast) & (X.condition == v[0]) & (X.patients == v[1])].shape[0] > 0)]):
        r = g[(g.condition == c) & (g.patients == p)].iloc[0]; col = HCC if p == 'CNV' else MUTED
        ax.plot([r.lo * factor, r.hi * factor], [i, i], color=col, lw=1.0); ax.scatter([r['mean'] * factor], [i], s=12, color=col, marker='^' if p == 'CNV' else 'o', lw=0, zorder=3); labs.append(f'{c} (n = {int(r.n_patients)})')
    ax.set_yticks(range(len(labs))); ax.set_yticklabels(labs, fontsize=5.5); ax.invert_yaxis(); ax.axvline(0, color=PALE, ls='--', lw=1.0, zorder=0)
cv.panel('a', 1, 1); cv.title(8, 2.2, 'sOUF, gain − neutral')
ax = cv.ax(46, 11, 42, 44); pan(ax, 'sOUF', 'gain-neutral'); ax.set_xlabel('%p (95% CI)', labelpad=1)
cv.panel('b', 94, 1); cv.title(101, 2.2, 'strict-mixed, gain − neutral')
ax = cv.ax(104, 11, 42, 44); pan(ax, 'strict-mixed', 'gain-neutral'); ax.set_xlabel('%p (95% CI)', labelpad=1); ax.set_yticklabels([])
cv.text(152, 20, 'Triangles, patients with\nCNV gain (n = 30; 14 for\nearlier controls).\nCircles, patients without\nCNV detected (n = 33).\nSegment-level, All markers.', 5.5, color=MUTED)
cv.panel('c', 1, 66); cv.title(8, 67.2, 'sOUF in gain and neutral segments separately')
ax = cv.ax(26, 78, 30, 36); cs = [('All markers (default)', 'CNV'), ('10,000-marker panel', 'CNV'), ('Strict: 4 CpG, 10 fragments', 'CNV')]
for i, (c, p) in enumerate(cs):
    for ct, off, col, mk in (('gain', -.15, HCC, '^'), ('neutral', .15, MUTED, 'o')):
        r = X[(X.state == 'sOUF') & (X.contrast == ct) & (X.condition == c)]
        if not len(r): continue
        r = r.iloc[0]; ax.plot([r.lo * 100, r.hi * 100], [i + off] * 2, color=col, lw=1.0); ax.scatter([r['mean'] * 100], [i + off], s=11, color=col, marker=mk, lw=0, zorder=3)
ax.set_xlim(1, 5.6); ax.set_yticks(range(3)); ax.set_yticklabels(['All markers\n(default)', '10,000-marker\npanel', 'Strict:\n4 CpG, 10 frag.'], fontsize=5.5); ax.invert_yaxis(); ax.set_xlabel('Healthy-referenced sOUF (%p)', labelpad=1)
ax.legend(handles=[Line2D([], [], color=HCC, marker='^', ms=3, label='Gain'), Line2D([], [], color=MUTED, marker='o', ms=3, label='Neutral')], frameon=False, fontsize=5.5, loc='center right', handlelength=1.2)
cv.panel('d', 92, 66); cv.title(99, 67.2, 'BFV before and after removing sOUF fragments (gain − neutral)')
ax = cv.ax(108, 78, 66, 36); st = ['BFV before removal', 'BFV after sOUF removal', 'BFV after random removal']; cols = ['#a0aec0', HCC, '#4a5568']
cn = [('10,000-marker panel', 'CNV'), ('Panel, balanced strata', 'CNV'), ('Panel, earlier matched controls', 'CNV'), ('Panel', 'No CNV')]
for i, (c, p) in enumerate(cn):
    for j, s_ in enumerate(st):
        r = X[(X.state == s_) & (X.contrast == 'gain-neutral') & (X.condition == c) & (X.patients == p)]
        if not len(r): continue
        r = r.iloc[0]; y = i + (j - 1) * .24; ax.plot([r.lo * 1e3, r.hi * 1e3], [y, y], color=cols[j], lw=1.0); ax.scatter([r['mean'] * 1e3], [y], s=9, color=cols[j], lw=0, zorder=3)
ax.set_xlim(-1, 7.5); ax.set_yticks(range(4)); ax.set_yticklabels(['Panel,\nvariable K (CNV)', 'Panel, balanced\nstrata (CNV)', 'Panel, earlier\ncontrols (CNV)', 'Panel\n(no CNV)'], fontsize=5.5); ax.invert_yaxis(); ax.axvline(0, color=PALE, ls='--', lw=1.0, zorder=0); ax.set_xlabel('Δ BFV, gain − neutral (×10$^{-3}$)', labelpad=1)
ax.legend(handles=[Line2D([], [], color=c, marker='o', ms=2.8, label=s_) for c, s_ in zip(cols, st)], frameon=False, fontsize=5.5, loc='lower right', handlelength=1.2, labelspacing=.25)
cv.save()
