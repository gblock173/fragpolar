import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
from sc import *
R = AN / 'mixture_and_length/results/'
m = pd.read_csv(R / 'len_length_metrics_per_sample.tsv', sep='\t')
m = m[(m.group != 'hepatitis') | (m.sid.isin(CFG.ids('hepatitis_included')))].copy()
m['group'] = m.group.map({'healthy': 'Healthy', 'liver': 'HCC', 'hepatitis': 'Hepatitis', 'pbmc': 'PBMC'})
m['frac300_pct'] = m.frac300 * 100; m['frac166_pct'] = m.frac166 * 100
wt(m[['sid', 'group', 'frac300_pct', 'median_len', 'frac166_pct']].rename(columns={'sid': 'sample', 'median_len': 'median_length_bp'}), 'S8abc.tsv')
a = pd.read_csv(R / 'len_adjust_retention_summary.tsv', sep='\t'); a = a[(a.set == 'all') & a.model.isin(['len1', 'len3', 'alpha', 'full'])]
wt(a[['model', 'component', 'n', 'AllAll_retention', 'ret_lo', 'ret_hi', 'median_retention_59']], 'S8d.tsv')
c = pd.read_csv(R / 'len_corr_summary.tsv', sep='\t'); c = c[c.subset.isin(['HCC', 'healthy']) & c.metric.isin(['frac300', 'frac166', 'median_len'])]
wt(c[['subset', 'component', 'metric', 'AllAll_rho', 'median_rho_59']], 'S8e.tsv')
A, D, E = rt('S8abc.tsv'), rt('S8d.tsv'), rt('S8e.tsv'); rng = np.random.default_rng(5)
G = ['Healthy', 'HCC', 'Hepatitis', 'PBMC']; GC_ = {'Healthy': HEALTHY, 'HCC': HCC, 'Hepatitis': HEP, 'PBMC': PBMC}; GM = {'Healthy': 'o', 'HCC': '^', 'Hepatitis': 'D', 'PBMC': 's'}
cv = SC(126, 8)
for k, (col, yl, ttl, lt) in enumerate((('frac300_pct', 'Fragments > 300 bp (%)', 'Long fragments', 'a'), ('median_length_bp', 'Median fragment length (bp)', 'Median length', 'b'), ('frac166_pct', 'Fragments 156–176 bp (%)', 'Mononucleosome mode', 'c'))):
    x0 = 1 + k * 60; cv.panel(lt, x0, 1); cv.title(x0 + 7, 2.2, ttl); ax = cv.ax(x0 + 12, 11, 42, 36)
    for i, g in enumerate(G):
        v = A[A.group == g][col].values; box(ax, [v], [i], [GC_[g]], w=.5)
        ax.scatter(i + rng.uniform(-.2, .2, len(v)), v, s=3 if len(v) > 100 else 7, marker=GM[g], color=GC_[g], alpha=.4 if len(v) > 100 else .85, lw=0, zorder=2)
    ax.set_xticks(range(4)); ax.set_xticklabels([f'{g}\n(n = {(A.group == g).sum():,})' for g in G], fontsize=5.5); ax.set_ylabel(yl, labelpad=1)
cv.panel('d', 1, 60); cv.title(8, 61.2, 'Effect retained after length adjustment (All markers)')
ax = cv.ax(24, 71, 44, 36); comps = ['strict-mixed', 'WFV', 'strictM', 'BFV', 'FPI', 'sOUF']; mods = [('len1', 'Length > 300 bp', '#a0aec0'), ('len3', '+ median, 166 bp mode', '#4a5568'), ('alpha', '+ PBMC-like score', '#d9ab00'), ('full', '+ age, sex, depth', '#1a1a1a')]
for j, (mod, nm, col) in enumerate(mods):
    for i, cp in enumerate(comps):
        r = D[(D.model == mod) & (D.component == cp)].iloc[0]; y = i + (j - 1.5) * .17
        ax.plot([r.ret_lo, r.ret_hi], [y, y], color=col, lw=1.0); ax.scatter([r.AllAll_retention], [y], s=7, color=col, lw=0, zorder=3)
ax.axvline(1, color=PALE, ls='--', lw=1.0, zorder=0); ax.set_yticks(range(6)); ax.set_yticklabels(comps); ax.invert_yaxis(); ax.set_xlabel('Retained effect (adjusted / unadjusted)', labelpad=1)
ax.legend(handles=[Line2D([], [], color=c, marker='o', ms=2.5, label=n) for _, n, c in mods], frameon=False, fontsize=5.5, loc='upper left', bbox_to_anchor=(1.02, 1.0), handlelength=1.2, labelspacing=.3)
cv.panel('e', 112, 60); cv.title(119, 61.2, 'Within-group Spearman ρ with length metrics')
cm = plt.get_cmap('RdBu_r'); norm = TwoSlopeNorm(vmin=-.6, vcenter=0, vmax=.6); cols = [(s, mt) for s in ('HCC', 'healthy') for mt in ('frac300', 'frac166', 'median_len')]
cl = {'frac300': '>300', 'frac166': '166 bp', 'median_len': 'Median'}; comps2 = ['strict-mixed', 'WFV', 'strictM', 'sOUF', 'BFV', 'FPI']
ov = cv.ov; hx, hy, cw, ch = 131, 74, 7, 5.5
for j, (s, mt) in enumerate(cols):
    ov.text(hx + j * cw + cw / 2, hy - 1, cl[mt], fontsize=5.5, ha='center', va='bottom', rotation=0)
    for i, cp in enumerate(comps2):
        r = E[(E.subset == s) & (E.metric == mt) & (E.component == cp)].iloc[0].AllAll_rho
        ov.add_patch(plt.Rectangle((hx + j * cw + (2 if j >= 3 else 0), hy + i * ch), cw, ch, facecolor=cm(norm(r)), edgecolor='white', lw=.5)); ov.text(hx + j * cw + (2 if j >= 3 else 0) + cw / 2, hy + i * ch + ch / 2, f'{r:.2f}', fontsize=5.5, ha='center', va='center')
for i, cp in enumerate(comps2): ov.text(hx - 1, hy + i * ch + ch / 2, cp, fontsize=5.5, ha='right', va='center')
ov.text(hx + 1.5 * cw, hy - 4.5, 'HCC', fontsize=5.6, ha='center', fontweight='bold'); ov.text(hx + 4.5 * cw + 2, hy - 4.5, 'Healthy', fontsize=5.6, ha='center', fontweight='bold')
cbar(cv, cm, norm, 138, 112, 26, 'Spearman ρ', [-.6, 0, .6])
cv.save()
