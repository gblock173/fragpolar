import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
from sc import *
V6 = CFG.path('prior_source_data_dir') + '/'
t = pd.read_csv(V6 + 'Figure1b_turnover.tsv', sep='\t'); wt(t, 'S2a.tsv')
c = pd.read_csv(V6 + 'Figure1d_tissue_composition.tsv', sep='\t').drop(columns=['sample_id']); c['tissue_group'] = c.tissue_group.map({'normal': 'Non-tumour', 'tumor': 'Tumour'}); wt(c, 'S2b.tsv')
r = pd.read_csv(V6 + 'Figure1d_residual_samples.tsv', sep='\t'); r['group'] = r.group.map({'healthy': 'Healthy cfDNA', 'HCC': 'HCC cfDNA', 'normal': 'Non-tumour tissue', 'tumour': 'Tumour tissue'}); wt(r, 'S2c.tsv')
g = pd.read_csv(V6 + 'Figure1g_CNV_samples.tsv', sep='\t'); wt(g, 'S2d.tsv')
T, C, Rr, G = rt('S2a.tsv'), rt('S2b.tsv'), rt('S2c.tsv'), rt('S2d.tsv'); rng = np.random.default_rng(4)
cv = SC(126, 2)
cv.panel('a', 1, 1); cv.title(8, 2.2, 'Cell-type DNA fractions from atlas deconvolution')
ax = cv.ax(24, 10, 56, 44); names = ['Blood-T', 'Blood-B', 'Blood-NK', 'Blood-Mono+Macro', 'Blood-Granul', 'Eryth-prog', 'Endothel', 'Liver-Hep']; T = T.set_index('cell_type').loc[names]
for j, (k, col, mk, nm) in enumerate((('pbmc', PBMC, 's', 'PBMC (30)'), ('cfdna_healthy', HEALTHY, 'o', 'Healthy (1,009)'), ('cfdna_hcc', HCC, '^', 'HCC (177)'))):
    y = np.arange(8) + (j - 1) * .26
    for yy, (m, a, b) in zip(y, zip(T['median_' + k], T['q25_' + k], T['q75_' + k])): ax.plot([a * 100, b * 100], [yy, yy], color=col, lw=1.0)
    ax.scatter(T['median_' + k] * 100, y, s=9, color=col, marker=mk, lw=0, zorder=3, label=nm)
ax.set_yticks(range(8)); ax.set_yticklabels(['T', 'B', 'NK', 'Mono/Macro', 'Granulocyte', 'Eryth-prog', 'Endothelial', 'Liver-Hep'], fontsize=5.5); ax.invert_yaxis(); ax.set_xlabel('DNA fraction (%), median and IQR', labelpad=1)
ax.legend(frameon=False, fontsize=5.5, loc='lower right', handletextpad=.2, borderaxespad=0)
cv.panel('b', 92, 1); cv.title(99, 2.2, 'Liver tissue composition (atlas)')
ax = cv.ax(106, 10, 70, 44); cats = ['Hepatocyte', 'Blood cells', 'Endothelial', 'Fibroblast / smooth muscle', 'Pancreas', 'Other']; cc = ['#E69F00', '#56B4E9', '#CC79A7', '#6A4C9C', '#F0E442', '#bdbdbd']
x0 = 0; ticks = []
for gname in ('Non-tumour', 'Tumour'):
    g_ = C[C.tissue_group == gname].sort_values('plot_order'); bot = np.zeros(len(g_)); xs = np.arange(len(g_)) + x0
    for k, col in zip(cats, cc): ax.bar(xs, g_[k].values * 100, bottom=bot, width=1, color=col, lw=0); bot += g_[k].values * 100
    ticks.append((x0 + len(g_) / 2 - .5, f'{gname} ({len(g_)})')); x0 += len(g_) + 3
ax.set_xticks([t[0] for t in ticks]); ax.set_xticklabels([t[1] for t in ticks]); ax.set_ylabel('DNA fraction (%)', labelpad=1); ax.set_xlim(-1, x0 - 2); ax.set_ylim(0, 100)
ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c_) for c_ in cc], labels=cats, frameon=False, fontsize=5.5, ncol=3, loc='upper left', bbox_to_anchor=(0, -.1), handlelength=.8, columnspacing=.8)
cv.panel('c', 1, 72); cv.title(8, 73.2, 'Atlas fit error, 1 − R$^2$')
ax = cv.ax(14, 84, 56, 28); gl = ['Non-tumour tissue', 'Tumour tissue', 'Healthy cfDNA', 'HCC cfDNA']; gcol = [NONTUM, TUMOUR, HEALTHY, HCC]; gmk = ['o', '^', 'o', '^']
for i, (gn, col, mk) in enumerate(zip(gl, gcol, gmk)):
    v = Rr[Rr.group == gn].one_minus_R2.values; box(ax, [v], [i], [col], w=.5)
    ax.scatter(i + rng.uniform(-.2, .2, len(v)), v, s=3 if len(v) > 100 else 7, marker=mk, color=col, alpha=.4 if len(v) > 100 else .85, lw=0, zorder=2)
ax.set_yscale('log'); ax.set_xticks(range(4)); ax.set_xticklabels([f'{g_.replace(" ", chr(10), 1)}\n(n = {(Rr.group == g_).sum():,})' for g_ in gl], fontsize=5.5); ax.set_ylabel('1 − R$^2$ (log)', labelpad=1)
cv.panel('d', 92, 72); cv.title(99, 73.2, 'Fit error and CNV tumour fraction (HCC, n = 46)')
ax = cv.ax(106, 84, 40, 28); ax.scatter(G.cnv_tumour_fraction * 100, G.one_minus_R2, s=9, marker='^', color=HCC, lw=0, alpha=.85); ax.set_xscale('log'); ax.set_yscale('log')
from scipy.stats import spearmanr
rho = spearmanr(G.cnv_tumour_fraction, G.one_minus_R2)[0]; ax.set_xlabel('CNV tumour fraction (%)', labelpad=1); ax.set_ylabel('1 − R$^2$ (log)', labelpad=1); ax.set_xticks([1, 3, 10, 30]); ax.set_xticklabels(['1', '3', '10', '30'])
ax.text(.04, .94, f'Spearman ρ = {rho:.2f}', transform=ax.transAxes, fontsize=5.5, va='top')
cv.save(); print(len(G), rho, Rr.groupby('group').one_minus_R2.median())
