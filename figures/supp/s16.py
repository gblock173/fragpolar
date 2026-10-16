"""Supplementary Figure 16: fragment length and CpG number (K) adjustment.
a, HCC - healthy difference of per-sample mean K and mean fragment length (all markers); b, effect retained after adjustment;
c, tumour fraction vs fragment length and <120 bp fraction vs sOUF (HCC with tumour fraction, n = 46)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
from sc import *
# ---- a
rows = []
for root, sp in (('train_rebuild', 'Train'), ('group_effects', 'Test')):
    g = pd.read_csv(AN / root / 'results/ge_sample_covariate_group_diffs.tsv', sep='\t')
    g = g[(g.split == sp) & (g.stratum == 'All|All') & g.metric.isin(['K', 'fraglen'])]
    for _, r in g.iterrows():
        rows.append(dict(split=sp, covariate={'K': 'mean CpGs per fragment (K)', 'fraglen': 'mean fragment length (bp)'}[r.metric], n_HCC=r.n_HCC, n_healthy=r.n_healthy,
                         mean_healthy=r.mean_healthy, mean_HCC=r.mean_HCC, delta=r.delta, ci_lo=r.ci_lo, ci_hi=r.ci_hi, p=r.p))
wt(pd.DataFrame(rows), 'S16a.tsv')
# ---- b
MN = {'opp': 'sOUF', 'B': 'BFV', 'R_ratio': 'FPI'}; rb = []
for root, sp in (('train_rebuild', 'Train'), ('group_effects', 'Test')):
    s = pd.read_csv(AN / root / 'results/ge_regression_summary.tsv', sep='\t')
    per = pd.read_csv(AN / root / 'results/ge_regression_adjusted_by_K_fraglen.tsv', sep='\t')
    for m, nm in MN.items():
        lab = 'Train(strict보유)' if (sp == 'Train' and m == 'opp') else sp
        x = s[(s.split == lab) & (s.metric == m)].iloc[0]; v = per[(per.split == lab) & (per.metric == m)].retained_pct
        assert len(v) == 59
        rb.append(dict(model='mean K + mean fragment length (59-strata median, IQR)', split=sp, metric=nm, retained_pct=x.median_retained_pct,
                       lo=np.percentile(v, 25), hi=np.percentile(v, 75)))
k = pd.read_csv(AN / 'fpi_k/results/fpik_regression_summary.tsv', sep='\t'); kb = pd.read_csv(AN / 'fpi_k/results/fpik_regression_by_stratum.tsv', sep='\t')
for sp in ('Train', 'Test'):
    x = k[(k.split == sp) & (k.metric == 'FPI_K')].iloc[0]; v = kb[(kb.split == sp) & (kb.metric == 'FPI_K')].retained_pct
    rb.append(dict(model='mean K + mean fragment length (59-strata median, IQR)', split=sp, metric='FPI_K', retained_pct=x.median_retained_pct, lo=np.percentile(v, 25), hi=np.percentile(v, 75)))
e = pd.read_csv(AN / 'mixture_and_length/results/len_adjust_retention_summary.tsv', sep='\t')
ML = {'len1': 'fraction >300 bp (all markers, 95% CI)', 'len3': '+ median length, 166-bp mode fraction (all markers, 95% CI)'}
for md, lab in ML.items():
    for m in ('sOUF', 'BFV', 'FPI'):
        x = e[(e.set == 'all') & (e.model == md) & (e.component == m)].iloc[0]
        rb.append(dict(model=lab, split='Train + Test', metric=m, retained_pct=x.AllAll_retention * 100, lo=x.ret_lo * 100, hi=x.ret_hi * 100))
wt(pd.DataFrame(rb), 'S16b.tsv')
# ---- c (correlations from pipelines/09_mixture_and_length/ml7_length_tumour_fraction.py)
F = pd.read_csv(AN / 'mixture_and_length/results/len_tumour_fraction_samples.tsv', sep='\t')
R = pd.read_csv(AN / 'mixture_and_length/results/len_tumour_fraction_spearman.tsv', sep='\t')
wt(F[['sample', 'split', 'f', 'median_len', 'frac_short120', 'sOUF']].rename(columns={'f': 'tumour_fraction'}), 'S16c.tsv')
wt(R, 'S16c_spearman.tsv')
A, Bt, C = rt('S16a.tsv'), rt('S16b.tsv'), rt('S16c.tsv'); print(Bt.round(1).to_string())

def fmt_p(p):
    if p >= .01: return f'p = {p:.2f}'
    m, e_ = f'{p:.1e}'.split('e'); return f'p = {m}×10$^{{{int(e_)}}}$'


cv = SC(150, 16); ov = cv.ov
cv.panel('a', 1, 1, 'Per-sample covariates, HCC − healthy (all markers)')
for j, (cov, xl) in enumerate((('mean CpGs per fragment (K)', 'Δ mean K (CpGs per fragment)'), ('mean fragment length (bp)', 'Δ mean fragment length (bp)'))):
    ax = cv.ax(20 + j * 46, 12, 32, 16)
    for i, sp in enumerate(('Train', 'Test')):
        r = A[(A.split == sp) & (A.covariate == cov)].iloc[0]
        ci_row(ax, r.delta, r.ci_lo, r.ci_hi, i, SPLIT_C[sp], SPLIT_MK[sp])
        ax.text(1.02, i, f'{r.delta:+.3f}'.replace('-', '−') if j == 0 else f'{r.delta:+.1f}'.replace('-', '−'), transform=ax.get_yaxis_transform(), va='center', fontsize=5.5)
        ax.text(1.02, i + .42, fmt_p(r.p), transform=ax.get_yaxis_transform(), va='center', fontsize=5.5, color=MUTED)
    refline(ax, 0, 'v'); ax.set_yticks([0, 1]); ax.set_yticklabels(['Train', 'Test'] if j == 0 else []); ax.set_ylim(1.6, -.6); ax.set_xlabel(xl, labelpad=1)
    if j == 0: ax.set_xlim(-.04, .06)
    else: ax.set_xlim(-9, 1)
cv.note(118, 12, 'Train 104 vs 673, Test 73 vs 336.\nMean over marker-covering fragments\nper sample; 95% bootstrap CI and\nMann–Whitney p.')
# ---- b
B0 = 40
cv.panel('b', 1, B0, 'Effect retained after adjustment (adjusted / unadjusted HCC coefficient)')
ax = cv.ax(42, B0 + 8, 76, 42); y = 0; yt = []; yl = []
mcol = {'mean K + mean fragment length (59-strata median, IQR)': None, ML['len1']: '#7a6a9a', ML['len3']: '#3a2a6a'}
for m in ('sOUF', 'BFV', 'FPI', 'FPI_K'):
    ax.text(-.44, y - .75, m.replace('FPI_K', 'FPI$_K$'), transform=ax.get_yaxis_transform(), fontsize=6, fontweight='bold', va='center', ha='left')
    for md in mcol:
        g = Bt[(Bt.metric == m) & (Bt.model == md)]
        if not len(g): continue
        if mcol[md] is None:
            for sp, off in (('Train', -.17), ('Test', .17)):
                r = g[g.split == sp].iloc[0]; ci_row(ax, r.retained_pct, r.lo, r.hi, y + off, SPLIT_C[sp], SPLIT_MK[sp], ms=3)
            lab = 'mean K + mean length'
            ax.text(1.01, y, '%.0f%% / %.0f%% (Train / Test)' % (g[g.split == 'Train'].retained_pct.iloc[0], g[g.split == 'Test'].retained_pct.iloc[0]), transform=ax.get_yaxis_transform(), va='center', fontsize=5.5)
        else:
            r = g.iloc[0]; ci_row(ax, r.retained_pct, r.lo, r.hi, y, mcol[md], 'D', ms=3)
            lab = '>300 bp fraction' if md == ML['len1'] else '+ median length, 166 bp'
            ax.text(1.01, y, f'{r.retained_pct:.0f}% [{r.lo:.0f}–{r.hi:.0f}]', transform=ax.get_yaxis_transform(), va='center', fontsize=5.5)
        yt.append(y); yl.append(lab); y += 1
    y += .9
ax.set_yticks(yt); ax.set_yticklabels(yl, fontsize=5.5); ax.set_ylim(y - .6, -1.3); ax.axvline(100, color=INK, lw=1.0, ls='--', zorder=0)
ax.set_xlim(40, 160); ax.set_xlabel('Effect retained (%)', labelpad=1); ax.tick_params(axis='y', length=0)
hs = [Line2D([], [], marker=SPLIT_MK[k], ls='-', color=SPLIT_C[k], ms=3, label=f'{k}: median [IQR], 59 strata') for k in ('Train', 'Test')]
hs += [Line2D([], [], marker='D', ls='-', color='#3a2a6a', ms=3, label='Length distribution [95% CI]')]
cv.fig.legend(handles=hs, loc='upper left', bbox_to_anchor=(143 / W, 1 - (B0 + 9) / cv.H), fontsize=5.5, frameon=False, handletextpad=.4, labelspacing=.6)
cv.note(143, B0 + 22, 'Mean K + mean length: depth, age,\nsex + mean K + mean fragment\nlength (Fig. 2d).\nLength distribution: fraction\n>300 bp; + median length and\n156–176 bp fraction (HCC 177 vs\nhealthy 1,009).')
# ---- c
C0 = 100
cv.panel('c', 1, C0, 'Fragment shortening tracks tumour fraction (HCC with tumour fraction, n = 46)')
specs = [('tumour_fraction', 'median_len', 'f', 'median_len', 'Tumour fraction (%)', 'Median fragment length (bp)', 100, 1),
         ('tumour_fraction', 'frac_short120', 'f', 'frac_short120', 'Tumour fraction (%)', 'Fragments <120 bp (%)', 100, 100),
         ('frac_short120', 'sOUF', 'frac_short120', 'sOUF', 'Fragments <120 bp (%)', 'sOUF (%)', 100, 100)]
for j, (xc, yc, rx, ry, xl, yl, fx, fy) in enumerate(specs):
    ax = cv.ax(16 + j * 58, C0 + 9, 40, 30)
    for sp in ('Train', 'Test'):
        g = C[C.split == sp]; ax.scatter(g[xc] * fx, g[yc] * fy, s=8, marker=SPLIT_MK[sp], color=SPLIT_C[sp], alpha=.8, lw=0, label=f'{sp} (n = {len(g)})')
    r = R[(R.x == rx) & (R.y == ry) & (R.subset == 'All46')].iloc[0]
    ax_title(ax, f'Spearman ρ = {r.rho:.2f} [{r.ci_lo:.2f}, {r.ci_hi:.2f}]'.replace('-', '−'), pad=2)
    ax.set_xlabel(xl, labelpad=1); ax.set_ylabel(yl, labelpad=1)
    if j == 0: ax.legend(loc='upper right', fontsize=5.5, frameon=False, handletextpad=.1, borderaxespad=0)
cv.save()
