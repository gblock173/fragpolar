"""Supplementary Figure 17: tumour fraction vs WFV, strict-mixed and FPI (HCC with tumour fraction, n = 46; all markers).
Values = sample - healthy Train median (same reference as Fig. 4c); rho and bootstrap CI from the tumour-fraction analysis table."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
from sc import *
CDR = AN / 'tumour_hepatitis_age/results'
f = pd.read_csv(CDR / 'tf_f_samples.tsv', sep='\t').set_index('AnalysisID')
st = pd.read_csv(CDR / 'sample_table.tsv', sep='\t')
L = pd.read_parquet(CDR / 'sample_strata_long.parquet'); L = L[L.stratum == 'All|All'].pivot(index='sid', columns='metric', values='value'); L['BV'] = L.B / L.V
ht = st[(st.group == 'healthy') & (st.split == 'Train')].AnalysisID; assert len(ht) == 673
EP = [('U', 'WFV', 1000, 'WFV − healthy median (×10$^{-3}$)'), ('mixed', 'strict-mixed', 100, 'Strict-mixed − healthy median (%p)'), ('BV', 'FPI', 1, 'FPI − healthy median')]
out = f[['f', 'split']].rename_axis('sample').reset_index()
for m, nm, sc_, _ in EP:
    out[f'{nm}_minus_healthy_median'] = L.loc[out['sample'], m].values - L.loc[ht, m].median()
wt(out.rename(columns={'f': 'tumour_fraction'}), 'S17.tsv')
dr = pd.read_csv(CDR / 'tf_dose_response.tsv', sep='\t'); dr = dr[(dr.stratum == 'All|All') & (dr.ref == 'median') & dr.endpoint.isin(['U', 'mixed', 'BV']) & dr.subset.isin(['All46', 'Train28', 'Test18'])]
dr = dr.assign(metric=dr.endpoint.map({'U': 'WFV', 'mixed': 'strict-mixed', 'BV': 'FPI'}))
wt(dr[['metric', 'subset', 'n', 'rho', 'rho_lo', 'rho_hi', 'slope_SD_per_0p1f', 'intercept_f0', 'icpt_lo', 'icpt_hi']], 'S17_spearman.tsv')
D, R = rt('S17.tsv'), rt('S17_spearman.tsv')
# check: recomputed rho equals the table
from scipy.stats import spearmanr
for m, nm, sc_, _ in EP:
    r = spearmanr(D.tumour_fraction, D[f'{nm}_minus_healthy_median'])[0]; t = R[(R.metric == nm) & (R.subset == 'All46')].rho.iloc[0]; assert abs(r - t) < 1e-5, (nm, r, t)
cv = SC(62, 17)
cv.panel('a', 1, 1, 'Tumour fraction and within-fragment mixing (HCC with tumour fraction, n = 46; all markers)')
for j, (m, nm, sc_, yl) in enumerate(EP):
    ax = cv.ax(16 + j * 58, 15, 40, 32)
    for sp in ('Train', 'Test'):
        g = D[D.split == sp]; ax.scatter(g.tumour_fraction * 100, g[f'{nm}_minus_healthy_median'] * sc_, s=8, marker=SPLIT_MK[sp], color=SPLIT_C[sp], alpha=.8, lw=0, label=f'{sp} (n = {len(g)})')
    zeroline(ax)
    a = R[(R.metric == nm) & (R.subset == 'All46')].iloc[0]; tr = R[(R.metric == nm) & (R.subset == 'Train28')].iloc[0]; te = R[(R.metric == nm) & (R.subset == 'Test18')].iloc[0]
    ax_title(ax, f'{nm[0].upper() + nm[1:]}: ' + f'ρ = {a.rho:.2f} [{a.rho_lo:.2f}, {a.rho_hi:.2f}]'.replace('-', '−'), pad=9)
    ax.text(0, 1.015, f'Train {tr.rho:.2f}, Test {te.rho:.2f}'.replace('-', '−'), transform=ax.transAxes, fontsize=5.5, color=MUTED, va='bottom')
    ax.set_xlabel('Tumour fraction (%)', labelpad=1); ax.set_ylabel(yl, labelpad=1)
    if j == 0: ax.legend(loc='lower right', fontsize=5.5, frameon=False, handletextpad=.1, borderaxespad=0)
cv.note(16, 56.5, 'Line at 0, healthy Train median (n = 673). Spearman ρ with HCC-sample bootstrap 95% CI (5,000 resamples).')
cv.save()
