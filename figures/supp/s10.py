import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
from sc import *
CO_DIR = AN / 'component_origin/results/'; THA_DIR = AN / 'tumour_hepatitis_age/results/'
cand = [('tumour_tissue_a0.1', 'Tumour tissue mixed in (10%)', 'tumour'), ('TISSUE_tumour_minus_adjacent(paired)', 'Tumour − adjacent non-tumour (22 pairs)', 'tumour'),
        ('nontumour_liver_tissue_a0.02', 'Non-tumour liver tissue mixed in (2%)', 'liver'), ('TISSUE_nontumour_minus_plasma', 'Non-tumour liver tissue − plasma', 'liver'),
        ('hepatocyte_a0.02', 'Hepatocyte reads mixed in (2%)', 'liver'), ('Colon-Left-Epithelial_a0.02', 'Colon epithelium mixed in (2%, control)', 'control'),
        ('WBC_bulk(CNVS)_a0.02', 'Leukocyte reads mixed in (2%)', 'blood'), ('Blood-Granulocytes_a0.02', 'Neutrophil reads mixed in (2%)', 'blood'),
        ('monocyte_up_neutrophil_down_d0.05', 'Monocyte up, neutrophil down (5%)', 'blood'), ('lymphocyte_up_neutrophil_down_d0.05', 'Lymphocyte up, neutrophil down (5%)', 'blood'),
        ('OBS_PBMC_minus_plasma', 'PBMC − plasma (observed)', 'blood'), ('OBS_Hepatitis_10', 'Hepatitis − healthy (observed, n = 10)', 'hepatitis')]
s = pd.read_csv(CO_DIR / 'similarity_full.tsv', sep='\t'); s = s[s.component.isin(['P', 'H'])]
rows = []
for cid, lab, cl in cand:
    for comp in ('P', 'H'):
        r = s[(s.candidate == cid) & (s.component == comp)]
        if len(r): r = r.iloc[0]; rows.append((cid, lab, cl, comp, r.cosine, r.cos_lo, r.cos_hi))
wt(pd.DataFrame(rows, columns=['candidate', 'label', 'origin_class', 'component', 'cosine', 'ci_lo', 'ci_hi']), 'S10a.tsv')
o = pd.read_csv(CO_DIR / 'observed_vs_observed_full.tsv', sep='\t'); o = o[o.component == 'H_nof131'][['vs', 'metrics', 'cosine', 'lo', 'hi']]; wt(o, 'S10c.tsv')
d = pd.read_csv(THA_DIR / 'hep_contrasts_included.tsv', sep='\t'); d = d[d.metric.isin(['mixed', 'U', 'strictM']) & d.stratum.isin(strata_order())]
dp = d.pivot_table(index=['metric', 'stratum'], columns='contrast', values='d').reset_index(); dp.columns.name = None; wt(dp.rename(columns={'HCC-H': 'HCC_minus_healthy', 'hep-H': 'hepatitis_minus_healthy', 'HCC-hep': 'HCC_minus_hepatitis'}), 'S10b.tsv')
A, O, Dd = rt('S10a.tsv'), rt('S10c.tsv'), rt('S10b.tsv')
A['candidate'] = A.candidate.str.replace('Hepatitis_10', 'Hepatitis_10')
CC = {'tumour': TUMOUR, 'liver': '#2e8b57', 'control': PALE, 'blood': PBMC, 'hepatitis': HEP}
cv = SC(124, 10)
cv.panel('a', 1, 1); cv.title(8, 2.2, 'Cosine of origin candidates with the two components')
ax = cv.ax(66, 12, 40, 52); n = len(cand)
for i, (cid, lab, cl) in enumerate(cand):
    for comp, off, filled in (('P', -.16, False), ('H', .16, True)):
        r = A[(A.candidate == cid) & (A.component == comp)]
        if not len(r): continue
        r = r.iloc[0]; ax.plot([r.ci_lo, r.ci_hi], [i + off] * 2, color=CC[cl], lw=1.0); ax.scatter([r.cosine], [i + off], s=12, color=CC[cl] if filled else 'white', edgecolor=CC[cl], linewidths=1.0, zorder=3)
ax.axvline(0, color=PALE, ls='--', lw=1.0, zorder=0); ax.set_yticks(range(n)); ax.set_yticklabels([c[1] for c in cand], fontsize=5.5); ax.invert_yaxis(); ax.set_xlim(-.7, 1.0); ax.set_xlabel('Cosine (bootstrap 95% CI)', labelpad=1)
cv.fig.legend(handles=[Line2D([], [], marker='o', color=MUTED, mfc='white', ls='', ms=3.5, label='P: sOUF, BFV, FPI vector'), Line2D([], [], marker='o', color=MUTED, ls='', ms=3.5, label='H: strict-mixed, WFV, uniform,\nfully methylated vector')], frameon=False, fontsize=5.5, loc='upper left', bbox_to_anchor=(123 / W, 1 - 50 / cv.H), ncol=1, labelspacing=.6)
cv.panel('b', 116, 1); cv.title(123, 2.2, 'H axis, observed groups')
ax = cv.ax(142, 12, 34, 24); nm = {'P_slope': 'P (tumour-fraction slope)', 'Hepatitis_10': 'Hepatitis (n = 10)', 'PBMC_minus_plasma': 'PBMC − plasma', 'H_HCCall': 'H, all HCC'}
O = O[O.vs.isin(nm)].reset_index(drop=True)
for i, r in O.iterrows(): ax.plot([r.lo, r.hi], [i] * 2, color=INK, lw=1.0); ax.scatter([r.cosine], [i], s=12, color=INK, lw=0, zorder=3)
ax.axvline(0, color=PALE, ls='--', lw=1.0, zorder=0); ax.set_yticks(range(len(O))); ax.set_yticklabels([nm[v] for v in O.vs], fontsize=5.5); ax.invert_yaxis(); ax.set_xlim(-.6, 1.05); ax.set_xlabel('Cosine with H (no-fraction HCC, n = 131)', labelpad=1, fontsize=5.5)
cv.panel('c', 1, 80); cv.title(8, 81.2, 'Per-stratum effects: hepatitis (n = 10) against HCC, both minus healthy')
ml = {'mixed': 'strict-mixed (%p)', 'U': 'WFV (×10$^{-2}$)', 'strictM': 'fully methylated fragments (%p)'}
for k, m in enumerate(['mixed', 'U', 'strictM']):
    ax = cv.ax(14 + k * 58, 92, 36, 26); g = Dd[Dd.metric == m]
    ax.scatter(g.HCC_minus_healthy * 100, g.hepatitis_minus_healthy * 100, s=6, color=HEP, marker='D', alpha=.7, lw=0)
    ax.axhline(0, color=GRID, lw=.8, zorder=0); ax.axvline(0, color=GRID, lw=.8, zorder=0)
    ns = ((g.HCC_minus_healthy * g.hepatitis_minus_healthy) > 0).sum()
    ax.set_title(f'{ml[m]}: same sign in {ns} of {len(g)}', fontsize=5.5, pad=2, loc='left'); ax.set_xlabel('HCC − healthy', labelpad=1, fontsize=5.5); ax.set_ylabel('Hepatitis − healthy' if k == 0 else '', labelpad=1, fontsize=5.5)
cv.save()
