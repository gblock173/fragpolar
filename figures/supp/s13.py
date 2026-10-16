"""Supplementary Fig. 13 (new, earlier design pass): marker strata (healthy-baseline band x genomic annotation).
Counts recomputed from <data_dir>/marker-matrix/marker annotation.parquet (one row per marker; band from Train-healthy AMF: U < 0.2, mid, M > 0.8;
annot / genomic_context exclusive classes; atlas cell-type marker flags) and checked against <data_dir>/strata_definition.tsv
(60 candidate strata, 'used' = >= 100 markers -> 59 analysed). Writes supp/source_data/S13.tsv. Design = earlier Fig. 2b."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys
sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / 'main'))
from style_v8 import *
from matplotlib.colors import Normalize

D3 = Path(CFG.DATA)
CS = pd.read_parquet(CFG.inp('marker_annotation')); SDf = pd.read_csv(CFG.inp('strata_definition'), sep='\t')
BASES = ['All', 'U', 'mid', 'M']
ANN = [('All', 'All markers', None), ('Island', 'CpG island', 'annot'), ('Shore', 'CpG shore', 'annot'), ('Enhancer', 'Enhancer', 'annot'),
       ('TSS', 'TSS', 'annot'), ('ReprPC', 'Polycomb (ReprPC)', 'annot'), ('promoter', 'Promoter', 'genomic_context'), ('exon', 'Exon', 'genomic_context'),
       ('intron', 'Intron', 'genomic_context'), ('5UTR', '5′UTR', 'genomic_context'), ('3UTR', '3′UTR', 'genomic_context'),
       ('intergenic', 'Intergenic', 'genomic_context'), ('Liver-Hep', 'Liver-Hep markers', 'atlas'), ('Immune5', 'Immune5 markers', 'atlas'),
       ('Eryth-prog', 'Eryth-prog markers', 'atlas')]
rows = []
for b in BASES:
    mb = np.ones(len(CS), bool) if b == 'All' else (CS.band == b).to_numpy()
    for a, lab, col in ANN:
        if col is None:
            ma = np.ones(len(CS), bool)
        elif col == 'atlas':
            ma = CS[f'atlas_{a}'].to_numpy(bool)
        else:
            ma = (CS[col] == a).to_numpy()
        n = int((mb & ma).sum()); rows.append(dict(stratum=f'{b}|{a}', baseline=b, annotation=a, label=lab, n_markers=n, analysed=n >= 100))
T = pd.DataFrame(rows)
chk = T.merge(SDf[['stratum', 'n_markers', 'used']], on='stratum', suffixes=('', '_def'))
assert len(chk) == 60 and (chk.n_markers == chk.n_markers_def).all() and (chk.analysed == chk.used).all(), chk[chk.n_markers != chk.n_markers_def]
assert (CS.band.eq('U') == (CS.baseline_amf < .2)).all() and (CS.band.eq('M') == (CS.baseline_amf > .8)).all()
T.to_csv(SSD / 'S13.tsv', sep='\t', index=False)
NM = T.pivot(index='annotation', columns='baseline', values='n_markers').loc[[a for a, _, _ in ANN], BASES]
MK = {b: int(NM.loc['All', b]) for b in BASES}; assert MK['U'] + MK['mid'] + MK['M'] == MK['All']
n_used = int(T.analysed.sum()); print(NM.to_string()); print('analysed', n_used, 'not analysed', T.loc[~T.analysed, ['stratum', 'n_markers']].values.tolist())

cv = Canvas(76, w=120); ov = cv.ov
X0 = 2.0
cv.text(X0, 1.0, f'{MK["All"]:,} markers, banded by Train-healthy median methylation (AMF)', FS['label'], 'bold')
BANDCOL = {'All': '#e8ecf1', 'U': '#ffffff', 'mid': '#a9a9a9', 'M': '#3a3a3a'}; BANDTXT = {'All': INK, 'U': INK, 'mid': INK, 'M': 'white'}
SX0, SWD = X0 + 24, 70.0; ytk = 5.6; ybar = ytk + 2.6; hb = 3.6
for v in (0, .2, .8, 1):
    xx = SX0 + SWD * v; ov.add_line(Line2D([xx, xx], [ybar - .7, ybar], color=MUTED, lw=1.0)); cv.text(xx, ytk, f'{v:g}', FS['tick'], ha='center', color=INK2)
for b, lo, hi, lab in (('U', 0, .2, 'U  < 0.2'), ('mid', .2, .8, 'mid  0.2–0.8'), ('M', .8, 1, 'M  > 0.8')):
    x0, x1 = SX0 + SWD * lo, SX0 + SWD * hi
    ov.add_patch(Rectangle((x0, ybar), x1 - x0, hb, fc=BANDCOL[b], ec=INK, lw=.6))
    cv.text((x0 + x1) / 2, ybar + hb / 2, lab, FS['label'], 'bold', ha='center', va='center_baseline', color=BANDTXT[b])
    cv.text((x0 + x1) / 2, ybar + hb + .6, f'{MK[b]:,}', FS['label'], ha='center')
yt0 = ybar + hb + 5.0; RH = 3.15; LWD = 24.0; CWD = 15.0; XT = X0 + LWD
cv.text(X0, yt0 + RH / 2, 'Annotation', FS['small'] + .2, 'bold', va='center_baseline', color=INK2)
for k, b in enumerate(BASES):
    ov.add_patch(Rectangle((XT + k * CWD + .25, yt0 + .2), CWD - .5, RH - .4, fc=BANDCOL[b], ec=INK if b == 'U' else 'none', lw=.6))
    cv.text(XT + (k + .5) * CWD, yt0 + RH / 2, b, FS['label'], 'bold', ha='center', va='center_baseline', color=BANDTXT[b])
cn = Normalize(2.0, np.log10(MK['All'])); ramp = plt.get_cmap('Greys')
for i, (a, lab, _) in enumerate(ANN):
    yy = yt0 + (i + 1) * RH
    cv.text(X0, yy + RH / 2, lab, FS['small'] + .2, va='center_baseline')
    for k, b in enumerate(BASES):
        n = int(NM.loc[a, b]); x_ = XT + k * CWD
        if n >= 100:
            ov.add_patch(Rectangle((x_ + .25, yy + .2), CWD - .5, RH - .4, fc=ramp(.04 + .30 * cn(np.log10(n))), ec='none'))
            cv.text(x_ + CWD / 2, yy + RH / 2, f'{n:,}', FS['small'] + .1, ha='center', va='center_baseline')
        else:
            ov.add_patch(Rectangle((x_ + .25, yy + .2), CWD - .5, RH - .4, fc='white', ec=INK2, lw=.5, hatch='/////'))
            cv.text(x_ + CWD / 2, yy + RH / 2, f'{n}', FS['small'] + .1, ha='center', va='center_baseline', color=INK2, bbox=dict(fc='white', ec='none', pad=.3))
for i_sep in (1, 6, 12):
    ys = yt0 + (i_sep + 1) * RH; ov.add_line(Line2D([X0, XT + 4 * CWD], [ys, ys], color=GRID, lw=.6, ls=(0, (2, 1.5))))
xn = XT + 4 * CWD + 2.0
for (i0, i1, lab) in ((1, 5, 'regulatory\nannotation'), (6, 11, 'gene\ncontext'), (12, 14, 'atlas\ncell-type\nmarker sets')):
    ya, yb = yt0 + (i0 + 1) * RH + .4, yt0 + (i1 + 2) * RH - .4
    ov.add_line(Line2D([xn - 1.0, xn - .5, xn - .5, xn - 1.0], [ya, ya, yb, yb], color=INK2, lw=.8))
    cv.note(xn + .5, (ya + yb) / 2, lab, va='center')
yend = yt0 + 16 * RH
cv.note(X0, yend + 1.6, f'{n_used} of 60 strata analysed (≥ 100 markers). Regulatory annotation and gene context classes are exclusive within each\n'
        'column; atlas marker sets overlap with them. Shading, log10 marker count; hatched, not analysed.')
cv.save(SUPP / 'Supplementary_Figure_S13', svg=False)
