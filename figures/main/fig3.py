import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys; sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from style_v8 import *
from matplotlib.colors import ListedColormap
cv = Canvas(198); ov = overlay(cv); rng = np.random.default_rng(20261015)
CL = {'opposite': HCC, 'mixed': '#a0aec0', 'concordant': HEALTHY}
MK = {'Train': 'o', 'Test': 's'}
def boxstrip(ax, vals, pos, color, marker, w=.3, s=3.5, ms_alpha=.5):
    v = np.asarray(vals, float); v = v[np.isfinite(v)]
    ax.boxplot(v, positions=[pos], widths=w, showfliers=False, whis=1.5, medianprops={'color': INK, 'linewidth': 1.15}, boxprops={'color': MUTED, 'linewidth': 1.}, whiskerprops={'color': MUTED, 'linewidth': 1.}, capprops={'color': MUTED, 'linewidth': 1.})
    ax.scatter(pos + rng.uniform(-.08, .08, len(v)), v, s=s, marker=marker, color=color, alpha=ms_alpha, lw=0, zorder=2)
# ---------------- a
cv.panel('a', 1, 1, 'Fragments at two representative markers (Train matched pair)')
R_ = rd('Fig3a.tsv'); S_ = rd('Fig3a_summary.tsv')
markers = list(R_.marker.unique()); cmap = ListedColormap(['#e2e8f0', '#ffffff', '#1a1a1a'])
rw, rh = 24, 36; x0 = 7; ytop = 20
for mi, mk in enumerate(markers):
    gx = x0 + mi * 56
    ov.text(gx + 25, 10.5, f'Marker {mi + 1}  ({mk})', fontsize=5.6, ha='center', va='bottom', color=MUTED)
    for gi, grp in enumerate(['HCC', 'Healthy']):
        sub = R_[(R_.marker == mk) & (R_.group == grp)]; sm = S_[(S_.marker == mk) & (S_.group == grp)].iloc[0]
        nr = int(sub.row.max()) + 1; nc = int(sub.groupby('row').size().max()); Y = np.zeros((nr, nc))
        for r_, gg in sub.groupby('row'):  # CpGs of one fragment, left-aligned in genomic order
            for c_, v in enumerate(gg.sort_values('pos').call): Y[int(r_), c_] = 2 if v == 1 else 1
        axr = cv.ax(gx + gi * 27 + 2.2, ytop, rw - 2.2, rh * nr / 45); axr.imshow(Y, cmap=cmap, vmin=0, vmax=2, aspect='auto', interpolation='nearest'); bare(axr)
        cls = sub.drop_duplicates('row').sort_values('row')['class'].tolist()
        for t, c in enumerate(cls): ov.add_patch(plt.Rectangle((gx + gi * 27, ytop + t * rh / 45), 1.6, rh / 45 * 1.0, color=CL[c], lw=0))
        ov.text(gx + gi * 27 + 2.2 + (rw - 2.2) / 2, ytop - 1.0, grp, fontsize=5.8, ha='center', va='bottom', fontweight='bold', color=HCC if grp == 'HCC' else HEALTHY)
        ov.text(gx + gi * 27 + 2.2 + (rw - 2.2) / 2, ytop + rh + 1.0, f'n = {sm.n_fragments}\nsOUF {sm.n_opposite / sm.n_fragments * 100:.0f}%\nBFV {sm.BFV:.3f}', fontsize=5.5, ha='center', va='top')
ov.text(7, ytop + rh + 11.0, 'Rows: up to 45 random fragments (≥ 4 CpGs), CpGs left-aligned in genomic order;\nblack, methylated; white, unmethylated; grey, beyond fragment end.', fontsize=5.5, color=MUTED, va='top')
for k, (nm, c) in enumerate((('opposite uniform', HCC), ('mixed', '#a0aec0'), ('same-direction uniform', HEALTHY))):
    xx = 7 + k * 25 + (6 if k == 2 else 0) * 0
    ov.add_patch(plt.Rectangle((7 + [0, 26, 40][k], ytop + rh + 17.4), 2, 2, color=c, lw=0)); ov.text(10 + [0, 26, 40][k], ytop + rh + 18.4, nm, fontsize=5.5, va='center')
# ---------------- b
cv.panel('b', 120, 1, 'Observed vs expected ΔBFV')
Bb = rd('Fig3b.tsv'); axb = cv.ax(130, 11, 46, 40); f = 1e3
for sp in ('Train', 'Test'):
    g = Bb[Bb.split == sp]
    axb.scatter(g.expected_dBFV_same_direction_mix * f, g.observed_dBFV * f, s=9, marker=MK[sp], facecolor='none', edgecolor=PALE, lw=.6, zorder=2)
    axb.scatter(g.expected_dBFV_opposite_mix * f, g.observed_dBFV * f, s=7, marker=MK[sp], color=SPLIT_C[sp], alpha=.8, lw=0, zorder=3)
lim = [-.5, Bb.observed_dBFV.max() * f * 1.08]; axb.plot(lim, lim, color=INK, lw=1.0, ls='--', zorder=1); axb.set_xlim(lim); axb.set_ylim(lim)
axb.set_xlabel('Expected ΔBFV, opposite-uniform mix ($\\times$10$^{-3}$)', labelpad=1, fontsize=5.6); axb.set_ylabel('Observed ΔBFV ($\\times$10$^{-3}$)', labelpad=1)
st = {sp: (g.observed_dBFV.sum() / g.expected_dBFV_opposite_mix.sum(), np.corrcoef(g.observed_dBFV, g.expected_dBFV_opposite_mix)[0, 1]) for sp, g in Bb.groupby('split')}
ov.text(120, 64.2, f"Observed/expected (sum over strata): {st['Train'][0]:.2f} (Train), {st['Test'][0]:.2f} (Test);\nPearson r: {st['Train'][1]:.2f}, {st['Test'][1]:.2f}", fontsize=5.5, va='top', color=MUTED)
h = [Line2D([], [], ls='', marker='o', ms=3, color=SPLIT_C['Train'], label='Train'), Line2D([], [], ls='', marker='s', ms=3, color=SPLIT_C['Test'], label='Test'),
     Line2D([], [], ls='', marker='o', ms=3, mfc='none', mec=PALE, label='Same-direction mix')]
axb.legend(handles=h, frameon=False, fontsize=5.5, loc='lower right', handletextpad=.1, borderaxespad=0)
ov.text(120, 60, 'Points, strata (n = 59 per split). Dashed line, observed = expected.', fontsize=5.5, color=MUTED, va='top')
# ---------------- c
R2 = 80
cv.panel('c', 1, R2 - 3, 'sOUF removal')
C_ = rd('Fig3c.tsv'); axc = cv.ax(11, R2 + 7, 34, 36)
for i, sp in enumerate(('Train', 'Test')):
    g = C_[C_.split == sp]; boxstrip(axc, g.reduction, i, SPLIT_C[sp], MK[sp], w=.45, s=4)
    a = g[g.stratum == 'All|All'].iloc[0]; axc.errorbar([i + .3], [a.reduction], yerr=[[a.reduction - a.red_lo], [a.red_hi - a.reduction]], fmt='D', ms=3.2, color=HCC, mfc=HCC, mec=HCC, elinewidth=1.0, capsize=0, zorder=5)
    ov.text(0, 0, '')
axc.axhline(.5, color=INK, lw=1., ls='--', zorder=1); axc.set_xlim(-.6, 1.8); axc.set_xticks([0, 1]); axc.set_xticklabels(['Train', 'Test']); axc.set_ylabel('Reduction of ΔBFV', labelpad=1)
nge = {sp: int((C_[C_.split == sp].reduction >= .5).sum()) for sp in ('Train', 'Test')}
ov.text(11, R2 + 46.5, f"≥ 0.5 in {nge['Train']} (Train), {nge['Test']} (Test) strata.\nDiamond: all markers, 95% CI.", fontsize=5.5, color=MUTED, va='top')
# ---------------- d
cv.panel('d', 52, R2 - 3, 'Same-number removal controls')
D_ = rd('Fig3d.tsv'); cats = ['sOUF', 'Random', 'Same-direction uniform', 'Strict-mixed']; axd = cv.ax(61, R2 + 7, 64, 36)
for i, c in enumerate(cats):
    for sp, off in (('Train', -.2), ('Test', .2)):
        g = D_[(D_.removed == c) & (D_.split == sp)]; boxstrip(axd, g.reduction, i + off, SPLIT_C[sp], MK[sp], w=.3, s=3)
        a = g[g.stratum == 'All|All'].iloc[0]; axd.scatter(i + off + .17, a.reduction, marker='D', s=9, color=HCC, lw=0, zorder=5)
        axd.text(i + off, 1.0, f'{a.reduction * 100:.0f}'.replace('-', '−'), transform=axd.get_xaxis_transform(), fontsize=5.5, ha='center', va='bottom', color=SPLIT_C[sp], fontweight='bold')
axd.axhline(0, color=GRID, lw=.8, zorder=0); axd.axhline(.5, color=INK, lw=1., ls='--', zorder=1)
axd.set_xlim(-.6, 3.6); axd.set_xticks(range(4)); axd.set_xticklabels(['sOUF', 'Random', 'Same-dir.\nuniform', 'Strict-\nmixed'], fontsize=5.5); axd.set_ylabel('Reduction of ΔBFV', labelpad=1)
ov.text(61, R2 + 49, 'Top numbers, all-markers reduction (%) in Train and Test; diamonds, all markers.', fontsize=5.5, color=MUTED, va='top')
# ---------------- e
cv.panel('e', 130, R2 - 3, 'Read-split check')
E_ = rd('Fig3e.tsv'); cats = ['Other half, sOUF', 'Same half, sOUF', 'Other half, same-direction']; axe = cv.ax(139, R2 + 7, 38, 36)
for i, c in enumerate(cats):
    for sp, off in (('Train', -.2), ('Test', .2)):
        g = E_[(E_.dose_source == c) & (E_.split == sp)]; boxstrip(axe, g.fraction_explained, i + off, SPLIT_C[sp], MK[sp], w=.3, s=3)
        a = g[g.stratum == 'All|All'].iloc[0]; axe.scatter(i + off + .17, a.fraction_explained, marker='D', s=9, color=HCC, lw=0, zorder=5)
        axe.text(i + off, 1.0, f'{a.fraction_explained * 100:.0f}'.replace('-', '−'), transform=axe.get_xaxis_transform(), fontsize=5.5, ha='center', va='bottom', color=SPLIT_C[sp], fontweight='bold')
axe.axhline(0, color=GRID, lw=.8, zorder=0); axe.set_xlim(-.6, 2.6); axe.set_xticks(range(3)); axe.set_xticklabels(['Other\nhalf\nsOUF', 'Same\nhalf\nsOUF', 'Other\nhalf\nsame-dir.'], fontsize=5.5)
axe.set_ylim(-.4, 1.2); axe.set_ylabel('ΔBFV explained', labelpad=1)
ov.text(130, R2 + 52, 'Top numbers, all-markers value (%).\nStrata outside the axis not shown.', fontsize=5.5, color=MUTED, va='top')
# ---------------- f
R3 = 144
cv.panel('f', 1, R3 - 3, 'Similarity to synthetic generators (f = 5%)')
F_ = rd('Fig3f.tsv'); axf = cv.ax(11, R3 + 7, 70, 34)
for i, gname in enumerate(('cos_G1', 'cos_G2', 'cos_G3')):
    for sp, off in (('Train', -.2), ('Test', .2)):
        g = F_[F_.split == sp]; boxstrip(axf, g[gname], i + off, SPLIT_C[sp], MK[sp], w=.3, s=3)
        axf.text(i + off, 1.0, f'{g[gname].median():.2f}'.replace('-', '−'), transform=axf.get_xaxis_transform(), fontsize=5.5, ha='center', va='bottom', color=SPLIT_C[sp], fontweight='bold')
axf.axhline(0, color=INK, lw=1., ls='--', zorder=1); axf.set_ylim(-1, 1); axf.set_xlim(-.6, 2.6); axf.set_xticks(range(3)); axf.set_xticklabels(['G1\nopposite uniform', 'G2\nwithin-fragment', 'G3\nindependent flips'], fontsize=5.5)
axf.set_ylabel('Cosine similarity', labelpad=1)
npol = F_.groupby('split').polarization.sum()
ov.text(11, R3 + 49, f"Medians on top. Strata meeting all polarisation criteria: {npol['Train']} (Train), {npol['Test']} (Test) of 59; within-fragment criterion: 0.", fontsize=5.5, color=MUTED, va='top')
# key
kx = 96; ky = R3 + 4
ov.text(kx, ky, 'Key', fontsize=6.0, fontweight='bold', va='top')
for k, (lab, m) in enumerate((('Train (HCC and matched healthy pairs, n = 104)', 'o'), ('Test (n = 73 pairs)', 's'))):
    ov.scatter([kx + 1.5], [ky + 6 + k * 4.5], s=14, marker=m, color=SPLIT_C['Train' if k == 0 else 'Test'], lw=0); ov.text(kx + 4.5, ky + 6 + k * 4.5, lab, fontsize=5.5, va='center')
ov.scatter([kx + 1.5], [ky + 15.5], s=14, marker='D', color=HCC, lw=0); ov.text(kx + 4.5, ky + 15.5, 'all-markers value (stratum All|All)', fontsize=5.5, va='center')
ov.text(kx, ky + 20, 'Boxes: median, quartiles, 1.5 × IQR whiskers over 59 strata.\nStrata overlap; counts describe extent, not independent tests.', fontsize=5.5, color=MUTED, va='top')
cv.save('Figure_3')
