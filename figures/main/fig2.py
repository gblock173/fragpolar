import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys; sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from style_v8 import *
import matplotlib as mpl
from matplotlib.colors import TwoSlopeNorm
cv = Canvas(185); ov = overlay(cv); rng = np.random.default_rng(20261015)
# ---------------- a
cv.panel('a', 1, 1, 'HCC − healthy effect on 12 summaries across 59 strata (Test, matched pairs)')
A = rd('Fig2a.tsv'); strata = [s for s in strata_order() if s in set(A.stratum)]; assert len(strata) == 59
hx, hy, cw, ch, gap = 20, 13, 2.5, 3.0, 1.6
cmap = plt.get_cmap('RdBu_r'); norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
heads = {'All': 'All markers', 'U': 'U baseline', 'mid': 'Mid baseline', 'M': 'M baseline'}
x = hx; xs = {}
last = None
for s in strata:
    b = s.split('|')[0]
    if last is not None and b != last: x += gap
    if b != last: ov.text(0, 0, ''); start = x; last = b; hstart = {'b': b, 'x': x}
    xs[s] = x; x += cw
    A_end = x
for b in heads:
    ss = [s for s in strata if s.startswith(b + '|')]; x0, x1 = xs[ss[0]], xs[ss[-1]] + cw
    ov.text((x0 + x1) / 2, hy - 2.4, heads[b], fontsize=5.8, ha='center', va='bottom', fontweight='bold')
    ov.plot([x0, x1], [hy - 1.7, hy - 1.7], color=MUTED, lw=1.0)
for mi, m in enumerate(METRICS):
    ov.text(hx - 3.5, hy + mi * ch + ch / 2, m[1], fontsize=5.6, ha='right', va='center')
    ov.add_patch(plt.Rectangle((hx - 2.7, hy + mi * ch + .2), 1.2, ch - .4, color=CAT[m[2]], lw=0))
    for s in strata:
        r = A[(A.stratum == s) & (A.metric == m[1])].iloc[0]
        ov.add_patch(plt.Rectangle((xs[s], hy + mi * ch), cw, ch, facecolor=cmap(norm(r.delta_z)), edgecolor='white', lw=.4))
        if r.q < .05: ov.scatter(xs[s] + cw / 2, hy + mi * ch + ch / 2, s=1.6, color=INK, lw=0, zorder=5)
for s in strata: ov.text(xs[s] + cw / 2, hy + 12 * ch + .8, ann_label(s.split('|')[1]), fontsize=5.5, rotation=90, ha='center', va='top')
cax = cv.fig.add_axes([(hx + 2) / W, (cv.H - 62.5) / cv.H, 26 / W, 1.8 / cv.H])
cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal'); cb.set_ticks([-1, 0, 1]); cb.outline.set_linewidth(.8)
cax.tick_params(labelsize=5.5, length=2, width=.8, pad=1); cb.set_label('Signed Δz (HCC − healthy)', fontsize=5.5, labelpad=1)
ov.scatter([hx + 36], [61.6], s=1.6, color=INK, lw=0); ov.text(hx + 38, 61.6, 'BH q < 0.05', fontsize=5.5, va='center')
ov.text(hx + 62, 61.6, 'Columns: genomic annotation within each healthy-baseline class', fontsize=5.5, va='center', color=MUTED)
# ---------------- b
B0 = 73
cv.panel('b', 1, B0 - 3, 'Train and Test effects agree')
Bd = rd('Fig2b.tsv'); units = {'sOUF': (100, 'Δ sOUF (%p)'), 'BFV': (1000, 'Δ BFV ($\\times$10$^{-3}$)'), 'FPI': (1, 'Δ FPI')}
for k, nm in enumerate(['sOUF', 'BFV', 'FPI']):
    ax = cv.ax(10 + k * 37, B0 + 8, 29, 29); g = Bd[Bd.metric == nm].dropna(subset=['delta_train', 'delta_test']); f = units[nm][0]
    xx, yy = g.delta_train * f, g.delta_test * f
    ax.scatter(xx, yy, s=6, color=INK, alpha=.55, lw=0, zorder=3)
    aa = g[g.stratum == 'All|All']; ax.scatter(aa.delta_train * f, aa.delta_test * f, s=16, facecolor='none', edgecolor=HCC, lw=1.0, zorder=4)
    lo, hi = min(xx.min(), yy.min(), 0), max(xx.max(), yy.max()); pad = (hi - lo) * .06
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=PALE, lw=1.0, ls='--', zorder=1)
    ax.axhline(0, color=GRID, lw=.8, zorder=0); ax.axvline(0, color=GRID, lw=.8, zorder=0)
    ax.set_xlim(lo - pad, hi + pad); ax.set_ylim(lo - pad, hi + pad)
    ax.set_xlabel(units[nm][1].replace('Δ', 'Δ Train — ', 1) if False else 'Train ' + units[nm][1], labelpad=1); ax.set_ylabel('Test ' + units[nm][1], labelpad=1)
    npos = [(g.delta_train > 0).sum(), (g.delta_test > 0).sum()]
    ax.set_title(f'{nm}: ↑ in {npos[0]} (Train), {npos[1]} (Test) of 59', fontsize=5.6, pad=2, loc='left')
ov.scatter([12], [B0 + 46.5], s=16, facecolor='none', edgecolor=HCC, lw=1.0); ov.text(14.5, B0 + 46.5, 'All markers', fontsize=5.5, va='center')
ov.text(36, B0 + 46.5, 'Dashed line, identity; points, strata', fontsize=5.5, va='center', color=MUTED)
# ---------------- c
cv.panel('c', 124, B0 - 3, 'ΔBFV decomposition')
Cd = rd('Fig2c.tsv'); axc = cv.ax(132, B0 + 8, 44, 29)
for sp, mk in (('Train', 'o'), ('Test', 's')):
    g = Cd[Cd.split == sp]
    axc.scatter(g.term_mean_shift * 1e3, g.term_FPI_change * 1e3, s=7, marker=mk, color=SPLIT_C[sp], alpha=.75, lw=0, label=sp, zorder=3)
axc.axhline(0, color=GRID, lw=.8, zorder=0); axc.axvline(0, color=GRID, lw=.8, zorder=0)
axc.set_xlabel('Mean-shift term ($\\times$10$^{-3}$)', labelpad=1); axc.set_ylabel('FPI term ($\\times$10$^{-3}$)', labelpad=1)
axc.legend(frameon=False, fontsize=5.5, loc='upper right', handletextpad=.1, borderaxespad=0, markerscale=1.3)
med = Cd[Cd.delta_BFV > 0].groupby('split')[['frac_mean_shift', 'frac_FPI_change']].median()  # same rule as A1_decomposition_summary.tsv (strata with positive delta-BFV)
ov.text(124, B0 + 45, f"Median share of \u0394BFV (strata with \u0394BFV > 0): FPI term\n{med.frac_FPI_change.min()*100:.0f}\u2013{med.frac_FPI_change.max()*100:.0f}%, mean-shift term " + f"{med.frac_mean_shift.min()*100:.0f} to {med.frac_mean_shift.max()*100:.0f}%".replace('-', '−') + " (Train, Test)", fontsize=5.5, va='top', color=MUTED)
# ---------------- d
D0 = 127
cv.panel('d', 1, D0 - 3, 'Effect retained after adjustment')
Dd = rd('Fig2d.tsv'); axd = cv.ax(12, D0 + 8, 74, 34); mets = ['sOUF', 'BFV', 'FPI', 'FPI_K']; YMAX = 200
out = {}
for i, m in enumerate(mets):
    for sp, off, mk in (('Train', -.2, 'o'), ('Test', .2, 's')):
        v = Dd[(Dd.metric == m) & (Dd.split == sp)].retained_pct.to_numpy(); out[(m, sp)] = int(((v < 0) | (v > YMAX)).sum())
        vc = np.clip(v, 0, YMAX)
        bp = axd.boxplot(v[(v >= 0) & (v <= YMAX)], positions=[i + off], widths=.3, vert=True, showfliers=False, patch_artist=True, whis=1.5,
                         medianprops={'color': INK, 'linewidth': 1.15}, boxprops={'edgecolor': MUTED, 'linewidth': 1., 'facecolor': 'none'}, whiskerprops={'color': MUTED, 'linewidth': 1.}, capprops={'color': MUTED, 'linewidth': 1.})
        axd.scatter(i + off + rng.uniform(-.08, .08, len(v)), vc, s=3.5, marker=mk, color=SPLIT_C[sp], alpha=.5, lw=0, zorder=2)
axd.axhline(100, color=INK, lw=1.0, ls='--', zorder=1); axd.set_ylim(0, YMAX); axd.set_xlim(-.6, 3.6)
axd.set_xticks(range(4)); axd.set_xticklabels(['sOUF', 'BFV', 'FPI', 'FPI$_{K}$']); axd.set_ylabel('Effect retained (%)', labelpad=1)
h = [Line2D([], [], ls='', marker='o', ms=3, color=SPLIT_C['Train'], lw=0, label='Train'), Line2D([], [], ls='', marker='s', ms=3, color=SPLIT_C['Test'], lw=0, label='Test')]
axd.legend(handles=h, frameon=False, fontsize=5.5, loc='lower right', bbox_to_anchor=(1, 1.0), ncol=2, handletextpad=.1, columnspacing=.8, borderaxespad=0)
nout = sum(out.values())
ov.text(12, D0 + 49, f'Each point, one stratum (n = 59); coefficient of model B (depth, age, sex + mean K,\nmean fragment length) / model A (depth, age, sex) \u00d7 100. {nout} points outside\n0–{YMAX}% shown at the axis edge.', fontsize=5.5, va='top', color=MUTED)
# ---------------- e
cv.panel('e', 94, D0 - 3, 'Per-patient polarisation (all markers)')
E_ = rd('Fig2e.tsv'); cats = ['Healthy', 'HCC, tumour fraction estimated', 'HCC, not estimated']
ht = E_[(E_.category == 'Healthy') & (E_.split == 'Train')]
for k, (m, lab, f) in enumerate((('sOUF', 'sOUF (%)', 100), ('BFV', 'BFV ($\\times$10$^{-3}$)', 1000))):
    ax = cv.ax(106 + k * 37, D0 + 8, 30, 34)
    for i, c in enumerate(cats):
        v = E_[E_.category == c][m].to_numpy() * f
        ax.boxplot(v, positions=[i], widths=.5, showfliers=False, whis=1.5, medianprops={'color': INK, 'linewidth': 1.15}, boxprops={'color': MUTED, 'linewidth': 1.}, whiskerprops={'color': MUTED, 'linewidth': 1.}, capprops={'color': MUTED, 'linewidth': 1.})
        if i == 0: ax.scatter(i + rng.uniform(-.2, .2, len(v)), v, s=2.5, color=HEALTHY, marker='o', alpha=.3, lw=0, zorder=2)
        elif i == 1: ax.scatter(i + rng.uniform(-.2, .2, len(v)), v, s=6, color=HCC, marker='^', alpha=.85, lw=0, zorder=2)
        else: ax.scatter(i + rng.uniform(-.2, .2, len(v)), v, s=4, facecolor='none', edgecolor=HCC, marker='^', alpha=.7, lw=.5, zorder=2)
        if i:
            z = (np.median(E_[E_.category == c][m]) - ht[m].median()) / ht[m].std()
            ax.text(i, 1.0, f'+{z:.1f} SD', transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=5.5, fontweight='bold')
    ax.set_xticks(range(3)); ax.set_xticklabels(['Healthy', 'HCC\nwith TF', 'HCC\nno TF'], fontsize=5.5); ax.set_ylabel(lab, labelpad=1); ax.set_xlim(-.6, 2.6)
ov.text(101, D0 + 49, 'TF, tumour fraction estimated from copy-number changes.\nSD, healthy Train SD.', fontsize=5.5, va='top', color=MUTED)
cv.save('Figure_2'); print(out)
