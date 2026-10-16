"""Shared helpers for Supplementary Figures (earlier). Reuses ../main/style_v8.py (180 mm, sans-serif)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys
sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / 'main'))
from style_v8 import *
import matplotlib as mpl
from matplotlib.colors import TwoSlopeNorm
DISP = {m[0]: m[1] for m in METRICS}
class SC(Canvas):
    def __init__(self, H, n):
        super().__init__(H); self.n = n; self.ov = overlay(self)
    def save(self):
        return Canvas._save_path(self, SUPP / f'Supplementary_Figure_S{self.n}', svg=False)
    def title(self, x, top, s): self.text(x, top, s, FS['title'], 'bold')
def wt(df, name): df.to_csv(SSD / name, sep='\t', index=False, float_format='%.6g'); return df
def rt(name): return pd.read_csv(SSD / name, sep='\t')
def heat(cv, df, x0, y0, cw=2.5, ch=3.0, gap=1.6, vmax=1.0, label=True, vcol='delta_z', qcol='q', metrics=None, cmap='RdBu_r', catbars=True):
    """df: stratum, metric(display), value, q. Returns (cmap, norm)."""
    ov = cv.ov; strata = [s for s in strata_order() if s in set(df.stratum)]
    mets = metrics or [m[1] for m in METRICS]; cat = {m[1]: m[2] for m in METRICS}
    cm = plt.get_cmap(cmap); norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    x = x0; xs = {}; last = None
    for s in strata:
        b = s.split('|')[0]
        if last is not None and b != last: x += gap
        last = b; xs[s] = x; x += cw
    heads = {'All': 'All markers', 'U': 'U baseline', 'mid': 'Mid baseline', 'M': 'M baseline'}
    if label:
        for b, t in heads.items():
            ss = [s for s in strata if s.startswith(b + '|')]
            if not ss: continue
            a, c = xs[ss[0]], xs[ss[-1]] + cw
            ov.text((a + c) / 2, y0 - 2.2, t, fontsize=FS['title'] - .5, ha='center', va='bottom', fontweight='bold'); ov.plot([a, c], [y0 - 1.5, y0 - 1.5], color=MUTED, lw=1.0)
    idx = df.set_index(['stratum', 'metric'])
    for mi, m in enumerate(mets):
        ov.text(x0 - 3.5, y0 + mi * ch + ch / 2, m, fontsize=5.5, ha='right', va='center')
        if catbars and m in cat: ov.add_patch(plt.Rectangle((x0 - 2.7, y0 + mi * ch + .2), 1.2, ch - .4, color=CAT[cat[m]], lw=0))
        for s in strata:
            if (s, m) not in idx.index: continue
            r = idx.loc[(s, m)]
            ov.add_patch(plt.Rectangle((xs[s], y0 + mi * ch), cw, ch, facecolor=cm(norm(r[vcol])), edgecolor='white', lw=.4))
            if qcol and r[qcol] < .05: ov.scatter(xs[s] + cw / 2, y0 + mi * ch + ch / 2, s=1.4, color=INK, lw=0, zorder=5)
    if label:
        for s in strata: ov.text(xs[s] + cw / 2, y0 + len(mets) * ch + .8, ann_label(s.split('|')[1]), fontsize=5.5, rotation=90, ha='center', va='top')
    return cm, norm, xs
def cbar(cv, cm, norm, x, top, w, label, ticks):
    cax = cv.fig.add_axes([x / W, (cv.H - top - 1.8) / cv.H, w / W, 1.8 / cv.H])
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cm, norm=norm, orientation='horizontal'); cb.set_ticks(ticks); cb.outline.set_linewidth(.8)
    cax.tick_params(labelsize=5.5, length=2, width=.8, pad=1); cb.set_label(label, fontsize=5.5, labelpad=1)
def box(ax, data, pos, colors, w=.5, jitter=None, rng=None, s=4, alpha=.5):
    ax.boxplot(data, positions=pos, widths=w, patch_artist=True, showfliers=False, whis=1.5, medianprops={'color': INK, 'linewidth': 1.1},
               boxprops={'edgecolor': MUTED, 'linewidth': 1., 'facecolor': 'none'}, whiskerprops={'color': MUTED, 'linewidth': 1.}, capprops={'color': MUTED, 'linewidth': 1.})
    if rng is not None:
        for v, x, c in zip(data, pos, colors): ax.scatter(x + rng.uniform(-.17, .17, len(v)), v, s=s, color=c, alpha=alpha, lw=0, zorder=2)
def forest(ax, labels, mean, lo, hi, colors, y0=0, sz=14, marker='o'):
    for i, (m, l, h, c) in enumerate(zip(mean, lo, hi, colors)):
        ax.plot([l, h], [y0 + i, y0 + i], color=c, lw=1.0); ax.scatter([m], [y0 + i], s=sz, color=c, marker=marker, zorder=3, lw=0)
    ax.set_yticks(range(y0, y0 + len(labels))); ax.set_yticklabels(labels); ax.invert_yaxis(); ax.axvline(0, color=PALE, lw=1.0, ls='--', zorder=0)
