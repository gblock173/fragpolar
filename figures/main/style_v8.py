"""style_v8 — common design module for earlier main and supplementary figures (Nature Communications, 180 mm).

Design reference: v4_natcomm_KO Figure 2 (ms3_fig2.py + ms_feature_gallery.py): definition cards with category-coloured
border and tinted header, category strip with direct label above contiguous cards, Low/High 4-fragment examples with
computed values in bold and the derivation in grey, heatmap cells with printed values and outlined cells of interest.

Rules (NATURE_COMMUNICATIONS_FIGURE_GUIDE.md, FIGURE_DESIGN_GUIDELINES.md):
  - width 180 mm, sans-serif (Liberation Sans = Arial metrics), text 5.5-8.5 pt, panel letters 8.5 pt bold
  - axis / whisker / reference lines 1.0 pt, median / emphasis 1.15 pt; card border 1.0 pt
  - PDF/SVG vector, editable text (pdf.fonttype 42, svg.fonttype none); PNG 300 dpi for preview
  - category colours = earlier Fig. 2a CATCOL (CVD-checked, ms3_categories.py); group colours fixed below
Coordinates: mm from the top-left corner of the canvas (y downward) for every helper.
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
from pathlib import Path
import itertools
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Circle, FancyBboxPatch, Patch

AN = Path(CFG.WORK); FIG = Path(CFG.FIG); SD = FIG / 'source_data'; SUPP = FIG / 'supp'; SSD = SUPP / 'source_data'
W = 180.0; MM = 1 / 25.4; PT2MM = 25.4 / 72

# ---------------------------------------------------------------- typography and lines
FONT = 'Liberation Sans'
FS = dict(letter=8.5, title=6.5, label=6.0, tick=5.5, text=5.5, small=5.5, card_abbr=7.0, card_name=5.8, card_def=5.5, card_val=6.3)
LW = dict(axis=1.0, ref=1.0, median=1.15, card=1.0, cell=.5)

# ---------------------------------------------------------------- colours
INK = '#1a1a1a'; INK2 = '#555555'; MUTED = '#4a5568'; PALE = '#a0aec0'; GRID = '#cbd5e0'; METH = '#262626'
HEALTHY = '#2C7DA0'; HCC = '#C76B43'; HEP = '#7B5EA7'
TUMOUR = '#9b2c2c'; NONTUM = '#56B4E9'; PBMC = '#718096'
GROUP_COL = {'Healthy': HEALTHY, 'HCC': HCC, 'Hepatitis': HEP}
GROUP_MK = {'Healthy': 'o', 'HCC': '^', 'Hepatitis': 'D'}
SPLIT_C = {'Train': '#4a5568', 'Test': '#E69F00'}; SPLIT_MK = {'Train': 'o', 'Test': 's'}
# six feature categories: colours identical to earlier Fig. 2a (ms3_categories.CATCOL); names = earlier terminology
CAT_ORDER = ['Average level', 'Within-fragment mixing', 'Continuity and arrangement', 'Reference deviation', 'Pattern diversity',
             'Between-fragment difference']
CAT = {'Average level': '#4a3aa7', 'Within-fragment mixing': '#1baf7a', 'Continuity and arrangement': '#aa4499',
       'Reference deviation': '#d9ab00', 'Pattern diversity': '#e377c2', 'Between-fragment difference': '#008300'}
INTERP = 'Fragment-state indices'; CAT_INTERP_COL = '#3d4852'
# 12 summaries: (internal name, display name, category)
METRICS = [('amf', 'AMF', 'Average level'), ('hof', 'HOF', 'Reference deviation'), ('mhl', 'MHL', 'Continuity and arrangement'),
           ('uhl', 'UHL', 'Continuity and arrangement'), ('arrangement_excess_per_link', 'J', 'Continuity and arrangement'),
           ('pdr', 'PDR', 'Within-fragment mixing'), ('within_read_U', 'WFV', 'Within-fragment mixing'),
           ('between_read_B_emp', 'BFV', 'Between-fragment difference'), ('fdrp', 'FDRP', 'Between-fragment difference'),
           ('qfdrp', 'qFDRP', 'Between-fragment difference'), ('pm', 'PM', 'Pattern diversity'), ('me', 'ME', 'Pattern diversity')]
# order = Figure 1 (card order of panel b); used by every 12-summary figure
METRIC_CAT = {m[1]: m[2] for m in METRICS}
# fragment-state colours (sOUF / strict-mixed / healthy-direction uniform), shared by Fig. 3-5
STATE = {'sOUF': '#882255', 'strict-mixed': '#44AA99', 'healthy-direction uniform': '#DDCC77'}
DIV = {'up_ci': '#B2182B', 'up': '#E8A09A', 'dn': '#A6CBE3', 'dn_ci': '#2166AC'}


def tint(hexcol, a=.18):
    h = hexcol.lstrip('#'); rgb = [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    return tuple(1 - a * (1 - c) for c in rgb)


def catcol(metric):
    return CAT[METRIC_CAT[metric]] if metric in METRIC_CAT else CAT_INTERP_COL


def setup():
    plt.rcParams.update({'font.family': FONT, 'font.size': FS['text'], 'axes.titlesize': FS['title'], 'axes.labelsize': FS['label'],
                         'xtick.labelsize': FS['tick'], 'ytick.labelsize': FS['tick'], 'legend.fontsize': FS['text'],
                         'axes.linewidth': LW['axis'], 'axes.edgecolor': MUTED, 'axes.labelcolor': INK, 'text.color': INK,
                         'xtick.color': MUTED, 'ytick.color': MUTED, 'xtick.labelcolor': INK2, 'ytick.labelcolor': INK2,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'xtick.major.width': LW['axis'], 'ytick.major.width': LW['axis'], 'xtick.major.size': 2., 'ytick.major.size': 2.,
                         'xtick.major.pad': 1.5, 'ytick.major.pad': 1.5, 'axes.labelpad': 1.5,
                         'lines.linewidth': 1.0, 'legend.frameon': False, 'legend.handletextpad': .4, 'legend.borderaxespad': 0,
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white',
                         'axes.unicode_minus': True, 'mathtext.default': 'regular'})


setup()


# ---------------------------------------------------------------- canvas
class Canvas:
    """180-mm canvas with mm coordinates (top-left origin) and a full-size overlay axes `ov` for schematics."""

    def __init__(self, H, w=W):
        setup(); self.H = H; self.W = w
        self.fig = plt.figure(figsize=(w * MM, H * MM), facecolor='white')
        self.ov = self.fig.add_axes([0, 0, 1, 1]); self.ov.set_xlim(0, w); self.ov.set_ylim(H, 0); self.ov.axis('off')
        self.ov.patch.set_alpha(0); self.ov.set_zorder(5)

    def ax(self, x, top, w, h, zorder=1):
        a = self.fig.add_axes([x / self.W, (self.H - top - h) / self.H, w / self.W, h / self.H]); style(a); a.set_zorder(zorder); return a

    def text(self, x, top, s, size=FS['text'], weight='normal', **kw):
        kw.setdefault('va', 'top'); kw.setdefault('color', INK)
        return self.ov.text(x, top, s, fontsize=size, fontweight=weight, **kw)

    def note(self, x, top, s, **kw):
        kw.setdefault('color', INK2); kw.setdefault('linespacing', 1.25); return self.text(x, top, s, FS['small'], **kw)

    def panel(self, letter, x, top, title=None):
        """Panel letter (8.5 pt bold) and an optional short title on the same top line."""
        self.text(x, top - .4, letter, FS['letter'], 'bold')
        if title:
            self.text(x + 4.6, top + .25, title, FS['title'], 'bold')

    def save(self, path_stem, svg=True):
        p = Path(path_stem); p.parent.mkdir(parents=True, exist_ok=True)
        self.fig.savefig(str(p) + '.pdf')
        if svg:
            self.fig.savefig(str(p) + '.svg')
        self.fig.savefig(str(p) + '.png', dpi=300)
        q = qa(self.fig)
        print(f'saved {p.name}: {self.W:.0f} x {self.H:.0f} mm; QA {q}')
        plt.close(self.fig)
        return q


def style(ax):
    ax.minorticks_off(); ax.tick_params(direction='out', length=2, width=LW['axis'], pad=1.5)
    for sp in ax.spines.values():
        sp.set_linewidth(LW['axis']); sp.set_color(MUTED)


def bare(ax):
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])


def ax_title(ax, s, pad=3):
    ax.set_title(s, fontsize=FS['title'] - .3, fontweight='bold', pad=pad, loc='left')


def refline(ax, v, axis='h', **kw):
    kw = {**dict(color=MUTED, lw=LW['ref'], ls=(0, (3, 2)), zorder=1), **kw}
    return ax.axhline(v, **kw) if axis == 'h' else ax.axvline(v, **kw)


def zeroline(ax, axis='h'):
    return refline(ax, 0, axis, color=GRID, ls='-', lw=.8, zorder=0)


def boxstrip(ax, vals, pos, color, marker='o', rng=None, w=.32, s=4, alpha=.5, jit=None, face=True):
    v = np.asarray(vals, float); v = v[np.isfinite(v)]
    ax.boxplot(v, positions=[pos], widths=w, showfliers=False, whis=1.5, patch_artist=True,
               medianprops={'color': INK, 'linewidth': LW['median']},
               boxprops={'edgecolor': MUTED, 'linewidth': 1., 'facecolor': 'none'},
               whiskerprops={'color': MUTED, 'linewidth': 1.}, capprops={'color': MUTED, 'linewidth': 1.})
    if rng is not None:
        j = jit if jit is not None else w * .3
        x = pos + rng.uniform(-j, j, len(v))
        if face:
            ax.scatter(x, v, s=s, marker=marker, color=color, alpha=alpha, lw=0, zorder=2)
        else:
            ax.scatter(x, v, s=s, marker=marker, facecolor='none', edgecolor=color, alpha=alpha, lw=.5, zorder=2)


def ci_row(ax, est, lo, hi, y, color, marker='o', ms=3.4, lw=1.2):
    ax.plot([lo, hi], [y, y], color=color, lw=lw, solid_capstyle='butt', zorder=2)
    ax.plot(est, y, marker, color=color, ms=ms, mew=0, zorder=3)


def hcolorbar(cv, cmap, norm, x, top, w, label, ticks, ticklabels=None, h=1.6):
    cax = cv.fig.add_axes([x / cv.W, (cv.H - top - h) / cv.H, w / cv.W, h / cv.H])
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal'); cb.set_ticks(ticks)
    if ticklabels:
        cb.ax.set_xticklabels(ticklabels)
    cb.outline.set_linewidth(.8); cax.tick_params(labelsize=FS['tick'], length=1.8, width=.8, pad=1)
    cv.text(x, top - .5, label, FS['small'], va='bottom', color=INK2)
    return cb


def legend_handles(items):
    """items: (label, dict(marker=, color=, mfc=, ls=)) -> Line2D handles."""
    out = []
    for lab, kw in items:
        kw = {**dict(marker='o', ls='', ms=3.4, mew=0, color=MUTED), **kw}; out.append(Line2D([], [], label=lab, **kw))
    return out


# ---------------------------------------------------------------- cards (earlier Fig. 2a style)
def lh(size):
    return size * 1.17 * PT2MM


def card(ov, x, y, w, h, col, head_h, dashed=False):
    ov.add_patch(Rectangle((x, y), w, head_h, fc=tint(col, .16), ec='none', zorder=0))
    ov.add_patch(Rectangle((x, y), w, h, fc='none', ec=col, lw=LW['card'], zorder=5, ls=(0, (3, 1.6)) if dashed else '-'))


def cat_strip(ov, x0, x1, y, label, col, size=6.2):
    ov.text(x0 + .2, y, label, fontsize=size, fontweight='bold', va='top', color=INK)
    ov.add_patch(Rectangle((x0, y + lh(size) + .25), x1 - x0, .75, fc=col, ec='none', zorder=2))
    return lh(size) + 1.3


def fragments(ov, x, y, reads, dx=2.6, dy=2.2, rad=.78, col=None, hl=None, hl_rows=None, line_col='#8c8c8c'):
    """Block of fragments (rows) x CpGs as lollipops, top-left (x, y) mm. hl: 'rows' (tint whole rows in hl_rows),
    'run1'/'run0' (tint runs >= 2), 'switch' (bars between switching CpGs), ('pair', a, b) bracket."""
    n = len(reads[0]); span = (n - 1) * dx; tint_c = tint(col or MUTED, .42)
    for i, r in enumerate(reads):
        yc = y + rad + i * dy
        if hl == 'rows' and hl_rows is not None and i in hl_rows:
            ov.add_patch(FancyBboxPatch((x - rad - .35, yc - rad - .28), span + 2 * rad + .7, 2 * rad + .56,
                                        boxstyle='round,pad=0,rounding_size=%.2f' % (rad + .28), fc=tint_c, ec='none', zorder=1))
        if hl in ('run1', 'run0'):
            st = '1' if hl == 'run1' else '0'; j = 0
            while j < n:
                if r[j] == st:
                    k = j
                    while k + 1 < n and r[k + 1] == st:
                        k += 1
                    if k > j:
                        ov.add_patch(FancyBboxPatch((x + j * dx - rad - .3, yc - rad - .28), (k - j) * dx + 2 * rad + .6, 2 * rad + .56,
                                                    boxstyle='round,pad=0,rounding_size=%.2f' % (rad + .28), fc=tint_c, ec='none', zorder=1))
                    j = k + 1
                else:
                    j += 1
        ov.add_line(Line2D([x - rad - .5, x + span + rad + .5], [yc, yc], color=line_col, lw=.6, zorder=2, solid_capstyle='butt'))
        for j, c in enumerate(r):
            if c in '01':
                ov.add_patch(Circle((x + j * dx, yc), rad, fc=METH if c == '1' else 'white', ec=METH, lw=.5, zorder=3))
        if hl == 'switch':
            for j in range(n - 1):
                if r[j] != r[j + 1]:
                    xm = x + (j + .5) * dx
                    ov.add_line(Line2D([xm, xm], [yc - rad - .25, yc + rad + .25], color=col, lw=1.0, zorder=4, solid_capstyle='butt'))
    if isinstance(hl, tuple) and hl[0] == 'pair':
        a, b = hl[1], hl[2]; xb = x + span + rad + .9; ya = y + rad + a * dy; yb = y + rad + b * dy
        ov.add_line(Line2D([xb - .45, xb, xb, xb - .45], [ya, ya, yb, yb], color=col, lw=1.0, zorder=4, solid_joinstyle='miter'))
    return span + 2 * rad, (len(reads) - 1) * dy + 2 * rad


# ---------------------------------------------------------------- QA
def qa(fig):
    """Text QA on the rendered figure: min font size, text-text overlaps (> 0.4 px), text outside the canvas."""
    fig.canvas.draw(); r = fig.canvas.get_renderer(); fb = fig.bbox; items = []
    for ax in fig.axes:
        ts = list(ax.texts)
        if ax.axison:
            x0_, x1_ = sorted(ax.get_xlim()); y0_, y1_ = sorted(ax.get_ylim()); e_ = 1e-9
            ts += [t for t, v in zip(ax.get_xticklabels(), ax.get_xticks()) if t.get_visible() and x0_ - e_ <= v <= x1_ + e_]
            ts += [t for t, v in zip(ax.get_yticklabels(), ax.get_yticks()) if t.get_visible() and y0_ - e_ <= v <= y1_ + e_]
            ts += [ax.xaxis.label, ax.yaxis.label, ax.title]
        lg = ax.get_legend()
        if lg is not None:
            ts += lg.get_texts()
        for t in ts:
            if t.get_text().strip() and t.get_visible():
                items.append((t, t.get_window_extent(r)))
    for t in fig.texts:
        if t.get_text().strip():
            items.append((t, t.get_window_extent(r)))
    ov = []
    for (t1, b1), (t2, b2) in itertools.combinations(items, 2):
        if min(b1.x1, b2.x1) - max(b1.x0, b2.x0) > .4 and min(b1.y1, b2.y1) - max(b1.y0, b2.y0) > .4:
            ov.append((t1.get_text()[:24], t2.get_text()[:24]))
    out = [t.get_text()[:24] for t, b in items if b.x0 < fb.x0 - .5 or b.x1 > fb.x1 + .5 or b.y0 < fb.y0 - .5 or b.y1 > fb.y1 + .5]
    sizes = [t.get_fontsize() for t, _ in items]
    return dict(n_text=len(items), min_pt=round(min(sizes), 2) if sizes else None, overlaps=ov[:12], n_overlaps=len(ov), outside=out[:8])


def rd(name, supp=False):
    return pd.read_csv((SSD if supp else SD) / name, sep='\t')


def strata_order():
    base = ['All', 'U', 'mid', 'M']
    ann = ['All', 'Island', 'Shore', 'Enhancer', 'TSS', 'ReprPC', 'promoter', 'exon', 'intron', '5UTR', '3UTR', 'intergenic', 'Liver-Hep', 'Immune5',
           'Eryth-prog']
    return [f'{b}|{a}' for b in base for a in ann]


BASE_HEAD = {'All': 'All markers', 'U': 'U baseline', 'mid': 'Mid baseline', 'M': 'M baseline'}
# display names of genomic annotations (internal stratum keys such as 'U|Island' are unchanged)
ANN_LABEL = {'Island': 'CpG island', 'nonIsland': 'non-CpG island'}
ANN_LONG = {'Island': 'CpG island', 'nonIsland': 'non-CpG island', 'Shore': 'CpG shore', 'ReprPC': 'Polycomb-repressed', 'TSS': 'TSS',
            'Enhancer': 'Enhancer', 'promoter': 'Promoter', 'intergenic': 'Intergenic', 'exon': 'Exon', 'intron': 'Intron', 'All': 'All'}


def ann_label(a, long=False):
    """Display label of an annotation key ('Island' -> 'CpG island')."""
    return (ANN_LONG if long else ANN_LABEL).get(a, a)


def stratum_label(s, long=False):
    """'U|Island' -> 'U|CpG island' (display only)."""
    b, a = s.split('|'); return f'{b}|{ann_label(a, long)}'


# ---------------------------------------------------------------- compatibility with the earlier earlier figure scripts
def _save_compat(self, name):
    """Canvas.save accepts a full path stem, a figure number (int -> Figure_<n>) or a bare name (-> figures/<name>)."""
    if isinstance(name, int):
        name = f'Figure_{name}'
    p = Path(name)
    return Canvas._save_path(self, p if p.is_absolute() else FIG / p)


Canvas._save_path = Canvas.save
Canvas.save = _save_compat
Canvas.label = lambda self, l, x, top, title=None: self.panel(l, x, top, title)


def overlay(cv):
    return cv.ov


Fig = Canvas
title = ax_title
MK = GROUP_MK; COL = GROUP_COL


def boxjitter(ax, groups, positions, colors, markers, rng, s=3.5, alpha=.4, widths=.5):
    vals = [np.asarray(v, float) for v in groups]
    ax.boxplot(vals, positions=positions, widths=widths, vert=True, patch_artist=True, showfliers=False, whis=1.5,
               medianprops={'color': INK, 'linewidth': LW['median']}, boxprops={'edgecolor': MUTED, 'linewidth': 1., 'facecolor': 'none'},
               whiskerprops={'color': MUTED, 'linewidth': 1.}, capprops={'color': MUTED, 'linewidth': 1.})
    for v, x, c, m in zip(vals, positions, colors, markers):
        j = rng.uniform(-.2, .2, len(v))
        ax.scatter(x + j, v, s=s if len(v) > 50 else s * 2.2, color=c, marker=m, alpha=alpha if len(v) > 50 else .85, linewidths=0, zorder=2)


def strata_x(strata, bases=('All', 'U', 'mid', 'M'), gap=1.5, start=0.0):
    """x positions (data units) of strata grouped by healthy baseline in the order `bases`, annotation order of strata_order().
    Returns (ordered strata, dict stratum -> x, dict base -> (x_first, x_last))."""
    order = strata_order(); out = []; xs = {}; span = {}; x = start
    for b in bases:
        ss = [s for s in order if s.startswith(b + '|') and s in set(strata)]
        if not ss: continue
        if out: x += gap
        x0 = x
        for s in ss: xs[s] = x; out.append(s); x += 1
        span[b] = (x0, x - 1)
    return out, xs, span
