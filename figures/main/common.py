"""Shared style for earlier Figures 1-3 (Nature Communications, 180 mm wide)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
AN = Path(CFG.WORK); FIG = Path(CFG.FIG); SD = FIG / 'source_data'
SD.mkdir(parents=True, exist_ok=True)
W = 180.; MM = 1 / 25.4
HEALTHY = '#2C7DA0'; HCC = '#C76B43'
INK = '#1a1a1a'; MUTED = '#4a5568'; PALE = '#a0aec0'; GRID = '#cbd5e0'
CAT = {'Average level': '#4a3aa7', 'Within-fragment mixing': '#1baf7a', 'Continuity and arrangement': '#aa4499',
       'Reference deviation': '#d9ab00', 'Pattern diversity': '#e377c2', 'Between-fragment difference': '#008300'}
# 12 summaries: (internal name used in generator step tables, display name, category)
METRICS = [('amf', 'AMF', 'Average level'), ('pdr', 'PDR', 'Within-fragment mixing'), ('within_read_U', 'WFV', 'Within-fragment mixing'),
           ('mhl', 'MHL', 'Continuity and arrangement'), ('uhl', 'UHL', 'Continuity and arrangement'),
           ('arrangement_excess_per_link', 'J', 'Continuity and arrangement'), ('hof', 'HOF', 'Reference deviation'),
           ('pm', 'PM', 'Pattern diversity'), ('me', 'ME', 'Pattern diversity'),
           ('between_read_B_emp', 'BFV', 'Between-fragment difference'), ('fdrp', 'FDRP', 'Between-fragment difference'),
           ('qfdrp', 'qFDRP', 'Between-fragment difference')]
SPLIT_C = {'Train': '#4a5568', 'Test': '#E69F00'}
def setup():
    plt.rcParams.update({'font.family': 'Liberation Sans', 'font.size': 6, 'axes.titlesize': 6.5, 'axes.labelsize': 6.0,
        'xtick.labelsize': 5.6, 'ytick.labelsize': 5.6, 'axes.linewidth': 1.0, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
        'text.color': INK, 'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.spines.top': False, 'axes.spines.right': False,
        'xtick.major.width': 1.0, 'ytick.major.width': 1.0, 'xtick.major.size': 2., 'ytick.major.size': 2., 'lines.linewidth': 1.0,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none', 'savefig.facecolor': 'white', 'axes.unicode_minus': True})
class Canvas:
    def __init__(self, H):
        setup(); self.H = H
        self.fig = plt.figure(figsize=(W * MM, H * MM), facecolor='white')
    def ax(self, x, top, w, h):
        a = self.fig.add_axes([x / W, (self.H - top - h) / self.H, w / W, h / self.H]); style(a); return a
    def text(self, x, top, s, size=5.8, weight='normal', **kw):
        kw.setdefault('va', 'top'); return self.fig.text(x / W, 1 - top / self.H, s, fontsize=size, fontweight=weight, **kw)
    def panel(self, letter, x, top): self.text(x, top, letter, 8.5, 'bold')
    def save(self, name):
        for ext in ('pdf', 'svg'): self.fig.savefig(FIG / f'{name}.{ext}')
        self.fig.savefig(FIG / f'{name}.png', dpi=300)
def style(ax):
    ax.minorticks_off(); ax.tick_params(direction='out', length=2, width=1, pad=1.5)
    for sp in ax.spines.values(): sp.set_linewidth(1); sp.set_color(MUTED)
def bare(ax):
    for s in ax.spines.values(): s.set_visible(False)
    ax.set_xticks([]); ax.set_yticks([])
def rd(name): return pd.read_csv(SD / name, sep='\t')
def strata_order():
    base = ['All', 'U', 'mid', 'M']
    ann = ['All', 'Island', 'Shore', 'Enhancer', 'TSS', 'ReprPC', 'promoter', 'exon', 'intron', '5UTR', '3UTR', 'intergenic', 'Liver-Hep', 'Immune5', 'Eryth-prog']
    return [f'{b}|{a}' for b in base for a in ann]
def overlay(cv):
    """Full-canvas axes in mm (y downward) for schematics and text."""
    a = cv.fig.add_axes([0, 0, 1, 1]); a.set_xlim(0, W); a.set_ylim(cv.H, 0); a.axis('off'); a.patch.set_alpha(0); return a
def rbox(a, x, y, w, h, ec=MUTED, fc='none', lw=1.0, r=1.2, ls='-'):
    from matplotlib.patches import FancyBboxPatch
    p = FancyBboxPatch((x + r, y + r), w - 2 * r, h - 2 * r, boxstyle=f'round,pad={r},rounding_size={r}', ec=ec, fc=fc, lw=lw, ls=ls); a.add_patch(p); return p
def arrow(a, x0, y0, x1, y1, color=MUTED, lw=1.0):
    a.annotate('', xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle='-|>', color=color, lw=lw, shrinkA=0, shrinkB=0, mutation_scale=6))
