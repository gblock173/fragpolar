"""Figure 1 (earlier design pass): fragment-state framework.
a  same mean (mu = 0.5), different fragment states: three constructed sets, WFV/BFV split and FPI (source_data/Fig1a.tsv)
b  definition cards for the 12 summaries (6 categories) + 3 interpretive indices, Low/High examples with values
   recomputed from the definitions (prep_fig1b_cards.py -> source_data/Fig1b.tsv); card design = earlier Fig. 2a
c  CpG-order reversal counterexample (source_data/Fig1c.tsv, Fig1c_positions.tsv)
d  generator response, f = 5 %, Train seed 20 (source_data/Fig1d.tsv), earlier Fig. 2c heatmap style
Former panel e (study design) moved to Supplementary Fig. 1c.
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys
sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from style_v8 import *
from matplotlib.colors import TwoSlopeNorm
from matplotlib.font_manager import FontProperties

H = 228.5
cv = Canvas(H); ov = cv.ov
R_ = cv.fig.canvas.get_renderer()
ORDER1 = ['AMF', 'HOF', 'MHL', 'UHL', 'J', 'PDR', 'WFV', 'BFV', 'FDRP', 'qFDRP', 'PM', 'ME']   # = card order of panel b
_md = {m[1]: m for m in METRICS}; METRICS1 = [_md[k] for k in ORDER1]


def tw(s, size, weight='normal'):
    w, _, _ = R_.get_text_width_height_descent(s, FontProperties(family=FONT, size=size, weight=weight), ismath=False)
    return w / cv.fig.dpi * 25.4


def wrap(s, width, size, weight='normal'):
    out, cur = [], ''
    for wd in s.split(' '):
        c = (cur + ' ' + wd).strip()
        if cur and tw(c, size, weight) > width:
            out.append(cur); cur = wd
        else:
            cur = c
    return out + ([cur] if cur else [])


# ================================================================ a
YA = 1.0
cv.panel('a', 1, YA, 'Same mean, different fragment states (μ = 0.5; 12 fragments × 4 CpGs)')
A = rd('Fig1a.tsv'); sets = list(A.set.unique()); S = A.drop_duplicates('set').set_index('set')
SETLAB = {'Opposite uniform': 'Opposite uniform', 'Within-fragment mixed': 'Within-fragment mixed', 'Combined': 'Combined'}
CW_, CB_ = CAT['Within-fragment mixing'], CAT['Between-fragment difference']
BW = 30.0                                                         # bar length for mu(1 - mu) = 0.25
for i, s in enumerate(sets):
    x0 = 6 + i * 60; yb = YA + 9.6
    g = A[A.set == s].sort_values('constructed_fragment_index'); pats = g.binary_pattern.astype(str).str.zfill(4).tolist()
    cv.text(x0, YA + 6.0, SETLAB[s], FS['label'], 'bold')
    fragments(ov, x0 + .7, yb, pats, dx=2.3, dy=1.72, rad=.62)
    # stacked horizontal bar WFV | BFV (length = share of mu(1 - mu))
    xb = x0 + 14.5; yc = yb + 8.2; hb = 4.2; w_, b_ = float(S.loc[s, 'WFV']), float(S.loc[s, 'BFV'])
    ov.add_patch(Rectangle((xb, yc), BW * w_ / .25, hb, fc=CW_, ec='none'))
    ov.add_patch(Rectangle((xb + BW * w_ / .25, yc), BW * b_ / .25, hb, fc=CB_, ec='none'))
    for v_, x_, c_ in ((w_, xb + BW * w_ / .25 / 2, 'white'), (b_, xb + BW * (w_ + b_ / 2) / .25, 'white')):
        if v_ > .02:
            ov.text(x_, yc + hb / 2, f'{v_:.3f}', fontsize=FS['small'], ha='center', va='center_baseline', color=c_)
    for t in (0, .125, .25):
        xt = xb + BW * t / .25; ov.add_line(Line2D([xt, xt], [yc + hb, yc + hb + .9], color=MUTED, lw=LW['axis']))
        ov.text(xt, yc + hb + 1.2, f'{t:g}', fontsize=FS['tick'], ha='center', va='top', color=INK2)
    ov.add_line(Line2D([xb, xb + BW], [yc + hb, yc + hb], color=MUTED, lw=LW['axis']))
    ov.text(xb, yc - 1.0, f'FPI = {S.loc[s, "FPI"]:.1f}', fontsize=FS['label'] + .3, fontweight='bold', va='bottom')
# key a
yk = YA + 32.6; xk = 6
ov.add_patch(Circle((xk + .7, yk), .62, fc=METH, ec=METH, lw=.5)); cv.text(xk + 2, yk, 'methylated CpG', va='center_baseline')
ov.add_patch(Circle((xk + 21.5, yk), .62, fc='white', ec=METH, lw=.5)); cv.text(xk + 22.8, yk, 'unmethylated CpG; one row, one fragment', va='center_baseline')
xk2 = 80.5
for k, (lab, c) in enumerate((('WFV', CW_), ('BFV', CB_))):
    ov.add_patch(Rectangle((xk2 + k * 13, yk - .9), 2.4, 1.8, fc=c, ec='none')); cv.text(xk2 + k * 13 + 3.1, yk, lab, va='center_baseline')
cv.text(xk2 + 27, yk, 'bar length: share of μ(1 − μ) = 0.25 (WFV + BFV);  FPI = BFV/μ(1 − μ)', va='center_baseline', color=INK2)

# ================================================================ b: definition cards
YB = 39.0
cv.panel('b', 1, YB, 'Twelve summaries and three fragment-state indices, with Low and High examples')
C = rd('Fig1b.tsv')
ROWS = [['AMF', 'HOF', 'MHL', 'UHL', 'J'], ['PDR', 'WFV', 'BFV', 'FDRP', 'qFDRP'], ['PM', 'ME', 'sOUF', 'FPI', 'FPI_K']]
GXC = 1.6                                                         # gap between categories (within a category: GX)
def cardx(row):
    cats_ = [info[f].category for f in row]; ng = sum(cats_[k] != cats_[k - 1] for k in range(1, len(row)))
    gc = GX + (GXC - GX) * 2 / max(ng, 1) if ng else GX                 # every row spans the full width
    xs, x = [], X0
    for k in range(len(row)):
        if k:
            x += CWD + (gc if cats_[k] != cats_[k - 1] else GX)
        xs.append(x)
    return xs
X0, WB = 1.5, 177.0; GX = 1.6; NC = 5; CWD = (WB - (NC - 1) * GX) / NC; PAD = 1.1; IW = CWD - 2 * PAD
GAPB = 2.0; BWD = (IW - GAPB) / 2
DX, DY, RAD = 2.5, 2.05, .72
info = {r.card: r for r in C.drop_duplicates('card').itertuples()}
vals = {(r.card, r.example): r for r in C.itertuples()}
HL = {'PDR': 'mixed', 'WFV': 'mixed', 'MHL': 'run1', 'UHL': 'run0', 'J': 'switch', 'HOF': 'tail', 'sOUF': 'p1',
      'FDRP': ('pair', 0, 1), 'qFDRP': ('pair', 1, 2)}
DISPNAME = {'FPI_K': 'FPI$_K$'}
y = YB + 5.4
for ri, row in enumerate(ROWS):
    # text wrapping per row
    txt = {f: (wrap(info[f].full_name, IW, FS['card_name']), wrap(info[f].definition, IW, FS['card_def'])) for f in row}
    nn = max(len(t[0]) for t in txt.values()); nd = max(len(t[1]) for t in txt.values())
    head = PAD * .8 + lh(FS['card_abbr']) + nn * lh(FS['card_name']) + .5
    blk_h = 3 * DY + 2 * RAD
    body = .5 + nd * lh(FS['card_def']) + .5 + lh(FS['small']) + .3 + blk_h + .5 + lh(FS['card_val']) + lh(FS['small']) + PAD * .6
    ch = head + body
    # category strips over contiguous groups
    cats = [info[f].category for f in row]; XS = cardx(row)
    for c, grp in itertools.groupby(enumerate(row), key=lambda t: cats[t[0]]):
        grp = list(grp); i0, i1 = grp[0][0], grp[-1][0]
        xs_ = XS[i0]; xe_ = XS[i1] + CWD
        col = CAT.get(c, CAT_INTERP_COL)
        sh = cat_strip(ov, xs_, xe_, y, c, col)
    y += sh
    for ci, f in enumerate(row):
        x = XS[ci]; col = catcol(f) if f in METRIC_CAT else CAT_INTERP_COL
        card(ov, x, y, CWD, ch, col, head, dashed=f not in METRIC_CAT)
        name, dfn = txt[f]; yy = y + PAD * .8
        if f == 'FPI_K':
            t_ = cv.text(x + PAD, yy, 'FPI', FS['card_abbr'], 'bold'); cv.text(x + PAD + tw('FPI', FS['card_abbr'], 'bold') + .15, yy + 1.15, 'K', FS['small'], 'bold')
        else:
            cv.text(x + PAD, yy, f, FS['card_abbr'], 'bold')
        yy += lh(FS['card_abbr'])
        cv.text(x + PAD, yy, '\n'.join(name), FS['card_name'], linespacing=1.17); yy = y + head + .5
        cv.text(x + PAD, yy, '\n'.join(dfn).replace('₂', '$_2$'), FS['card_def'], color=INK2, linespacing=1.17); yy += nd * lh(FS['card_def']) + .5
        for bi, ex in enumerate(('Low', 'High')):
            bx = x + PAD + bi * (BWD + GAPB); cx = bx + BWD / 2; r = vals[(f, ex)]
            reads = r.fragments.split(); span = 3 * DX
            cv.text(cx, yy, ex, FS['small'], color=INK2, ha='center')
            ry = yy + lh(FS['small']) + .3; fx = cx - span / 2 - (.45 if isinstance(HL.get(f), tuple) else 0)
            hl = HL.get(f); rows_ = None; kind = hl
            ps = [s_.count('1') / 4 for s_ in reads]
            if hl == 'mixed':
                kind, rows_ = 'rows', [k for k, s_ in enumerate(reads) if '0' in s_ and '1' in s_]
            elif hl == 'tail':
                kind, rows_ = 'rows', [k for k, p in enumerate(ps) if p >= .7]
            elif hl == 'p1':
                kind, rows_ = 'rows', [k for k, p in enumerate(ps) if p == 1]
            fragments(ov, fx, ry, reads, dx=DX, dy=DY, rad=RAD, col=col, hl=kind, hl_rows=rows_)
            vy = ry + blk_h + .5
            cv.text(cx, vy, r.value_text, FS['card_val'], 'bold', ha='center')
            cv.text(cx, vy + lh(FS['card_val']), r.derivation, FS['small'], color=INK2, ha='center')
        xm = x + PAD + BWD + GAPB / 2
        ov.add_line(Line2D([xm, xm], [yy + .3, yy + lh(FS['small']) + .3 + blk_h + .5 + lh(FS['card_val']) + lh(FS['small']) - .3], color='#d9d9d9', lw=.5))
    y += ch + (1.5 if ri < 2 else 0)
yk = y + 2.4
ov.add_patch(FancyBboxPatch((X0, yk - 1.0), 4.0, 2.0, boxstyle='round,pad=0,rounding_size=1.0', fc='#d4d4d4', ec='none'))
cv.text(X0 + 4.8, yk, 'what the summary counts;  p, methylated fraction of a fragment;  values: formula applied to the 4 fragments shown '
        '(no depth gates);  HOF, sOUF: healthy reference unmethylated', va='center_baseline', color=INK2)
YB_END = yk + 2.0
print(f'b: card width {CWD:.1f} mm; cards end {y:.1f} mm; key {yk:.1f}')

# ================================================================ c
YC = YB_END + 4.0
cv.panel('c', 1, YC, 'CpG-order reversal leaves all 12 summaries unchanged')
cv.note(6, YC + 5.2, 'Patterns 0001, 0011, 0111, 1111 × 10, 20, 30, 40 fragments;\nreversed set: each pattern read backwards')
P = rd('Fig1c_positions.tsv'); axp = cv.ax(12, YC + 15.5, 27, 33)
for o, mk, col, ls in (('Original', 'o', MUTED, '-'), ('Reversed', 's', '#E69F00', (0, (3, 1.5)))):
    g = P[P.order == o]; axp.plot(g.cpg_position, g.methylated_fraction, ls=ls, marker=mk, ms=3.2, color=col, label=o, lw=1.1, mew=0)
axp.set_xticks([1, 2, 3, 4]); axp.set_ylim(.3, 1.05); axp.set_xlim(.6, 4.4)
axp.set_xlabel('CpG position'); axp.set_ylabel('Methylated fraction')
axp.legend(fontsize=FS['small'], loc='lower center', ncol=2, handlelength=2.0, bbox_to_anchor=(.5, 1.02), columnspacing=.9)
F_ = rd('Fig1c.tsv'); F_ = F_[F_.healthy_reference_direction == 0].reset_index(drop=True)
axf = cv.ax(56, YC + 15.5, 26, 33)
for i, m in enumerate(METRICS1):
    r = F_[F_.feature == {'WFV': 'U', 'BFV': 'B', 'HOF': 'A'}.get(m[1], m[1])].iloc[0]; col = CAT[m[2]]
    axf.scatter(r.original, i, s=13, color=col, zorder=3, lw=0)
    axf.scatter(r.reversed_order, i, s=32, facecolor='none', edgecolor=INK, lw=.8, zorder=4)
axf.set_ylim(11.6, -.6); axf.set_yticks(range(12)); axf.set_yticklabels([m[1] for m in METRICS1]); axf.set_xlim(-.04, .86)
axf.set_xticks([0, .25, .5, .75]); axf.set_xticklabels(['0', '0.25', '0.5', '0.75'])
for t, m in zip(axf.get_yticklabels(), METRICS1):
    t.set_color(INK)
axf.set_xlabel('Summary value'); axf.tick_params(axis='y', length=0, pad=6.5)
for i, m in enumerate(METRICS1):
    ov.add_patch(Rectangle((56 - 1.9, YC + 15.5 + (i + .6) * 33 / 12.2 - 1.05), 1.1, 2.1, fc=CAT[m[2]], ec='none', zorder=6))
axf.xaxis.grid(True, color=GRID, lw=.6); axf.set_axisbelow(True)
axf.legend(handles=[Line2D([], [], ls='', marker='o', ms=3.3, color=MUTED, mew=0, label='original'),
                    Line2D([], [], ls='', marker='o', ms=5.0, mfc='none', mec=INK, mew=.8, label='reversed')],
           loc='lower center', bbox_to_anchor=(.42, 1.02), ncol=2, fontsize=FS['small'], handletextpad=.2, columnspacing=.9)
cv.note(46, YC + 55.4, 'max |difference| = 0 (both healthy directions;\nHOF shown for healthy μ < 0.5)')

# ================================================================ d
XD = 88.0
cv.panel('d', XD, YC, 'Generator responses in healthy cfDNA (f = 5%)')
D = rd('Fig1d.tsv'); strata = ['All|All', 'U|All', 'mid|All', 'M|All']; gens = ['G1', 'G2', 'G3']
cv.note(XD + 4.6, YC + 5.2, 'G1, healthy-opposite uniform fragments added; G2, within-fragment\ndiscordant fragments added; G3, independent CpG flips')
hx = XD + 13.0; gapg = 1.4; cw = (178.5 - hx - 3 * gapg) / 12; chh = 3.15
hy = YC + 21.0
cmap = plt.get_cmap('RdBu_r'); norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
OUT = {('G1', 'PDR'), ('G1', 'WFV'), ('G2', 'HOF')}           # component not added by the generator (dilution only)
idx = D.set_index(['stratum', 'generator', 'metric'])
for si, s in enumerate(strata):
    gx = hx + si * (3 * cw + gapg)
    cv.text(gx + 1.5 * cw, hy - 6.0, BASE_HEAD[s.split('|')[0]], FS['small'], 'bold', ha='center')
    ov.add_line(Line2D([gx + .3, gx + 3 * cw - .3], [hy - 3.3, hy - 3.3], color=MUTED, lw=LW['axis']))
    for gi, g in enumerate(gens):
        cv.text(gx + (gi + .5) * cw, hy - 2.6, g, FS['small'], ha='center')
        for mi, m in enumerate(METRICS1):
            r = idx.loc[(s, g, m[1])]; v = float(r.scaled); z = float(r.mean_delta_z)
            xx, yy = gx + gi * cw, hy + mi * chh
            ov.add_patch(Rectangle((xx + .12, yy + .12), cw - .24, chh - .24, fc=cmap(norm(v)), ec='none'))
            zt = ('−' if z <= -.05 else '') + (f'{abs(z):.1f}' if abs(z) < 9.95 else f'{abs(z):.0f}')
            cv.text(xx + cw / 2, yy + chh / 2, zt, FS['small'], ha='center', va='center_baseline',
                    color='white' if abs(v) > .62 else INK)
            if (g, m[1]) in OUT:
                ov.add_patch(Rectangle((xx + .12, yy + .12), cw - .24, chh - .24, fc='none', ec=INK, lw=1.0, zorder=8))
for mi, m in enumerate(METRICS1):
    ov.add_patch(Rectangle((hx - 2.0, hy + mi * chh + .15), 1.2, chh - .3, fc=CAT[m[2]], ec='none'))
    cv.text(hx - 2.7, hy + mi * chh + chh / 2, m[1], FS['tick'] + .3, ha='right', va='center_baseline')
yend = hy + 12 * chh
hcolorbar(cv, cmap, norm, hx, yend + 4.6, 26, 'Colour: mean Δz scaled per summary', [-1, 0, 1], ['−1', '0', '1'])
xr = hx + 38.5
cv.note(xr, yend + 1.6, 'Numbers: mean Δz (n = 20 Train healthy)')
ov.add_patch(Rectangle((xr, yend + 5.2), 3.0, 2.0, fc='none', ec=INK, lw=1.0))
cv.note(xr + 4.0, yend + 4.9, 'component not added by the generator\n(change by dilution only)')
print('d end', yend + 10)
q = cv.save(FIG / 'Figure_1')
