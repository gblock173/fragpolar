"""Supplementary Fig. 1 (earlier design pass): sample flow (a), comparison samples (b) and study design (c; former main Fig. 1e).
Source data: supp/source_data/S1.tsv (a, b; unchanged) and S1c.tsv (= former source_data/Fig1e.tsv)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys
sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / 'main'))
from style_v8 import *

S = rd('S1.tsv', supp=True); E = rd('S1c.tsv', supp=True).set_index('group')
cv = Canvas(150); ov = cv.ov


def rbox(x, y, w, h, ec=MUTED, lw=1.0, ls='-', fc='white'):
    ov.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=1.0', ec=ec, fc=fc, lw=lw, ls=ls, zorder=2))


def arrow(x0, y0, x1, y1, color=MUTED):
    ov.annotate('', xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle='-|>', color=color, lw=1.0, shrinkA=0, shrinkB=0, mutation_scale=6))


def box_(x, y, w, h, title, lines, ec=MUTED, tc=INK, fc='white'):
    rbox(x, y, w, h, ec=ec, fc=fc)
    ov.add_patch(Rectangle((x + .5, y + .5), w - 1, 4.4, fc=tint(ec, .14), ec='none', zorder=2))
    cv.text(x + w / 2, y + 1.3, title, FS['label'], 'bold', ha='center', color=tc, zorder=3)
    for i, (t, c) in enumerate(lines):
        cv.text(x + w / 2, y + 6.4 + i * 3.3, t, FS['text'] + .2, ha='center', color=c, zorder=3)


def n_(node, grp):
    r = S[(S.node == node) & (S.group.str.startswith(grp))].iloc[0]; return f'{int(r.n):,}'


def e_(g, c='n_total'):
    return f'{int(E.loc[g, c]):,}'


# ---------------- a
cv.panel('a', 1, 1, 'Plasma cfDNA cohort')
box_(3, 9, 40, 17, 'Clinical sheet', [(f'Healthy {n_("sheet", "Healthy")}', HEALTHY), (f'HCC {n_("sheet", "HCC")}', HCC)])
box_(57, 9, 46, 17, 'Contamination QC evaluated', [(f'Healthy {n_("qc", "Healthy")}', HEALTHY), (f'HCC {n_("qc", "HCC")}', HCC)])
box_(117, 9, 60, 17, 'Analysis set', [(f'Healthy {n_("analysis", "Healthy")} (Train 673, Test 336)', HEALTHY),
                                      (f'HCC {n_("analysis", "HCC")} (Train 104, Test 73)', HCC)], ec=INK)
arrow(43.3, 17.5, 56.7, 17.5); arrow(103.3, 17.5, 116.7, 17.5)
box_(3, 31, 46, 20, 'Removed before QC', [(f'Healthy {n_("removed_1", "Healthy")}: outside the', HEALTHY), ('predefined split (59 used as null set)', INK2),
                                          (f'HCC {n_("removed_1", "HCC")}: stage unknown', HCC)], ec=PALE, tc=MUTED)
box_(57, 31, 46, 20, 'Removed by QC', [(f'Healthy {n_("removed_2", "Healthy")} (Train 2, Test 1)', HEALTHY), (f'HCC {n_("removed_2", "HCC")} (Train 10, Test 4)', HCC)],
     ec=PALE, tc=MUTED)
box_(111, 31, 66, 20, 'Subsets of the analysis set', [('Matched HCC–healthy pairs (year, depth): Train 104, Test 73', INK),
                                                      ('HCC with tumour fraction (copy number): 46 (Train 28, Test 18)', HCC),
                                                      ('HCC without tumour-fraction estimate: 131', HCC)], ec=INK)
arrow(50, 17.5, 50, 30.7, PALE); arrow(80, 26.3, 80, 30.7, PALE); arrow(147, 26.3, 147, 30.7, INK)
ov.add_line(Line2D([43.3, 50], [17.5, 17.5], color=PALE, lw=1.0, zorder=1))

# ---------------- b
YB = 57
cv.panel('b', 1, YB, 'Comparison samples (no Train/Test split)')
bx = [(3, 40, 'Hepatitis plasma', [('n = 10', HEP), ('internal comparison set', INK2), ('no age or sex data', INK2)], HEP),
      (47, 48, 'Tissue fragments', [('Tumour 26, non-tumour 27', TUMOUR), ('22 donor pairs', INK2), ('atlas set: 49 tumour, 32 non-tumour', INK2)], TUMOUR),
      (99, 37, 'PBMC genomic DNA', [('30 libraries', PBMC), ('14 plasma-linked pairs', INK2), ('(Train 8, Test 6)', INK2)], PBMC),
      (140, 37, 'Generator base', [('20 Train healthy', HEALTHY), ('fixed by sample-ID rule', INK2), ('independent of results', INK2)], HEALTHY)]
for x, w, t, ln, ec in bx:
    box_(x, YB + 7.5, w, 18, t, ln, ec=ec)

# ---------------- c: study design (former Fig. 1e)
YC = 90
cv.panel('c', 1, YC, 'Study design')
groups = [('Healthy cfDNA', HEALTHY, 'o', f"n = {e_('Healthy cfDNA')}"), ('HCC cfDNA', HCC, '^', f"n = {e_('HCC cfDNA')}"),
          ('Hepatitis cfDNA', HEP, 'D', f"n = {e_('Hepatitis cfDNA')}"),
          ('Liver tissue libraries', TUMOUR, 's', f"tumour {e_('Tumour tissue libraries')}, non-tumour {e_('Adjacent non-tumour tissue libraries')} ({e_('Tissue donor pairs')} pairs)"),
          ('PBMC', PBMC, 'P', f"n = {e_('PBMC')}")]
gx, gw, gh = 3, 50, 8.6; ys = [YC + 7 + k * 10.4 for k in range(5)]
for (nm, col, mk, sub), y in zip(groups, ys):
    rbox(gx, y, gw, gh, ec=col); ov.scatter([gx + 3.6], [y + gh / 2], s=16, marker=mk, color=col, lw=0, zorder=4)
    cv.text(gx + 7, y + 1.1, nm, FS['label'], 'bold', zorder=4); cv.text(gx + 7, y + 4.6, sub, FS['text'], color=INK2, zorder=4)
sx, sw = 70, 50
rbox(sx, YC + 7, sw, 15.5, ec=INK)
cv.text(sx + 2, YC + 8.2, 'Train', FS['label'], 'bold', zorder=4)
cv.text(sx + 2, YC + 11.8, f"healthy {e_('Healthy cfDNA', 'n_train')}, HCC {e_('HCC cfDNA', 'n_train')}; "
        f"{e_('Matched HCC-healthy pairs (year, depth)', 'n_train')} matched pairs\nbaseline, direction and covariates fixed", FS['text'], color=INK2,
        linespacing=1.3, zorder=4)
rbox(sx, YC + 27.5, sw, 15.5, ec=INK)
cv.text(sx + 2, YC + 28.7, 'Test', FS['label'], 'bold', zorder=4)
cv.text(sx + 2, YC + 32.3, f"healthy {e_('Healthy cfDNA', 'n_test')}, HCC {e_('HCC cfDNA', 'n_test')}; "
        f"{e_('Matched HCC-healthy pairs (year, depth)', 'n_test')} matched pairs\nsame definitions; internal re-evaluation", FS['text'], color=INK2,
        linespacing=1.3, zorder=4)
rbox(sx, YC + 48, sw, 9.5, ec=PBMC, ls=(0, (3, 2)))
cv.text(sx + 2, YC + 49.2, 'Comparison data', FS['label'], 'bold', zorder=4)
cv.text(sx + 2, YC + 52.8, 'hepatitis, tissue, PBMC (no split)', FS['text'], color=INK2, zorder=4)
# connectors groups -> splits
xj = 62
for y in (ys[0] + gh / 2, ys[1] + gh / 2):
    ov.add_line(Line2D([gx + gw + .3, xj], [y, y], color=MUTED, lw=1.0))
ov.add_line(Line2D([xj, xj], [ys[0] + gh / 2, YC + 35.2], color=MUTED, lw=1.0))
arrow(xj, YC + 14.7, sx - .3, YC + 14.7); arrow(xj, YC + 35.2, sx - .3, YC + 35.2)
xj2 = 58
for y in (ys[2] + gh / 2, ys[3] + gh / 2, ys[4] + gh / 2):
    ov.add_line(Line2D([gx + gw + .3, xj2], [y, y], color=PALE, lw=1.0))
ov.add_line(Line2D([xj2, xj2], [ys[2] + gh / 2, ys[4] + gh / 2], color=PALE, lw=1.0)); arrow(xj2, YC + 52.7, sx - .3, YC + 52.7, PALE)
# analyses
ax_, aw = 135, 43
items = [('Fig. 2', 'HCC polarisation', 'Train and Test, 59 strata'), ('Fig. 3', 'Opposite-uniform fragments', 'matched pairs, Train and Test'),
         ('Fig. 4', 'Tumour link', f"tumour fraction n = {e_('HCC with tumour-fraction estimate')}; tissue; PBMC"),
         ('Fig. 5', 'Hepatitis', f"n = {e_('Hepatitis cfDNA')}"), ('Fig. 6', 'Ageing', 'healthy cohort by age')]
cys = []
for k, (a1, a2, a3) in enumerate(items):
    y = YC + 7 + k * 10.4; cys.append(y + 4.3); rbox(ax_, y, aw, 8.6, ec=GRID)
    cv.text(ax_ + 2, y + 1.1, f'{a1}  {a2}', FS['text'] + .2, 'bold', zorder=4); cv.text(ax_ + 2, y + 4.6, a3, FS['text'], color=INK2, zorder=4)
xj3 = 128
for y in (YC + 14.7, YC + 35.2, YC + 52.7):
    ov.add_line(Line2D([sx + sw + .3, xj3], [y, y], color=MUTED, lw=1.0))
ov.add_line(Line2D([xj3, xj3], [cys[0], cys[-1]], color=MUTED, lw=1.0))
for y in cys:
    arrow(xj3, y, ax_ - .3, y)
cv.save(SUPP / 'Supplementary_Figure_S1', svg=False)
