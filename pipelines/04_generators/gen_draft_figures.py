"""generator step figure drafts: Fig 1d (generator 12-feature response heatmaps, Train seed, f = 0.05, with the observed Train matched HCC-healthy vector) and Fig 3 supplementary (per-stratum cosine distributions)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
from pathlib import Path
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
GEN_DIR = Path(CFG.work('generators')); RES = GEN_DIR / 'results'; FIGD = GEN_DIR / 'figures'
PF = ['amf', 'pdr', 'pm', 'me', 'mhl', 'uhl', 'fdrp', 'qfdrp', 'within_read_U', 'between_read_B_emp', 'hof', 'arrangement_excess_per_link']
LAB = ['AMF', 'PDR', 'PM', 'ME', 'MHL', 'UHL', 'FDRP', 'qFDRP', 'U', 'B', 'A', 'J']
plt.rcParams.update({'font.size': 8, 'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False})
RT = pd.read_csv(RES / 'gen_response_TrainSeed20.tsv', sep='\t'); OB = pd.read_csv(RES / 'gen_observed_cells_deltaZ.tsv', sep='\t')
SN = list(dict.fromkeys(RT.stratum)); NAMES = {'G1': 'G1 opposite-uniform', 'G2': 'G2 within-fragment discordant', 'G3': 'G3 independent flips'}
def mat(df, col): return df.pivot(index='stratum', columns='feature', values=col).reindex(index=SN, columns=PF)
# ---- Fig 1d: heatmaps, 59 strata x 12 features; each panel scaled by its own 98th percentile |value| (units differ by orders of magnitude)
panels = [(f'Observed: Train HCC - matched healthy (n=104 pairs)', mat(OB[OB.contrast == 'Train_matched'], 'delta_z'))] + [(NAMES[g] + ' (f=0.05, 20 Train seeds)', mat(RT[RT.condition == f'{g}_0.05'], 'mean')) for g in ('G1', 'G2', 'G3')]
fig, ax = plt.subplots(1, 4, figsize=(11.5, 9.5), sharey=True)
for a, (t, M) in zip(ax, panels):
    v = np.nanpercentile(np.abs(M.to_numpy()), 98); im = a.imshow(M.to_numpy(), aspect='auto', cmap='RdBu_r', vmin=-v, vmax=v)
    a.set_xticks(range(12)); a.set_xticklabels(LAB, rotation=90); a.set_title(t, fontsize=7.5); cb = fig.colorbar(im, ax=a, fraction=.03, pad=.12, orientation='horizontal'); cb.set_label('mean delta z (own scale)', fontsize=7)
ax[0].set_yticks(range(len(SN))); ax[0].set_yticklabels(SN, fontsize=5.5)
fig.suptitle('Fig. 1d draft - generator response of 12 metrics (Train healthy seeds, fixed f = 0.05) vs observed Train vector', fontsize=9); fig.tight_layout(); fig.savefig(FIGD / 'Fig1d_generator_response_heatmap_TrainSeed.png', dpi=200); fig.savefig(FIGD / 'Fig1d_generator_response_heatmap_TrainSeed.pdf'); plt.close(fig)
# ---- Fig 3 supplementary: cosine distributions
CS = pd.read_csv(RES / 'gen_cosine_by_stratum.tsv', sep='\t'); col = {'G1': '#0072B2', 'G2': '#D55E00', 'G3': '#009E73'}
sets = [('TrainSeed20', 'Train_matched', 'Train obs (matched)'), ('TrainSeed20', 'Test_paired', 'Test obs (73 pairs)'), ('TestSeed20_orig', 'Test_paired', 'orig: Test seed, Test obs')]
fig, ax = plt.subplots(1, 2, figsize=(9, 3.8), sharey=True); rng = np.random.default_rng(1)
for a, how, ttl in [(ax[0], 'fixed', 'a  fixed f = 0.05 (primary)'), (ax[1], 'best', 'b  best f per generator and stratum (original rule; inflated)')]:
    for i, (sd, k, lab) in enumerate(sets):
        g = CS[(CS.seedset == sd) & (CS.contrast == k)]
        for j, G in enumerate(['G1', 'G2', 'G3']):
            x = g[f'{G}_{how}'].dropna().to_numpy(); pos = i + (j - 1) * .26
            bp = a.boxplot([x], positions=[pos], widths=.2, showfliers=False, patch_artist=True, medianprops=dict(color='k')); bp['boxes'][0].set(facecolor=col[G], alpha=.55, edgecolor=col[G])
            a.scatter(pos + rng.uniform(-.07, .07, len(x)), x, s=3, color=col[G], alpha=.6, lw=0)
    a.axhline(0, color='grey', lw=.6); a.axhline(.8, color='grey', lw=.6, ls=':'); a.set_xticks(range(3)); a.set_xticklabels([s[2] for s in sets], fontsize=7); a.set_title(ttl, fontsize=8, loc='left')
ax[0].set_ylabel('cosine (observed vs generator), 59 strata'); [ax[1].scatter([], [], color=col[G], label=NAMES[G], s=10) for G in col]; ax[1].legend(fontsize=6.5, frameon=False, loc='lower right')
fig.tight_layout(); fig.savefig(FIGD / 'Fig3_suppl_cosine_distribution.png', dpi=200); fig.savefig(FIGD / 'Fig3_suppl_cosine_distribution.pdf'); plt.close(fig)
print('figures done')
