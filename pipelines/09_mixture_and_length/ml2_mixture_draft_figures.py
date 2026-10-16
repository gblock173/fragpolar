import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
E = CFG.work('mixture_and_length')
P = {'Train': pd.read_csv(E + '/results/mix_strata_point.tsv', sep='\t'), 'Test': pd.read_csv(E + '/results/mix_strata_point_Test.tsv', sep='\t')}
rows = []
for sp, d in P.items():
    ok = d.exp_M3_dU.notna()
    rows.append(dict(split=sp, sumratio_dU_M3=d.obs_dU[ok].sum() / d.exp_M3_dU[ok].sum(), AllAll_dU_obs=d.obs_dU[0], AllAll_dU_M1=d.exp_M1_dU[0], AllAll_dU_M3=d.exp_M3_dU[0],
                     AllAll_dFPI_obs=d.obs_dFPI[0], AllAll_dFPI_M1=d.exp_M1_dFPI[0], n_strata_obsB_pos=int((d.obs_dB > 0).sum()), n_strata_expM1_pos=int((d.exp_M1 > 0).sum()),
                     n_strata_ctrl_neg=int((d.exp_M1c < 0).sum()), median_resid_B=(d.obs_dB - d.exp_M1).median(), AllAll_resid_B=d.obs_dB[0] - d.exp_M1[0],
                     n_strata_obs_gt_M1=int((d.obs_dB > d.exp_M1).sum()), n_strata=int(d.obs_dB.notna().sum())))
pd.DataFrame(rows).to_csv(E + '/results/mix_extra_summary.tsv', sep='\t', index=False); print(pd.DataFrame(rows).T)
fig, ax = plt.subplots(1, 3, figsize=(13, 4))
for sp, c in [('Train', '#3b6ea5'), ('Test', '#c0392b')]:
    d = P[sp]; ax[0].scatter(d.exp_M1 * 100, d.obs_dB * 100, s=14, c=c, label=sp, alpha=.8)
    ax[1].scatter(d.exp_M1c * 100, d.obs_dB * 100, s=14, c=c, alpha=.8)
    ax[2].scatter(d.exp_M1_dU * 100, d.obs_dU * 100, s=14, c=c, alpha=.8)
m = max(P['Train'].obs_dB.max(), P['Train'].exp_M1.max()) * 100
ax[0].plot([0, m], [0, m], 'k--', lw=.8); ax[0].set_xlabel('expected ΔBFV from ΔsOUF (%p x100)'); ax[0].set_ylabel('observed ΔBFV (x100)'); ax[0].set_title('opposite-uniform mixing (59 strata)'); ax[0].legend()
ax[1].set_xlabel('expected ΔBFV, same δ as healthy-direction uniform (x100)'); ax[1].set_title('control'); ax[1].axvline(0, c='gray', lw=.5)
ax[2].plot([-.8, 0], [-.8, 0], 'k--', lw=.8); ax[2].set_xlabel('expected ΔWFV from ΔsOUF (x100)'); ax[2].set_ylabel('observed ΔWFV (x100)'); ax[2].set_title('WFV is not explained by sOUF')
plt.tight_layout(); plt.savefig(E + '/figures/mix_obs_vs_expected.png', dpi=150)
