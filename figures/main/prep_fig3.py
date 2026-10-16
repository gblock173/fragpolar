import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys; sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from common import *
# b: observed vs expected dBFV (mixture step)
e = pd.concat([pd.read_csv(AN / 'mixture_and_length/results/mix_strata_point.tsv', sep='\t'), pd.read_csv(AN / 'mixture_and_length/results/mix_strata_point_Test.tsv', sep='\t')])
e[['split', 'stratum', 'delta_sOUF', 'obs_dB', 'exp_M1', 'exp_M1c']].rename(columns={'obs_dB': 'observed_dBFV', 'exp_M1': 'expected_dBFV_opposite_mix', 'exp_M1c': 'expected_dBFV_same_direction_mix'}).to_csv(SD / 'Fig3b.tsv', sep='\t', index=False)
for sp, g in e.groupby('split'): print(sp, 'sum ratio', g.obs_dB.sum() / g.exp_M1.sum(), 'pearson', np.corrcoef(g.obs_dB, g.exp_M1)[0, 1])
# c: native reduction
n = pd.read_csv(AN / 'direct_fragments/results_genpairs/df_native_reduction_by_stratum.tsv', sep='\t')
n[['split', 'stratum', 'n_pairs', 'dB', 'dB_opp_removed', 'reduction', 'red_lo', 'red_hi']].rename(columns={'dB': 'dBFV', 'dB_opp_removed': 'dBFV_after_sOUF_removal'}).to_csv(SD / 'Fig3c.tsv', sep='\t', index=False)
# d: matched-dose controls
m = pd.read_csv(AN / 'direct_fragments/results_genpairs/df_matched_dose_by_stratum.tsv', sep='\t')
name = {'rand': 'Random', 'conc': 'Same-direction uniform', 'mix': 'Strict-mixed'}
rows = []
for _, r in n.iterrows(): rows.append((r.split, 'sOUF', r.stratum, r.reduction, r.red_lo, r.red_hi, np.nan, np.nan, np.nan))
for _, r in m.iterrows(): rows.append((r.split, name[r.control], r.stratum, r.red_ctl, r.red_ctl_lo, r.red_ctl_hi, r.red_diff, r.red_diff_lo, r.red_diff_hi))
pd.DataFrame(rows, columns=['split', 'removed', 'stratum', 'reduction', 'red_lo', 'red_hi', 'diff_vs_sOUF', 'diff_lo', 'diff_hi']).to_csv(SD / 'Fig3d.tsv', sep='\t', index=False)
# e: read split
b = pd.read_csv(AN / 'direct_fragments/results_genpairs/df_readsplit_deletion.tsv', sep='\t'); rows = []
for _, r in b[(b.B_half == 'A')].iterrows():
    if r.dose_half == 'B': rows.append((r.split, 'Other half, sOUF', r.stratum, r.red_opp, r.red_opp_lo, r.red_opp_hi)); rows.append((r.split, 'Other half, same-direction', r.stratum, r.red_conc, np.nan, np.nan))
    else: rows.append((r.split, 'Same half, sOUF', r.stratum, r.red_opp, r.red_opp_lo, r.red_opp_hi))
pd.DataFrame(rows, columns=['split', 'dose_source', 'stratum', 'fraction_explained', 'lo', 'hi']).to_csv(SD / 'Fig3e.tsv', sep='\t', index=False)
# f: cosine
c = pd.read_csv(AN / 'generators/results/gen_cosine_by_stratum_rebuilt.tsv', sep='\t'); c = c[(c.seedset == 'TrainSeed20') & c.contrast.isin(['Train_matched', 'Test_paired'])]
pol = pd.read_csv(AN / 'generators/results/gen_polarization_by_stratum_rebuilt.tsv', sep='\t'); pol = pol[(pol.seedset == 'TrainSeed20') & (pol.f_rule == 'fixed')]
o = c[['contrast', 'stratum', 'G1_0.05', 'G2_0.05', 'G3_0.05']].merge(pol[['contrast', 'stratum', 'polarization', 'within_fragment']], on=['contrast', 'stratum'])
o['split'] = o.contrast.map({'Train_matched': 'Train', 'Test_paired': 'Test'}); o = o.rename(columns={'G1_0.05': 'cos_G1', 'G2_0.05': 'cos_G2', 'G3_0.05': 'cos_G3'})
o[['split', 'stratum', 'cos_G1', 'cos_G2', 'cos_G3', 'polarization', 'within_fragment']].to_csv(SD / 'Fig3f.tsv', sep='\t', index=False)
print(o.groupby('split')[['cos_G1', 'cos_G2', 'cos_G3', 'polarization', 'within_fragment']].agg({'cos_G1': 'median', 'cos_G2': 'median', 'cos_G3': 'median', 'polarization': 'sum', 'within_fragment': 'sum'}))
print(pd.DataFrame(rows).groupby([0, 1])[3].median())
