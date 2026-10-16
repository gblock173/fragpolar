import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys; sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from common import *
disp = {m[0]: m[1] for m in METRICS}
# a: Test matched pairs, signed z
d = pd.read_csv(AN / 'generators/results/gen_observed_cells_deltaZ_rebuilt.tsv', sep='\t')
d = d[d.contrast == 'Test_paired'].copy(); d['metric'] = d.feature.map(disp)
d[['stratum', 'metric', 'delta_z', 'p', 'q', 'n']].rename(columns={'n': 'n_pairs'}).to_csv(SD / 'Fig2a.tsv', sep='\t', index=False)
# b: strata effects Train (train-rebuild step) vs Test (A)
def eff(root):
    e = pd.read_csv(AN / root / 'results/ge_strata_effects_all_metrics.tsv', sep='\t')
    if root == 'train_rebuild':  # train-rebuild step labels its all-777 strict-metric rows 'Train(strict보유)'
        x = e[e.split == 'Train(strict보유)']; assert (x.n_healthy == 673).all() and (x.n_HCC == 104).all(); e = e.copy(); e.loc[e.split == 'Train(strict보유)', 'split'] = 'Train'
    return e[e.split.isin(['Train', 'Test'])]
At, Ae = eff('train_rebuild'), eff('group_effects')
rows = []
for m, nm in (('opp', 'sOUF'), ('B', 'BFV'), ('R_ratio', 'FPI')):
    tr = At[(At.split == 'Train') & (At.metric == m)].set_index('stratum'); te = Ae[(Ae.split == 'Test') & (Ae.metric == m)].set_index('stratum')
    for s in tr.index: rows.append((nm, s, tr.loc[s, 'delta'], tr.loc[s, 'q_BH'], te.loc[s, 'delta'], te.loc[s, 'q_BH']))
pd.DataFrame(rows, columns=['metric', 'stratum', 'delta_train', 'q_train', 'delta_test', 'q_test']).to_csv(SD / 'Fig2b.tsv', sep='\t', index=False)
# c: delta-BFV decomposition
dec = pd.concat([pd.read_csv(AN / 'train_rebuild/results/ge_deltaB_decomposition.tsv', sep='\t').query("split=='Train'"),
                 pd.read_csv(AN / 'group_effects/results/ge_deltaB_decomposition.tsv', sep='\t').query("split=='Test'")])
dec[['split', 'stratum', 'delta_B', 'term_V_shift', 'term_R_change', 'term_interaction_resid', 'frac_V_shift', 'frac_R_change']].rename(columns={'delta_B': 'delta_BFV', 'term_V_shift': 'term_mean_shift', 'term_R_change': 'term_FPI_change', 'frac_V_shift': 'frac_mean_shift', 'frac_R_change': 'frac_FPI_change'}).to_csv(SD / 'Fig2c.tsv', sep='\t', index=False)
# d: retention after covariate adjustment
rt = []
for root, sp in (('train_rebuild', 'Train'), ('group_effects', 'Test')):
    r = pd.read_csv(AN / root / 'results/ge_regression_adjusted_by_K_fraglen.tsv', sep='\t'); r = r[r.split == sp] if root == 'group_effects' else r[r.split.isin(['Train', 'Train(strict보유)'])]
    for m, nm in (('opp', 'sOUF'), ('B', 'BFV'), ('R_ratio', 'FPI')):
        for _, x in r[r.metric == m].iterrows(): rt.append((sp, nm, x.stratum, x.coef_A, x.coef_B, x.retained_pct, x.B_q))
k = pd.read_csv(AN / 'fpi_k/results/fpik_regression_by_stratum.tsv', sep='\t')
for sp in ('Train', 'Test'):
    for _, x in k[(k.split == sp) & (k.metric == 'FPI_K')].iterrows(): rt.append((sp, 'FPI_K', x.stratum, x.coef_A, x.coef_B, x.retained_pct, x.B_q))
pd.DataFrame(rt, columns=['split', 'metric', 'stratum', 'coef_A', 'coef_B', 'retained_pct', 'q_B']).to_csv(SD / 'Fig2d.tsv', sep='\t', index=False)
# e: per-patient All|All sOUF and BFV
L = pd.read_parquet(AN / 'tumour_hepatitis_age/results/sample_strata_long.parquet'); L = L[L.stratum == 'All|All']
W_ = L.pivot(index='sid', columns='metric', values='value')[['opp', 'B']].rename(columns={'opp': 'sOUF', 'B': 'BFV'}).reset_index()
T = pd.read_csv(AN / 'tumour_hepatitis_age/results/sample_table.tsv', sep='\t')[['AnalysisID', 'group', 'split']].rename(columns={'AnalysisID': 'sid'})
fs = pd.read_csv(AN / 'tumour_hepatitis_age/results/tf_f_samples.tsv', sep='\t')[['AnalysisID', 'f']].rename(columns={'AnalysisID': 'sid', 'f': 'tumour_fraction'})
X = W_.merge(T, on='sid').merge(fs, on='sid', how='left')
X = X[X.group.isin(['healthy', 'liver'])].copy()
X['category'] = np.where(X.group == 'healthy', 'Healthy', np.where(X.tumour_fraction.notna(), 'HCC, tumour fraction estimated', 'HCC, not estimated'))
X.to_csv(SD / 'Fig2e.tsv', sep='\t', index=False)
print(X.category.value_counts())
ht = X[(X.category == 'Healthy') & (X.split == 'Train')]
for c in ['HCC, tumour fraction estimated', 'HCC, not estimated']:
    g = X[X.category == c]
    print(c, {m: round((g[m].median() - ht[m].median()) / ht[m].std(), 2) for m in ['sOUF', 'BFV']}, {m: round((g[m].median() - ht[m].median()) / (1.4826 * (ht[m] - ht[m].median()).abs().median()), 2) for m in ['sOUF', 'BFV']})
print(rt and pd.DataFrame(rt, columns=['s','m','st','r','c','ret','q']).groupby(['s','m']).ret.describe())
print(d.delta_z.describe())
