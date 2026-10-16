"""Source data for Fig. 6d (redesigned 2026-10-15): per-stratum effects of three sources on sOUF, BFV and strict-mixed in healthy Train SD units.
HCC and hepatitis: group mean difference vs healthy (HCC 177, hepatitis 10 [listed samples], healthy 1,009; group-wise bootstrap 95% CI) from the
public-code hepatitis contrast table; ageing: Test healthy age slope per 10 years (sex, log depth adjusted). All effects / healthy Train SD (n = 673)
of the same stratum and metric, so each source keeps its own scale (difference vs per-10-year slope)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys; sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from style_v8 import *
CD = AN / 'tumour_hepatitis_age/results'
HEPC = CD / 'hep_contrasts_included.tsv'   # hepatitis: listed samples only (sample_lists.hepatitis_included)
MET = {'opp': 'sOUF', 'B': 'BFV', 'mixed': 'strict-mixed'}
st = pd.read_csv(CD / 'sample_table.tsv', sep='\t'); ht = st[(st.group == 'healthy') & (st.split == 'Train')].AnalysisID; assert len(ht) == 673
L = pd.read_parquet(CD / 'sample_strata_long.parquet'); L = L[L.metric.isin(MET) & L.sid.isin(ht) & L.stratum.isin(strata_order())]
sd = L.groupby(['stratum', 'metric']).value.std(ddof=1).rename('healthy_train_SD')
D = pd.read_csv(HEPC, sep='\t'); D = D[D.metric.isin(MET) & D.contrast.isin(['HCC-H', 'hep-H']) & D.stratum.isin(strata_order())]
D = D.rename(columns={'d': 'effect', 'lo': 'ci_lo', 'hi': 'ci_hi'}); D['source'] = D.contrast.map({'HCC-H': 'HCC − healthy', 'hep-H': 'Hepatitis − healthy'})
AG = pd.read_csv(CD / 'age_and_HCC_effects_all_strata.tsv', sep='\t')
# check: the age table's Train rows carry the same healthy Train SD
chk = AG[(AG.kind == 'age_per10y') & (AG.subset == 'Train') & AG.metric.isin(MET) & AG.stratum.isin(strata_order())].set_index(['stratum', 'metric']).healthy_SD
assert np.allclose(chk.values, sd.loc[chk.index].values, rtol=1e-6), 'SD mismatch'
AG = AG[(AG.kind == 'age_per10y') & (AG.subset == 'Test') & AG.metric.isin(MET) & AG.stratum.isin(strata_order())]
AG = AG.rename(columns={'coef': 'effect', 'lo': 'ci_lo', 'hi': 'ci_hi'}).assign(source='Age (per 10 years, Test healthy)')
X = pd.concat([D[['source', 'stratum', 'metric', 'effect', 'ci_lo', 'ci_hi']], AG[['source', 'stratum', 'metric', 'effect', 'ci_lo', 'ci_hi']]])
X = X.join(sd, on=['stratum', 'metric'])
for c in ('effect', 'ci_lo', 'ci_hi'): X[c + '_SD'] = X[c] / X.healthy_train_SD
X['direction'] = np.select([X.ci_lo > 0, X.ci_hi < 0], ['increase (CI > 0)', 'decrease (CI < 0)'], 'CI includes 0')
X['metric'] = X.metric.map(MET); X['baseline'] = X.stratum.str.split('|').str[0]; X['stratum_label'] = X.stratum.map(stratum_label)
assert X.groupby(['source', 'metric']).size().eq(59).all()
X[['source', 'metric', 'baseline', 'stratum', 'stratum_label', 'effect', 'ci_lo', 'ci_hi', 'healthy_train_SD', 'effect_SD', 'ci_lo_SD', 'ci_hi_SD', 'direction']].to_csv(
    SD / 'Fig6d.tsv', sep='\t', index=False, float_format='%.6g')
print(X.groupby(['source', 'metric']).direction.value_counts().unstack().fillna(0).astype(int))
