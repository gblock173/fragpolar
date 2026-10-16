"""ml7: fragment length vs tumour fraction (HCC with CNV tumour fraction, n = 46) — Spearman rho with HCC-sample bootstrap 95% CI.
Inputs (read only): <work>/mixture_and_length/results/len_length_metrics_per_sample.tsv (ml5), <work>/tumour_hepatitis_age/results/tf_f_samples.tsv (tha2),
<work>/tumour_hepatitis_age/results/sample_strata_long.parquet (tha1; All|All). Output: <work>/mixture_and_length/results/len_tumour_fraction_spearman.tsv,
len_tumour_fraction_samples.tsv. Seed 20261015, 5,000 bootstrap resamples."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
AN = Path(CFG.WORK); OUT = AN / 'mixture_and_length/results'
E = pd.read_csv(AN / 'mixture_and_length/results/len_length_metrics_per_sample.tsv', sep='\t').set_index('sid')
f = pd.read_csv(AN / 'tumour_hepatitis_age/results/tf_f_samples.tsv', sep='\t').set_index('AnalysisID')
L = pd.read_parquet(AN / 'tumour_hepatitis_age/results/sample_strata_long.parquet'); L = L[L.stratum == 'All|All'].pivot(index='sid', columns='metric', values='value')
L['FPI'] = L.B / L.V
d = f[['f', 'split']].join(E[['median_len', 'frac_short120', 'frac166', 'frac300']]).join(L[['opp', 'B', 'mixed', 'U', 'FPI']])
assert len(d) == 46 and d.notna().all().all()
d = d.rename(columns={'opp': 'sOUF', 'B': 'BFV', 'mixed': 'strict_mixed', 'U': 'WFV'})
d.rename_axis('sample').reset_index().to_csv(OUT / 'len_tumour_fraction_samples.tsv', sep='\t', index=False, float_format='%.6g')
rng = np.random.default_rng(20261015); NB = 5000
pairs = [('f', 'median_len'), ('f', 'frac_short120'), ('f', 'frac166'), ('f', 'frac300'), ('f', 'sOUF'), ('f', 'BFV'),
         ('frac_short120', 'sOUF'), ('median_len', 'sOUF'), ('frac_short120', 'BFV'), ('median_len', 'BFV')]
rows = []
for x, y in pairs:
    for sub, m in (('All46', np.ones(len(d), bool)), ('Train', (d.split == 'Train').to_numpy()), ('Test', (d.split == 'Test').to_numpy())):
        a, b = d[x].to_numpy()[m], d[y].to_numpy()[m]; r, p = spearmanr(a, b)
        bs = np.array([spearmanr(a[i], b[i])[0] for i in rng.integers(0, len(a), (NB, len(a)))])
        rows.append(dict(x=x, y=y, subset=sub, n=int(m.sum()), rho=r, ci_lo=np.nanpercentile(bs, 2.5), ci_hi=np.nanpercentile(bs, 97.5), p=p))
R = pd.DataFrame(rows); R.to_csv(OUT / 'len_tumour_fraction_spearman.tsv', sep='\t', index=False, float_format='%.6g'); print(R.round(3).to_string())
