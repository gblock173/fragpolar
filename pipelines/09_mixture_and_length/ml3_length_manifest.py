"""length step 표본 manifest: HCC 177 / 건강인 1009 / 간염(sample_lists.hepatitis_included) / PBMC 30 의 raw cpg_calls 경로."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import pandas as pd, os
V = CFG.WORK; OUT = V + '/mixture_and_length'
MAIN = CFG.CPG_CALLS + '/'
HEP = CFG.HEP_CPG_CALLS + '/'
PB = V + '/pbmc/data/cpg_calls/'
s = pd.read_csv(V + '/tumour_hepatitis_age/results/sample_table.tsv', sep='\t')
rows = []
s = s[(s.group != 'hepatitis') | s.AnalysisID.isin(CFG.ids('hepatitis_included'))]
for r in s.itertuples():
    fn = r.AnalysisID + '.cpg_calls.parquet'
    p = HEP + fn if (r.group == 'hepatitis' and os.path.exists(HEP + fn)) else MAIN + fn
    rows.append(dict(sid=r.AnalysisID, group=r.group, source=r.source, path=p, exists=os.path.exists(p)))
for f in sorted(os.listdir(PB)):
    rows.append(dict(sid=f.split('.')[0], group='pbmc', source='pbmc', path=PB + f, exists=True))
m = pd.DataFrame(rows); m['task'] = range(len(m))
print(m.groupby('group').exists.agg(['sum', 'size']))
assert m.exists.all()
m.to_csv(OUT + '/data/length_manifest.tsv', sep='\t', index=False)
