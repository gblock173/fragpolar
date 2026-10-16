"""direct-fragment step, whole-group comparison (--compare full): lists all analysis-freeze HCC and healthy cohort samples per split
(Train HCC vs Train healthy, Test HCC vs Test healthy) and the samples still to be extracted (no work/ml/<sid>.npz yet).
Writes results/df_groups.tsv (sid, grp, split), results/df_extra_samples.tsv (sid, task; input of df1/df2 via SAMPLES=...) and
results/df_union_samples.tsv (existing + new, input of df3 via SAMPLES=...). Existing samples are never re-extracted."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import os
import numpy as np,pandas as pd
from pathlib import Path
O=Path(CFG.work('direct_fragments'));(O/'results').mkdir(parents=True,exist_ok=True)
s=pd.read_csv(CFG.inp('samples_table'),sep='\t');s=s[s.set.eq('cohort')&s.group.isin(['liver','healthy'])].copy()
s['grp']=np.where(s.group.eq('liver'),'HCC','healthy')
g=s[['AnalysisID','grp','split']].rename(columns={'AnalysisID':'sid'});g.to_csv(O/'results/df_groups.tsv',sep='\t',index=False)
have={p.name[:-len('.QA.json')] for p in (O/'work/ml').glob('*.QA.json')}
new=g[~g.sid.isin(have)].sort_values('sid');new[['sid']].assign(task=np.arange(len(new))+1).to_csv(O/'results/df_extra_samples.tsv',sep='\t',index=False)
u=sorted(have&set(g.sid))+list(new.sid);pd.DataFrame(dict(sid=u,task=np.arange(len(u))+1)).to_csv(O/'results/df_union_samples.tsv',sep='\t',index=False)
print(g.groupby(['grp','split']).size().to_string());print('already extracted',len(have&set(g.sid)),'to extract',len(new))
