"""direct-fragment step: HCC-healthy 1:1 pairs (nearest log_depth within era, no replacement, HCC sorted by AnalysisID; same rule as the one-pass engine pairing step).
Applied separately to Train (104 HCC; healthy pool Train 673) and Test (73 HCC; pool Test 336). analysis-freeze sample list."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import pandas as pd,numpy as np
from pathlib import Path
R=Path(CFG.DATA);O=Path(CFG.work('direct_fragments'))
s=pd.read_csv(CFG.inp('samples_table'),sep='\t');s=s[s.set.eq('cohort')&s.group.isin(['liver','healthy'])].copy();s['era']=s.run.astype(str).str[:2]
rows=[]
for sp in ['Train','Test']:
 L=s[s.group.eq('liver')&s.split.eq(sp)].sort_values('AnalysisID');H=s[s.group.eq('healthy')&s.split.eq(sp)];used=set()
 for r in L.itertuples():
  c=H[~H.AnalysisID.isin(used)&H.era.eq(r.era)]
  if not len(c):c=H[~H.AnalysisID.isin(used)]
  j=(c.log_depth-r.log_depth).abs().idxmin();used.add(H.AnalysisID[j])
  rows.append(dict(split=sp,pair_id=f'{sp}_{r.AnalysisID}',hcc=r.AnalysisID,healthy=H.AnalysisID[j],era=r.era,ld_hcc=r.log_depth,ld_healthy=H.log_depth[j],era_match=H.era[j]==r.era,cv_fold_hcc=r.cv_fold))
P=pd.DataFrame(rows);P.to_csv(O/'results/df_pairs.tsv',sep='\t',index=False)
long=pd.concat([P[['split','pair_id','hcc']].rename(columns={'hcc':'sid'}).assign(group='HCC'),P[['split','pair_id','healthy']].rename(columns={'healthy':'sid'}).assign(group='healthy')])
long['task']=np.arange(len(long))+1;long.to_csv(O/'results/df_samples.tsv',sep='\t',index=False)
print(P.groupby('split').agg(n=('hcc','size'),era_ok=('era_match','mean'),mad_ld=('ld_hcc',lambda x:0)).to_string());print((P.ld_hcc-P.ld_healthy).abs().groupby(P.split).describe())
# compare with legacy earlier matching for Test
m=pd.read_csv(CFG.inp('prior_engine_manifest'),sep='\t');m=m[m.role.eq('healthy_matched')];leg=set(zip(m.matched_to,m.sid))
now=set(zip(P[P.split.eq('Test')].hcc,P[P.split.eq('Test')].healthy));print('Test pairs overlapping legacy v5 pairs:',len(leg&now),'of',len(now),'legacy',len(leg))
print(len(long),long.sid.is_unique)
