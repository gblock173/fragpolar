"""Fig 4a pair-level values (22 tumour/adjacent pairs) re-aggregated from component-origin step tissue sufficient statistics
(analysis/component_origin/data/tissue_suff) with the component-origin step stratum rules. Writes source_data/Fig4a_pairs.tsv."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import numpy as np, pandas as pd, scipy.sparse as sp, warnings; warnings.filterwarnings('ignore')
D=CFG.DATA; E=CFG.work('component_origin'); SD=CFG.FIG+'/source_data/'
st=pd.read_parquet(CFG.inp('marker_annotation')); NM=len(st)
sd=pd.read_csv(CFG.inp('strata_definition'),sep='\t'); sd=sd[sd.used]
b0=st.baseline_fullTrain.to_numpy(np.float64); bsel=np.where(b0<.5,1,np.where(b0>.5,0,-1))
masks={'All markers':np.ones(NM,bool),'U-baseline':(st.band=='U').to_numpy(),'M-baseline':(st.band=='M').to_numpy()}
# sanity: the three strata exist in the definition
assert {'All|All','U|All','M|All'} <= set(sd.stratum)
tt=pd.read_csv(f'{E}/data/tissue_tasks.tsv',sep='\t',dtype=str,keep_default_na=False)
def strata(sid):
    d=pd.read_parquet(f'{E}/data/tissue_suff/{sid}.parquet').set_index('marker_id').reindex(st.marker_id)
    n=d.n_reads.to_numpy(float)
    with np.errstate(invalid='ignore',divide='ignore'):
        sM=d.n_strictM.to_numpy()/n; sU=d.n_strictU.to_numpy()/n
    ok=n>=5; sM=np.where(ok,sM,np.nan); sU=np.where(ok,sU,np.nan)
    opp=np.where(bsel==1,sM,np.where(bsel==0,sU,np.nan)); conc=np.where(bsel==1,sU,np.where(bsel==0,sM,np.nan)); mixed=1-sM-sU
    out={}
    for k,m in masks.items():
        v=m&np.isfinite(sM)&np.isfinite(sU); vo=v&np.isfinite(opp)
        out[(k,'sOUF')]=opp[vo].mean() if vo.sum()>=20 else np.nan
        out[(k,'healthy-direction uniform')]=conc[vo].mean() if vo.sum()>=20 else np.nan
        out[(k,'strict-mixed')]=mixed[v].mean() if v.sum()>=20 else np.nan
    return out
import os
tt=tt[[os.path.exists(f'{E}/data/tissue_suff/{s}.parquet') for s in tt.sid]]
V={s:strata(s) for s in tt.sid}; rows=[]
for dn in sorted(set(tt.donor[tt.donor!=''])):
    ts=tt[(tt.donor==dn)&(tt.role=='tumour')].sid; ns=tt[(tt.donor==dn)&(tt.role=='normal')].sid
    if len(ts) and len(ns):
        for key in V[ts.iloc[0]]:
            rows.append(dict(donor=dn,markers=key[0],state=key[1],diff_pp=100*(np.nanmean([V[s][key] for s in ts])-np.nanmean([V[s][key] for s in ns]))))
P=pd.DataFrame(rows); P.to_csv(SD+'Fig4a_pairs.tsv',sep='\t',index=False,float_format='%.6g')
print(P.donor.nunique(),'pairs')
def wilcoxon(v):
    """two-sided exact signed-rank test (continuous data, no ties; scipy unusable in this env)"""
    v=v[v!=0]; n=len(v); r=np.argsort(np.argsort(np.abs(v)))+1; W=int(r[v>0].sum()); tot=n*(n+1)//2
    cnt=np.zeros(tot+1); cnt[0]=1
    for k in range(1,n+1):
        new=cnt.copy(); new[k:]+=cnt[:tot+1-k]; cnt=new
    cnt/=cnt.sum(); lo=min(W,tot-W)
    return type('R',(),{'pvalue':min(1.0,2*cnt[:lo+1].sum())})
rows=[]
for (m,stt),g in P.groupby(['markers','state']):
    v=g.diff_pp.to_numpy(); rows.append(dict(markers=m,state=stt,mean_diff_pp=v.mean(),median_diff_pp=np.median(v),q25_pp=np.percentile(v,25),q75_pp=np.percentile(v,75),n_pairs=len(v),p=wilcoxon(v).pvalue))
A=pd.DataFrame(rows); o=np.argsort(A.p.to_numpy()); p=A.p.to_numpy()[o]; n=len(p)
q=np.minimum.accumulate((p*n/np.arange(1,n+1))[::-1])[::-1]; A['BH_q']=np.nan; A.loc[A.index[o],'BH_q']=np.minimum(q,1)
A.to_csv(SD+'Fig4a.tsv',sep='\t',index=False,float_format='%.6g'); print(A.to_string())
