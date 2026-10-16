"""Fig 5a source data: panel AUC recomputed from panel scores + samples table.tsv. bootstrap 2000, seed 20261015."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import pandas as pd, numpy as np
from scipy.stats import rankdata
OUT=CFG.FIG+'/source_data/Fig5a.tsv'
s=pd.read_table(CFG.inp('panel_scores')); m=pd.read_table(CFG.inp('samples_table'))
d=m.merge(s,on='AnalysisID')
H=d[(d.group=='healthy')&(d.split=='Test')]; C=d[(d.group=='liver')&(d.split=='Test')]
E=d[(d.group=='hepatitis')&d.AnalysisID.isin(CFG.ids('hepatitis_included'))]
print(len(H),len(C),len(E))
def auc(a,b):  # P(a>b)+0.5 ties
    r=rankdata(np.r_[a,b]); n1=len(a); return (r[:n1].sum()-n1*(n1+1)/2)/(n1*len(b))
rng=np.random.default_rng(20261015)
rows=[]
for sc,lab in [('panel_PDR_30','PDR'),('panel_FDRP_30','FDRP'),('panel_B_30','BFV')]:
    for name,pos,neg,pl,nl in [('HCC vs healthy',C,H,'HCC','Healthy'),('Hepatitis vs healthy',E,H,'Hepatitis','Healthy'),('HCC vs hepatitis',C,E,'HCC','Hepatitis')]:
        a=pos[sc].values; b=neg[sc].values; est=auc(a,b)
        bs=[auc(rng.choice(a,len(a)),rng.choice(b,len(b))) for _ in range(2000)]
        lo,hi=np.percentile(bs,[2.5,97.5])
        rows.append(dict(panel=lab,contrast=name,n_pos=len(a),n_neg=len(b),auc=est,ci_lo=lo,ci_hi=hi))
r=pd.DataFrame(rows); r.to_csv(OUT,sep='\t',index=False,float_format='%.4f'); print(r)
