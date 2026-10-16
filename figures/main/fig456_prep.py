"""Prepare source-data TSVs for Fig 4-6 (earlier). Only reads existing result tables; writes figures/source_data/Fig{4,5,6}*.tsv."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import pandas as pd, numpy as np
from scipy.stats import spearmanr
THA_DIR=CFG.work('tumour_hepatitis_age','results')+'/'; CO_DIR=CFG.work('component_origin','results')+'/'
SD=CFG.FIG+'/source_data/'
rd=lambda f: pd.read_table(f)
def w(df,name): df.to_csv(SD+name,sep='\t',index=False,float_format='%.6g'); print(name,df.shape)
rng=np.random.default_rng(20261015)

# ---- Fig 4a: tissue pairs (22)
g=rd(CFG.inp('prior_tissue_group_effects'))
g=g[(g.contrast=='Tumour_vs_ownNormal_trim')&g.stratum.isin(['All|All','U|All','M|All'])&g.metric.isin(['hopp_strict','mixed_strict','hcon_strict'])]
lab={'hopp_strict':'sOUF','mixed_strict':'strict-mixed','hcon_strict':'healthy-direction uniform'}
sl={'All|All':'All markers','U|All':'U-baseline','M|All':'M-baseline'}
a=pd.DataFrame({'markers':g.stratum.map(sl),'state':g.metric.map(lab),'mean_diff_pp':g.mean_diff*100,'median_diff_pp':g.median_diff*100,'n_pairs':g.n,'BH_q':g.q})
w(a,'Fig4a.tsv')
# ---- Fig 4b: CNV gain-neutral
e=rd(CFG.inp('prior_cnv_effect_summary'))
e=e[(e.contrast=='gain-neutral')&(e.level=='segment')&(e.stratum=='All')&e.scenario.isin(['full_primary','fixed_depth'])&e.metric.isin(['strict_opposite','strict_mixed'])]
b=pd.DataFrame({'condition':e.scenario.map({'full_primary':'default','fixed_depth':'strict (4 CpG, 10 fragments)'}),'state':e.metric.map({'strict_opposite':'sOUF','strict_mixed':'strict-mixed'}),
  'n_patients':e.n_patients,'n_healthy_matched':e.n_healthy,'mean_pp':e['mean']*100,'ci_lo_pp':e.lo*100,'ci_hi_pp':e.hi*100})
w(b,'Fig4b.tsv')
# ---- samples
st=rd(THA_DIR+'sample_table.tsv'); L=pd.read_parquet(THA_DIR+'sample_strata_long.parquet')
def vals(metric,stratum): return L[(L.metric==metric)&(L.stratum==stratum)].set_index('sid').value
opp=vals('opp','All|All'); B=vals('B','All|All'); mx=vals('mixed','All|All')
htr=st[(st.group=='healthy')&(st.split=='Train')].AnalysisID
base_opp=opp.loc[htr.values].median(); base_B=B.loc[htr.values].median(); print('base',base_opp,base_B)
f=rd(THA_DIR+'tf_f_samples.tsv').rename(columns={'AnalysisID':'sample'})
c=f[['sample','f','split']].copy(); c['sOUF_diff_pp']=(opp.loc[f['sample']].values-base_opp)*100; c['BFV_diff_x1e2']=(B.loc[f['sample']].values-base_B)*100
w(c,'Fig4c.tsv')
rows=[]
def sp(x,y,sub):
    r=spearmanr(x,y)[0]; bs=[]
    for _ in range(2000):
        i=rng.integers(0,len(x),len(x)); bs.append(spearmanr(x[i],y[i])[0])
    return dict(subset=sub,n=len(x),rho=r,ci_lo=np.nanpercentile(bs,2.5),ci_hi=np.nanpercentile(bs,97.5))
for yname in ['sOUF_diff_pp','BFV_diff_x1e2']:
    for sub,m in [('all',c.split.notna()),('Train',c.split=='Train'),('Test',c.split=='Test')]:
        d=c[m]; r=sp(d.f.values,d[yname].values,sub); r['y']=yname; rows.append(r)
w(pd.DataFrame(rows),'Fig4c_stats.tsv'); print(pd.DataFrame(rows))
# ---- Fig 4d cosines
s=rd(CO_DIR+'similarity_full.tsv'); s=s[s.component=='P'].set_index('candidate')
pick=[('tumour_tissue_a0.1','Tumour tissue DNA mixed in (10%)','tumour'),('TISSUE_tumour_minus_adjacent(paired)','Tumour minus adjacent non-tumour (22 pairs)','tumour'),
 ('nontumour_liver_tissue_a0.02','Non-tumour liver tissue mixed in (2%)','liver'),('hepatocyte_a0.02','Hepatocyte reads mixed in (2%)','liver'),('Colon-Left-Epithelial_a0.02','Colon epithelium mixed in (2%, control)','control'),
 ('WBC_bulk(CNVS)_a0.02','Leukocyte reads mixed in (2%)','blood'),('Blood-Granulocytes_a0.02','Neutrophil reads mixed in (2%)','blood'),('OBS_PBMC_minus_plasma','PBMC minus plasma (observed)','blood')]
dd=pd.DataFrame([dict(candidate_id=k,label=l,origin_class=cl,cosine=s.loc[k,'cosine'],ci_lo=s.loc[k,'cos_lo'],ci_hi=s.loc[k,'cos_hi']) for k,l,cl in pick])
# reference: HCC Train vs Test effect-vector cosine (adjusted z, opp/B/BV, 59 strata)
A=rd(THA_DIR+'age_and_HCC_effects_all_strata.tsv'); A=A[(A.kind=='HCC_vs_healthy_adj')&A.metric.isin(['opp','B','BV'])&~A.stratum.str.endswith('|nonIsland')]
tr=A[A.subset=='Train'].set_index(['stratum','metric']).z; te=A[A.subset=='Test'].set_index(['stratum','metric']).z; te=te.loc[tr.index]
cs=float(tr@te/np.linalg.norm(tr)/np.linalg.norm(te)); print('HCC Train-Test cosine',cs,len(tr))
dd=pd.concat([dd,pd.DataFrame([dict(candidate_id='HCC_TrainTest',label='HCC Train vs Test effect vectors (reference)',origin_class='reference',cosine=cs,ci_lo=np.nan,ci_hi=np.nan)])],ignore_index=True)
w(dd[dd.origin_class!='reference'],'Fig4d.tsv'); w(dd[dd.origin_class=='reference'],'Fig4d_reference_check.tsv')

# ---- Fig 5b: per-sample distributions
sel=st[((st.group=='healthy')|(st.group=='liver')|((st.group=='hepatitis')&st.AnalysisID.isin(CFG.ids('hepatitis_included'))))].copy()
sel['grp']=np.where(sel.group=='healthy','Healthy',np.where(sel.group=='liver','HCC','Hepatitis'))
print(sel.grp.value_counts())
b5=pd.DataFrame({'sample':sel.AnalysisID.values,'group':sel.grp.values,'split':sel.split.values,'sOUF_pct':opp.loc[sel.AnalysisID].values*100,'BFV_x1e2':B.loc[sel.AnalysisID].values*100,'strict_mixed_pct':mx.loc[sel.AnalysisID].values*100})
w(b5,'Fig5b.tsv')
# ---- Fig 5c: 59-strata contrasts
D=rd(THA_DIR+'hep_contrasts_included.tsv'); D=D[~D.stratum.str.endswith('|nonIsland')]
assert D.stratum.nunique()==59
rowsc=[]
for stt,gg in D.groupby('stratum'):
    r={'stratum':stt,'baseline':stt.split('|')[0]}
    for met,nm,sc in [('opp','sOUF',100),('B','BFV',100),('mixed','strict_mixed',100)]:
        for con,cn in [('HCC-H','HCC'),('hep-H','hep')]:
            x=gg[(gg.metric==met)&(gg.contrast==con)].iloc[0]
            r[f'{nm}_{cn}_d']=x.d*sc; r[f'{nm}_{cn}_lo']=x.lo*sc; r[f'{nm}_{cn}_hi']=x.hi*sc
    rowsc.append(r)
w(pd.DataFrame(rowsc),'Fig5c.tsv')
# Fig 5 summary contrasts (all markers) for annotation
w(D[D.stratum=='All|All'].assign(d_pp=lambda x:x.d*100,lo_pp=lambda x:x.lo*100,hi_pp=lambda x:x.hi*100),'Fig5b_contrasts_AllMarkers.tsv')

# ---- Fig 6a
S=rd(THA_DIR+'age_direction_summary.tsv'); S=S[(S.strata_set=='all59')&S.metric.isin(['opp','B'])&S.subset.isin(['Train_strictAvail','Test'])]
S=S.assign(subset=S.subset.replace({'Train_strictAvail':'Train'}),metric=S.metric.map({'opp':'sOUF','B':'BFV'}))
w(S[['metric','subset','n_strata','age_up','age_down','age_up_CI','age_down_CI','hcc_up','hcc_down','hcc_up_CI','hcc_down_CI']],'Fig6a.tsv')
# ---- Fig 6b/c
hs=st[(st.group=='healthy')&st.split.isin(['Train','Test'])].dropna(subset=['age','sex','log_depth'])
print('healthy with covariates',hs.split.value_counts().to_dict())
Ui=vals('mixed','U|Island'); Un=vals('mixed','U|nonIsland')
out=[];stats=[]
for sp_ in ['Train','Test']:
    d=hs[hs.split==sp_]; X=np.c_[np.ones(len(d)),d.age/10,d.sex==2,d.log_depth]
    for nm,ser in [('Island',Ui),('nonIsland',Un)]:
        y=ser.loc[d.AnalysisID].values*100
        beta=np.linalg.lstsq(X,y,rcond=None)[0]
        # bootstrap CI for age slope
        bs=[]
        for _ in range(1000):
            i=rng.integers(0,len(y),len(y)); bs.append(np.linalg.lstsq(X[i],y[i],rcond=None)[0][1])
        stats.append(dict(split=sp_,series=nm,slope_pp_per10y=beta[1],ci_lo=np.percentile(bs,2.5),ci_hi=np.percentile(bs,97.5),intercept=beta[0],n=len(y)))
        if nm=='Island':
            adj=y-beta[2]*(X[:,2]-X[:,2].mean())-beta[3]*(X[:,3]-X[:,3].mean())
            out.append(pd.DataFrame({'sample':d.AnalysisID.values,'split':sp_,'age':d.age.values,'cpg_island_mixed_adj_pp':adj,'cpg_island_mixed_raw_pp':y}))
w(pd.concat(out),'Fig6b.tsv'); print(pd.DataFrame(stats))
I=rd(THA_DIR+'age_island_vs_nonisland.tsv'); I=I[(I.metric=='mixed')&(I.comparator=='U|nonIsland')&I.subset.isin(['Train','Test'])]
I=I.assign(slope_pp_per10y=I.coef_per10y*100,ci_lo=I.lo*100,ci_hi=I.hi*100,series=I.series.map({'Island':'CpG island (U-baseline)','comparator':'Non-CpG island (U-baseline)','Island_minus_comparator':'CpG island minus non-CpG island'}))
w(I[['subset','series','slope_pp_per10y','ci_lo','ci_hi','n']],'Fig6c.tsv')
# ---- Fig 6d: redesigned 2026-10-15 (per-stratum effects in healthy Train SD units) -> prep_fig6d.py
