import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
E=CFG.work('component_origin')
d=pd.read_csv(f'{E}/results/similarity_full.tsv',sep='\t'); v=pd.read_csv(f'{E}/results/effect_vectors_z_full.tsv',sep='\t')
rows=[('종양 조직 10% 혼합','tumour_tissue_a0.1'),('종양-인접 비종양 (쌍 22, 직접)','TISSUE_tumour_minus_adjacent(paired)'),('종양-비종양 10% 혼합','tumour_minus_nontumour_a0.10'),
 ('비종양 간 조직-혈장 (직접)','TISSUE_nontumour_minus_plasma'),('비종양 간 조직 2% 혼합','nontumour_liver_tissue_a0.02'),('간세포 2% 혼합','hepatocyte_a0.02'),
 ('혈구: 호중구 2% 혼합','Blood-Granulocytes_a0.02'),('혈구: 백혈구 전체 2% 혼합','WBC_bulk(CNVS)_a0.02'),('혈구: 림프구↑ 호중구↓ 5%','lymphocyte_up_neutrophil_down_d0.05'),
 ('혈구: 단핵구↑ 호중구↓ 5%','monocyte_up_neutrophil_down_d0.05'),('대장 상피 2% 혼합 (대조)','Colon-Left-Epithelial_a0.02'),('PBMC-혈장 (관측)','OBS_PBMC_minus_plasma'),('간염 10-건강인 (관측)','OBS_Hepatitis_10')]
plt.rcParams['font.family']='NanumGothic' if any('Nanum' in f.name for f in matplotlib.font_manager.fontManager.ttflist) else 'DejaVu Sans'
fig,ax=plt.subplots(1,2,figsize=(13,5.5),sharey=True)
for a,(comp,tt) in zip(ax,[('P','P 양극화 성분 (종양분율 기울기; sOUF, BFV, FPI)'),('H','H 균질화 성분 (종양분율 미측정 131; mixed, WFV, 건강방향 균일, strictM)')]):
    x=d[d.component==comp].set_index('candidate')
    for i,(lab,k) in enumerate(rows):
        r=x.loc[k]; a.barh(i,r.cosine,color='#c0504d' if r.cosine<0 else '#4f81bd'); a.plot([r.cos_lo,r.cos_hi],[i,i],'k-',lw=1.2)
    a.set_yticks(range(len(rows))); a.set_yticklabels([r[0] for r in rows],fontsize=8); a.invert_yaxis(); a.axvline(0,color='k',lw=.6); a.set_xlim(-1,1); a.set_xlabel('층x지표 z 벡터 코사인 (95% bootstrap CI)'); a.set_title(tt,fontsize=9)
plt.tight_layout(); plt.savefig(f'{E}/figures/co_cosine_by_origin.png',dpi=160); plt.close()
fig,ax=plt.subplots(2,3,figsize=(13,8))
pairs=[('P_slope',['OBS_P_slope'],'tumour_tissue_a0.1',['opp','B','BV']),('P_slope',[],'TISSUE_tumour_minus_adjacent(paired)',['opp','B','BV']),('P_slope',[],'hepatocyte_a0.02',['opp','B','BV']),
       ('H_nof131',[],'TISSUE_nontumour_minus_plasma',['mixed','U','conc','strictM']),('H_nof131',[],'hepatocyte_a0.02',['mixed','U','conc','strictM']),('H_nof131',[],'OBS_PBMC_minus_plasma',['mixed','U','conc','strictM'])]
cols=dict(opp='C0',B='C1',BV='C2',mixed='C3',U='C4',conc='C5',strictM='C6')
for a,(o,_,s,mets) in zip(ax.ravel(),pairs):
    vo=v[v.vector==o].set_index('stratum'); vs=v[v.vector==s].set_index('stratum')
    for m in mets: a.scatter(vs[m],vo[m],s=12,c=cols[m],label=m)
    a.axhline(0,color='gray',lw=.5); a.axvline(0,color='gray',lw=.5); a.set_xlabel(s+' (z)',fontsize=8); a.set_ylabel(o+' (z)',fontsize=8); a.legend(fontsize=7)
plt.tight_layout(); plt.savefig(f'{E}/figures/co_strata_scatter.png',dpi=140)
