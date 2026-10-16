"""Representative raster (rule fixed before looking at effects): Train, first HCC by sorted AnalysisID and its matched healthy.
Markers: first two markers (index order) in the HCC sample with n>=40, >=3 opposite, concordant and mixed fragment each.
Up to 45 K>=4 fragments per sample drawn with fixed seed; rows ordered opposite / mixed / concordant."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import numpy as np,pandas as pd,pyarrow.parquet as pq,matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from pathlib import Path
R=Path(CFG.DATA);O=Path(CFG.work('direct_fragments'))
RAW=Path(CFG.CPG_CALLS)
P=pd.read_csv(O/'results/df_pairs.tsv',sep='\t');pr=P[P.split.eq('Train')].sort_values('hcc').iloc[0];hcc,hl=pr.hcc,pr.healthy
st=pd.read_parquet(CFG.inp('marker_annotation'),columns=['marker_id',CFG.col('healthy_direction'),'baseline_fullTrain']);z=np.load(O/'work/ml'/(hcc+'.npz'))
cand=np.flatnonzero((z['n']>=40)&(z['r']>=3)&(z['c_con']>=1)&(z['c_mix']>=1))[:2];print('markers',cand,st.marker_id.iloc[cand].tolist())
def frags(sid,mid):
 d=pq.read_table(RAW/(sid+'.cpg_calls.parquet'),columns=['marker_id','fragment_id','pos','call'],filters=[('marker_id','=',mid)]).to_pandas();d=d[d.call.isin([0,1])].drop_duplicates(['fragment_id','pos'])
 bad=d.groupby(['fragment_id','pos']).call.nunique()
 return d
fig,ax=plt.subplots(2,len(cand),figsize=(4.2*len(cand),7.5),squeeze=False);cm=ListedColormap(['#bbbbbb','#ffffff','#222222'])
rec=[]
for c,mi in enumerate(cand):
 mid=st.marker_id.iloc[mi];dirn=int(st[CFG.col('healthy_direction')].iloc[mi])
 for rr,(sid,lab) in enumerate(((hcc,'HCC'),(hl,'healthy (matched)'))):
  d=frags(sid,mid);pos=np.sort(d.pos.unique());g=d.groupby('fragment_id');k=g.call.size();m=g.call.sum();keep=k[k>=4].index;d=d[d.fragment_id.isin(keep)]
  kk=k[keep];mm=m[keep];opp=(mm==kk) if dirn==1 else (mm==0);con=(mm==0) if dirn==1 else (mm==kk);cls=np.where(opp,0,np.where(con,2,1))
  p=(mm/kk);Bm=p.var(ddof=0);nfr=len(kk);rec.append(dict(sample=sid,group=lab,marker=mid,direction=dirn,n=nfr,opposite=int(opp.sum()),mixed=int((cls==1).sum()),concordant=int(con.sum()),B=float(Bm)))
  rng=np.random.default_rng(20261014);sel=rng.choice(nfr,min(45,nfr),replace=False);ids=np.array(keep)[sel];cl=cls[sel];pp=p.to_numpy()[sel];o=np.lexsort((pp,cl));ids=ids[o];cl=cl[o]
  idx={f:i for i,f in enumerate(ids)};pi={q:j for j,q in enumerate(pos)}
  sub=d[d.fragment_id.isin(ids)]
  Y=np.zeros((len(ids),len(pos)));
  for f,q,v in zip(sub.fragment_id,sub.pos,sub.call):Y[idx[f],pi[q]]=2 if v==1 else 1
  a=ax[rr,c];a.imshow(Y,cmap=cm,vmin=0,vmax=2,aspect='auto',interpolation='nearest')
  for t,cc_ in enumerate(cl):a.add_patch(plt.Rectangle((-1.2,t-.5),.6,1,color=['#c0392b','#e0a030','#2a7ab0'][cc_],clip_on=False))
  a.set_xticks([]);a.set_yticks([]);a.set_title(f'{lab} {sid}\n{mid}  n={nfr}, opp={int(opp.sum())}, B={Bm:.3f}',fontsize=8)
fig.text(.5,.01,'rows: red=opposite strict-uniform, orange=mixed, blue=concordant strict-uniform; black=methylated, white=unmethylated, grey=not covered; <=45 fragments',ha='center',fontsize=7)
fig.tight_layout(rect=(0,.03,1,1));fig.savefig(O/'figures/df_raster_Train_example.png',dpi=200);pd.DataFrame(rec).to_csv(O/'results/df_raster_Train_example.tsv',sep='\t',index=False);print(pd.DataFrame(rec))
