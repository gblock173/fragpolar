"""direct-fragment step: per-sample x 59-stratum summaries from marker-level npz (stratum value = unweighted mean over markers with finite value; >=20 markers, as in strata step).
Output work/stratum_values.pkl: dict name -> array (n_samples, 59) plus nm_* marker counts."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys,os,pickle,time
import numpy as np,pandas as pd
from pathlib import Path
R=Path(CFG.DATA);O=Path(CFG.work('direct_fragments'))
st=pd.read_parquet(CFG.inp('marker_annotation'));dr=st[CFG.col('healthy_direction')].to_numpy(float)
sd=pd.read_csv(CFG.inp('strata_definition'),sep='\t');sd=sd[sd.used];STR=list(sd.stratum)
def mask(r):
 bm=np.ones(len(st),bool) if r.baseline=='All' else (st.band==r.baseline).to_numpy()
 am=np.ones(len(st),bool) if r.axis=='all' else (st.annot==r.annotation).to_numpy() if r.axis=='annot' else (st.genomic_context==r.annotation).to_numpy() if r.axis=='context' else st['atlas_'+r.annotation].to_numpy()
 return bm&am
M=np.stack([mask(r) for r in sd.itertuples()]).astype(np.float32);assert M.shape[0]==59
def smean(v,sup=None):
 f=np.isfinite(v) if sup is None else (np.isfinite(v)&sup);x=np.where(f,v,0).astype(np.float32);cnt=M@f.astype(np.float32);s=M@x
 with np.errstate(invalid='ignore',divide='ignore'):out=np.where(cnt>=20,s/cnt,np.nan)
 return out,cnt
def var(S1,S2,n):
 with np.errstate(invalid='ignore',divide='ignore'):
  mu=S1/n;return np.maximum(S2/n-mu*mu,0)
def sample(z):
 o={};d=lambda k:z[k].astype(np.float64)
 n=d('n');r=d('r');Bf=d('B_full');Bo=d('B_opp');Bc=d('B_conc');Br=(d('B_rand_a')+d('B_rand_b'))/2;Bm=(d('B_mix_a')+d('B_mix_b'))/2
 o['full_native'],o['nm_full_native']=smean(Bf);o['opp_native'],o['nm_opp_native']=smean(Bo)
 for nm,Bx,ba,bb in (('rand',Br,d('B_rand_a'),d('B_rand_b')),('conc',Bc,Bc,Bc),('mix',Bm,d('B_mix_a'),d('B_mix_b'))):
  sup=np.isfinite(Bx)&np.isfinite(Bo)&np.isfinite(Bf)
  o[f'full_{nm}'],o[f'nm_{nm}']=smean(Bf,sup);o[f'opp_{nm}'],_=smean(Bo,sup);o[f'ctl_{nm}'],_=smean(Bx,sup)
  if nm!='conc':o[f'ctlA_{nm}'],_=smean(ba,sup);o[f'ctlB_{nm}'],_=smean(bb,sup)
  rf=np.where(sup&(n>0),r/np.maximum(n,1),np.nan);o[f'rfrac_{nm}'],_=smean(rf,sup)
 # direct-fragment step half quantities
 nA,nB=d('n_A'),d('n_B');S1={'A':d('S1_A'),'B':d('S1_B')};S2={'A':d('S2_A'),'B':d('S2_B')};nn={'A':nA,'B':nB};op={'A':d('opp_A'),'B':d('opp_B')};cn={'A':d('con_A'),'B':d('con_B')}
 V=(nA>=5)&(nB>=5);Bh={};mu={};rho={};rc={}
 for h in 'AB':
  Bh[h]=var(S1[h],S2[h],nn[h]);mu[h]=S1[h]/np.maximum(nn[h],1);rho[h]=op[h]/np.maximum(nn[h],1);rc[h]=cn[h]/np.maximum(nn[h],1)
  o[f'h_B_{h}'],o['nm_V']=smean(Bh[h],V);o[f'h_mu_{h}'],_=smean(mu[h],V);o[f'h_rho_{h}'],_=smean(rho[h],V);o[f'h_rc_{h}'],_=smean(rc[h],V)
  den=mu[h]*(1-mu[h]);Rh=np.where(den>1e-3,Bh[h]/np.maximum(den,1e-3),np.nan);o[f'h_R_{h}'],_=smean(Rh,V)
 o['h_mu_full'],_=smean((S1['A']+S1['B'])/np.maximum(nA+nB,1),V)
 for X,Y in (('A','A'),('B','B'),('A','B'),('B','A')):
  c_o=rho[Y]*nn[X];c_c=rc[Y]*nn[X];n2o=nn[X]-c_o;n2c=nn[X]-c_c;sup=V&(n2o>=5)&(n2c>=5)
  def dele(c,n2,p):
   with np.errstate(invalid='ignore',divide='ignore'):
    m2=(S1[X]-c*p)/n2;v2=(S2[X]-c*p)/n2;return np.maximum(v2-m2*m2,0)
  o[f'x{X}{Y}_full'],o[f'nm_x{X}{Y}']=smean(Bh[X],sup);o[f'x{X}{Y}_opp'],_=smean(dele(c_o,n2o,dr),sup);o[f'x{X}{Y}_con'],_=smean(dele(c_c,n2c,1-dr),sup);o[f'x{X}{Y}_rho'],_=smean(rho[Y],sup)
 return o
if __name__=='__main__':
 t=time.time();S=pd.read_csv(os.environ.get('SAMPLES',O/'results/df_samples.tsv'),sep='\t');vals={}
 for i,sid in enumerate(S.sid):
  z=np.load(O/'work/ml'/(sid+'.npz'));o=sample(z)
  for k,v in o.items():vals.setdefault(k,[]).append(v)
  if i%50==0:print(i,sid,round(time.time()-t),flush=True)
 vals={k:np.array(v) for k,v in vals.items()};vals['strata']=STR;vals['sid']=list(S.sid)
 pickle.dump(vals,open(os.environ.get('PKL',O/'work/stratum_values.pkl'),'wb'));print('done',round(time.time()-t))
