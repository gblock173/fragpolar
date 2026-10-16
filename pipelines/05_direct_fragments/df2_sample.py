"""direct-fragment step/direct-fragment step per-sample marker-level quantities from the (marker,half,K,M) histogram.
Whole sample (both halves): n, sums of p and p^2, opposite/concordant/mixed counts; r = #opposite strict-uniform fragments (p0 = healthy_direction: 1 -> opposite = all methylated).
Matched-dose removals of r fragments (n-r>=5 required; control NA if fewer than r eligible fragments):
  opp (all r opposite), rand (r random of all, REPS reps, mean), conc (r concordant strict-uniform; all identical so deterministic), mix (r random strict-mixed, REPS reps, mean).
Half arrays (A/B read-split) for direct-fragment step."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import os,sys,time,json,hashlib
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np,pandas as pd
from pathlib import Path
R=Path(CFG.DATA);O=Path(CFG.work('direct_fragments'))
REPS=20;NM=190218
def var(S1,S2,n):
 with np.errstate(invalid='ignore',divide='ignore'):
  mu=S1/n;return np.maximum(S2/n-mu*mu,0)
def run(sid,reps=REPS,hist=None):
 t=time.time();st=pd.read_parquet(CFG.inp('marker_annotation'),columns=[CFG.col('healthy_direction')]);dr=st[CFG.col('healthy_direction')].to_numpy(np.int8);assert len(dr)==NM
 H=pd.read_parquet(hist or O/'work/hist'/(sid+'.parquet'));mk=H.marker.to_numpy(np.int64);hf=H.half.to_numpy(np.int8);k=H.k.to_numpy(np.int64);m=H.m.to_numpy(np.int64);c=H.cnt.to_numpy(np.int64)
 p=m/k;p0=dr[mk];opp=np.where(p0==1,m==k,m==0);con=np.where(p0==1,m==0,m==k);mix=~(opp|con)
 bc=lambda w,sel=None:np.bincount(mk if sel is None else mk[sel],weights=w if sel is None else w[sel],minlength=NM)
 out={}
 def agg(sel):
  w=c*sel;return dict(n=bc(w),S1=bc(w*p),S2=bc(w*p*p),opp=bc(w*opp),con=bc(w*con),mix=bc(w*mix))
 W=agg(np.ones(len(H),bool));A=agg(hf==0);Bh=agg(hf==1)
 n=W['n'];r=W['opp'];ok=n>=5;n2=n-r;p0f=dr.astype(float)
 Bfull=np.where(ok,var(W['S1'],W['S2'],n),np.nan)
 okd=ok&(n2>=5)
 Bopp=np.where(okd,var(W['S1']-r*p0f,W['S2']-r*p0f,n2),np.nan)           # p0^2=p0
 pc=1-p0f;conc_ok=okd&(W['con']>=r)
 Bconc=np.where(conc_ok,var(W['S1']-r*pc,W['S2']-r*pc,n2),np.nan)
 # random removal (all fragments) and strict-mixed removal, only where r>0 (r=0 is a no-op identical to Bfull)
 rng=np.random.default_rng(int(hashlib.sha256(('df_matched_'+sid).encode()).hexdigest()[:8],16))
 def sim(sel_rows,elig,key):
  # sel_rows: histogram rows eligible to be removed; elig: marker mask where control is defined and r>0
  rows=np.flatnonzero(sel_rows&elig[mk]);cnt=c[rows];mkr=np.repeat(mk[rows],cnt);pp=np.repeat(p[rows],cnt)
  uq,inv=np.unique(mkr,return_inverse=True);start=np.r_[0,np.cumsum(np.bincount(inv))[:-1]];rr=r[uq].astype(np.int64)
  S1=np.bincount(mk,weights=c*p*sel_rows,minlength=NM) if key=='mix' else W['S1'];S2=np.bincount(mk,weights=c*p*p*sel_rows,minlength=NM) if key=='mix' else W['S2']
  res=[]
  for rep in range(reps):
   o=np.argsort(inv+rng.random(len(inv)));pos=np.empty(len(inv),np.int64);pos[o]=np.arange(len(inv))-start[inv[o]]  # random rank within marker
   rem=pos<rr[inv];R1=np.bincount(inv[rem],weights=pp[rem],minlength=len(uq));R2=np.bincount(inv[rem],weights=pp[rem]**2,minlength=len(uq))
   res.append((uq,R1,R2))
  return res
 def finish(res,elig,Sall1,Sall2):
  a=np.full(NM,np.nan);b=np.full(NM,np.nan)
  if not res:return a,b
  uq=res[0][0];vals=[]
  for _,R1,R2 in res:
   vals.append(var(Sall1[uq]-R1,Sall2[uq]-R2,n2[uq]))
  vals=np.array(vals);h=len(vals)//2;a[uq]=vals[:h].mean(0);b[uq]=vals[h:].mean(0);return a,b
 el=okd&(r>0)
 res=sim(np.ones(len(H),bool),el,'rand');Br_a,Br_b=finish(res,el,W['S1'],W['S2']);del res
 elm=el&(W['mix']>=r);res=sim(mix,elm,'mix');Sm1=bc(c*p*mix);Sm2=bc(c*p*p*mix)
 Bm_a,Bm_b=finish(res,elm,W['S1'],W['S2']);del res
 for a_ in (Br_a,Br_b):a_[(r==0)&okd]=Bfull[(r==0)&okd]
 for a_ in (Bm_a,Bm_b):a_[(r==0)&okd]=Bfull[(r==0)&okd]
 out.update(n=n.astype(np.int32),r=r.astype(np.int32),c_con=W['con'].astype(np.int32),c_mix=W['mix'].astype(np.int32),S1=W['S1'].astype(np.float32),
  B_full=Bfull.astype(np.float32),B_opp=Bopp.astype(np.float32),B_conc=Bconc.astype(np.float32),B_rand_a=Br_a.astype(np.float32),B_rand_b=Br_b.astype(np.float32),B_mix_a=Bm_a.astype(np.float32),B_mix_b=Bm_b.astype(np.float32))
 for nm,X in (('A',A),('B',Bh)):
  for key in ('n','opp','con'):out[f'{key}_{nm}']=X[key].astype(np.int32)
  out[f'S1_{nm}']=X['S1'].astype(np.float32);out[f'S2_{nm}']=X['S2'].astype(np.float64)
 out['S2']=W['S2']
 d=O/'work/ml';d.mkdir(parents=True,exist_ok=True);np.savez_compressed(d/(sid+'.npz'),**out)
 qa=dict(sid=sid,seconds=time.time()-t,n_ge5=int(ok.sum()),r_pos=int(((r>0)&ok).sum()),Bfull_mean=float(np.nanmean(Bfull)),Bopp_mean=float(np.nanmean(Bopp)),Brand_mean=float(np.nanmean(Br_a)),Bconc_mean=float(np.nanmean(Bconc)),Bmix_mean=float(np.nanmean(Bm_a)),
  sim_n_rand=int(el.sum()),sim_n_mix=int(elm.sum()),frag_total=int(n.sum()))
 (d/(sid+'.QA.json')).write_text(json.dumps(qa,indent=1));print(json.dumps(qa),flush=True)
if __name__=='__main__':
 sid=sys.argv[1] if not sys.argv[1].isdigit() else pd.read_csv(Path(os.environ.get('SAMPLES',O/'results/df_samples.tsv')),sep='\t').set_index('task').sid[int(sys.argv[1])]
 run(sid)
