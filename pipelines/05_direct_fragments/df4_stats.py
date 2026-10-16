"""direct-fragment step: group effects (HCC - matched healthy, paired), matched-dose reductions, read-split analyses; bootstrap over pairs (B=2000, seed 20261014)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import os,pickle
import numpy as np,pandas as pd
from scipy.stats import rankdata
from pathlib import Path
O=Path(CFG.work('direct_fragments'));BOOT=2000;SEED=20261014
V=pickle.load(open(os.environ.get('PKL',O/'work/stratum_values.pkl'),'rb'));STR=V['strata'];sidx={s:i for i,s in enumerate(V['sid'])}
# comparison mode: --compare matched|full (or FRAGPOLAR_COMPARE). matched = paired HCC - matched healthy (pairs file); full = HCC mean - healthy mean over the whole group.
MODE=(_sys.argv[_sys.argv.index('--compare')+1] if '--compare' in _sys.argv else os.environ.get('FRAGPOLAR_COMPARE','matched'));assert MODE in ('matched','full')
P=pd.read_csv(os.environ.get('PAIRS',O/'results/df_pairs.tsv'),sep='\t') if MODE=='matched' else None
GR=pd.read_csv(os.environ.get('GROUPS',O/'results/df_groups.tsv'),sep='\t') if MODE=='full' else None  # columns sid, grp (HCC|healthy), split
def nz(x):return np.where(np.isfinite(x),x,0)
def idx(sp):
 if MODE=='full':
  g=GR[GR.split.eq(sp)];return np.array([sidx[s] for s in g[g.grp.eq('HCC')].sid]),np.array([sidx[s] for s in g[g.grp.eq('healthy')].sid])
 return pair_idx(sp)
def pair_idx(sp):
 q=P[P.split.eq(sp)];return np.array([sidx[s] for s in q.hcc]),np.array([sidx[s] for s in q.healthy])
def diffs(name,ih,ic):return V[name][ih]-V[name][ic]
def bmean(D,W):
 ok=np.isfinite(D);return (W@nz(D))/np.maximum(W@ok.astype(float),1)*np.where((W@ok.astype(float))>0,1,np.nan)
def ci(x,ax=0):return np.nanpercentile(x,2.5,axis=ax),np.nanpercentile(x,97.5,axis=ax)
def point(D):
 ok=np.isfinite(D);return np.where(ok.sum(0)>0,nz(D).sum(0)/np.maximum(ok.sum(0),1),np.nan),ok.sum(0)
def run(sp):
 ih,ic=idx(sp);rng=np.random.default_rng(SEED);npair=len(ih);W=rng.multinomial(npair,np.ones(npair)/npair,size=BOOT).astype(float)
 if MODE=='full':  # unpaired: group means, each group resampled separately
  W0=rng.multinomial(len(ic),np.ones(len(ic))/len(ic),size=BOOT).astype(float)
  def PB(name):
   X1,X0=V[name][ih],V[name][ic];p1,n1=point(X1);p0,_=point(X0);return p1-p0,n1,bmean(X1,W)-bmean(X0,W0)
 else:
  def PB(name):
   D=diffs(name,ih,ic);p,n=point(D);return p,n,bmean(D,W)
 rows=[];nat=[]
 # --- direct-fragment step native (existing definition)
 (pf,npf,bf),(po,_,bo)=PB('full_native'),PB('opp_native')
 red=1-po/pf;bred=1-bo/bf;bad=bf<=0;bred=np.where(bad,np.nan,bred);lo,hi=ci(bred)
 for j,s in enumerate(STR):nat.append(dict(split=sp,stratum=s,n_pairs=int(npf[j]),dB=pf[j],dB_lo=ci(bf)[0][j],dB_hi=ci(bf)[1][j],dB_opp_removed=po[j],reduction=red[j],red_lo=lo[j],red_hi=hi[j],frac_boot_nonpos_dB=float(bad[:,j].mean())))
 # --- direct-fragment step matched-dose
 for c in ('rand','conc','mix'):
  (pf,npf,bf),(po,_,bo),(pc,_,bc)=PB(f'full_{c}'),PB(f'opp_{c}'),PB(f'ctl_{c}');bad=bf<=0
  rv=lambda a,b:np.where(bad,np.nan,1-a/b)
  ro,rc_=1-po/pf,1-pc/pf;diff=ro-rc_;bdiff=np.where(bad,np.nan,(bc-bo)/bf);babs=bc-bo;l1,h1=ci(bdiff);l2,h2=ci(babs);lo_,ho_=ci(rv(bo,bf));lc,hc=ci(rv(bc,bf))
  rfr=np.nanmean(V[f'rfrac_{c}'][np.r_[ih,ic]],axis=0);nm=np.nanmean(V[f'nm_{c}'][np.r_[ih,ic]],axis=0)
  for j,s in enumerate(STR):rows.append(dict(split=sp,control=c,stratum=s,n_pairs=int(npf[j]),mean_markers_support=nm[j],mean_removed_fraction=rfr[j],dB_full=pf[j],dB_after_opp=po[j],dB_after_ctl=pc[j],
   red_opp=ro[j],red_opp_lo=lo_[j],red_opp_hi=ho_[j],red_ctl=rc_[j],red_ctl_lo=lc[j],red_ctl_hi=hc[j],red_diff=diff[j],red_diff_lo=l1[j],red_diff_hi=h1[j],absdiff_ctl_minus_opp=(pc-po)[j],absdiff_lo=l2[j],absdiff_hi=h2[j],frac_boot_nonpos_dB=float(bad[:,j].mean())))
 b1=pd.DataFrame(rows);nt=pd.DataFrame(nat)
 # --- direct-fragment step (a) read-split effects
 a=[]
 qs=['h_B_A','h_B_B','h_R_A','h_R_B','h_rho_A','h_rho_B','h_rc_A','h_rc_B','h_mu_A','h_mu_B']
 for q in qs:
  p,n_,b=PB(q);l,h=ci(b)
  for j,s in enumerate(STR):a.append(dict(split=sp,stratum=s,quantity=q,n_pairs=int(n_[j]),delta=p[j],lo=l[j],hi=h[j]))
 b2a=pd.DataFrame(a)
 # (c) deletion with dose from X's own half (XX) or the other half (XY)
 cc=[]
 for X,Y in (('A','B'),('B','A'),('A','A'),('B','B')):
  nmx=f'x{X}{Y}';(pf,npf,bf),(po,_,bo),(pc,_,bc)=PB(nmx+'_full'),PB(nmx+'_opp'),PB(nmx+'_con');bad=bf<=0
  ro,rcn=1-po/pf,1-pc/pf;bro=np.where(bad,np.nan,1-bo/bf);brc=np.where(bad,np.nan,1-bc/bf);bd=np.where(bad,np.nan,(bc-bo)/bf);lo,hi=ci(bro);l3,h3=ci(bd)
  for j,s in enumerate(STR):cc.append(dict(split=sp,B_half=X,dose_half=Y,stratum=s,n_pairs=int(npf[j]),dB_full=pf[j],dB_after_opp=po[j],dB_after_conc=pc[j],red_opp=ro[j],red_opp_lo=lo[j],red_opp_hi=hi[j],red_conc=rcn[j],red_diff=(ro-rcn)[j],red_diff_lo=l3[j],red_diff_hi=h3[j],frac_boot_nonpos_dB=float(bad[:,j].mean())))
 b2c=pd.DataFrame(cc)
 # (b) correlations within group, across samples
 cr=[]
 def rk(x):return rankdata(x)
 def spear(x,y):
  ok=np.isfinite(x)&np.isfinite(y)
  if ok.sum()<10:return np.nan
  return float(np.corrcoef(rk(x[ok]),rk(y[ok]))[0,1])
 def partial(x,y,z):
  ok=np.isfinite(x)&np.isfinite(y)&np.isfinite(z)
  if ok.sum()<10:return np.nan
  rx,ry,rz=rk(x[ok]),rk(y[ok]),rk(z[ok]);res=lambda u:u-np.polyval(np.polyfit(rz,u,1),rz)
  return float(np.corrcoef(res(rx),res(ry))[0,1])
 for g,ii in (('HCC',ih),('healthy',ic)):
  for j,s in enumerate(STR):
   z=V['h_mu_full'][ii,j]
   for nm,(xb,yr) in {'same_A':('h_B_A','h_rho_A'),'same_B':('h_B_B','h_rho_B'),'cross_AB':('h_B_A','h_rho_B'),'cross_BA':('h_B_B','h_rho_A'),'cross_AB_R':('h_R_A','h_rho_B'),'same_A_R':('h_R_A','h_rho_A')}.items():
    x,y=V[xb][ii,j],V[yr][ii,j];cr.append(dict(split=sp,group=g,stratum=s,pair=nm,n=int(np.sum(np.isfinite(x)&np.isfinite(y))),spearman=spear(x,y),partial_mu=partial(x,y,z)))
 return b1,nt,b2a,b2c,pd.DataFrame(cr)
if __name__=='__main__':
 out={k:[] for k in ('b1','nat','b2a','b2c','b2b')}
 for sp in ('Train','Test'):
  for k,d in zip(out,run(sp)):out[k].append(d)
 names={'b1':'df_matched_dose_by_stratum','nat':'df_native_reduction_by_stratum','b2a':'df_readsplit_group_effects','b2c':'df_readsplit_deletion','b2b':'df_readsplit_correlations'}
 for k,l in out.items():pd.concat(l).to_csv(Path(os.environ.get('OUTD',O/'results'))/(names[k]+'.tsv'),sep='\t',index=False,float_format='%.6g')
 print('written')
