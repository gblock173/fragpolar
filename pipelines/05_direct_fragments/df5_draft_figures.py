import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import numpy as np,pandas as pd,matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import os;O=Path(CFG.work('direct_fragments'));Rz=Path(os.environ.get('OUTD',O/'results'));F=Path(os.environ.get('FIGD',O/'figures'))
b1=pd.read_csv(Rz/'df_matched_dose_by_stratum.tsv',sep='\t');nat=pd.read_csv(Rz/'df_native_reduction_by_stratum.tsv',sep='\t')
a=pd.read_csv(Rz/'df_readsplit_group_effects.tsv',sep='\t');c=pd.read_csv(Rz/'df_readsplit_deletion.tsv',sep='\t');cr=pd.read_csv(Rz/'df_readsplit_correlations.tsv',sep='\t')
COL={'rand':'#7f7f7f','conc':'#2a7ab0','mix':'#e0a030'};LAB={'rand':'random','conc':'concordant strict-uniform','mix':'strict-mixed'}
plt.rcParams.update({'font.size':8})
# Fig B1a: All|All reductions
fig,ax=plt.subplots(1,2,figsize=(9,3.6),sharey=True)
for k,sp in enumerate(('Train','Test')):
 g=b1[(b1.split==sp)&(b1.stratum=='All|All')].set_index('control');x=[];lab=[];cols=[]
 o=g.iloc[0];xs=[('opposite (all r)',o.red_opp,o.red_opp_lo,o.red_opp_hi,'#c0392b')]+[(f'{LAB[cc]} (r matched)',g.loc[cc,'red_ctl'],g.loc[cc,'red_ctl_lo'],g.loc[cc,'red_ctl_hi'],COL[cc]) for cc in ('rand','conc','mix')]
 for i,(l,v,lo,hi,cl) in enumerate(xs):ax[k].barh(i,v,color=cl,xerr=[[v-lo],[hi-v]],capsize=2)
 ax[k].set_yticks(range(4));ax[k].set_yticklabels([t[0] for t in xs]);ax[k].invert_yaxis();ax[k].set_xlabel('fraction of HCC-healthy dB removed');ax[k].set_title(f'{sp}: All markers (n pairs={int(o.n_pairs)})');ax[k].axvline(0,color='k',lw=.5)
fig.tight_layout();fig.savefig(F/'df_AllMarkers_reduction.png',dpi=200);plt.close()
# Fig direct-fragment step: forest of (red_opp - red_ctl) over 59 strata
for sp in ('Train','Test'):
 sub=b1[b1.split==sp];order=sub[sub.control=='conc'].sort_values('red_diff').stratum.tolist() if False else list(dict.fromkeys(sub.stratum))
 fig,axs=plt.subplots(1,3,figsize=(11,11),sharey=True)
 for k,cc in enumerate(('rand','conc','mix')):
  g=sub[sub.control==cc].set_index('stratum').loc[order];y=np.arange(len(order))
  bad=g.frac_boot_nonpos_dB>0.05;axs[k].errorbar(g.red_diff,y,xerr=[np.nan_to_num(g.red_diff-g.red_diff_lo),np.nan_to_num(g.red_diff_hi-g.red_diff)],fmt='o',ms=3,color=COL[cc],ecolor='#999999',lw=.8)
  axs[k].scatter(g.red_diff[bad],y[bad],s=14,facecolor='none',edgecolor='k',zorder=3);axs[k].axvline(0,color='k',lw=.6);axs[k].set_title(f'opp-removal minus {LAB[cc]} removal\nreduction difference ({sp})');axs[k].set_xlabel('difference in fraction of dB removed')
 axs[0].set_yticks(range(len(order)));axs[0].set_yticklabels(order,fontsize=6);axs[0].invert_yaxis();fig.text(.01,.005,'open circle: bootstrap dB(HCC-healthy) not clearly positive (>5% draws <=0); ratio unstable',fontsize=6)
 fig.tight_layout();fig.savefig(F/f'df_forest_{sp}.png',dpi=170);plt.close()
# Fig direct-fragment step
fig,ax=plt.subplots(1,3,figsize=(13,4))
for sp,mk in (('Train','o'),('Test','s')):
 x=a[(a.split==sp)&(a.quantity=='h_B_A')].set_index('stratum').delta;y=a[(a.split==sp)&(a.quantity=='h_rho_B')].set_index('stratum').delta
 ax[0].scatter(x,y.loc[x.index],s=12,marker=mk,label=sp,alpha=.8)
ax[0].axhline(0,color='k',lw=.5);ax[0].axvline(0,color='k',lw=.5);ax[0].set_xlabel('dB in half A (HCC-healthy)');ax[0].set_ylabel('d opposite fraction in half B');ax[0].legend();ax[0].set_title('(a) 59 strata: independent halves')
data=[];lab=[]
for sp in ('Train','Test'):
 for g in ('HCC','healthy'):
  for nm,l in (('same_A','same half'),('cross_AB','cross half')):
   v=cr[(cr.split==sp)&(cr.group==g)&(cr.pair==nm)].spearman.dropna();data.append(v.values);lab.append(f'{sp}\n{g}\n{l}')
ax[1].boxplot(data,labels=lab,showfliers=False);ax[1].tick_params(axis='x',labelsize=5);ax[1].set_ylabel('across-sample Spearman, B vs opposite fraction\n(per stratum)');ax[1].set_title('(b) same-half vs cross-half correlation');ax[1].axhline(0,color='k',lw=.5)
data=[];lab=[]
for sp in ('Train','Test'):
 for X,Y,l in (('A','A','A->A (same)'),('A','B','A->B (indep)'),('B','A','B->A (indep)')):
  v=c[(c.split==sp)&(c.B_half==X)&(c.dose_half==Y)].red_opp.dropna();data.append(v.values);lab.append(f'{sp}\n{l}')
 v=c[(c.split==sp)&(c.B_half=='A')&(c.dose_half=='B')].red_conc.dropna();data.append(v.values);lab.append(f'{sp}\nA->B conc')
ax[2].boxplot(data,labels=lab,showfliers=False);ax[2].tick_params(axis='x',labelsize=5);ax[2].set_ylabel('fraction of dB(A) explained by dose-matched removal');ax[2].set_title('(c) dose from other half');ax[2].axhline(0,color='k',lw=.5)
fig.tight_layout();fig.savefig(F/'df_readsplit.png',dpi=200)
print('figs done')
