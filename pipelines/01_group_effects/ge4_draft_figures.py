"""그림 초안 (영문 라벨, matplotlib)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
R=CFG.work('group_effects'); T=R+'/results'; F=R+'/figures'
E=pd.read_csv(T+'/ge_strata_effects_all_metrics.tsv',sep='\t'); D=pd.read_csv(T+'/ge_deltaB_decomposition.tsv',sep='\t')
RG=pd.read_csv(T+'/ge_regression_adjusted_by_K_fraglen.tsv',sep='\t'); TS=pd.read_csv(T+'/ge_binstandardized_strata_effects.tsv',sep='\t'); M=pd.read_csv(T+'/ge_margin_summary.tsv',sep='\t')
def save(fig,n): fig.tight_layout(); fig.savefig(f'{F}/{n}.png',dpi=160); fig.savefig(f'{F}/{n}.pdf'); plt.close(fig)
# Fig group-effect step: Delta B and Delta R_B per stratum (Train vs Test)
fig,ax=plt.subplots(1,3,figsize=(13,4.2))
for a,(m,t) in zip(ax,[('B','Delta B (HCC - healthy)'),('R_ratio','Delta R_B = B/[mu(1-mu)]'),('U','Delta U')]):
    x=E[(E.split=='Train')&(E.metric==m)].set_index('stratum').delta; y=E[(E.split=='Test')&(E.metric==m)].set_index('stratum').delta
    a.scatter(x,y.reindex(x.index),s=14); lim=[min(x.min(),y.min()),max(x.max(),y.max())]; a.plot(lim,lim,'k:',lw=.8); a.axhline(0,c='gray',lw=.5); a.axvline(0,c='gray',lw=.5)
    a.set_xlabel('Train (104 vs 673)'); a.set_ylabel('Test (73 vs 336)'); a.set_title(f'{t}\nTrain>0: {(x>0).sum()}/59, Test>0: {(y>0).sum()}/59',fontsize=9)
save(fig,'Fig_ge1_delta_B_R_U_by_stratum')
# decomposition
fig,ax=plt.subplots(1,2,figsize=(11,4.2),sharey=True)
for a,sp in zip(ax,['Train','Test']):
    d=D[D.split==sp].sort_values('delta_B'); y=np.arange(len(d))
    a.barh(y,d.term_V_shift,label='mean-shift term  dV*R',color='#c0504d'); a.barh(y,d.term_R_change,left=d.term_V_shift.clip(upper=0)*0+np.where(d.term_V_shift>0,d.term_V_shift,0),label='R_B change term  dR*V',color='#4f81bd',alpha=.8)
    a.plot(d.delta_B,y,'k.',ms=3,label='total delta B'); a.set_title(sp); a.set_xlabel('contribution to delta B'); a.set_yticks([])
ax[0].legend(fontsize=7); save(fig,'Fig_ge1_decomposition')
# group-effect step retention
fig,ax=plt.subplots(1,4,figsize=(15,3.8))
for a,(sp,mt) in zip(ax,[('Train','B'),('Test','B'),('Train','R_ratio'),('Test','R_ratio')]):
    g=RG[(RG.split==sp)&(RG.metric==mt)]; a.scatter(g.coef_A,g.coef_B,s=14); lim=[min(g.coef_A.min(),g.coef_B.min()),max(g.coef_A.max(),g.coef_B.max())]; a.plot(lim,lim,'k:',lw=.8)
    a.set_xlabel('group coef: depth+age+sex'); a.set_ylabel('+ mean_K + mean_fraglen'); a.set_title(f'{mt} {sp}',fontsize=9)
save(fig,'Fig_ge2_regression_retention')
# group-effect step
fig,ax=plt.subplots(1,3,figsize=(12,3.8))
for bn,ls in [('full','-'),('trH','--')]:
    for sp,c in [('Test','C0'),('Train(strict보유)','C1')]:
        g=M[(M.baseline==bn)&(M.split==sp)].sort_values('margin')
        ax[0].plot(g.margin,g.median_dOpp,ls,c=c,marker='o',label=f'{bn}/{sp.replace("(strict보유)"," strict-subset")}'); ax[1].plot(g.margin,g.median_reduction_pct,ls,c=c,marker='o'); ax[2].plot(g.margin,g.n_dOpp_q05_pos,ls,c=c,marker='o')
ax[0].set_ylabel('median delta opposite-uniform'); ax[1].set_ylabel('median % delta B removed by opposite fragments'); ax[2].set_ylabel('strata (q<.05, delta opp>0)')
for a in ax: a.set_xlabel('margin m: |baseline-0.5|>=m')
ax[0].legend(fontsize=7); save(fig,'Fig_ge3_margin')
print('ok')
