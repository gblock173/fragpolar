"""Figure 5: hepatitis. Reads only figures/source_data/Fig5*.tsv."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys; sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from style_v8 import *
import pandas as pd
SD=CFG.FIG+'/source_data/'
rd=lambda n: pd.read_csv(SD+n,sep='\t')
rng=np.random.default_rng(20261015)
F=Fig(122)
# ---- a
F.label('a',1,1,'Marker-panel AUC'); ax=F.ax(28,12,52,42); A=rd('Fig5a.tsv')
PC={'PDR':('#7F7F7F','s'),'FDRP':('#CC79A7','D'),'BFV':('#009E73','o')}
cons=['HCC vs healthy','Hepatitis vs healthy','HCC vs hepatitis']; off={'PDR':.22,'FDRP':0,'BFV':-.22}
for ci,c in enumerate(cons):
    y0=2-ci
    for p in ['PDR','FDRP','BFV']:
        r=A[(A.panel==p)&(A.contrast==c)].iloc[0]; y=y0+off[p]; col,m=PC[p]
        ax.plot([r.ci_lo,r.ci_hi],[y,y],color=col,lw=1.2,solid_capstyle='butt'); ax.plot(r.auc,y,m,color=col,ms=3.3,mew=0)
        ax.text(1.02,y,f'{r.auc:.2f}',transform=ax.get_yaxis_transform(),fontsize=5.5,va='center')
ax.axvline(.5,color=MUTED,lw=1,ls=(0,(3,2))); ax.set_xlim(.2,1.0); ax.set_ylim(-.5,2.5)
ax.set_yticks([2,1,0]); ax.set_yticklabels(['HCC vs\nhealthy\n(73 vs 336)','Hepatitis vs\nhealthy\n(10 vs 336)','HCC vs\nhepatitis\n(73 vs 10)'],fontsize=5.5); ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False)
ax.set_xlabel('AUC (95% CI)'); ax.set_xticks([.25,.5,.75,1])
ax.legend(handles=[Line2D([],[],marker=m,ls='-',color=c,ms=3.3,mew=0,label=p) for p,(c,m) in PC.items()],loc='lower center',bbox_to_anchor=(.5,1.0),ncol=3,fontsize=5.5,handlelength=1.1,handletextpad=.3,columnspacing=.8)
# ---- b
F.label('b',96,1,'Per-sample values (all markers)'); Bd=rd('Fig5b.tsv'); Bc=rd('Fig5b_contrasts_AllMarkers.tsv')
order=['Healthy','HCC','Hepatitis']
specs=[('sOUF_pct','sOUF (%)','sOUF','opp'),('BFV_x1e2','BFV (×0.01)','BFV','B'),('strict_mixed_pct','Strict-mixed (%)','Strict-mixed','mixed')]
for k,(col,yl,tt,_) in enumerate(specs):
    ax=F.ax(110+k*24.5,12,17,42)
    boxjitter(ax,[Bd[Bd.group==g][col] for g in order],[0,1,2],[COL[g] for g in order],[MK[g] for g in order],rng,s=2.2)
    ax.set_xticks([0,1,2]); ax.set_xticklabels(['H','HCC','Hep'],fontsize=5.5); ax.set_ylabel(yl,fontsize=5.5,labelpad=0.5); title(ax,tt,pad=3)
    ax.set_xlim(-.6,2.6)
F.text(104,62,'H, healthy (n = 1,009); HCC (n = 177); Hep, hepatitis (n = 10)',5.5,color=MUTED)
# ---- c
F.label('c',1,69,'HCC and hepatitis effects across 59 strata'); C=rd('Fig5c.tsv')
bc={'U':('o','#0072B2'),'mid':('s','#7F7F7F'),'M':('^','#D55E00'),'All':('D','#111111')}
for k,(nm,lab,u) in enumerate([('sOUF','sOUF','%p'),('BFV','BFV','×0.01')]):
    ax=F.ax(14+k*92,81,62,34)
    for b,(m,cl) in bc.items():
        d=C[C.baseline==b]; ax.scatter(d[f'{nm}_HCC_d'],d[f'{nm}_hep_d'],s=9 if b!='All' else 11,marker=m,color=cl,alpha=.8,linewidths=0,label={'U':'U-baseline','mid':'Mid-baseline','M':'M-baseline','All':'All-marker baseline'}[b])
    ax.axhline(0,color=MUTED,lw=1,ls=(0,(3,2)),zorder=1); ax.axvline(0,color=MUTED,lw=1,ls=(0,(3,2)),zorder=1)
    ax.set_xlabel(f'HCC − healthy ({u})'); ax.set_ylabel(f'Hepatitis − healthy ({u})'); title(ax,lab,pad=3)
    nh=(C[f'{nm}_hep_hi']<0).sum(); nu=(C[f'{nm}_hep_lo']>0).sum()
    ax.text(.98,.03,f'HCC > 0: {(C[f"{nm}_HCC_d"]>0).sum()}/59\nHepatitis < 0: {(C[f"{nm}_hep_d"]<0).sum()}/59\n(CI excludes 0: {nh})',transform=ax.transAxes,ha='right',va='bottom',fontsize=5.5)
    if k==0: ax.legend(loc='lower left',bbox_to_anchor=(.25,1.02),ncol=4,fontsize=5.5,handletextpad=.1,columnspacing=1.0,borderaxespad=0)
F.save(5)
