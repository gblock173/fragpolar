"""Figure 4: tumour association. Reads only figures/source_data/Fig4*.tsv."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys; sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from style_v8 import *
import pandas as pd
SD=CFG.FIG+'/source_data/'
rd=lambda n: pd.read_csv(SD+n,sep='\t')
STC={'sOUF':'#882255','strict-mixed':'#44AA99','healthy-direction uniform':'#DDCC77'}
F=Fig(140)
# ---- a
F.label('a',1,1,'Tissue fragment states (22 pairs)'); ax=F.ax(13,13,70,44); A=rd('Fig4a.tsv')
mk=['All markers','U-baseline','M-baseline']; states=['sOUF','strict-mixed','healthy-direction uniform']
bw=.25
for i,s in enumerate(states):
    d=A[A.state==s].set_index('markers').loc[mk]
    x=np.arange(3)+(i-1)*bw
    ax.bar(x,d.mean_diff_pp,bw*.9,color=STC[s],label={'sOUF':'sOUF','strict-mixed':'Strict-mixed','healthy-direction uniform':'Healthy-direction uniform'}[s],linewidth=0)
    for xx,v,q in zip(x,d.mean_diff_pp,d.BH_q):
        pass
Pp=rd('Fig4a_pairs.tsv'); rj=np.random.default_rng(20261015)
for i,s in enumerate(states):
    for j,m in enumerate(mk):
        v=Pp[(Pp.state==s)&(Pp.markers==m)].diff_pp.to_numpy(); xc=j+(i-1)*bw
        ax.scatter(xc+rj.uniform(-.06,.06,len(v)),v,s=2.2,color='k',alpha=.55,lw=0,zorder=4)
        qv=A[(A.state==s)&(A.markers==m)].BH_q.iloc[0]; mv=A[(A.state==s)&(A.markers==m)].mean_diff_pp.iloc[0]
        if qv<0.05:
            ys_=v.max()+.8 if mv>0 else v.min()-.8
            ax.text(xc,ys_,'*',ha='center',va='bottom' if mv>0 else 'top',fontsize=8,fontweight='bold',color=INK,zorder=6)
        q1,md,q3=np.percentile(v,[25,50,75]); ax.plot([xc,xc],[q1,q3],color='k',lw=.7,zorder=5); ax.plot([xc-.06,xc+.06],[md,md],color='k',lw=.9,zorder=5)
ax.axhline(0,color=MUTED,lw=1); ax.set_xticks(range(3)); ax.set_xticklabels(['All\nmarkers','U-baseline\nmarkers','M-baseline\nmarkers'])
ax.set_ylabel('Tumour − adjacent non-tumour (%p)'); ax.set_ylim(-35,26); ax.set_yticks([-30,-20,-10,0,10,20])
ax.legend(loc='lower left',bbox_to_anchor=(0,1.0),ncol=3,fontsize=5.5,handlelength=.9,handletextpad=.4,columnspacing=1.0,borderaxespad=0)
F.text(13,64,'Bars, mean; dots, pairs (all shown); black line, median and IQR; *, BH q < 0.05',5.5,color=MUTED)
# ---- b
F.label('b',96,1,'CNV gain vs neutral (30 HCC)'); ax=F.ax(124,13,34,44); B=rd('Fig4b.tsv')
rows=[('default','sOUF'),('default','strict-mixed'),('strict (4 CpG, 10 fragments)','sOUF'),('strict (4 CpG, 10 fragments)','strict-mixed')]
ys=[3.35,2.65,1.35,.65]
for (c,s),y in zip(rows,ys):
    r=B[(B.condition==c)&(B.state==s)].iloc[0]
    ax.plot([r.ci_lo_pp,r.ci_hi_pp],[y,y],color=STC[s],lw=1.2,solid_capstyle='butt'); ax.plot(r.mean_pp,y,'o',color=STC[s],ms=3.6,mew=0)
    ax.text(1.04,y,transform=ax.get_yaxis_transform(),s=f'{r.mean_pp:+.2f} [{r.ci_lo_pp:.2f}, {r.ci_hi_pp:.2f}]'.replace('-','−'),fontsize=5.5,va='center',ha='left',color=INK,clip_on=False)
ax.axvline(0,color=MUTED,lw=1,ls=(0,(3,2))); ax.set_xlim(-0.4,1.6); ax.set_ylim(0,4.8)
ax.set_yticks([3,1]); ax.set_yticklabels(['Default\ncondition','Matched\n(4 CpG,\n10 fragments)'],fontsize=5.5); ax.tick_params(axis='y',length=0)
ax.spines['left'].set_visible(False); ax.set_xlabel('Gain − neutral (%p)')
ax.text(-.4,3.95,'',fontsize=5)
ax.legend(handles=[Line2D([],[],marker='o',ls='-',color=STC[s],ms=3.6,mew=0,label=l) for s,l in [('sOUF','sOUF'),('strict-mixed','Strict-mixed')]],loc='upper right',bbox_to_anchor=(1.0,1.02),fontsize=5.5,handlelength=1.2,handletextpad=.4)
# ---- c
F.label('c',1,71,'Association with tumour fraction'); C=rd('Fig4c.tsv'); S=rd('Fig4c_stats.tsv')
for k,(col,yl,tt,sx) in enumerate([('sOUF_diff_pp','sOUF − healthy median (%p)','sOUF','sOUF_diff_pp'),('BFV_diff_x1e2','BFV − healthy median (×0.01)','BFV','BFV_diff_x1e2')]):
    ax=F.ax(13+k*42,80,28,43)
    for sp,fc in [('Train',HCC),('Test','none')]:
        d=C[C.split==sp]; ax.scatter(d.f,d[col],s=9,marker='^',facecolors=fc,edgecolors=HCC,linewidths=.8,zorder=3,label=f'{sp} (n = {len(d)})')
    ax.axhline(0,color=MUTED,lw=1,ls=(0,(3,2)),zorder=1); ax.set_xlabel('Tumour fraction f'); ax.set_ylabel(yl); ax.set_xlim(0,.5); mx=C[col].max(); ax.set_ylim(-0.32*mx,1.32*mx)
    s=S[S.y==sx].set_index('subset')
    ax.text(.03,.97,f"ρ = {s.loc['all','rho']:.2f} [{s.loc['all','ci_lo']:.2f}, {s.loc['all','ci_hi']:.2f}]\nTrain {s.loc['Train','rho']:.2f}, Test {s.loc['Test','rho']:.2f}",transform=ax.transAxes,va='top',fontsize=5.5)
    title(ax,tt,pad=3)
    if k==0: ax.legend(loc='lower right',bbox_to_anchor=(1.03,.05),fontsize=5.5,handletextpad=.1,borderaxespad=0)
# ---- d
F.label('d',96,71,'Candidate origins (59 strata)'); ax=F.ax(132,80,32,50); D=rd('Fig4d.tsv')
cc={'tumour':TUMOUR,'liver':'#E69F00','control':PALE,'blood':NONTUM}
short={'tumour_tissue_a0.1':'Tumour tissue mixed in (10%)','TISSUE_tumour_minus_adjacent(paired)':'Tumour − adjacent\nnon-tumour (22 pairs)','nontumour_liver_tissue_a0.02':'Non-tumour liver\ntissue mixed in (2%)',
 'hepatocyte_a0.02':'Hepatocyte reads\nmixed in (2%)','Colon-Left-Epithelial_a0.02':'Colon epithelium\nmixed in (2%, control)','WBC_bulk(CNVS)_a0.02':'Leukocyte reads\nmixed in (2%)','Blood-Granulocytes_a0.02':'Neutrophil reads\nmixed in (2%)','OBS_PBMC_minus_plasma':'PBMC − plasma\n(observed)'}
short['tumour_tissue_a0.1']='Tumour tissue\nmixed in (10%)'
n=len(D)
for i,r in D.iterrows():
    y=n-1-i; ax.plot([r.ci_lo,r.ci_hi],[y,y],color=cc[r.origin_class],lw=1.2,solid_capstyle='butt'); ax.plot(r.cosine,y,'o',color=cc[r.origin_class],ms=3.8,mew=0)
    ax.text(1.05,y,f'{r.cosine:.2f}'.replace('-','−'),fontsize=5.5,va='center',ha='left',clip_on=False)
ax.set_yticks(range(n)); ax.set_yticklabels([short[c] for c in D.candidate_id][::-1],fontsize=5.5)
ax.axvline(0,color=MUTED,lw=1,ls=(0,(3,2))); ax.set_xlim(-.6,1); ax.set_ylim(-.6,n-.4); ax.set_xlabel('Cosine with polarisation vector'); ax.tick_params(axis='y',length=0)
ax.set_xticks([-.5,0,.5,1]); ax.spines['left'].set_visible(False)
F.save(4)
