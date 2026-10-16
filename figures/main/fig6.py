"""Figure 6: ageing. Reads only figures/source_data/Fig6*.tsv."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys; sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from style_v8 import *
import pandas as pd
from matplotlib.patches import Patch
SD=CFG.FIG+'/source_data/'
rd=lambda n: pd.read_csv(SD+n,sep='\t')
UPC,UPL,DNL,DNC='#B2182B','#E8A09A','#A6CBE3','#2166AC'   # increase (CI), increase, decrease, decrease (CI)
def divbar(ax,y,upCI,upns,dnns,dnCI,h=.62):
    ax.barh(y,-dnCI,h,color=DNC,linewidth=0); ax.barh(y,-dnns,h,left=-dnCI,color=DNL,linewidth=0)
    ax.barh(y,upns,h,color=UPL,linewidth=0); ax.barh(y,upCI,h,left=upns,color=UPC,linewidth=0)
def leg(ax,**kw):
    ax.legend(handles=[Patch(color=DNC,label='Decrease, CI excludes 0'),Patch(color=DNL,label='Decrease, CI includes 0'),Patch(color=UPL,label='Increase, CI includes 0'),Patch(color=UPC,label='Increase, CI excludes 0')],fontsize=5.5,handlelength=.9,handletextpad=.4,ncol=4,columnspacing=.9,**kw)
F=Fig(190); rng=np.random.default_rng(20261015)
# ---- a
F.label('a',1,1,'Age and HCC effects'); ax=F.ax(32,11,54,40); A=rd('Fig6a.tsv')
rows=[];y=0;yt=[];yl=[]
for met in ['sOUF','BFV']:
    for who,sub in [('Age, Train','Train'),('Age, Test','Test'),('HCC, Train','Train'),('HCC, Test','Test')]:
        r=A[(A.metric==met)&(A['subset']==sub)].iloc[0]
        if who.startswith('Age'): divbar(ax,y,r.age_up_CI,r.age_up-r.age_up_CI,r.age_down-r.age_down_CI,r.age_down_CI)
        else: divbar(ax,y,r.hcc_up_CI,r.hcc_up-r.hcc_up_CI,r.hcc_down-r.hcc_down_CI,r.hcc_down_CI)
        yt.append(y); yl.append(who); y-=1
    y-=.6
ax.set_yticks(yt); ax.set_yticklabels(yl,fontsize=5.5); ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False)
ax.axvline(0,color=INK,lw=1); ax.set_xlim(-59,59); ax.set_xticks([-59,-40,-20,0,20,40,59]); ax.set_xticklabels(['59','40','20','0','20','40','59'])
ax.set_xlabel('Number of strata (of 59)')
for lab,yy in [('sOUF',yt[0]+.9),('BFV',yt[4]+.9)]: ax.text(-86,yy,lab,fontsize=6,fontweight='bold',ha='left',va='center',clip_on=False)
# ---- b
F.label('b',96,1,'CpG island strict-mixed fraction with age'); Bd=rd('Fig6b.tsv'); C=rd('Fig6c.tsv')
for k,sp in enumerate(['Train','Test']):
    ax=F.ax(114+k*33,11,28,40); d=Bd[Bd.split==sp]; ctr=d.cpg_island_mixed_adj_pp.mean()
    ax.scatter(d.age,d.cpg_island_mixed_adj_pp-ctr,s=3.2,color=HEALTHY,alpha=.35,linewidths=0,marker='o')
    xs=np.array([20,80]); sl=C[(C['subset']==sp)&C.series.str.startswith('CpG island (')].iloc[0].slope_pp_per10y
    ax.plot(xs,sl*(xs-d.age.mean())/10,color=INK,lw=1.2)
    nl=C[(C['subset']==sp)&C.series.str.startswith('Non-CpG')].iloc[0].slope_pp_per10y
    ax.plot(xs,nl*(xs-d.age.mean())/10,color=PALE,lw=1.2,ls=(0,(3,2)))
    ax.set_xlabel('Age (years)'); ax.set_ylim(-2.0,2.5); ax.set_xlim(18,85)
    if k==0: ax.set_ylabel('CpG island strict-mixed,\nadjusted (%p, centred)')
    else: ax.set_yticklabels([])
    title(ax,f'{sp} (n = {len(d)})',pad=3)
    r=C[(C['subset']==sp)&C.series.str.startswith('CpG island (')].iloc[0]
    ax.text(.04,.96,f'{r.slope_pp_per10y:+.2f} %p / 10 y',transform=ax.transAxes,va='top',fontsize=5.5)
F.text(104,59,'Solid line, CpG island; dashed line, non-CpG island (both U-baseline markers)',5.5,color=MUTED)
# ---- c
F.label('c',1,68,'Strict-mixed slope with age'); ax=F.ax(34,80,46,30)
sers=['CpG island (U-baseline)','Non-CpG island (U-baseline)','CpG island minus non-CpG island']; cols=[HEALTHY,PALE,INK]
for si,(s_,c_) in enumerate(zip(sers,cols)):
    for ki,sp in enumerate(['Train','Test']):
        r=C[(C['subset']==sp)&(C.series==s_)].iloc[0]; yy=(2-si)+(.17 if ki==0 else -.17); mk='o' if sp=='Train' else 's'
        ax.plot([r.ci_lo,r.ci_hi],[yy,yy],color=c_,lw=1.2,solid_capstyle='butt'); ax.plot(r.slope_pp_per10y,yy,mk,color=c_,ms=3.4,mew=0)
        ax.text(1.03,yy,f'{r.slope_pp_per10y:+.2f} [{r.ci_lo:.2f}, {r.ci_hi:.2f}]'.replace('-','−'),transform=ax.get_yaxis_transform(),fontsize=5.5,va='center')
ax.axvline(0,color=MUTED,lw=1,ls=(0,(3,2))); ax.set_ylim(-.6,2.6); ax.set_yticks([2,1,0]); ax.set_yticklabels(['CpG island','Non-CpG\nisland','CpG island −\nnon-CpG island'],fontsize=5.5)
ax.tick_params(axis='y',length=0); ax.spines['left'].set_visible(False); ax.set_xlim(-.2,.3); ax.set_xlabel('Slope (%p per 10 years)')
ax.legend(handles=[Line2D([],[],marker='o',ls='',color=MUTED,ms=3.4,mew=0,label='Train'),Line2D([],[],marker='s',ls='',color=MUTED,ms=3.4,mew=0,label='Test')],loc='lower right',bbox_to_anchor=(1.0,1.0),ncol=2,fontsize=5.5,handletextpad=.1)
leg(F.fig,loc='upper left',bbox_to_anchor=(1/W,1-62.5/F.H))
# ---- d (redesigned 2026-10-15): per-stratum effects, each source in healthy Train SD units
F.label('d',1,120,'Three sources compared stratum by stratum (effect and 95% CI, healthy Train SD units)'); D=rd('Fig6d.tsv')
MET=['sOUF','BFV','strict-mixed']; BASES=('U','mid','M','All'); GAPB=5.0
pos={}; blocks={}; x0=0.
for m in MET:
    st_,xs_,sp_=strata_x(D.stratum.unique(),bases=BASES,gap=1.6,start=x0)
    for s_ in st_: pos[(m,s_)]=xs_[s_]
    blocks[m]=(x0,max(xs_.values()),sp_); x0=max(xs_.values())+GAPB
SRC=[('HCC − healthy','HCC − healthy\n(SD)',(-3,6.5)),('Hepatitis − healthy','Hepatitis −\nhealthy (SD)',(-2.7,2.3)),('Age (per 10 years, Test healthy)','Age\n(SD per 10 y)',(-.42,.36))]
COLD={'increase (CI > 0)':DIV['up_ci'],'decrease (CI < 0)':DIV['dn_ci'],'CI includes 0':PALE}
Y0,HH,GG=133,13.5,3.0
for i,(src,yl,yl_) in enumerate(SRC):
    ax=F.ax(24,Y0+i*(HH+GG),152,HH); g=D[D.source==src]
    for m in MET:
        gm=g[g.metric==m]; xx=np.array([pos[(m,s_)] for s_ in gm.stratum])
        for dirn,col in COLD.items():
            k=(gm.direction==dirn).to_numpy()
            ax.vlines(xx[k],gm.ci_lo_SD.to_numpy()[k],gm.ci_hi_SD.to_numpy()[k],color=col,lw=.55,zorder=2)
            if dirn=='CI includes 0': ax.scatter(xx[k],gm.effect_SD.to_numpy()[k],s=4,facecolor='white',edgecolor=MUTED,lw=.45,zorder=3)
            else: ax.scatter(xx[k],gm.effect_SD.to_numpy()[k],s=4,color=col,lw=0,zorder=3)
        a0,a1,sp_=blocks[m]
        for b,(b0,b1) in sp_.items():
            if b!=BASES[0]: ax.axvline(b0-.8,color=GRID,lw=.5,zorder=0)
        nu=int((gm.direction=='increase (CI > 0)').sum()); nd=int((gm.direction=='decrease (CI < 0)').sum())
        ax.text(a0,1.0,f'↑{nu}  ↓{nd}',transform=ax.get_xaxis_transform(),fontsize=5.5,va='bottom',ha='left',color=INK2)
        if i==0: ax.text((a0+a1)/2,1.0+ 7.0/HH,m,transform=ax.get_xaxis_transform(),fontsize=6.2,fontweight='bold',ha='center',va='bottom')
    ax.axhline(0,color=INK,lw=.8,zorder=1); ax.set_ylim(*yl_); ax.set_xlim(-1.5,x0-GAPB+1.5); ax.set_ylabel(yl,labelpad=1,fontsize=5.8)
    ax.spines['bottom'].set_visible(False); ax.set_xticks([]); 
    for m in MET[1:]:
        ax.axvline(blocks[m][0]-GAPB/2-.4,color=MUTED,lw=.8,zorder=0)
for m in MET:
    for b,(b0,b1) in blocks[m][2].items():
        ax.text((b0+b1)/2,-.04,BASE_HEAD[b].replace(' baseline','').replace('All markers','All'),transform=ax.get_xaxis_transform(),fontsize=5.5,ha='center',va='top',color=INK2)
F.text(24,Y0+3*(HH+GG)+2.5,'Healthy baseline of the strata (annotation order within each group as in Fig. 2a)',5.5,color=MUTED)
from matplotlib.lines import Line2D as _L2
hd=[_L2([],[],marker='o',ls='',ms=3,color=DIV['up_ci'],mew=0,label='Increase, CI excludes 0'),_L2([],[],marker='o',ls='',ms=3,color=DIV['dn_ci'],mew=0,label='Decrease, CI excludes 0'),
    _L2([],[],marker='o',ls='',ms=3,mfc='white',mec=MUTED,mew=.6,label='CI includes 0')]
F.fig.legend(handles=hd,loc='upper left',bbox_to_anchor=(122/W,1-72/F.H),fontsize=5.5,frameon=False,handletextpad=.3,labelspacing=.6)
F.text(122,68.5,'Key for panel d',5.8,'bold')
F.text(122,88,'Each row on its own scale, effect / healthy\nTrain SD of the stratum. HCC (n = 177) and\nhepatitis (n = 10) minus healthy (n = 1,009),\ngroup bootstrap 95% CI; age, Test healthy\nslope per 10 years (sex, depth adjusted).\n↑/↓, strata with CI excluding 0.',5.5,color=MUTED)
F.save(6)
