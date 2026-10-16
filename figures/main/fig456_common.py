"""Shared style for Fig 4-6 (earlier). Palette from earlier palettes (accessible variants), 180 mm width."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
W=180.; MM=1/25.4
INK='#1a1a1a'; MUTED='#4a5568'; GRID='#cbd5e0'; PALE='#a0aec0'
HEALTHY='#2C7DA0'; HCC='#C76B43'; HEP='#7B5EA7'; TUMOUR='#9b2c2c'; NONTUM='#56B4E9'; PBMC='#718096'
MK={'Healthy':'o','HCC':'^','Hepatitis':'D'}; COL={'Healthy':HEALTHY,'HCC':HCC,'Hepatitis':HEP}
plt.rcParams.update({'font.family':'Liberation Sans','font.size':5.8,'axes.titlesize':6.5,'axes.labelsize':6.0,'xtick.labelsize':5.5,'ytick.labelsize':5.5,
 'axes.linewidth':1.0,'axes.edgecolor':MUTED,'axes.labelcolor':INK,'text.color':INK,'xtick.color':MUTED,'ytick.color':MUTED,'axes.spines.top':False,'axes.spines.right':False,
 'xtick.major.width':1.0,'ytick.major.width':1.0,'xtick.major.size':2.,'ytick.major.size':2.,'lines.linewidth':1.0,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white','legend.frameon':False})
class Fig:
    def __init__(self,h): self.H=h; self.fig=plt.figure(figsize=(W*MM,h*MM),facecolor='white')
    def ax(self,x,top,w,h):
        a=self.fig.add_axes([x/W,(self.H-top-h)/self.H,w/W,h/self.H]); style(a); return a
    def text(self,x,top,s,size=5.8,weight='normal',**kw): return self.fig.text(x/W,1-top/self.H,s,fontsize=size,fontweight=weight,va='top',**kw)
    def label(self,l,x,top): self.text(x,top,l,8.5,'bold')
    def save(self,n):
        import os
        b=f'{CFG.FIG}/Figure_{n}'
        for e in ['pdf','svg']: self.fig.savefig(f'{b}.{e}')
        self.fig.savefig(f'{b}.png',dpi=300)
def style(ax):
    ax.minorticks_off(); ax.tick_params(direction='out',length=2,width=1,pad=1.5)
    for s in ax.spines.values(): s.set_linewidth(1); s.set_color(MUTED)
def title(ax,s,pad=3): ax.set_title(s,fontsize=6.3,fontweight='bold',pad=pad,loc='left')
def boxjitter(ax,groups,positions,colors,markers,rng,s=3.5,alpha=.4,widths=.5):
    vals=[np.asarray(v,float) for v in groups]
    ax.boxplot(vals,positions=positions,widths=widths,vert=True,patch_artist=True,showfliers=False,whis=1.5,medianprops={'color':INK,'linewidth':1.15},
        boxprops={'edgecolor':MUTED,'linewidth':1.,'facecolor':'none'},whiskerprops={'color':MUTED,'linewidth':1.},capprops={'color':MUTED,'linewidth':1.})
    for v,x,c,m in zip(vals,positions,colors,markers):
        j=rng.uniform(-.2,.2,len(v)); ax.scatter(x+j,v,s=s if len(v)>50 else s*2.2,color=c,marker=m,alpha=alpha if len(v)>50 else .85,linewidths=0,zorder=2)
