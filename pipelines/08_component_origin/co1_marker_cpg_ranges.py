"""marker(190,218) -> wgbstools hg38 CpG index 범위 [a,b] (전역 1-based index). marker 창은 extract_cpg_calls와 같이 start <= pos0 < end (pos0 = CpG.bed 2열 - 1)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import numpy as np, pandas as pd, gzip
E=CFG.work('component_origin')
st=pd.read_parquet(CFG.inp('marker_annotation'))
cg=pd.read_csv(CFG.path('cpg_bed'),sep='\t',header=None,names=['chrom','p1','idx'])
assert (np.diff(cg.idx.to_numpy())==1).all() and cg.idx.iloc[0]==1
pos0=cg.p1.to_numpy()-1; ch=cg.chrom.to_numpy(); out=np.zeros((len(st),3),dtype=np.int64)
for c in st.chrom.unique():
    m=np.flatnonzero(ch==c); p=pos0[m]; ix=np.flatnonzero((st.chrom==c).to_numpy())
    a=np.searchsorted(p,st.start.to_numpy()[ix],'left'); b=np.searchsorted(p,st.end.to_numpy()[ix],'left')  # [a,b)
    out[ix,0]=m[0]+1+a; out[ix,1]=m[0]+1+b-1; out[ix,2]=b-a
d=pd.DataFrame(dict(marker_row=np.arange(len(st)),marker_id=st.marker_id,chrom=st.chrom,a=out[:,0],b=out[:,1],ncpg=out[:,2]))
d.to_parquet(E+'/data/marker_cpg_ranges.parquet',index=False)
print(d.ncpg.describe()); print('empty',(d.ncpg==0).sum())
s=d.sort_values('a'); print('order same as row order:',(s.marker_row.to_numpy()==np.arange(len(d))).all(),' overlaps:',(s.a.to_numpy()[1:]<=s.b.to_numpy()[:-1]).sum(), 'chroms',d.chrom.unique()[:30])
