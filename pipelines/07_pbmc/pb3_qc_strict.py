"""재집계 suff vs core 행렬 mu/U/B/n 일치 확인(재집계 표본 중 일부)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import numpy as np, pandas as pd, glob, os, json
D=CFG.DATA; PB_DIR=CFG.work('pbmc')
cs=pd.read_csv(CFG.inp('marker_matrix_samples'),sep='\t',low_memory=False,usecols=['AnalysisID']); cs['col']=np.arange(len(cs)); cs=cs.set_index('AnalysisID')
st=pd.read_parquet(CFG.inp('marker_annotation')); mm={f:np.memmap(f'{CFG.inp("marker_matrices")}/{f}.f32',dtype=np.float32,mode='r',shape=(2044,190218)) for f in ['n','mu','B']}
rows=[]
for f in sorted(glob.glob(PB_DIR+'/data/plasma_strict_suff/*.parquet'))[::6]:
    sid=os.path.basename(f)[:-8]; s=pd.read_parquet(f).set_index('marker_id').reindex(st.marker_id); n=s.n_reads.to_numpy(float)
    mu=s.sum_mean_meth.to_numpy()/n; B=np.maximum(s.sum_mean_meth_sq.to_numpy()/n-mu**2,0); i=cs.loc[sid,'col']
    cn,cmu,cB=[np.asarray(mm[k][i],float) for k in ['n','mu','B']]; ok=np.isfinite(mu)&np.isfinite(cmu)
    rows.append(dict(sid=sid,n_compared=int(ok.sum()),max_abs_dn=float(np.nanmax(abs(np.nan_to_num(n)-np.nan_to_num(cn)))),max_abs_dmu=float(np.max(abs(mu[ok]-cmu[ok]))),mean_abs_dmu=float(np.mean(abs(mu[ok]-cmu[ok]))),max_abs_dB=float(np.nanmax(abs(B[ok]-cB[ok])))))
d=pd.DataFrame(rows); d.to_csv(PB_DIR+'/results/plasma_recompute_vs_core_qc.tsv',sep='\t',index=False); print(d.describe().to_string())
