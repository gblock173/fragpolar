"""train-rebuild step producer 출력 vs 종양·간염·연령 단계 strict_recompute 출력 (공통 완료 표본) 일치 확인."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import os, numpy as np, pandas as pd
TR_DIR=CFG.WORK
cd=TR_DIR+'/tumour_hepatitis_age/results/strict_recompute'; me=TR_DIR+'/train_rebuild/data/suff'
ids=[f[:-8] for f in os.listdir(me) if f.endswith('.parquet') and os.path.exists(f'{cd}/{f}')]
rows=[]
for s in sorted(ids):
    a=pd.read_parquet(f'{me}/{s}.parquet').set_index('marker_id'); b=pd.read_parquet(f'{cd}/{s}.parquet').set_index('marker_id')
    j=a.index.intersection(b.index); r=dict(sample=s,n_me=len(a),n_cd=len(b),n_common=len(j))
    for c in ['n_reads','sum_mean_meth','sum_mean_meth_sq','n_strictM','n_strictU']:
        x=a.loc[j,c].to_numpy(float); y=b.loc[j,c].to_numpy(float); r[c+'_maxabs']=float(np.nanmax(np.abs(x-y))) if len(j) else np.nan
    rows.append(r)
T=pd.DataFrame(rows); T.to_csv(TR_DIR+'/train_rebuild/results/parity_vs_strict_recompute.tsv',sep='\t',index=False); print(T.to_string())
