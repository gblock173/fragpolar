"""Seed-set consistency of generator response: Train-seed (generator step) vs original Test-seed (earlier analysis) mean response, per generator and f (cosine per stratum, all-cell Pearson)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import numpy as np, pandas as pd
A=CFG.work('generators','results')+'/'
T=pd.read_csv(A+'gen_response_TrainSeed20.tsv',sep='\t'); O=CFG.unmap_features(pd.read_csv(CFG.inp('prior_perturbation_response'),sep='\t')); O=O[O.condition.str.match(r'^G[123]_')]
m=T.merge(O,on=['condition','stratum','feature'],suffixes=('_tr','_te')); rows=[]
for c,g in m.groupby('condition'):
    cs=[]
    for s,h in g.groupby('stratum'):
        a=h.mean_tr.to_numpy(); b=h.mean_te.to_numpy(); ok=np.isfinite(a)&np.isfinite(b)
        if ok.sum()>=4: cs.append(a[ok]@b[ok]/np.linalg.norm(a[ok])/np.linalg.norm(b[ok]))
    ok=g.mean_tr.notna()&g.mean_te.notna(); rows.append(dict(condition=c,median_cos_across_strata=np.median(cs),min_cos=np.min(cs),pearson_all_cells=np.corrcoef(g.mean_tr[ok],g.mean_te[ok])[0,1],median_abs_diff=float(np.nanmedian(np.abs(g.mean_tr-g.mean_te)))))
r=pd.DataFrame(rows); r.to_csv(A+'gen_seed_consistency_Train_vs_Test.tsv',sep='\t',index=False); print(r.round(4).to_string())
