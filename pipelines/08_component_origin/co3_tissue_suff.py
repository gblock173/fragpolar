"""조직 cpg_calls(종양 26 / 등록 normal 27) -> marker 충분통계. PBMC step/혈장과 같은 fragpolar.extract.raw_extract(K>=4, 중복/충돌 처리) 재사용. usage: co3_tissue_suff.py SID"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, os, json, time, pandas as pd
from fragpolar.extract import raw_extract, check_metrics
E = CFG.work('component_origin')
sid = sys.argv[1]; tk = pd.read_csv(E + '/data/tissue_tasks.tsv', sep='\t', dtype=str, keep_default_na=False).set_index('sid')
u = pd.read_parquet(CFG.inp('marker_regions')); u = u[u.primary_eligible]
regions = pd.read_parquet(CFG.inp('extraction_regions')); regions = regions[regions.marker_id.isin(set(u.marker_id))]
out = f'{E}/data/tissue_suff/{sid}.parquet'; t = time.time()
d, qa = raw_extract(tk.loc[sid, 'path'], regions); qa.update(check_metrics(d)); d.to_parquet(out + '.tmp', index=False, compression='zstd'); os.replace(out + '.tmp', out)
print(json.dumps(dict(id=sid, sec=time.time() - t, markers=len(d), qa=qa), default=str), flush=True)
