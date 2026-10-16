"""core suff가 legacy(재사용)라 strictM/strictU가 NaN인 혈장 Train 422개를, 기존 cpg_calls에서 동일 규칙(fragpolar.extract.raw_extract)으로 재집계.
usage: pb2_plasma_strict.py CHUNK NCHUNK
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, os, json, time, pandas as pd
from fragpolar.extract import raw_extract, check_metrics
PB_DIR = CFG.work('pbmc')
ch, nch = int(sys.argv[1]), int(sys.argv[2])
miss = pd.read_csv(PB_DIR + '/data/plasma_missing_strict.tsv', sep='\t').iloc[ch::nch]
u = pd.read_parquet(CFG.inp('marker_regions')); u = u[u.primary_eligible]
regions = pd.read_parquet(CFG.inp('extraction_regions')); regions = regions[regions.marker_id.isin(set(u.marker_id))]
for sid in miss.AnalysisID:
    out = f'{PB_DIR}/data/plasma_strict_suff/{sid}.parquet'
    if os.path.exists(out): continue
    t = time.time(); p = CFG.cpg_calls(sid)
    d, qa = raw_extract(p, regions); qa.update(check_metrics(d)); d.to_parquet(out + '.tmp', index=False, compression='zstd'); os.replace(out + '.tmp', out)
    print(json.dumps(dict(id=sid, sec=time.time() - t, markers=len(d))), flush=True)
