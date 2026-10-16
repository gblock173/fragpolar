"""BAM -> cpg_calls(fragpolar/vendor/extract_cpg_calls.py, 혈장과 동일 인자) -> suff(fragpolar.extract.raw_extract, 동일 규칙).
usage: pb1_extract.py --id ID --bam BAM --outdir DIR [--keep-slim]
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, os, json, time, argparse, subprocess, resource
from pathlib import Path
import pandas as pd
VEN = _pl.Path(__file__).resolve().parents[2] / 'src' / 'fragpolar' / 'vendor'
MT = CFG.path('marker_table_tsv')
REF = CFG.path('reference_fasta')
ap = argparse.ArgumentParser(); ap.add_argument('--id', required=True); ap.add_argument('--bam', required=True)
ap.add_argument('--outdir', required=True); ap.add_argument('--slim', action='store_true', help='cpg_calls를 분석 필요 7열로 재저장(zstd)')
a = ap.parse_args(); od = Path(a.outdir); (od / 'cpg_calls').mkdir(parents=True, exist_ok=True); (od / 'suff').mkdir(exist_ok=True)
cp = od / 'cpg_calls' / f'{a.id}.cpg_calls.parquet'; sp = od / 'suff' / f'{a.id}.parquet'; st = time.time()
if not cp.exists():
    tmp = od / 'cpg_calls' / f'{a.id}.cpg_calls.tmp.parquet'
    subprocess.run([sys.executable, '-u', str(VEN / 'extract_cpg_calls.py'), '--bam', a.bam, '--ref', REF, '--markers', MT, '--out', str(tmp), '--threads', '4'], check=True, cwd=VEN)
    os.replace(tmp, cp)
t_ext = time.time() - st
from fragpolar.extract import raw_extract, check_metrics
u = pd.read_parquet(CFG.inp('marker_regions')); u = u[u.primary_eligible]
regions = pd.read_parquet(CFG.inp('extraction_regions')); regions = regions[regions.marker_id.isin(set(u.marker_id))]
d, qa = raw_extract(str(cp), regions); qa.update(check_metrics(d))
tmp = sp.with_suffix('.tmp'); d.to_parquet(tmp, index=False, compression='zstd'); os.replace(tmp, sp)
rec = dict(id=a.id, bam=a.bam, seconds_extract=t_ext, seconds_total=time.time() - st, cpg_bytes=cp.stat().st_size, markers=len(d), QA=qa,
           job=os.getenv('JOB_ID'),
           peak_RSS_GiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20, peak_RSS_children_GiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 2**20)
(od / 'suff' / f'{a.id}.json').write_text(json.dumps(rec, indent=2, default=str)); print(json.dumps(rec, default=str), flush=True)
