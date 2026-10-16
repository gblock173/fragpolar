"""strict 통계(n_strictM/n_strictU)가 없는 표본에 대해 raw cpg_calls를 한 번 순회하여 fragpolar.engine(= one-pass engine).aggregate 와 동일 규칙으로
K>=4 fragment 합계(SUMS)를 다시 만든다. 기존 코드 import 만 하고 수정하지 않는다. gold/panel 은 계산하지 않음. 결과: results/strict_recompute/{sid}.parquet
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, json, time
import numpy as np, pandas as pd, pyarrow.parquet as pq
from fragpolar import engine as E
sid = sys.argv[1]; t0 = time.time()
OUT = CFG.work('tumour_hepatitis_age', 'results', 'strict_recompute')
CC = CFG.cpg_calls(sid)
u = E.marker_table(); ids = u.marker_id.to_numpy(); nm = len(u); pos_of = pd.Index(ids); mkcode = lambda s: pos_of.get_indexer(s.to_numpy())
mstart = u.start.to_numpy(np.int64); mend = u.end.to_numpy(np.int64)
acc = np.zeros((nm, len(E.SUMS))); kh = np.zeros(64); sh = np.zeros(100)
pf = pq.ParquetFile(CC); nrg = pf.metadata.num_row_groups; carry = None
def process(d):
    d = d[d.marker_id.isin(pos_of)]
    if not len(d): return
    c, mk, fr, pos, call, _, _ = E.clean_sort(d, mkcode)
    span = (c.fragment_end.to_numpy() - c.fragment_start.to_numpy()).astype(np.int64); comp = c.is_complete_pair.to_numpy(bool)
    inreg = (pos >= mstart[mk]) & (pos < mend[mk])
    global acc; acc += E.aggregate(mk, fr, pos, call, span, comp, inreg, nm, kh, sh)
for j in range(nrg):
    d = pf.read_row_group(j, columns=E.FIELDS, use_threads=False).to_pandas()
    if carry is not None: d = pd.concat([carry, d], ignore_index=True)
    if not len(d): carry = None; continue
    last = d.marker_id.iloc[-1]; end = d.marker_id.eq(last).to_numpy(); carry = d[end]; d = d[~end]
    if len(d): process(d)
    if j % 50 == 0: print(j, nrg, round(time.time() - t0), flush=True)
if carry is not None and len(carry): process(carry)
s = pd.DataFrame(acc, columns=E.SUMS)[['n_reads', 'sum_mean_meth', 'sum_mean_meth_sq', 'n_strictM', 'n_strictU']]; s.insert(0, 'marker_id', ids); s.to_parquet(f'{OUT}/{sid}.parquet', index=False, compression='zstd')
print('done', sid, time.time() - t0, flush=True)
