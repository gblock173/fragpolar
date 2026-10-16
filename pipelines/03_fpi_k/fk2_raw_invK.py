"""FPI_K step: raw cpg_calls에서 marker별 정확한 E[1/K] (K>=4 fragment)와 1/mean_K를 비교하기 위한 합 계산.
fragpolar.engine(= one-pass engine, 수정 없음)의 marker_table/clean_sort/FIELDS와 같은 규칙(+row-group carry-over)을 그대로 쓴다.
출력: results/raw_invK_<sid>.npz  (marker별 n, sum_k, sum_invk, sum_p, sum_p2, K 히스토그램)
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, time, json
import numpy as np, pandas as pd, pyarrow.parquet as pq
from fragpolar import engine as E
sid = sys.argv[1]
path = CFG.cpg_calls(sid)
u = E.marker_table(); ids = u.marker_id.to_numpy(); nm = len(u); pos_of = pd.Index(ids)
mkcode = lambda s: pos_of.get_indexer(s.to_numpy())
acc = np.zeros((nm, 5)); kh = np.zeros(64); t0 = time.time()


def process(d):
    d = d[d.marker_id.isin(pos_of)]
    if not len(d): return
    c, mk, fr, pos, call, _, _ = E.clean_sort(d, mkcode)
    if not len(mk): return
    new = np.r_[True, (mk[1:] != mk[:-1]) | (fr[1:] != fr[:-1])]; idx = np.flatnonzero(new); rid = np.cumsum(new) - 1
    k = np.bincount(rid); m = np.bincount(rid, weights=call.astype(np.float64))
    keep = k >= 4; g = mk[idx][keep]; k = k[keep]; m = m[keep]; p = m / k
    kh[:] += np.bincount(np.minimum(k, 63), minlength=64)
    for j, v in enumerate([np.ones(len(k)), k, 1.0 / k, p, p * p]):
        acc[:, j] += np.bincount(g, weights=np.asarray(v, np.float64), minlength=nm)


pf = pq.ParquetFile(path); nrg = pf.metadata.num_row_groups; carry = None
for j in range(nrg):
    d = pf.read_row_group(j, columns=E.FIELDS, use_threads=False).to_pandas()
    if carry is not None: d = pd.concat([carry, d], ignore_index=True)
    if not len(d): carry = None; continue
    last = d.marker_id.iloc[-1]; end = d.marker_id.eq(last).to_numpy(); carry = d[end]; d = d[~end]
    if len(d): process(d)
    if j % 25 == 0: print(j, nrg, round(time.time() - t0), flush=True)
if carry is not None and len(carry): process(carry)
np.savez_compressed(CFG.work('fpi_k', 'results', f'raw_invK_{sid}.npz'), acc=acc, kh=kh, marker_id=ids)
print('done', sid, time.time() - t0)
