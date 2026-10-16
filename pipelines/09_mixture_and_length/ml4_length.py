"""length step: raw cpg_calls에서 표본별 fragment 길이 분포(1 bp 히스토그램)를 만든다. 읽기 전용.
- fragment_id 로 중복 제거(같은 fragment 의 여러 CpG 행, 여러 marker 행 -> 1개).
- 길이 = fragment_end - fragment_start (파이프라인의 span 정의와 동일).
- 4개 히스토그램: 전체 / 완전한 쌍(is_complete_pair)만 / K>=4(고유 CpG 위치 4개 이상) / 완전한 쌍 & K>=4.
  K 는 row group 안에서 센 값이다(row group 경계에서 약간 과소할 수 있음).
- 길이는 0..1999 bp, 2000 이상은 마지막 bin 에 합친다.
usage: python ml4_length.py <task_index> [limit_rowgroups]
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, time, os
import numpy as np, pandas as pd, pyarrow.parquet as pq
OUT = CFG.work('mixture_and_length')
t = int(sys.argv[1]); lim = int(sys.argv[2]) if len(sys.argv) > 2 else None
m = pd.read_csv(OUT + '/data/length_manifest.tsv', sep='\t'); r = m[m.task == t].iloc[0]
fn = f'{OUT}/results/length_hist/{r.sid}.npz'
if os.path.exists(fn) and lim is None: print('skip', r.sid); sys.exit(0)
t0 = time.time(); pf = pq.ParquetFile(r.path); nrg = pf.metadata.num_row_groups; nrg = nrg if lim is None else min(lim, nrg)
H, S, C, KK, rows = [], [], [], [], 0
for j in range(nrg):
    d = pf.read_row_group(j, columns=['fragment_id', 'pos', 'fragment_start', 'fragment_end', 'is_complete_pair'], use_threads=False).to_pandas()
    rows += len(d)
    d = d.drop_duplicates(['fragment_id', 'pos'])
    fi, uq = pd.factorize(d.fragment_id)
    K = np.bincount(fi); first = np.unique(fi, return_index=True)[1]
    sp = (d.fragment_end.to_numpy(np.int64) - d.fragment_start.to_numpy(np.int64))[first]; cp = d.is_complete_pair.to_numpy(bool)[first]
    H.append(pd.util.hash_pandas_object(pd.Series(uq), index=False).to_numpy(np.uint64)); S.append(sp.astype(np.int32)); C.append(cp); KK.append(K.astype(np.int16))
    if j % 40 == 0: print(r.sid, j, nrg, round(time.time() - t0), flush=True)
H = np.concatenate(H); S = np.concatenate(S); C = np.concatenate(C); KK = np.concatenate(KK)
_, ix = np.unique(H, return_index=True); S, C, KK = S[ix], C[ix], KK[ix]
NB = 2000
def hist(mask): return np.bincount(np.clip(S[mask], 0, NB - 1), minlength=NB).astype(np.int64)
allm = np.ones(len(S), bool)
np.savez_compressed(fn if lim is None else fn.replace('.npz', '_test.npz'), h_all=hist(allm), h_comp=hist(C), h_k4=hist(KK >= 4), h_comp_k4=hist(C & (KK >= 4)),
                    n_rows=rows, n_frag=len(S), n_frag_dup_removed_across_rg=int(len(H) - len(S)), n_neg_span=int((S < 0).sum()))
print('done', r.sid, len(S), round(time.time() - t0), flush=True)
