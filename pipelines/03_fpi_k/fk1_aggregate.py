"""FPI_K step: K 보정 FPI용 추가 집계. 군 비교 단계의 agg_sums.npy(ΣB, ΣV, cnt, ΣK, Σfl)는 재사용하고,
여기서는 h=1/mean_K 가 들어가는 합(Σ h·V 등)만 marker_matrices에서 새로 계산한다. 읽기 전용(기존 파일 수정 없음).
마스크는 군 비교 단계의 'n5'(n>=5 & mu,U,B 유한)와 동일하고, 집계 cnt가 agg_sums와 일치하는지 assert로 확인한다.
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, json, time
import numpy as np
A = CFG.WORK
sys.path.insert(0, CFG.pipeline_dir('group_effects'))
import ge1_aggregate as a1
from multiprocessing import Pool
OUT = A + '/fpi_k/results'
NAMES = ['n5|cnt', 'n5|hV', 'n5|h', 'n5|hV_floor', 'n10|cnt', 'n10|hV', 'rmk|cnt', 'rmk|FPIK_mk', 'rmk|h']


def work(i):
    row = a1.row
    n = row('n', i); mu = row('mu', i); U = row('U', i); B = row('B', i); K = row('mean_K', i)
    fin = np.isfinite(mu) & np.isfinite(U) & np.isfinite(B)
    with np.errstate(all='ignore'):
        V = mu * (1 - mu); h = 1.0 / K
        kok = np.isfinite(h) & (K >= 1)
        cols = []
        for thr in (5, 10):
            m = (n >= thr) & fin
            cols.append(m.astype(float))
            cols.append(np.where(m & kok, h * V, 0.0))
            if thr == 5:
                cols.append(np.where(m & kok, h, 0.0))
                cols.append((m & ~kok).astype(float))          # K 결측으로 h를 못 구한 marker 수 (0이어야 정상)
        m = (n >= 5) & fin & (V >= 0.01) & kok
        fk = (B - h * V) / ((1 - h) * V)
        cols.append(m.astype(float)); cols.append(np.where(m, fk, 0.0)); cols.append(np.where(m, h, 0.0))
    X = np.column_stack(cols); X[~np.isfinite(X)] = 0.0
    return i, (a1.G['S'] @ X).astype(np.float64)


if __name__ == '__main__':
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    a1.setup(); cols = a1.G['sel'].col.to_numpy(); t0 = time.time(); res = {}
    with Pool(nproc) as p:
        for k, (i, o) in enumerate(p.imap_unordered(work, cols, chunksize=4)):
            res[i] = o
            if (k + 1) % 100 == 0: print('done', k + 1, time.time() - t0, flush=True)
    arr = np.stack([res[i] for i in cols])
    # 컬럼 배치: [n5cnt, n5hV, n5h, n5_noK, n10cnt, n10hV, rmkcnt, rmkFPIK, rmkh]
    names = ['n5|cnt', 'n5|hV', 'n5|h', 'n5|noK', 'n10|cnt', 'n10|hV', 'rmkK|cnt', 'rmkK|FPIK_mk', 'rmkK|h']
    assert arr.shape[2] == len(names), arr.shape
    # 군 비교 단계 집계와 같은 marker 집합인지 확인
    old = np.load(A + '/group_effects/results/agg_sums.npy', mmap_mode='r'); om = json.load(open(A + '/group_effects/results/agg_meta.json'))
    oc = om['names'].index('n5|cnt'); assert om['samples'] == a1.G['sel'].AnalysisID.tolist()
    diff = np.abs(np.asarray(old[:, :, oc]) - arr[:, :, 0]).max(); print('max |cnt - agg_sums n5|cnt| =', diff, ' max noK =', arr[:, :, 3].max())
    assert diff == 0
    np.save(OUT + '/fpik_sums.npy', arr)
    json.dump(dict(names=names, strata=a1.G['strata'], samples=a1.G['sel'].AnalysisID.tolist(), max_cnt_diff=float(diff), max_noK=float(arr[:, :, 3].max())), open(OUT + '/fpik_meta.json', 'w'))
    print('saved', arr.shape, time.time() - t0)
