"""Loyfer 2023 PAT(hg38) -> marker(190,218)별 fragment 상태 충분통계 (우리 정의와 동일: 마커 창 안 관측 CpG K>=4인 read 1개 = fragment 1개,
p = 메틸 C 수 / 관측 CpG 수, strictM = p==1, strictU = p==0).
PAT 한 줄 = (chr, 전역 CpG index, 패턴[C 메틸, T 비메틸, . 미관측], count). 마커 창 밖 CpG는 잘라내고(cpg_calls와 같음) 관측 CpG만 센다.
두 변형을 동시에 계산한다.
  full  : read 전체(PAT는 read pair가 합쳐져 있어 cfDNA보다 길 수 있음)
  t167  : read가 167 bp를 넘으면 read 안에서 무작위 167 bp 창(균등, seed 고정)만 남김 -> cfDNA 길이에 맞춘 민감도.
usage: co2_pat_extract.py PAT.gz OUT.parquet [--max-lines N]
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, time, json, argparse
import numpy as np, pandas as pd
from numba import njit

E = CFG.work('component_origin')
CPG = CFG.path('cpg_bed')
WIN = 167
NACC = 8  # n, sumK, sumP, sumP2, nM, nU, (K hist handled separately)


@njit(cache=True)
def kernel(starts, cnts, buf, offs, A, B, pos0, acc, acc_t, khist, khist_t, span_stats, win, seed):
    np.random.seed(seed)
    nm = len(A)
    for i in range(len(starts)):
        s = starts[i]
        o0 = offs[i]; L = offs[i + 1] - o0
        e = s + L - 1
        j = np.searchsorted(B, s)  # 첫 b >= s
        while j < nm and A[j] <= e:
            lo = max(A[j], s); hi = min(B[j], e)
            K = 0; m = 0; first = -1; last = -1
            for c in range(lo, hi + 1):
                ch = buf[o0 + c - s]
                if ch == 67 or ch == 84:  # C or T
                    K += 1
                    if ch == 67:
                        m += 1
                    if first < 0:
                        first = pos0[c]
                    last = pos0[c]
            cn = cnts[i]
            if K >= 4:
                p = m / K
                acc[j, 0] += cn; acc[j, 1] += cn * K; acc[j, 2] += cn * p; acc[j, 3] += cn * p * p
                if m == K: acc[j, 4] += cn
                if m == 0: acc[j, 5] += cn
                khist[min(K, 40)] += cn
            if K >= 4 and last - first >= win:
                span_stats[0] += cn
                for _ in range(cn):  # 복사본마다 독립 창
                    w = first + np.random.randint(0, last - first - win + 2)
                    K2 = 0; m2 = 0
                    for c in range(lo, hi + 1):
                        ch = buf[o0 + c - s]
                        if (ch == 67 or ch == 84) and pos0[c] >= w and pos0[c] < w + win:
                            K2 += 1
                            if ch == 67:
                                m2 += 1
                    if K2 >= 4:
                        p = m2 / K2
                        acc_t[j, 0] += 1; acc_t[j, 1] += K2; acc_t[j, 2] += p; acc_t[j, 3] += p * p
                        if m2 == K2: acc_t[j, 4] += 1
                        if m2 == 0: acc_t[j, 5] += 1
                        khist_t[min(K2, 40)] += 1
            elif K >= 4:
                p = m / K
                acc_t[j, 0] += cn; acc_t[j, 1] += cn * K; acc_t[j, 2] += cn * p; acc_t[j, 3] += cn * p * p
                if m == K: acc_t[j, 4] += cn
                if m == 0: acc_t[j, 5] += cn
                khist_t[min(K, 40)] += cn
            j += 1


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('pat'); ap.add_argument('out'); ap.add_argument('--max-lines', type=int, default=0)
    a = ap.parse_args(); t0 = time.time()
    mr = pd.read_parquet(E + '/data/marker_cpg_ranges.parquet').sort_values('a').reset_index(drop=True)
    A = mr.a.to_numpy(np.int64); B = mr.b.to_numpy(np.int64); bmax = int(B.max())
    cg = pd.read_csv(CPG, sep='\t', header=None, names=['chrom', 'p1', 'idx'], usecols=['p1']); pos0 = np.r_[0, cg.p1.to_numpy(np.int64) - 1]  # pos0[idx]
    nm = len(A); acc = np.zeros((nm, 8)); acc_t = np.zeros((nm, 8)); kh = np.zeros(41, np.int64); kht = np.zeros(41, np.int64); ss = np.zeros(2, np.int64)
    nl = 0; nl_used = 0
    rd = pd.read_csv(a.pat, sep='\t', header=None, names=['chr', 'idx', 'pat', 'cnt'], usecols=[1, 2, 3], dtype={1: np.int64, 2: str, 3: np.int64}, chunksize=4_000_000)
    for ch in rd:
        nl += len(ch)
        if ch.idx.iloc[0] > bmax: break
        ch = ch[ch.idx <= bmax]
        L = ch.pat.str.len().to_numpy(np.int64); offs = np.r_[0, np.cumsum(L)]
        buf = np.frombuffer(''.join(ch.pat.tolist()).encode('ascii'), dtype=np.uint8)
        kernel(ch.idx.to_numpy(np.int64), ch.cnt.to_numpy(np.int64), buf, offs, A, B, pos0, acc, acc_t, kh, kht, ss, WIN, 20261014)
        nl_used += len(ch)
        if a.max_lines and nl >= a.max_lines: break
        print(json.dumps(dict(lines=nl, sec=round(time.time() - t0))), flush=True)
    cols = ['n_reads', 'sum_ncpg', 'sum_mean_meth', 'sum_mean_meth_sq', 'n_strictM', 'n_strictU']
    d = pd.DataFrame({'marker_row': mr.marker_row.to_numpy(), 'marker_id': mr.marker_id.to_numpy()})
    for k, c in enumerate(cols):
        d[c] = acc[:, k]; d[c + '_t'] = acc_t[:, k]
    d = d.sort_values('marker_row').reset_index(drop=True); d.to_parquet(a.out, index=False, compression='zstd')
    meta = dict(pat=a.pat, lines=nl, lines_used=nl_used, sec=time.time() - t0, reads_full=float(acc[:, 0].sum()), reads_t167=float(acc_t[:, 0].sum()),
                reads_span_gt166=int(ss[0]), K_hist_full=kh.tolist(), K_hist_t167=kht.tolist(), mean_K_full=float(acc[:, 1].sum() / max(acc[:, 0].sum(), 1)),
                mean_K_t167=float(acc_t[:, 1].sum() / max(acc_t[:, 0].sum(), 1)))
    open(a.out.replace('.parquet', '.json'), 'w').write(json.dumps(meta)); print(json.dumps({k: v for k, v in meta.items() if 'hist' not in k}))


if __name__ == '__main__':
    main()
