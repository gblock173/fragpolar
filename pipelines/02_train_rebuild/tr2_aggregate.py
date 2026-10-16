"""군 비교 단계 1단계: marker_matrices(표본 x marker, float32, 행=표본)에서 analysis-freeze 표본(healthy/HCC)만 읽어
표본 x 59층 집계량을 만든다. (기존 코드/결과 수정 없음, 읽기 전용)

층 정의: strata_definition.tsv(used) = band(marker annotation.band) x annotation (panel toolkit.stratum_masks와 동일 논리).
유효 marker: n>=5 (earlier analysis와 동일, 민감도 n>=10). 층 평균은 표본마다 유효 marker의 단순 평균(층 내 marker 수>=20 규칙은 2단계에서 적용).
방향: baseline<0.5 이면 건강 기준 반대 균일 = strictM, 아니면 strictU. baseline 변형 2개:
  full  = marker annotation.baseline_fullTrain (기존)
  trainH= 이 연구 Train healthy(673)의 marker별 평균 mu 로 새로 계산 (Test 정보 미사용)
bin 경계(K x fraglen 4x4 사분위)는 무작위 Train 표본 100개(고정 seed)에서만 정한다.
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, json, time
import numpy as np, pandas as pd
import scipy.sparse as sp
from multiprocessing import Pool

D = CFG.DATA; TR_DIR = CFG.work('train_rebuild'); OUT = TR_DIR + '/results'
NS, NM = 2044, 190218
MATS = ['mu', 'U', 'B', 'strictM', 'strictU', 'n', 'mean_K', 'mean_fraglen', 'mean_complete_span']
MM = {f: np.memmap(f'{CFG.inp("marker_matrices")}/{f}.f32', dtype=np.float32, mode='r', shape=(NS, NM)) for f in MATS}
MARGINS = [0.0, 0.1, 0.2, 0.3]
G = {}


# train-rebuild step: legacy 422 Train 표본은 raw 재생성 sufficient parquet에서 s32_marker-matrix_core.rederive와 동일 식으로 행을 만든다(float32 저장까지 동일).
LEG = set(open(TR_DIR + '/data/legacy422_ids.txt').read().split()); REB = {}; CUR = {}
SC = ['marker_id', 'n_reads', 'sum_ncpg', 'sum_fraglen', 'n_complete', 'sum_complete_span', 'sum_mean_meth', 'sum_mean_meth_sq', 'n_strictM', 'n_strictU']


def dv(x, y):
    return np.divide(x, y, out=np.full(np.shape(x), np.nan), where=np.asarray(y) > 0)


def rebuilt_rows(sid):
    if sid in REB: return REB[sid]
    d = pd.read_parquet(f'{TR_DIR}/data/suff/{sid}.parquet', columns=SC); assert d.marker_id.is_unique
    d = d.set_index('marker_id').reindex(G['strata_ids'])
    n = d.n_reads.to_numpy(float); Mm = d.sum_mean_meth.to_numpy(float); M2 = d.sum_mean_meth_sq.to_numpy(float)
    mu = dv(Mm, n); p2 = dv(M2, n); U = mu - p2; B = p2 - mu * mu
    assert not (np.any(U[np.isfinite(U)] < -1e-10) or np.any(B[np.isfinite(B)] < -1e-10))
    out = dict(mu=mu, U=np.maximum(U, 0), B=np.maximum(B, 0), strictM=dv(d.n_strictM.to_numpy(float), n), strictU=dv(d.n_strictU.to_numpy(float), n),
               n=n, mean_K=dv(d.sum_ncpg.to_numpy(float), n), mean_fraglen=dv(d.sum_fraglen.to_numpy(float), n),
               mean_complete_span=dv(d.sum_complete_span.to_numpy(float), d.n_complete.to_numpy(float)))
    for k in out:
        v = out[k].copy(); v[n < 1] = np.nan; out[k] = v.astype(np.float32)
    REB.clear(); REB[sid] = out
    return out


def row(f, i):
    sid = G['id_by_col'][i]
    if sid in LEG: return np.asarray(rebuilt_rows(sid)[f], dtype=np.float64)
    return np.asarray(MM[f][i], dtype=np.float64)


def setup(test=False):
    cs = pd.read_csv(CFG.inp('marker_matrix_samples'), sep='\t', low_memory=False, usecols=['AnalysisID', 'group', 'split'])
    ms = pd.read_csv(CFG.inp('samples_table'), sep='\t')
    ms = ms[ms.group.isin(['healthy', 'liver']) & ms.split.isin(['Train', 'Test'])]
    cs['col'] = np.arange(len(cs))
    sel = cs.merge(ms[['AnalysisID', 'age', 'sex', 'log_depth', 'cv_fold']], on='AnalysisID', how='inner', suffixes=('', '_ms')).sort_values('col').reset_index(drop=True)
    assert len(sel) == 1186, len(sel)
    st = pd.read_parquet(CFG.inp('marker_annotation'))
    sd = pd.read_csv(CFG.inp('strata_definition'), sep='\t'); sd = sd[sd.used]
    masks = []
    for r in sd.itertuples():
        bm = np.ones(NM, bool) if r.baseline == 'All' else (st.band == r.baseline).to_numpy()
        am = np.ones(NM, bool) if r.axis == 'all' else (st.annot == r.annotation).to_numpy() if r.axis == 'annot' else \
            (st.genomic_context == r.annotation).to_numpy() if r.axis == 'context' else st['atlas_' + r.annotation].to_numpy()
        masks.append(bm & am)
    masks = np.stack(masks)
    assert len(sd) == 59 and (masks.sum(1) == sd.n_markers.to_numpy()).all(), 'strata marker count mismatch'
    G.update(strata_ids=st.marker_id.to_numpy(), id_by_col=dict(zip(sel.col, sel.AnalysisID)), sel=sel, strata=sd.stratum.tolist(), S=sp.csr_matrix(masks.astype(np.float64)), b_full=st.baseline_fullTrain.to_numpy(np.float64))
    G['masks'] = masks


def pass1():
    """baseline_trainH, bin edges, 진단용 Test-healthy 평균 mu."""
    sel = G['sel']; t0 = time.time()
    trH = sel[(sel.group == 'healthy') & (sel.split == 'Train')].col.to_numpy()
    teH = sel[(sel.group == 'healthy') & (sel.split == 'Test')].col.to_numpy()
    def meanmu(cols):
        s = np.zeros(NM); c = np.zeros(NM)
        for i in cols:
            n = row('n', i); mu = row('mu', i); ok = (n >= 5) & np.isfinite(mu); s += np.where(ok, mu, 0); c += ok
        return np.where(c >= 20, s / np.maximum(c, 1), np.nan), c
    bT, cT = meanmu(trH); bE, cE = meanmu(teH)
    tr = sel[sel.split == 'Train'].col.to_numpy(); rng = np.random.default_rng(20261014)
    pick = rng.choice(tr, 100, replace=False); Ks = []; Ls = []
    for i in pick:
        n = row('n', i); K = row('mean_K', i); L = row('mean_fraglen', i); ok = (n >= 5) & np.isfinite(K) & np.isfinite(L)
        Ks.append(K[ok]); Ls.append(L[ok])
    Ks = np.concatenate(Ks); Ls = np.concatenate(Ls); q = [.25, .5, .75]
    ke = np.quantile(Ks, q); le = np.quantile(Ls, q)
    G.update(b_trainH=bT, kedge=ke, ledge=le)
    pd.DataFrame(dict(marker=np.arange(NM), baseline_fullTrain=G['b_full'], baseline_trainH=bT, n_trainH=cT, mean_mu_TestHealthy_diag=bE, n_TestH=cE)).to_parquet(f'{OUT}/marker_baselines.parquet')
    json.dump(dict(kedge=ke.tolist(), ledge=le.tolist(), bin_pick_samples=sel.loc[sel.col.isin(pick), 'AnalysisID'].tolist(), seconds=time.time() - t0), open(f'{OUT}/bin_edges.json', 'w'), indent=1)
    print('pass1 done', time.time() - t0, 'edges', ke, le, flush=True)


def work(i):
    n = row('n', i); mu = row('mu', i); U = row('U', i); B = row('B', i); sM = row('strictM', i); sU = row('strictU', i)
    K = row('mean_K', i); fl = row('mean_fraglen', i); span = row('mean_complete_span', i)
    fin = np.isfinite(mu) & np.isfinite(U) & np.isfinite(B)           # 핵심(strict 불필요)
    fins = fin & np.isfinite(sM) & np.isfinite(sU)                     # strict 통계가 있는 표본만 (Train 일부 표본은 suff에 n_strictM 결측)
    V = mu * (1 - mu); mixed = 1 - sM - sU
    groups = []   # (tag, mask, {name: values})
    def add(tag, mask, vals):
        groups.append((tag, mask, vals))
    for thr in (5, 10):
        vm = (n >= thr) & fin; tg = f'n{thr}'
        d = dict(mu=mu, U=U, B=B, V=V)
        if thr == 5:
            d.update(K=K, fl=fl)
        add(tg, vm, d)
        if thr == 5:
            add('span', vm & np.isfinite(span), dict(span=span))
            vr = vm & (V >= 0.01)
            add('rmk', vr, dict(Rmk=np.where(vr, B / np.where(V > 0, V, 1), 0)))
            vs = (n >= 5) & fins
            add('n5s', vs, dict(sM=sM, sU=sU, mixed=mixed))
    vm = (n >= 5) & fin; vs = (n >= 5) & fins
    for bn, b in (('full', G['b_full']), ('trH', G['b_trainH'])):
        opp_base = np.where(b < .5, sM, sU); conc_base = np.where(b < .5, sU, sM); p0 = np.where(b < .5, 1., 0.)
        for m in MARGINS:
            mk = vs & np.isfinite(b) & (np.abs(b - .5) >= m)
            f = opp_base
            with np.errstate(invalid='ignore', divide='ignore'):
                mu2 = (mu - f * p0) / (1 - f); m2 = (B + mu * mu - f * p0) / (1 - f)
                Bdel = np.maximum(m2 - mu2 ** 2, 0)
            mkd = mk & (n * (1 - f) >= 5) & (f < 1) & np.isfinite(Bdel)
            tag = f'{bn}_m{int(round(m*10))}'
            add(tag, mk, dict(opp=opp_base, conc=conc_base, B=B, U=U, V=V, mixed=mixed))
            add(tag + 'd', mkd, dict(Bdel=Bdel, Bsub=B))
    # bin (K x fraglen 4x4), baseline full, margin 0
    kb = np.digitize(K, G['kedge']); lb = np.digitize(fl, G['ledge']); bid = kb * 4 + lb
    opp0 = np.where(G['b_full'] < .5, sM, sU)
    okc = vm & np.isfinite(K) & np.isfinite(fl); oks = vs & np.isfinite(K) & np.isfinite(fl)
    for k in range(16):
        add(f'bin{k}', okc & (bid == k), dict(B=B, V=V, U=U, mu=mu))
        add(f'bins{k}', oks & (bid == k), dict(mixed=mixed, opp=opp0))
    cols = []; names = []
    for tag, mask, vals in groups:
        cols.append(mask.astype(np.float64)); names.append(f'{tag}|cnt')
        for nm, v in vals.items():
            cols.append(np.where(mask, v, 0.0)); names.append(f'{tag}|{nm}')
    X = np.column_stack(cols)
    X[~np.isfinite(X)] = 0.0
    out = G['S'] @ X          # 59 x ncols (합)
    return i, out.astype(np.float64), names


def main():
    test = len(sys.argv) > 1 and sys.argv[1] == 'test'
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    setup(); pass1()
    cols = G['sel'].col.to_numpy()
    if test:
        cols = cols[:3]
    t0 = time.time(); res = {}; names = None
    with Pool(nproc) as p:
        for k, (i, out, nm) in enumerate(p.imap_unordered(work, cols, chunksize=4)):
            res[i] = out; names = nm
            if (k + 1) % 100 == 0:
                print('done', k + 1, time.time() - t0, flush=True)
    arr = np.stack([res[i] for i in cols])      # 표본 x 59 x C  (합)
    tag = '_test' if test else ''
    np.save(f'{OUT}/agg_sums{tag}.npy', arr)
    json.dump(dict(names=names, strata=G['strata'], samples=G['sel'].AnalysisID.tolist()), open(f'{OUT}/agg_meta{tag}.json', 'w'))
    G['sel'].to_csv(f'{OUT}/agg_samples{tag}.tsv', sep='\t', index=False)
    print('saved', arr.shape, time.time() - t0, flush=True)


if __name__ == '__main__':
    main()
