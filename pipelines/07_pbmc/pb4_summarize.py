"""표본 x 59층 요약값 (PBMC suff + 혈장 marker_matrices). 규칙은 tumour_hepatitis_age/tha1_extract.py와 동일:
유효 marker n>=5, 층 내 유효 marker>=20, marker 동일 가중 평균. opposite = baseline_fullTrain<0.5 -> strictM, >0.5 -> strictU (0.5는 제외).
strict-mixed = 1 - strictM - strictU. 추가: PBMC/혈장 공통 marker 민감도(common), mean_K, mean_fraglen.
mu/U/B는 core rederive와 동일: U=max(mu-p2,0), B=max(p2-mu^2,0).
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import numpy as np, pandas as pd, scipy.sparse as sp, time, os
D = CFG.DATA; PB_DIR = CFG.work('pbmc')
NS, NM = 2044, 190218
MM = {f: np.memmap(f'{CFG.inp("marker_matrices")}/{f}.f32', dtype=np.float32, mode='r', shape=(NS, NM)) for f in ['mu', 'U', 'B', 'strictM', 'strictU', 'n', 'mean_K', 'mean_fraglen']}
cs = pd.read_csv(CFG.inp('marker_matrix_samples'), sep='\t', low_memory=False, usecols=['AnalysisID']); cs['col'] = np.arange(len(cs))
ms = pd.read_csv(CFG.inp('samples_table'), sep='\t')
ms = ms[ms.group.isin(['healthy', 'liver']) & ms.split.isin(['Train', 'Test'])].merge(cs, on='AnalysisID', how='left')
assert len(ms) == 1186 and ms.col.notna().all()
st = pd.read_parquet(CFG.inp('marker_annotation'))
sd = pd.read_csv(CFG.inp('strata_definition'), sep='\t'); sd = sd[sd.used]
masks = []
for r in sd.itertuples():
    bm = np.ones(NM, bool) if r.baseline == 'All' else (st.band == r.baseline).to_numpy()
    am = np.ones(NM, bool) if r.axis == 'all' else (st.annot == r.annotation).to_numpy() if r.axis == 'annot' else \
        (st.genomic_context == r.annotation).to_numpy() if r.axis == 'context' else st['atlas_' + r.annotation].to_numpy()
    masks.append(bm & am)
names = sd.stratum.tolist(); masks = np.stack(masks); assert len(names) == 59 and (masks.sum(1) == sd.n_markers.to_numpy()).all()
S = sp.csr_matrix(masks.astype(np.float64))
b0 = st.baseline_fullTrain.to_numpy(np.float64); bsel = np.where(b0 < .5, 1, np.where(b0 > .5, 0, -1))

def plasma(i):
    g = lambda f: np.asarray(MM[f][i], dtype=np.float64)
    return g('n'), g('mu'), g('U'), g('B'), g('strictM'), g('strictU'), g('mean_K'), g('mean_fraglen')

def pbmc(path):
    s = pd.read_parquet(path).set_index('marker_id').reindex(st.marker_id); n = s.n_reads.to_numpy(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        mu = s.sum_mean_meth.to_numpy() / n; p2 = s.sum_mean_meth_sq.to_numpy() / n
        U = np.maximum(mu - p2, 0); B = np.maximum(p2 - mu * mu, 0)
        sM = s.n_strictM.to_numpy(float) / n; sU = s.n_strictU.to_numpy(float) / n
        K = s.sum_ncpg.to_numpy(float) / n; FL = s.sum_fraglen.to_numpy(float) / n
    n = np.where(n >= 1, n, np.nan)
    return n, mu, U, B, sM, sU, K, FL

def summarize(sid, arrs, common=None):
    n, mu, U, B, sM, sU, K, FL = arrs
    ok = np.isfinite(n) & (n >= 5) & np.isfinite(mu) & np.isfinite(U) & np.isfinite(B)
    if common is not None: ok &= common
    oks = ok & np.isfinite(sM) & np.isfinite(sU)
    opp = np.where(bsel == 1, sM, np.where(bsel == 0, sU, np.nan)); conc = np.where(bsel == 1, sU, np.where(bsel == 0, sM, np.nan)); oko = oks & np.isfinite(opp)
    out = []
    for name, v, vm in [('mu', mu, ok), ('U', U, ok), ('B', B, ok), ('V', U + B, ok), ('strictM', sM, oks), ('strictU', sU, oks), ('mixed', 1 - sM - sU, oks), ('opp', opp, oko), ('conc', conc, oko), ('meanK', K, ok), ('fraglen', FL, ok)]:
        v = np.nan_to_num(v); cnt = S @ vm.astype(np.float64); sm = S @ np.where(vm, v, 0.0)
        out.append(pd.DataFrame(dict(sid=sid, stratum=names, metric=name, value=np.where(cnt >= 20, sm / np.maximum(cnt, 1), np.nan), n_markers=cnt.astype(int))))
    return pd.concat(out)

if __name__ == '__main__':
    t0 = time.time()
    tk = pd.read_csv(f'{PB_DIR}/data/pbmc_tasks.tsv', sep='\t')
    P = {r.pbmc_lab_id: pbmc(f'{PB_DIR}/data/suff/{r.pbmc_lab_id}.parquet') for r in tk.itertuples()}
    # 공통 marker: PBMC 30 중 >=90% 에서 n>=5 이고 healthy Train 중 >=90% 에서 n>=5
    cov_p = np.mean([np.isfinite(a[0]) & (a[0] >= 5) for a in P.values()], axis=0)
    trh = ms[(ms.group == 'healthy') & (ms.split == 'Train')].col.to_numpy().astype(int)
    cov_h = np.zeros(NM)
    for i in trh: cov_h += (np.asarray(MM['n'][i], dtype=np.float64) >= 5)
    cov_h /= len(trh); common = (cov_p >= .9) & (cov_h >= .9)
    print('common markers', common.sum(), 'pbmc cov mean', cov_p.mean(), 'healthy cov mean', cov_h.mean(), flush=True)
    pd.DataFrame(dict(marker_id=st.marker_id, pbmc_cov=cov_p, healthy_train_cov=cov_h, common=common)).to_parquet(f'{PB_DIR}/results/marker_coverage.parquet', index=False)
    parts = []; partsC = []
    for k, (sid, a) in enumerate(P.items()):
        d = summarize(sid, a); d['kind'] = 'PBMC'; parts.append(d); d = summarize(sid, a, common); d['kind'] = 'PBMC'; partsC.append(d)
    for k, r in enumerate(ms.itertuples()):
        sp_ = f'{PB_DIR}/data/plasma_strict_suff/{r.AnalysisID}.parquet'
        a = pbmc(sp_) if os.path.exists(sp_) else plasma(int(r.col)); d = summarize(r.AnalysisID, a); d['kind'] = 'plasma'; parts.append(d)
        d = summarize(r.AnalysisID, a, common); d['kind'] = 'plasma'; partsC.append(d)
        if (k + 1) % 200 == 0: print(k + 1, time.time() - t0, flush=True)
    pd.concat(parts).to_parquet(f'{PB_DIR}/results/sample_strata_long.parquet', index=False)
    pd.concat(partsC).to_parquet(f'{PB_DIR}/results/sample_strata_long_commonmarkers.parquet', index=False)
    ms.assign(strict_from_recompute=[os.path.exists(f'{PB_DIR}/data/plasma_strict_suff/{a}.parquet') for a in ms.AnalysisID]).to_csv(f'{PB_DIR}/results/plasma_sample_table.tsv', sep='\t', index=False)
    print('done', time.time() - t0)
