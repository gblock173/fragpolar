"""종양·간염·연령 단계 1단계: 표본 x 층(59 + 비-island 보조 2개) 요약값을 만든다. (읽기 전용 입력)
입력: marker_matrices(표본 x marker float32, 행=표본), 간염(sample_lists.hepatitis_included)은 marker 행렬에 없으므로 suff parquet에서 같은 공식으로 재유도.
유효 marker: n>=5 (earlier analysis와 동일), 층 내 유효 marker >= 20 (panel toolkit.stratum_means와 동일). marker 동일 가중 평균.
opposite: baseline_fullTrain<0.5 이면 strictM, >0.5 이면 strictU (=0.5 인 165개는 opposite 계산에서 제외).
strict-mixed = 1 - strictM - strictU (0<p<1). 0.3<p<0.7 부분 메틸은 marker 행렬/suff에 없다.
지표: opp, conc, strictM, strictU, mixed, B, U, mu, V(=U+B). B/mu(1-mu) 는 후처리에서 mean(B)/mean(V) 로 계산.
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, time, json, os
import numpy as np, pandas as pd, scipy.sparse as sp

D = CFG.DATA; OUT = CFG.work('tumour_hepatitis_age', 'results')
NS, NM = 2044, 190218
MM = {f: np.memmap(f'{CFG.inp("marker_matrices")}/{f}.f32', dtype=np.float32, mode='r', shape=(NS, NM)) for f in ['mu', 'U', 'B', 'strictM', 'strictU', 'n']}
test = len(sys.argv) > 1 and sys.argv[1] == 'test'

cs = pd.read_csv(CFG.inp('marker_matrix_samples'), sep='\t', low_memory=False, usecols=['AnalysisID'])
cs['col'] = np.arange(len(cs))
ms = pd.read_csv(CFG.inp('samples_table'), sep='\t')
HEP = set(CFG.ids('hepatitis_included'))   # 간염은 시료 품질 기준을 충족한 표본만 사용
ms = ms[ms.group.isin(['healthy', 'liver']) | ((ms.group == 'hepatitis') & ms.AnalysisID.isin(HEP))].copy()
hep = pd.read_csv(CFG.inp('hepatitis_samples_table'), sep='\t')[['AnalysisID', 'suff']]
ms = ms.merge(cs, on='AnalysisID', how='left').merge(hep, on='AnalysisID', how='left')
assert ms.col.notna().sum() == 1186 and ms.suff.notna().sum() == len(HEP), (ms.col.notna().sum(), ms.suff.notna().sum())
if test: ms = pd.concat([ms[ms.group == 'healthy'].head(2), ms[ms.group == 'liver'].head(1), ms[ms.group == 'hepatitis'].head(2)])

st = pd.read_parquet(CFG.inp('marker_annotation'))
sd = pd.read_csv(CFG.inp('strata_definition'), sep='\t'); sd = sd[sd.used]
masks, names = [], []
for r in sd.itertuples():
    bm = np.ones(NM, bool) if r.baseline == 'All' else (st.band == r.baseline).to_numpy()
    am = np.ones(NM, bool) if r.axis == 'all' else (st.annot == r.annotation).to_numpy() if r.axis == 'annot' else \
        (st.genomic_context == r.annotation).to_numpy() if r.axis == 'context' else st['atlas_' + r.annotation].to_numpy()
    masks.append(bm & am); names.append(r.stratum)
assert len(names) == 59 and (np.array([m.sum() for m in masks]) == sd.n_markers.to_numpy()).all()
for b in ['U', 'M', 'mid']:                                   # 보조: 비-island 합산 (연령 분석용)
    masks.append((st.band == b).to_numpy() & (st.annot != 'Island').to_numpy()); names.append(f'{b}|nonIsland')
masks.append((st.annot != 'Island').to_numpy()); names.append('All|nonIsland')
S = sp.csr_matrix(np.stack(masks).astype(np.float64))
b0 = st.baseline_fullTrain.to_numpy(np.float64)
bsel = np.where(b0 < 0.5, 1, np.where(b0 > 0.5, 0, -1))      # 1: opposite=strictM, 0: opposite=strictU

def rows_core(i):
    g = lambda f: np.asarray(MM[f][i], dtype=np.float64)
    return g('n'), g('mu'), g('U'), g('B'), g('strictM'), g('strictU')

def rows_suff(path):
    s = pd.read_parquet(path).set_index('marker_id').reindex(st.marker_id)
    n = s.n_reads.to_numpy(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        mu = s.sum_mean_meth.to_numpy() / n; p2 = s.sum_mean_meth_sq.to_numpy() / n
        U = np.maximum(mu - p2, 0); B = np.maximum(p2 - mu * mu, 0)
        sM = s.n_strictM.to_numpy(float) / n; sU = s.n_strictU.to_numpy(float) / n
    return n, mu, U, B, sM, sU

def summarize(sid, n, mu, U, B, sM, sU):
    ok = np.isfinite(n) & (n >= 5) & np.isfinite(mu) & np.isfinite(U) & np.isfinite(B)
    oks = ok & np.isfinite(sM) & np.isfinite(sU)
    mixed = 1 - sM - sU
    opp = np.where(bsel == 1, sM, np.where(bsel == 0, sU, np.nan)); conc = np.where(bsel == 1, sU, np.where(bsel == 0, sM, np.nan))
    oko = oks & np.isfinite(opp)
    out = []
    for name, v, vm in [('mu', mu, ok), ('U', U, ok), ('B', B, ok), ('V', U + B, ok), ('strictM', sM, oks), ('strictU', sU, oks), ('mixed', mixed, oks), ('opp', opp, oko), ('conc', conc, oko)]:
        cnt = S @ vm.astype(np.float64); sm = S @ np.where(vm, v, 0.0)
        val = np.where(cnt >= 20, sm / np.maximum(cnt, 1), np.nan)
        out.append(pd.DataFrame(dict(sid=sid, stratum=names, metric=name, value=val, n_markers=cnt.astype(int))))
    return pd.concat(out)

t0 = time.time(); parts = []; recomputed = []
for k, r in enumerate(ms.itertuples()):
    arrs = rows_core(int(r.col)) if np.isfinite(r.col) else rows_suff(r.suff)
    rp = f'{OUT}/strict_recompute/{r.AnalysisID}.parquet'
    if np.isfinite(r.col) and os.path.exists(rp) and not np.isfinite(arrs[4]).any():      # strict 결측 표본: raw 재계산 결과로 strict만 대체 (n·mu 패리티 확인)
        n2, mu2, _, _, sM2, sU2 = rows_suff(rp); ok = np.isfinite(arrs[0]) & (n2 > 0)
        assert np.nanmax(np.abs(n2[ok] - arrs[0][ok])) == 0 and np.nanmax(np.abs(mu2[ok] - arrs[1][ok])) < 1e-5, 'recompute parity fail ' + r.AnalysisID
        arrs = arrs[:4] + (sM2, sU2); recomputed.append(r.AnalysisID)
    parts.append(summarize(r.AnalysisID, *arrs))
    if (k + 1) % 100 == 0: print(k + 1, time.time() - t0, flush=True)
L = pd.concat(parts, ignore_index=True)
fn = 'sample_strata_long_TEST.parquet' if test else 'sample_strata_long.parquet'
L.to_parquet(f'{OUT}/{fn}', index=False)
if not test:
    ms.drop(columns=['suff']).to_csv(f'{OUT}/sample_table.tsv', sep='\t', index=False)
    # 패리티 점검용: 건강인 3명을 suff 경로로도 재계산
    chk = pd.read_csv(CFG.inp('marker_matrix_samples'), sep='\t', low_memory=False, usecols=['AnalysisID', 'suff']).merge(ms[['AnalysisID']], on='AnalysisID').head(3)
    pp = [summarize(r.AnalysisID, *rows_suff(r.suff)) for r in chk.itertuples()]
    pd.concat(pp).to_parquet(f'{OUT}/parity_suff_vs_core.parquet', index=False)
json.dump(recomputed, open(f'{OUT}/strict_recomputed_sids.json', 'w')); print('strict recomputed', len(recomputed)); print('done', len(L), time.time() - t0, flush=True)
