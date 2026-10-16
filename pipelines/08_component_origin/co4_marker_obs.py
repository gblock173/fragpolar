"""관측 혈장 marker 수준 집단 평균과 종양분율(f) 기울기. 표본 값 규칙은 tumour_hepatitis_age/tha1_extract.py와 동일
(marker_matrices + strict 결측 Train은 strict_recompute로 대체, 간염(목록)은 suff, PBMC는 PBMC step suff). 유효 표본 값: n>=5.
출력: data/marker_obs_groups.parquet (marker_row x 집단 x 지표 평균, 유효 표본 수), data/marker_obs_fslope.parquet
지표: mu, U, B, p2(=mu-U), sM, sU. (opp/conc/mixed는 marker 평균 sM,sU에서 계산: 선형이라 평균과 교환 가능)"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import numpy as np, pandas as pd, os, time
D = CFG.DATA; AN = CFG.WORK; E = AN + '/component_origin'; THA_DIR = AN + '/tumour_hepatitis_age'; PB_DIR = AN + '/pbmc'
NS, NM = 2044, 190218
MM = {f: np.memmap(f'{CFG.inp("marker_matrices")}/{f}.f32', dtype=np.float32, mode='r', shape=(NS, NM)) for f in ['mu', 'U', 'B', 'strictM', 'strictU', 'n']}
st = pd.read_parquet(CFG.inp('marker_annotation'))
cs = pd.read_csv(CFG.inp('marker_matrix_samples'), sep='\t', low_memory=False, usecols=['AnalysisID']); cs['col'] = np.arange(len(cs))
ms = pd.read_csv(CFG.inp('samples_table'), sep='\t')
HEP_IDS = set(CFG.ids('hepatitis_included'))   # 간염은 시료 품질 기준을 충족한 표본만 사용
ms = ms[(ms.group.isin(['healthy', 'liver']) | ((ms.group == 'hepatitis') & ms.AnalysisID.isin(HEP_IDS))) & ms.split.isin(['Train', 'Test', 'Hepatitis'])].merge(cs, on='AnalysisID', how='left')
hep = pd.read_csv(CFG.inp('hepatitis_samples_table'), sep='\t')[['AnalysisID', 'suff']]; ms = ms.merge(hep, on='AnalysisID', how='left')
assert len(ms) == 1186 + len(HEP_IDS)
fs = pd.read_csv(f'{THA_DIR}/results/tf_f_samples.tsv', sep='\t').set_index('AnalysisID')

def from_suff(path):
    s = pd.read_parquet(path).set_index('marker_id').reindex(st.marker_id); n = s.n_reads.to_numpy(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        mu = s.sum_mean_meth.to_numpy() / n; p2 = s.sum_mean_meth_sq.to_numpy() / n
        sM = s.n_strictM.to_numpy(float) / n; sU = s.n_strictU.to_numpy(float) / n
    return n, mu, np.maximum(mu - p2, 0), np.maximum(p2 - mu * mu, 0), sM, sU

def sample_arrays(r):
    if np.isfinite(r.col):
        g = lambda f: np.asarray(MM[f][int(r.col)], dtype=np.float64)
        a = (g('n'), g('mu'), g('U'), g('B'), g('strictM'), g('strictU'))
        rp = f'{THA_DIR}/results/strict_recompute/{r.AnalysisID}.parquet'
        if os.path.exists(rp) and not np.isfinite(a[4]).any():
            n2, mu2, _, _, sM2, sU2 = from_suff(rp); a = a[:4] + (sM2, sU2)
        return a
    return from_suff(r.suff)

GR = {'H_Train': lambda r: r.group == 'healthy' and r.split == 'Train', 'H_Test': lambda r: r.group == 'healthy' and r.split == 'Test',
      'HCC_all': lambda r: r.group == 'liver', 'HCC_Train': lambda r: r.group == 'liver' and r.split == 'Train', 'HCC_Test': lambda r: r.group == 'liver' and r.split == 'Test',
      'HCC_f46': lambda r: r.group == 'liver' and r.AnalysisID in fs.index, 'HCC_nof131': lambda r: r.group == 'liver' and r.AnalysisID not in fs.index,
      'Hepatitis': lambda r: r.group == 'hepatitis' and r.AnalysisID in HEP_IDS}  # 간염: 시료 품질 기준을 충족한 표본만 (config sample_lists.hepatitis_included)
MET = ['mu', 'U', 'B', 'sM', 'sU']
acc = {g: {m: np.zeros(NM) for m in MET} for g in GR}; cnt = {g: np.zeros(NM) for g in GR}
fl = []; fy = []; t0 = time.time()
for k, r in enumerate(ms.itertuples()):
    n, mu, U, B, sM, sU = sample_arrays(r); ok = np.isfinite(n) & (n >= 5) & np.isfinite(mu) & np.isfinite(U) & np.isfinite(B)
    oks = ok & np.isfinite(sM) & np.isfinite(sU)
    arr = dict(mu=mu, U=U, B=B, sM=sM, sU=sU)
    for g, fn in GR.items():
        if fn(r):
            for m in MET:
                v = oks
                acc[g][m] += np.where(v, np.nan_to_num(arr[m]), 0)
            cnt[g] += oks  # 지표 모두 유효한 표본만 센다 (strict 포함)
            if g == 'H_Train': pass
    if r.group == 'liver' and r.AnalysisID in fs.index:
        fl.append(fs.loc[r.AnalysisID, 'f']); fy.append(np.stack([np.where(oks, arr[m], np.nan) for m in MET]))
    if (k + 1) % 200 == 0: print(k + 1, time.time() - t0, flush=True)
# PBMC (PBMC step)
tk = pd.read_csv(f'{PB_DIR}/data/pbmc_tasks.tsv', sep='\t'); acc['PBMC'] = {m: np.zeros(NM) for m in MET}; cnt['PBMC'] = np.zeros(NM)
for r in tk.itertuples():
    n, mu, U, B, sM, sU = from_suff(f'{PB_DIR}/data/suff/{r.pbmc_lab_id}.parquet'); ok = np.isfinite(n) & (n >= 5) & np.isfinite(mu) & np.isfinite(sM) & np.isfinite(sU)
    for m, a in zip(MET, [mu, U, B, sM, sU]): acc['PBMC'][m] += np.where(ok, np.nan_to_num(a), 0)
    cnt['PBMC'] += ok
rows = {}
for g in acc:
    for m in MET:
        with np.errstate(invalid='ignore', divide='ignore'): rows[f'{g}__{m}'] = np.where(cnt[g] > 0, acc[g][m] / cnt[g], np.nan)
    rows[f'{g}__nsamp'] = cnt[g]
pd.DataFrame(rows).assign(marker_row=np.arange(NM)).to_parquet(f'{E}/data/marker_obs_groups.parquet', index=False)
# f 기울기 (46명): 표본이 해당 marker에서 유효할 때만. y는 원 지표 값. 기울기 단위: 지표/f
F = np.array(fl); Y = np.stack(fy)  # (46, 5, NM)
out = {}
for i, m in enumerate(MET):
    y = Y[:, i, :]; w = np.isfinite(y); nv = w.sum(0)
    fb = (w * F[:, None]).sum(0) / np.maximum(nv, 1); fc = (F[:, None] - fb) * w; yb = np.where(w, y, 0).sum(0) / np.maximum(nv, 1)
    sl = (fc * np.where(w, y - yb, 0)).sum(0) / np.maximum((fc ** 2).sum(0), 1e-12)
    out[f'slope_{m}'] = np.where(nv >= 20, sl, np.nan)
out['nvalid'] = nv
pd.DataFrame(out).assign(marker_row=np.arange(NM)).to_parquet(f'{E}/data/marker_obs_fslope.parquet', index=False)
print('done', time.time() - t0, 'f n', len(F), 'groups', {g: int(cnt[g].max()) for g in cnt})
