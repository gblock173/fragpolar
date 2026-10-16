"""PBMC step 분석: 군 수준 + 효과 비교(PBMC-healthy vs HCC-healthy) + 연결 쌍 + 그림. 입력 results/sample_strata_long*.parquet
usage: pb5_analyze.py [suffix]   suffix='' (기본, 각 표본의 유효 marker) | 'common' (공통 marker 민감도)
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, json, numpy as np, pandas as pd
from scipy.stats import pearsonr, spearmanr, rankdata
PB_DIR = CFG.work('pbmc')
variant = sys.argv[1] if len(sys.argv) > 1 else ''
fn = 'sample_strata_long_commonmarkers' if variant == 'common' else 'sample_strata_long'
tag = '_commonmarkers' if variant == 'common' else ''
L = pd.read_parquet(f'{PB_DIR}/results/{fn}.parquet')
ms = pd.read_csv(f'{PB_DIR}/results/plasma_sample_table.tsv', sep='\t')
tk = pd.read_csv(f'{PB_DIR}/data/pbmc_tasks.tsv', sep='\t')
grp = {'PBMC': tk.pbmc_lab_id.tolist(),
       'H_Train': ms[(ms.group == 'healthy') & (ms.split == 'Train')].AnalysisID.tolist(),
       'H_Test': ms[(ms.group == 'healthy') & (ms.split == 'Test')].AnalysisID.tolist(),
       'HCC_Train': ms[(ms.group == 'liver') & (ms.split == 'Train')].AnalysisID.tolist(),
       'HCC_Test': ms[(ms.group == 'liver') & (ms.split == 'Test')].AnalysisID.tolist()}
assert [len(v) for v in grp.values()] == [30, 673, 336, 104, 73]
strata = pd.Index(L.stratum.unique()); METR = ['opp', 'mixed', 'strictM', 'strictU', 'B', 'mu', 'U', 'V', 'conc', 'meanK', 'fraglen']
X = {}   # metric -> DataFrame sample x stratum
for m in METR:
    X[m] = L[L.metric == m].pivot(index='sid', columns='stratum', values='value').reindex(columns=strata)
B = 2000; rng = np.random.default_rng(20261014)
idx = {g: [X['opp'].index.get_loc(s) for s in v] for g, v in grp.items()}
draws = {g: rng.integers(0, len(v), size=(B, len(v))) for g, v in grp.items()}
def gmean(arr, g, d=None):
    a = arr[idx[g]] if d is None else arr[np.array(idx[g])[d]]
    with np.errstate(all='ignore'): return np.nanmean(a, axis=-2) if d is None else np.nanmean(a, axis=-2)
def ci(x): return np.nanpercentile(x, [2.5, 97.5], axis=0)
CONTR = {'PBMC_minus_HTrain': ('PBMC', 'H_Train'), 'HCCTrain_minus_HTrain': ('HCC_Train', 'H_Train'),
         'HCCTest_minus_HTrain': ('HCC_Test', 'H_Train'), 'HCCTest_minus_HTest': ('HCC_Test', 'H_Test'), 'PBMC_minus_HCCTrain': ('PBMC', 'HCC_Train')}
grows, erows = [], []; boot = {}
for m in METR:
    A = X[m].to_numpy(float); gm = {}; gb = {}
    for g in grp:
        gm[g] = np.nanmean(A[idx[g]], axis=0)
        gb[g] = np.stack([np.nanmean(A[np.array(idx[g])[draws[g][b]]], axis=0) for b in range(B)])
        lo, hi = ci(gb[g]); nval = np.sum(np.isfinite(A[idx[g]]), axis=0)
        for j, s in enumerate(strata): grows.append(dict(metric=m, group=g, stratum=s, n_samples_valid=int(nval[j]), mean=gm[g][j], lo=lo[j], hi=hi[j]))
    for c, (a, b) in CONTR.items():
        e = gm[a] - gm[b]; eb = gb[a] - gb[b]; lo, hi = ci(eb); boot[(m, c)] = (e, eb)
        for j, s in enumerate(strata): erows.append(dict(metric=m, contrast=c, stratum=s, effect=e[j], lo=lo[j], hi=hi[j], ci_excludes0=bool(np.isfinite(lo[j]) and (lo[j] > 0 or hi[j] < 0))))
pd.DataFrame(grows).to_csv(f'{PB_DIR}/results/group_means{tag}.tsv', sep='\t', index=False)
pd.DataFrame(erows).to_csv(f'{PB_DIR}/results/effects{tag}.tsv', sep='\t', index=False)
# 층 수준 효과 비교
def corr(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    return pearsonr(x[ok], y[ok])[0], spearmanr(x[ok], y[ok])[0], int(ok.sum())
crow = []
for m in ['opp', 'mixed', 'strictM', 'strictU', 'B', 'mu', 'U']:
    for hcc in ['HCCTrain_minus_HTrain', 'HCCTest_minus_HTrain']:
        e1, eb1 = boot[(m, 'PBMC_minus_HTrain')]; e2, eb2 = boot[(m, hcc)]
        r, rho, n = corr(e1, e2); rb = np.array([corr(eb1[b], eb2[b])[:2] for b in range(B)])
        ok = np.isfinite(e1) & np.isfinite(e2); same = int(((e1 > 0) == (e2 > 0))[ok].sum())
        # 기울기(PBMC효과 = a + b * HCC효과): OLS
        sl = np.polyfit(e2[ok], e1[ok], 1)[0]
        slb = np.array([np.polyfit(eb2[b][ok], eb1[b][ok], 1)[0] for b in range(B)])
        crow.append(dict(metric=m, hcc_contrast=hcc, n_strata=n, pearson=r, pearson_lo=np.nanpercentile(rb[:, 0], 2.5), pearson_hi=np.nanpercentile(rb[:, 0], 97.5),
                         spearman=rho, spearman_lo=np.nanpercentile(rb[:, 1], 2.5), spearman_hi=np.nanpercentile(rb[:, 1], 97.5),
                         slope_PBMCeff_on_HCCeff=sl, slope_lo=np.nanpercentile(slb, 2.5), slope_hi=np.nanpercentile(slb, 97.5), n_same_sign=same))
pd.DataFrame(crow).to_csv(f'{PB_DIR}/results/effect_correlation{tag}.tsv', sep='\t', index=False)
# 요약: HCC에서 opposite가 유의 증가한 층에서 PBMC도 증가?
e = pd.DataFrame(erows); w = e.pivot_table(index=['metric', 'stratum'], columns='contrast', values=['effect', 'lo', 'hi'])
srow = []
for m in ['opp', 'mixed', 'B']:
    ww = w.loc[m]; hu = (ww[('lo', 'HCCTrain_minus_HTrain')] > 0); hd = (ww[('hi', 'HCCTrain_minus_HTrain')] < 0)
    pu = (ww[('lo', 'PBMC_minus_HTrain')] > 0); pd_ = (ww[('hi', 'PBMC_minus_HTrain')] < 0)
    srow.append(dict(metric=m, n_strata=len(ww), HCC_up_ci=int(hu.sum()), HCC_down_ci=int(hd.sum()), PBMC_up_ci=int(pu.sum()), PBMC_down_ci=int(pd_.sum()),
                     HCC_up_and_PBMC_up=int((hu & pu).sum()), HCC_up_and_PBMC_down=int((hu & pd_).sum()), HCC_up_and_PBMC_ns=int((hu & ~pu & ~pd_).sum())))
pd.DataFrame(srow).to_csv(f'{PB_DIR}/results/stratum_direction_counts{tag}.tsv', sep='\t', index=False)
# 연결 쌍
pairs = tk[tk.allowed_healthy_pair.astype(str).eq('True')].copy(); prow = []
for r in pairs.itertuples():
    for m in ['opp', 'mixed', 'strictM', 'strictU', 'B', 'mu', 'U']:
        for s in strata:
            ph = X[m].loc[r.pbmc_lab_id, s]; hp = X[m].loc[r.healthy_cfdna_lab_id, s] if r.healthy_cfdna_lab_id in X[m].index else np.nan
            prow.append(dict(pbmc=r.pbmc_lab_id, healthy_plasma=r.healthy_cfdna_lab_id, split=r.original_split, metric=m, stratum=s, pbmc_value=ph, plasma_value=hp, diff=ph - hp))
P = pd.DataFrame(prow); P.to_csv(f'{PB_DIR}/results/paired_diffs{tag}.tsv', sep='\t', index=False)
ps = P.groupby(['split', 'metric', 'stratum']).agg(n_pairs_with_plasma=('plasma_value', lambda x: int(x.notna().sum())), mean_diff=('diff', 'mean'), median_diff=('diff', 'median'), min_diff=('diff', 'min'), max_diff=('diff', 'max'), n_pos=('diff', lambda x: int((x > 0).sum()))).reset_index()
ps.to_csv(f'{PB_DIR}/results/paired_summary{tag}.tsv', sep='\t', index=False)
json.dump(dict(groups={k: len(v) for k, v in grp.items()}, variant=variant or 'per_sample_valid_markers', bootstrap=B), open(f'{PB_DIR}/results/analysis_meta{tag}.json', 'w'))
print('done')
