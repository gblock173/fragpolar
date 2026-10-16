"""tumour-fraction step: 종양분율(f, tumour-fraction table, HCC 46명) 용량-반응. 사전 지정 층/지표 (결과 보기 전 고정).
차이 d = 표본값 - healthy Train(673) 중앙값 (주), 보조: healthy Train에서 적합한 y~age+sex+log_depth 기대값 대비 잔차.
Spearman rho(f,d) + 단순 선형 기울기 + f=0 절편; HCC 표본 bootstrap 5000 (seed 고정). Train 28 / Test 18 / 전체 46 모두 보고.
"""
import numpy as np, pandas as pd
from scipy import stats
from tha_common import *
W, nm, S = load()
rng = np.random.default_rng(SEED); NB = 5000
STRATA = ['All|All', 'U|All', 'M|All', 'mid|All', 'U|Island', 'All|Liver-Hep', 'All|Immune5']
ENDP = ['opp', 'mixed', 'B', 'BV', 'U', 'mu']
pa = pd.read_csv(CFG.inp('tumour_fraction_table'), sep='\t').set_index('AnalysisID')
hc = S[(S.group == 'healthy') & (S.split == 'Train')].index
hcc = S[S.group == 'liver']
f_ids = [i for i in pa.index if i in hcc.index]; assert len(f_ids) == 46
nof = [i for i in hcc.index if i not in f_ids]; assert len(nof) == 131
fsplit = pa.loc[f_ids, 'split']; assert (fsplit.values == S.loc[f_ids, 'split'].values).all()
cov = lambda ids: np.column_stack([S.loc[ids, 'age'] / 10, S.loc[ids, 'sex'] - 1, S.loc[ids, 'log_depth']])

def spear(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    return stats.spearmanr(x[ok], y[ok])[0] if ok.sum() > 4 else np.nan

def slope_icpt(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 5: return np.nan, np.nan
    b = np.polyfit(x[ok], y[ok], 1); return b[0], b[1]

rows = []; dist = []
for s in STRATA:
    for e in ENDP:
        c = (s, e); yh = W.loc[hc, c].to_numpy(float); med = np.nanmedian(yh); sd = np.nanstd(yh, ddof=1)
        # 공변량 기대값 모형 (healthy Train)
        Xh = cov(hc); beta, _ = ols(Xh, yh)
        exp = lambda ids: beta[0] + cov(ids) @ beta[1:]
        yfull = W.loc[f_ids, c].to_numpy(float); f = pa.loc[f_ids, 'f'].to_numpy(float)
        for dtype in ['median', 'covadj']:
            d = yfull - (med if dtype == 'median' else exp(f_ids))
            for sub, mask in [('All46', np.ones(46, bool)), ('Train28', (fsplit == 'Train').to_numpy()), ('Test18', (fsplit == 'Test').to_numpy())]:
                x = f[mask]; yy = d[mask]; ok = np.isfinite(yy); n = int(ok.sum())
                rho = spear(x, yy); sl, ic = slope_icpt(x, yy)
                bi = boot_idx(len(x), NB, rng); r_b = np.empty(NB); s_b = np.empty(NB); i_b = np.empty(NB)
                for k in range(NB):
                    r_b[k] = spear(x[bi[k]], yy[bi[k]]); s_b[k], i_b[k] = slope_icpt(x[bi[k]], yy[bi[k]])
                q = lambda v: tuple(np.nanpercentile(v, [2.5, 97.5]))
                rows.append(dict(stratum=s, endpoint=e, ref=dtype, subset=sub, n=n, healthy_ref=med if dtype == 'median' else np.nan, healthy_SD=sd, rho=rho, rho_lo=q(r_b)[0], rho_hi=q(r_b)[1],
                                 slope_per_unit_f=sl, slope_lo=q(s_b)[0], slope_hi=q(s_b)[1], intercept_f0=ic, icpt_lo=q(i_b)[0], icpt_hi=q(i_b)[1], f_range=f"{x.min():.3f}-{x.max():.3f}",
                                 slope_SD_per_0p1f=sl * 0.1 / sd, icpt_SD=ic / sd))
        # 보조: 공변량 보정 선형 기울기 (d_median ~ f + age + sex + depth), All46
        d = yfull - med; Xa = np.column_stack([f, cov(f_ids)]); b, _ = ols(Xa, d); bi = boot_idx(46, 2000, rng)
        bb = np.array([ols(Xa[i], d[i])[0][1] for i in bi])
        rows.append(dict(stratum=s, endpoint=e, ref='median', subset='All46_adjCov', n=46, slope_per_unit_f=b[1], slope_lo=np.percentile(bb, 2.5), slope_hi=np.percentile(bb, 97.5), healthy_SD=sd, healthy_ref=med, slope_SD_per_0p1f=b[1] * .1 / sd))
        # f 없는 131명 및 f 있는 46명, 전체 HCC 177, healthy Test 분포 (median 기준 차이)
        for lab, ids in [('HCC_f46', f_ids), ('HCC_noF131', nof), ('HCC_all177', list(hcc.index)), ('healthyTrain673', list(hc)), ('healthyTest336', list(S[(S.group == 'healthy') & (S.split == 'Test')].index))]:
            y = W.loc[ids, c].to_numpy(float) - med; y = y[np.isfinite(y)]
            dist.append(dict(stratum=s, endpoint=e, group=lab, n=len(y), median_diff=np.median(y), q25=np.percentile(y, 25), q75=np.percentile(y, 75), mean_diff=y.mean(), SD_units_median=np.median(y) / sd))
        a = W.loc[f_ids, c].to_numpy(float); b2 = W.loc[nof, c].to_numpy(float)
        a = a[np.isfinite(a)]; b2 = b2[np.isfinite(b2)]
        dist.append(dict(stratum=s, endpoint=e, group='MWU_f46_vs_noF131', n=len(a), median_diff=np.median(a) - np.median(b2), mean_diff=stats.mannwhitneyu(a, b2).pvalue))
        # 비 CNV 대리: 전체 177에서 1-R2 (atlas 미설명 분율), Liver-Hep 비율과의 Spearman
        ids = list(hcc.index); y = W.loc[ids, c].to_numpy(float) - med
        for prox in ['one_minus_R2', 'Liver-Hep']:
            dist.append(dict(stratum=s, endpoint=e, group=f'rho_vs_{prox}_HCC177', n=int(np.isfinite(y).sum()), median_diff=spear(S.loc[ids, prox].to_numpy(float), y)))
pd.DataFrame(rows).to_csv(RES + '/tf_dose_response.tsv', sep='\t', index=False)
pd.DataFrame(dist).to_csv(RES + '/tf_group_distributions.tsv', sep='\t', index=False)
pd.DataFrame(dict(AnalysisID=f_ids, f=pa.loc[f_ids, 'f'].values, split=fsplit.values, label=pa.loc[f_ids, 'label'].values)).to_csv(RES + '/tf_f_samples.tsv', sep='\t', index=False)
r = pd.DataFrame(rows); print(r[(r.subset == 'All46') & (r.ref == 'median')][['stratum', 'endpoint', 'n', 'rho', 'rho_lo', 'rho_hi', 'slope_SD_per_0p1f', 'icpt_SD', 'icpt_lo', 'icpt_hi']].round(3).to_string())
