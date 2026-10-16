"""ageing step: 연령 맥락 + 노화/HCC 효과의 층별 방향 비교.
모형(고정, Train 전에 지정): y ~ age/10 + (sex-1) + log_depth, 건강인 안에서 age/10 계수 = '10년당 기울기'.
Train에서 모형 확정 -> Test에 같은 식을 한 번 적용(재적합), 둘 다 보고. Train은 strict 통계가 있는 표본만 가능한 지표(opp/mixed/strictM/U)가 있어
  Train(지표별 가용 n) 과 Train_strictAvail(모든 지표 동일 표본) 둘을 제시한다.
HCC 효과: 같은 split 의 healthy+HCC 에서 y ~ HCC + age/10 + (sex-1) + log_depth 의 HCC 계수 (표본 분리 모형).
CI: 표본 복원추출 bootstrap (primary 5000, 전체 층 표 1000), seed 고정. 표준화 효과 = 계수 / 같은 split 건강인 y 의 SD.
"""
import sys, numpy as np, pandas as pd
from multiprocessing import Pool
from scipy import stats
from tha_common import *
W, nm, S = load()
STR59 = pd.read_csv(CFG.inp('strata_definition'), sep='\t').query('used').stratum.tolist()
EXTRA = ['U|nonIsland', 'M|nonIsland', 'mid|nonIsland', 'All|nonIsland']
ALLS = STR59 + EXTRA
has = W[('All|All', 'opp')].notna()
S['strict'] = has.reindex(S.index).fillna(False)

def design(ids, hcc=None):
    cols = [S.loc[ids, 'age'].to_numpy(float) / 10, S.loc[ids, 'sex'].to_numpy(float) - 1, S.loc[ids, 'log_depth'].to_numpy(float)]
    if hcc is not None: cols = [hcc] + cols
    return np.column_stack([np.ones(len(ids))] + cols)

def bols(X, y, B, seed, k):
    """point estimate + bootstrap of coefficient k (batched normal equations). NaN rows dropped."""
    ok = np.isfinite(y) & np.isfinite(X).all(1); X = X[ok]; y = y[ok]; n = len(y)
    if n < 30: return np.nan, np.nan, np.nan, n
    b0 = np.linalg.lstsq(X, y, rcond=None)[0][k]
    rng = np.random.default_rng(seed); idx = rng.integers(0, n, size=(B, n)); Xb = X[idx]; yb = y[idx]
    XtX = np.einsum('bnp,bnq->bpq', Xb, Xb); Xty = np.einsum('bnp,bn->bp', Xb, yb)
    bb = np.linalg.solve(XtX, Xty[..., None])[..., 0][:, k]
    return b0, np.percentile(bb, 2.5), np.percentile(bb, 97.5), n

def sets():
    g = S.group; sp = S.split
    return {'Train': S.index[(g == 'healthy') & (sp == 'Train')], 'Train_strictAvail': S.index[(g == 'healthy') & (sp == 'Train') & S.strict],
            'Test': S.index[(g == 'healthy') & (sp == 'Test')]}, {'Train': S.index[(g == 'liver') & (sp == 'Train')], 'Train_strictAvail': S.index[(g == 'liver') & (sp == 'Train') & S.strict],
            'Test': S.index[(g == 'liver') & (sp == 'Test')]}
HS, CS = sets()
HS['pooled'] = S.index[S.group == 'healthy']; CS['pooled'] = S.index[S.group == 'liver']

def job(args):
    col, B = args; out = []
    for sub in ['Train', 'Train_strictAvail', 'Test', 'pooled']:
        h = HS[sub]; c = CS[sub]
        y = W.loc[h, col].to_numpy(float); a = bols(design(h), y, B, SEED, 1)
        sd = np.nanstd(y, ddof=1)
        out.append(dict(kind='age_per10y', subset=sub, stratum=col[0], metric=col[1], coef=a[0], lo=a[1], hi=a[2], n=a[3], healthy_SD=sd))
        ids = list(h) + list(c); hc = np.r_[np.zeros(len(h)), np.ones(len(c))]
        y2 = W.loc[ids, col].to_numpy(float); b = bols(design(ids, hc), y2, B, SEED + 1, 1)
        out.append(dict(kind='HCC_vs_healthy_adj', subset=sub, stratum=col[0], metric=col[1], coef=b[0], lo=b[1], hi=b[2], n=b[3], healthy_SD=sd))
    return out

if __name__ == '__main__':
    cols = [(s, m) for s in ALLS for m in METRICS if (s, m) in W.columns]
    with Pool(8) as p: res = p.map(job, [(c, 1000) for c in cols], chunksize=4)
    T = pd.DataFrame([r for x in res for r in x]); T['z'] = T.coef / T.healthy_SD; T['zlo'] = T.lo / T.healthy_SD; T['zhi'] = T.hi / T.healthy_SD
    T.to_csv(RES + '/age_and_HCC_effects_all_strata.tsv', sep='\t', index=False)
    # ---- primary: U|Island vs U|nonIsland (strict-mixed), 표본 내 차이 Island - nonIsland ----
    prim = []
    comp = ['U|nonIsland', 'U|Shore', 'U|Enhancer', 'U|TSS', 'U|ReprPC', 'U|promoter', 'U|intergenic']
    for metric in ['mixed', 'B', 'U', 'opp', 'BV']:
        for ref in comp:
            for band in (['U'] if ref.startswith('U') else []):
                for sub in ['Train', 'Train_strictAvail', 'Test', 'pooled']:
                    h = HS[sub]; X = design(h)
                    yi = W.loc[h, ('U|Island', metric)].to_numpy(float); yr = W.loc[h, (ref, metric)].to_numpy(float)
                    for lab, y in [('Island', yi), ('comparator', yr), ('Island_minus_comparator', yi - yr)]:
                        r = bols(X, y, 5000, SEED + 7, 1)
                        prim.append(dict(metric=metric, comparator=ref, subset=sub, series=lab, coef_per10y=r[0], lo=r[1], hi=r[2], n=r[3], mean_level=float(np.nanmean(y)), healthy_SD=float(np.nanstd(y, ddof=1))))
    # 특이성 대조: M, mid 기저 Island vs nonIsland
    for band in ['M', 'mid']:
        for metric in ['mixed', 'B', 'U']:
            for sub in ['Train_strictAvail', 'Test']:
                h = HS[sub]; X = design(h)
                yi = W.loc[h, (f'{band}|Island', metric)].to_numpy(float); yr = W.loc[h, (f'{band}|nonIsland', metric)].to_numpy(float)
                for lab, y in [('Island', yi), ('comparator', yr), ('Island_minus_comparator', yi - yr)]:
                    r = bols(X, y, 5000, SEED + 7, 1)
                    prim.append(dict(metric=metric, comparator=f'{band}|nonIsland', subset=sub, series=lab, coef_per10y=r[0], lo=r[1], hi=r[2], n=r[3], mean_level=float(np.nanmean(y)), healthy_SD=float(np.nanstd(y, ddof=1))))
    # 0.3<p<0.7 부분 메틸 (U|Island 만 samples table.UIsland_mixed, 전 표본)
    for sub in ['Train', 'Train_strictAvail', 'Test', 'pooled']:
        h = HS[sub]; y = S.loc[h, 'UIsland_mixed'].to_numpy(float); r = bols(design(h), y, 5000, SEED + 7, 1)
        prim.append(dict(metric='partial_0.3-0.7(UIsland_mixed)', comparator='-', subset=sub, series='Island', coef_per10y=r[0], lo=r[1], hi=r[2], n=r[3], mean_level=float(np.nanmean(y)), healthy_SD=float(np.nanstd(y, ddof=1))))
    pd.DataFrame(prim).to_csv(RES + '/age_island_vs_nonisland.tsv', sep='\t', index=False)
    # UIsland_mixed(0.3<p<0.7) vs core strict-mixed(U|Island) 일치도
    h = HS['Train_strictAvail']; q = pd.DataFrame(dict(a=S.loc[h, 'UIsland_mixed'], b=W.loc[h, ('U|Island', 'mixed')])).dropna()
    pd.DataFrame([dict(n=len(q), spearman_partial_vs_strictmixed=stats.spearmanr(q.a, q.b)[0], mean_partial=q.a.mean(), mean_strictmixed=q.b.mean())]).to_csv(RES + '/age_partial_vs_strict_mixed_check.tsv', sep='\t', index=False)
    # ---- 방향 비교: 59층, 지표별 ----
    rows = []
    for sub in ['Train_strictAvail', 'Test', 'pooled']:
        for m in ['opp', 'mixed', 'B', 'U', 'mu', 'BV']:
            a = T[(T.kind == 'age_per10y') & (T.subset == sub) & (T.metric == m) & T.stratum.isin(STR59)].set_index('stratum')
            h = T[(T.kind == 'HCC_vs_healthy_adj') & (T.subset == sub) & (T.metric == m) & T.stratum.isin(STR59)].set_index('stratum')
            j = a.join(h, lsuffix='_age', rsuffix='_hcc', how='inner').dropna(subset=['z_age', 'z_hcc'])
            sa = lambda x, s: (x['lo' + s] > 0) | (x['hi' + s] < 0)
            sg_a = np.sign(j.coef_age); sg_h = np.sign(j.coef_hcc)
            for band_lab, mask in [('all59', np.ones(len(j), bool))] + [(f'baseline_{b}', j.index.str.startswith(b + '|')) for b in ['All', 'U', 'mid', 'M']]:
                jj = j[mask]
                if len(jj) < 4: continue
                rho = stats.spearmanr(jj.z_age, jj.z_hcc)[0]
                sga = np.sign(jj.coef_age); sgh = np.sign(jj.coef_hcc); sca = sa(jj, '_age'); sch = sa(jj, '_hcc')
                rows.append(dict(subset=sub, metric=m, strata_set=band_lab, n_strata=len(jj), age_up=int((sga > 0).sum()), age_down=int((sga < 0).sum()), hcc_up=int((sgh > 0).sum()), hcc_down=int((sgh < 0).sum()),
                                 age_up_CI=int((sca & (sga > 0)).sum()), age_down_CI=int((sca & (sga < 0)).sum()), hcc_up_CI=int((sch & (sgh > 0)).sum()), hcc_down_CI=int((sch & (sgh < 0)).sum()),
                                 same_dir=int((sga == sgh).sum()), opposite_dir=int((sga != sgh).sum()), same_dir_both_CI=int((sca & sch & (sga == sgh)).sum()), opp_dir_both_CI=int((sca & sch & (sga != sgh)).sum()),
                                 spearman_z_age_vs_z_hcc=rho, median_z_age=jj.z_age.median(), median_z_hcc=jj.z_hcc.median()))
            j.reset_index().to_csv(RES + f'/age_direction_by_stratum_{sub}_{m}.tsv', sep='\t', index=False)
    pd.DataFrame(rows).to_csv(RES + '/age_direction_summary.tsv', sep='\t', index=False)
    print(pd.DataFrame(rows).query("strata_set=='all59'").round(2).to_string())
    P = pd.DataFrame(prim); print(P[(P.metric == 'mixed') & (P.comparator == 'U|nonIsland')].round(4).to_string())
