"""length step: 고분자 DNA(길이) 지표와 균질화/양극화 성분의 관계, 보정 유지율, 같은 기관 비교. 입력은 length_hist/*.npz 와 종양·간염·연령 단계 표본x층 표."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import json, glob, os, sys
import numpy as np, pandas as pd
from scipy import stats
V = CFG.WORK; E = V + '/mixture_and_length'; RES = E + ('/results_dbg' if os.environ.get('LEN_DBG') else '/results')
rng = np.random.default_rng(20261014); NB = int(os.environ.get('LEN_NB', 1000))
man = pd.read_csv(E + '/data/length_manifest.tsv', sep='\t')
# ---- 길이 지표
rows = []
for r in man.itertuples():
    f = f'{E}/results/length_hist/{r.sid}.npz'
    if not os.path.exists(f): continue
    z = np.load(f); d = {'sid': r.sid, 'group': r.group, 'source': r.source}
    for tag, key in [('', 'h_all'), ('_k4', 'h_k4')]:
        h = z[key].astype(float); N = h.sum(); cs = np.cumsum(h); L = np.arange(len(h))
        d['n_frag' + tag] = N; d['frac300' + tag] = h[301:].sum() / N; d['frac400' + tag] = h[401:].sum() / N
        d['median_len' + tag] = L[np.searchsorted(cs, N / 2)]; d['frac166' + tag] = h[156:177].sum() / N
        d['frac_short120' + tag] = h[:120].sum() / N; d['mean_len' + tag] = (h * L).sum() / N
        d['frac_cap500' + tag] = h[500:].sum() / N
    d['frac_complete'] = z['h_comp'].sum() / z['h_all'].sum(); d['frac_k4'] = z['h_k4'].sum() / z['h_all'].sum()
    rows.append(d)
LM = pd.DataFrame(rows); LM.to_csv(RES + '/len_length_metrics_per_sample.tsv', sep='\t', index=False)
print('length metrics', len(LM), LM.group.value_counts().to_dict(), flush=True)
S = pd.read_csv(V + '/tumour_hepatitis_age/results/sample_table.tsv', sep='\t')
S['grp'] = S.group.replace({'liver': 'HCC'})
S = S[(S.group != 'hepatitis') | S.AnalysisID.isin(CFG.ids('hepatitis_included'))].copy()   # 간염은 목록의 표본만
Lg = pd.read_parquet(V + '/tumour_hepatitis_age/results/sample_strata_long.parquet')
W = Lg.pivot_table(index='sid', columns=['stratum', 'metric'], values='value')
strata = [s for s in W.columns.get_level_values(0).unique() if '|nonIsland' not in s]
COMP = ['mixed', 'strictM', 'opp', 'B', 'BV', 'U']
NAME = {'mixed': 'strict-mixed', 'strictM': 'strictM', 'opp': 'sOUF', 'B': 'BFV', 'BV': 'FPI', 'U': 'WFV'}
Y = {}
for s in strata:
    for c in COMP:
        Y[(s, c)] = (W[(s, 'B')] / W[(s, 'V')]) if c == 'BV' else W[(s, c)]
Y = pd.DataFrame(Y)
D = S.set_index('AnalysisID').join(LM.set_index('sid')[[c for c in LM.columns if c not in ('sid', 'group', 'source')]], how='left')
# PBMC-유사 gDNA QC 점수 (deconvolution QC): healthy 는 5-fold OOF, cancer 는 5-fold consensus
Q = CFG.path('pbmc_qc_dir') + '/'
a1 = pd.read_csv(Q + 'healthy_5fold_oof_pbmc_scores.tsv', sep='\t', usecols=['sample_id', 'alpha_blood6_raw', 'cosine_blood6'])
a2 = pd.read_csv(Q + 'cancer_5fold_consensus_pbmc_scores.tsv', sep='\t', usecols=['sample_id', 'alpha_blood6_raw', 'cosine_blood6'])
Al = pd.concat([a1, a2]).drop_duplicates('sample_id').set_index('sample_id').rename(columns={'alpha_blood6_raw': 'pbmc_alpha', 'cosine_blood6': 'pbmc_cos'})
D = D.join(Al, how='left')
print('pbmc alpha available', D.groupby('grp').pbmc_alpha.apply(lambda x: x.notna().sum()).to_dict(), flush=True)
D.to_csv(RES + '/len_sample_master.tsv', sep='\t')
MET = ['frac300', 'frac400', 'median_len', 'frac166', 'frac300_k4', 'pbmc_alpha']
MET_ADJ = ['frac300', 'median_len', 'frac166']

# ---- A. 군별 길이 지표 요약
out = []
for g, d in D.groupby('grp'):
    for m in ['frac300', 'frac400', 'frac_cap500', 'median_len', 'mean_len', 'frac166', 'frac_short120', 'frac300_k4', 'median_len_k4', 'frac_complete', 'pbmc_alpha', 'n_frag']:
        x = d[m].dropna()
        if len(x): out.append(dict(group=g, metric=m, n=len(x), median=x.median(), q25=x.quantile(.25), q75=x.quantile(.75), mean=x.mean(), sd=x.std()))
LS = pd.DataFrame(out); LS.to_csv(RES + '/len_length_group_summary.tsv', sep='\t', index=False)
# 군 간 비교 (HCC vs healthy 등), Cliff delta
def cliff(a, b):
    u = stats.mannwhitneyu(a, b); return 2 * u.statistic / (len(a) * len(b)) - 1, u.pvalue
cmp = []
H = D[D.grp == 'healthy']
for g in ['HCC', 'hepatitis', 'pbmc']:
    X = D[D.grp == g]
    for m in ['frac300', 'frac400', 'median_len', 'frac166', 'frac300_k4', 'pbmc_alpha']:
        a, b = X[m].dropna(), H[m].dropna()
        if len(a) > 2 and len(b) > 2:
            cd, p = cliff(a, b); cmp.append(dict(group=g, ref='healthy', metric=m, n=len(a), median_diff=a.median() - b.median(), mean_diff=a.mean() - b.mean(), smd=(a.mean() - b.mean()) / b.std(), cliff_delta=cd, p_mwu=p))
pd.DataFrame(cmp).to_csv(RES + '/len_length_group_compare.tsv', sep='\t', index=False)
# ---- B. 군 내 상관
def zs(x): return (x - x.mean()) / x.std()
def resid(y, X):
    X1 = np.column_stack([np.ones(len(y)), X]); return y - X1 @ np.linalg.lstsq(X1, y, rcond=None)[0]
cor = []
for gname, sub in [('healthy', D[D.grp == 'healthy']), ('HCC', D[D.grp == 'HCC']), ('healthy_Train', D[(D.grp == 'healthy') & (D.split == 'Train')]), ('healthy_Test', D[(D.grp == 'healthy') & (D.split == 'Test')]),
                   ('HCC_Train', D[(D.grp == 'HCC') & (D.split == 'Train')]), ('HCC_Test', D[(D.grp == 'HCC') & (D.split == 'Test')])]:
    for s in strata:
        for c in COMP:
            y = Y[(s, c)].reindex(sub.index)
            for m in MET:
                ok = y.notna() & sub[m].notna()
                if ok.sum() < 20: continue
                r_, p_ = stats.spearmanr(sub.loc[ok, m], y[ok])
                cv = sub.loc[ok, ['log_depth']].assign(age=sub.loc[ok, 'age'], sex=sub.loc[ok, 'sex']).dropna()
                pr = np.nan
                if len(cv) > 20:
                    k = cv.index; pr = stats.spearmanr(resid(stats.rankdata(sub.loc[k, m]), cv.apply(stats.rankdata).to_numpy()), resid(stats.rankdata(y[k]), cv.apply(stats.rankdata).to_numpy()))[0]
                cor.append(dict(subset=gname, stratum=s, component=NAME[c], metric=m, n=int(ok.sum()), spearman=r_, p=p_, partial_spearman_depth_age_sex=pr))
COR = pd.DataFrame(cor); COR.to_csv(RES + '/len_corr_all_strata.tsv', sep='\t', index=False)
cs_ = COR.groupby(['subset', 'component', 'metric']).agg(AllAll_rho=('spearman', lambda x: np.nan), median_rho_59=('spearman', 'median'), n_pos=('spearman', lambda x: (x > 0).sum()), n_neg=('spearman', lambda x: (x < 0).sum()), n_strata=('spearman', 'size')).reset_index()
aa = COR[COR.stratum == 'All|All'].set_index(['subset', 'component', 'metric'])
cs_['AllAll_rho'] = [aa.spearman.get((a, b, c), np.nan) for a, b, c in zip(cs_.subset, cs_.component, cs_.metric)]
cs_['AllAll_p'] = [aa.p.get((a, b, c), np.nan) for a, b, c in zip(cs_.subset, cs_.component, cs_.metric)]
cs_['AllAll_partial_rho'] = [aa.partial_spearman_depth_age_sex.get((a, b, c), np.nan) for a, b, c in zip(cs_.subset, cs_.component, cs_.metric)]
cs_.to_csv(RES + '/len_corr_summary.tsv', sep='\t', index=False)
# ---- C. 보정 유지율
def design(sub, kind):
    cols = {'len1': ['frac300'], 'len3': MET_ADJ, 'len3q': MET_ADJ, 'alpha': MET_ADJ + ['pbmc_alpha'], 'full': MET_ADJ + ['pbmc_alpha', 'age', 'sex', 'log_depth']}[kind]
    X = sub[cols].apply(zs)
    if kind == 'len3q':
        for c in MET_ADJ: X[c + '_sq'] = X[c] ** 2
    return X
def effect(y, h, X):
    """y: 표본 벡터, h: HCC 더미. 반환 (조정 전, 조정 후) HCC 계수."""
    ok = np.isfinite(y) & np.isfinite(X).all(1); y, h, X = y[ok], h[ok], X[ok]
    b0 = np.linalg.lstsq(np.column_stack([np.ones(len(y)), h]), y, rcond=None)[0][1]
    b1 = np.linalg.lstsq(np.column_stack([np.ones(len(y)), h, X]), y, rcond=None)[0][1]
    return b0, b1, ok.sum()
ret = []; main = D[D.grp.isin(['HCC', 'healthy'])].copy(); main['h'] = (main.grp == 'HCC').astype(float)
sets = {'all': main, 'Train': main[main.split == 'Train'], 'Test': main[main.split == 'Test']}
for sn, sub in sets.items():
    for kind in ['len1', 'len3', 'len3q', 'alpha', 'full']:
        X = design(sub, kind).to_numpy(); h = sub.h.to_numpy()
        for s in strata:
            for c in COMP:
                b0, b1, n = effect(Y[(s, c)].reindex(sub.index).to_numpy(float), h, X)
                ret.append(dict(set=sn, model=kind, stratum=s, component=NAME[c], n=n, effect_raw=b0, effect_adj=b1, retention=b1 / b0 if b0 != 0 else np.nan))
RET = pd.DataFrame(ret); RET.to_csv(RES + '/len_adjust_retention_all_strata.tsv', sep='\t', index=False)
# 요약 + All|All bootstrap
bs_rows = []
for sn, sub in sets.items():
    for kind in ['len1', 'len3', 'len3q', 'alpha', 'full']:
        X = design(sub, kind).to_numpy(); h = sub.h.to_numpy()
        gi = [np.where(h == 0)[0], np.where(h == 1)[0]]
        for c in COMP:
            y = Y[('All|All', c)].reindex(sub.index).to_numpy(float)
            bt = []
            for _ in range(NB):
                ix = np.concatenate([rng.choice(g_, len(g_)) for g_ in gi])
                try:
                    b0, b1, _n = effect(y[ix], h[ix], X[ix]); bt.append(b1 / b0)
                except Exception: pass
            b0, b1, n = effect(y, h, X)
            q = RET[(RET.set == sn) & (RET.model == kind) & (RET.component == NAME[c])]
            bs_rows.append(dict(set=sn, model=kind, component=NAME[c], n=n, AllAll_effect_raw=b0, AllAll_effect_adj=b1, AllAll_retention=b1 / b0, ret_lo=np.nanpercentile(bt, 2.5), ret_hi=np.nanpercentile(bt, 97.5),
                                median_retention_59=q.retention.median(), min_retention_59=q.retention.min(), max_retention_59=q.retention.max()))
pd.DataFrame(bs_rows).to_csv(RES + '/len_adjust_retention_summary.tsv', sep='\t', index=False)
# 분위 표준화 / 공통 지지 범위 (frac300), All|All
sup = []
for c in COMP:
    y = Y[('All|All', c)].reindex(main.index); x = main.frac300
    ok = y.notna() & x.notna(); dd = pd.DataFrame(dict(y=y[ok], x=x[ok], h=main.h[ok]))
    def stdz(dd):
        dd = dd.assign(q=pd.qcut(dd.x, 5, labels=False, duplicates='drop')); tot = 0; wsum = 0
        for q_, g_ in dd.groupby('q'):
            a, b = g_[g_.h == 1], g_[g_.h == 0]
            if len(a) >= 3 and len(b) >= 3: tot += len(a) * (a.y.mean() - b.y.mean()); wsum += len(a)
        return tot / wsum
    def trim(dd):
        hh, cc = dd[dd.h == 0].x, dd[dd.h == 1].x; lo = max(hh.quantile(.05), cc.quantile(.05)); hi = min(hh.quantile(.95), cc.quantile(.95)); g_ = dd[(dd.x >= lo) & (dd.x <= hi)]
        return g_[g_.h == 1].y.mean() - g_[g_.h == 0].y.mean(), len(g_[g_.h == 1]), len(g_[g_.h == 0])
    raw = dd[dd.h == 1].y.mean() - dd[dd.h == 0].y.mean(); st_ = stdz(dd); tr = trim(dd)
    bt1, bt2 = [], []
    gi = [dd[dd.h == 0].index.to_numpy(), dd[dd.h == 1].index.to_numpy()]
    for _ in range(NB):
        ix = np.concatenate([rng.choice(g_, len(g_)) for g_ in gi]); b_ = dd.loc[ix].reset_index(drop=True)
        r0 = b_[b_.h == 1].y.mean() - b_[b_.h == 0].y.mean(); bt1.append(stdz(b_) / r0); bt2.append(trim(b_)[0] / r0)
    sup.append(dict(component=NAME[c], effect_raw=raw, effect_quintile_standardized=st_, retention_quintile=st_ / raw, q_lo=np.percentile(bt1, 2.5), q_hi=np.percentile(bt1, 97.5),
                    effect_common_support=tr[0], n_hcc_support=tr[1], n_healthy_support=tr[2], retention_support=tr[0] / raw, s_lo=np.percentile(bt2, 2.5), s_hi=np.percentile(bt2, 97.5)))
pd.DataFrame(sup).to_csv(RES + '/len_frac300_quintile_support_AllAll.tsv', sep='\t', index=False)
print('retention done', flush=True)
# ---- D. 같은 기관 비교 (내부 참고용; 원고에는 기관 정보를 쓰지 않는다)
def boot_diff(a, b, nb=int(os.environ.get('LEN_NB2', 2000))):
    """a, b: (n x p) 배열. 평균/중앙값 차이와 bootstrap CI."""
    ia = rng.integers(0, len(a), (nb, len(a))); ib = rng.integers(0, len(b), (nb, len(b)))
    dm = np.nanmean(a, 0) - np.nanmean(b, 0); dmed = np.nanmedian(a, 0) - np.nanmedian(b, 0)
    bm = np.stack([np.nanmean(a[i], 0) - np.nanmean(b[j], 0) for i, j in zip(ia, ib)]); bmed = np.stack([np.nanmedian(a[i], 0) - np.nanmedian(b[j], 0) for i, j in zip(ia, ib)])
    return dm, np.percentile(bm, 2.5, 0), np.percentile(bm, 97.5, 0), dmed, np.percentile(bmed, 2.5, 0), np.percentile(bmed, 97.5, 0)
def contrasts(ga, gb, label):
    cols = [(s, c) for s in strata for c in COMP] + [('LEN', m) for m in ['frac300', 'frac400', 'median_len', 'frac166', 'pbmc_alpha']]
    def mat(g): return np.column_stack([Y[(s, c)].reindex(g.index).to_numpy(float) if s != 'LEN' else g[c].to_numpy(float) for s, c in cols])
    a, b = mat(ga), mat(gb); dm, l1, h1, dmed, l2, h2 = boot_diff(a, b)
    sdh = np.nanstd(b, 0, ddof=1)
    return pd.DataFrame(dict(comparison=label, stratum=[s for s, _ in cols], component=[NAME.get(c, c) for _, c in cols], n_A=len(ga), n_B=len(gb), mean_A=np.nanmean(a, 0), mean_B=np.nanmean(b, 0),
                             diff_mean=dm, ci_lo=l1, ci_hi=h1, diff_median=dmed, med_ci_lo=l2, med_ci_hi=h2, smd_vs_B_sd=dm / sdh))
hcc = D[D.grp == 'HCC']; hl = D[D.grp == 'healthy']
ERA = CFG.param('same_site_era')   # 수집 연도 코드 (config)
snu = hcc[hcc.index.isin(CFG.ids('same_site_hcc')) & (hcc.era == ERA)]; snuh = hl[hl.index.isin(CFG.ids('same_site_healthy')) & (hl.era == ERA)]
C = pd.concat([contrasts(hcc, hl, 'FULL: HCC(all) vs healthy(all)'), contrasts(snu, snuh, 'SAME-INST: HCC vs healthy (same collection site, same era)'),
               contrasts(hcc[hcc.era == ERA], hl[hl.era == ERA], 'ERA: HCC vs healthy (same era, all sites)')])
C.to_csv(RES + '/len_same_institution_contrasts_all_strata.tsv', sep='\t', index=False)
aa = C[C.stratum == 'All|All'].copy(); aa.to_csv(RES + '/len_same_institution_AllAll.tsv', sep='\t', index=False)
# 59층 방향 요약
sm = []
for lab, g in C.groupby('comparison'):
    for comp in [NAME[c] for c in COMP]:
        q = g[(g.component == comp) & (g.stratum.isin(strata))]
        sm.append(dict(comparison=lab, component=comp, n_strata=len(q), n_pos=(q.diff_mean > 0).sum(), n_pos_CI=(q.ci_lo > 0).sum(), n_neg=(q.diff_mean < 0).sum(), n_neg_CI=(q.ci_hi < 0).sum(), median_diff=q.diff_mean.median()))
pd.DataFrame(sm).to_csv(RES + '/len_same_institution_strata_signcounts.tsv', sep='\t', index=False)
# 같은 기관 효과 / 전체 효과 비율 (All|All)
f = aa[aa.comparison.str.startswith('FULL')].set_index('component'); s_ = aa[aa.comparison.str.startswith('SAME')].set_index('component')
pd.DataFrame(dict(full=f.diff_mean, same_inst=s_.diff_mean, ratio=s_.diff_mean / f.diff_mean)).to_csv(RES + '/len_same_institution_vs_full_ratio.tsv', sep='\t')
# 출처별 중앙값 (strict-mixed, sOUF, FPI)
src = []
for (g, sc), d in D[D.grp.isin(['HCC', 'healthy'])].groupby(['grp', 'source']):
    r = dict(group=g, source=sc, n=len(d), eras=','.join(map(str, sorted(d.era.dropna().unique().astype(int)))))
    for c in ['mixed', 'opp', 'BV', 'B', 'strictM', 'U']: r[NAME[c] + '_median'] = Y[('All|All', c)].reindex(d.index).median()
    for m in ['frac300', 'median_len', 'pbmc_alpha']: r[m + '_median'] = d[m].median()
    src.append(r)
pd.DataFrame(src).to_csv(RES + '/len_by_source_AllAll_medians.tsv', sep='\t', index=False)
# 같은 기관 세트에서 길이 지표 보정 후 효과 유지율
sub = pd.concat([snu, snuh]).copy(); sub['h'] = (sub.group == 'liver').astype(float); rr = []
for kind in ['len1', 'len3', 'alpha', 'full']:
    X = design(sub, kind).to_numpy()
    for c in COMP:
        b0, b1, n = effect(Y[('All|All', c)].reindex(sub.index).to_numpy(float), sub.h.to_numpy(), X); rr.append(dict(model=kind, component=NAME[c], n=n, effect_raw=b0, effect_adj=b1, retention=b1 / b0))
pd.DataFrame(rr).to_csv(RES + '/len_same_institution_adjusted_AllAll.tsv', sep='\t', index=False)
print('done', flush=True)
