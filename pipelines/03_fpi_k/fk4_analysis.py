"""FPI_K step: K 보정 FPI (FPI_K) 59층 분석. 군 비교 단계 집계(agg_sums.npy: ΣB, ΣV, cnt, ΣK, Σfl)와 FPI_K step 집계(a6_sums.npy: Σh·V 등)를 합쳐 쓴다.
정의(층 s, 표본마다): FPI   = ΣB/ΣV                         (군 비교 단계의 R_ratio, ratio-of-means)
                     h_w   = Σ(h·V)/ΣV,  h = 1/mean_K       (V 가중 binomial 바닥)
                     FPI_K = Σ(B − c·h·V) / Σ((1 − c·h)·V)  (c=1이 주 방법; c>1은 Jensen 보정 민감도)
 marker 평균형: FPI_K_mk = mean_marker[(B − hV)/((1−h)V)]  (V>=0.01)  -- 정의 선택은 Train healthy robust CV로 (군 비교 단계와 같은 규칙).
검정/CI/BH/회귀는 군 비교 단계 ge2_analysis.py와 동일(군별 표본 bootstrap 2000, Mann-Whitney, HC3 OLS)이다. 층당 marker>=20.
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import json, os, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings('ignore')
A = CFG.WORK; RA = A + '/group_effects/results'; RES = A + '/fpi_k/results'
NB = 2000; MINM = 20
sums = np.load(f'{RA}/agg_sums.npy', mmap_mode='r'); meta = json.load(open(f'{RA}/agg_meta.json')); names = meta['names']; strata = meta['strata']; col = {n: i for i, n in enumerate(names)}
s6 = np.load(f'{RES}/fpik_sums.npy'); m6 = json.load(open(f'{RES}/fpik_meta.json')); col6 = {n: i for i, n in enumerate(m6['names'])}
assert m6['samples'] == meta['samples']
smp = pd.read_csv(f'{RA}/agg_samples.tsv', sep='\t'); smp['hcc'] = (smp.group == 'liver').astype(int)
N, S = len(smp), len(strata)
G = lambda nm: np.asarray(sums[:, :, col[nm]], float)
G6 = lambda nm: s6[:, :, col6[nm]]
cnt = G('n5|cnt'); ok = cnt >= MINM
SB, SV, SK, SF = G('n5|B'), G('n5|V'), G('n5|K'), G('n5|fl'); ShV = G6('n5|hV')
mean = lambda x: np.where(ok, x / np.where(cnt > 0, cnt, np.nan), np.nan)
K = mean(SK); FL = mean(SF)
cg = None
jf = f'{RES}/fpik_raw_jensen_summary.json'
if os.path.exists(jf): cg = json.load(open(jf))['c_marker_Vweighted_AllAll_mean']; print('c_global (Train healthy raw, marker-level V-weighted):', cg)
c_s = pd.read_csv(f'{RES}/fpik_khist_correction_by_sample.tsv', sep='\t').set_index('AnalysisID').c


def fpik(c):
    with np.errstate(all='ignore'):
        return np.where(ok, (SB - c * ShV) / (SV - c * ShV), np.nan)


FPI = np.where(ok, SB / SV, np.nan); FPIK = fpik(1.0); HW = np.where(ok, ShV / SV, np.nan)
rm = G6('rmkK|cnt'); FPIK_mk = np.where(rm >= MINM, G6('rmkK|FPIK_mk') / np.where(rm > 0, rm, np.nan), np.nan)
FPI_mk = np.where(G('rmk|cnt') >= MINM, G('rmk|Rmk') / np.where(G('rmk|cnt') > 0, G('rmk|cnt'), np.nan), np.nan)
ARR = dict(FPI=FPI, FPI_K=FPIK, h_floor=HW, FPI_markermean=FPI_mk, FPI_K_markermean=FPIK_mk)
if cg: ARR['FPI_K_cglobal'] = fpik(cg)
cs = smp.AnalysisID.map(c_s).to_numpy(float)                     # hist 보유 표본(Test 173 + Train HCC 7)만 유한
with np.errstate(all='ignore'):
    ARR['FPI_K_csample'] = np.where(np.isfinite(cs)[:, None] & ok, (SB - cs[:, None] * ShV) / (SV - cs[:, None] * ShV), np.nan)


def bh(p):
    p = np.asarray(p, float); q = np.full(len(p), np.nan); k_ = np.isfinite(p); pp = p[k_]
    if not len(pp): return q
    o = np.argsort(pp); r = pp[o] * len(pp) / (np.arange(len(pp)) + 1); r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(len(pp)); out[o] = np.minimum(r, 1); q[k_] = out; return q


WC = {}
def boot_w(n, key):
    if (n, key) not in WC: WC[(n, key)] = np.random.default_rng(1000 + n * 7 + key).multinomial(n, np.ones(n) / n, size=NB).astype(float)
    return WC[(n, key)]
def bmean(X, W):
    o = np.isfinite(X); return (W @ np.where(o, X, 0)) / np.where(W @ o > 0, W @ o, np.nan)
def contrast(Xc, Xh, key):
    P = Xc.shape[1]; bd = bmean(Xc, boot_w(Xc.shape[0], key)) - bmean(Xh, boot_w(Xh.shape[0], key + 1))
    d = np.nanmean(Xc, 0) - np.nanmean(Xh, 0); md = np.nanmedian(Xc, 0) - np.nanmedian(Xh, 0); lo, hi = np.nanpercentile(bd, [2.5, 97.5], axis=0); p = np.full(P, np.nan)
    for j in range(P):
        a = Xc[:, j]; b = Xh[:, j]; a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
        if len(a) >= 5 and len(b) >= 5: p[j] = stats.mannwhitneyu(a, b, alternative='two-sided').pvalue
    return d, md, lo, hi, p, np.isfinite(Xc).sum(0), np.isfinite(Xh).sum(0), np.nanmean(Xc, 0), np.nanmean(Xh, 0)


def effect_table(arrs, extra_mask=None):
    rows = []
    for nm, X in arrs.items():
        for sp, key in (('Train', 0), ('Test', 10)):
            m = (smp.split == sp).to_numpy()
            if nm == 'FPI_K_csample': m = m & np.isfinite(cs)
            h = X[m & (smp.hcc == 0).to_numpy()]; c = X[m & (smp.hcc == 1).to_numpy()]
            if len(h) < 5 or len(c) < 5: continue
            d, md, lo, hi, p, nc, nh, mc_, mh_ = contrast(c, h, key); q = bh(p)
            for j, s in enumerate(strata):
                rows.append(dict(split=sp, metric=nm, stratum=s, n_healthy=nh[j], n_HCC=nc[j], mean_healthy=mh_[j], mean_HCC=mc_[j], delta=d[j], median_delta=md[j], ci_lo=lo[j], ci_hi=hi[j], p=p[j], q_BH=q[j]))
    return pd.DataFrame(rows)


T = effect_table(ARR); T.to_csv(f'{RES}/fpik_strata_effects.tsv', sep='\t', index=False)
out = []
for (sp, mt), g in T.groupby(['split', 'metric']):
    g = g[np.isfinite(g.delta)]
    out.append(dict(split=sp, metric=mt, n_strata=len(g), n_delta_pos=int((g.delta > 0).sum()), n_delta_neg=int((g.delta < 0).sum()), n_q05_pos=int(((g.delta > 0) & (g.q_BH < .05)).sum()), n_q05_neg=int(((g.delta < 0) & (g.q_BH < .05)).sum()),
                    n_CI_pos=int((g.ci_lo > 0).sum()), n_CI_neg=int((g.ci_hi < 0).sum()), median_delta=g.delta.median(), n_median_delta_pos=int((g.median_delta > 0).sum()),
                    allall_healthy=g.loc[g.stratum == 'All|All', 'mean_healthy'].squeeze(), allall_HCC=g.loc[g.stratum == 'All|All', 'mean_HCC'].squeeze(), allall_delta=g.loc[g.stratum == 'All|All', 'delta'].squeeze(),
                    allall_ci_lo=g.loc[g.stratum == 'All|All', 'ci_lo'].squeeze(), allall_ci_hi=g.loc[g.stratum == 'All|All', 'ci_hi'].squeeze(), allall_q=g.loc[g.stratum == 'All|All', 'q_BH'].squeeze()))
SUM = pd.DataFrame(out); SUM.to_csv(f'{RES}/fpik_effect_summary.tsv', sep='\t', index=False); print(SUM.to_string())

# ---- 정의 결정 (군 비교 단계와 같은 규칙: Train healthy robust CV 가 더 작은 정의를 59층 중 >=30에서 보이면 ratio-of-means)
trh = ((smp.split == 'Train') & (smp.hcc == 0)).to_numpy()
def rcv(X):
    x = X[trh]; med = np.nanmedian(x, 0); mad = np.nanmedian(np.abs(x - med), 0) * 1.4826; return mad / np.abs(med)
cr, cm = rcv(FPIK), rcv(FPIK_mk); nb = int((cr < cm).sum())
dec = dict(n_strata_ratio_lower_CV=nb, n_strata=S, median_CV_ratio_of_means=float(np.nanmedian(cr)), median_CV_marker_mean=float(np.nanmedian(cm)), primary='FPI_K (ratio-of-means)' if nb >= 30 else 'FPI_K_markermean',
           rule='Train healthy 표본 간 robust CV(MAD*1.4826/|중앙값|)가 더 작은 정의를 59층 중 >=30에서 보이면 ratio-of-means(Σ(B−hV)/Σ((1−h)V))를 1차로 채택')
json.dump(dec, open(f'{RES}/fpik_definition_decision.json', 'w'), ensure_ascii=False, indent=1); print(dec)
pd.DataFrame(dict(stratum=strata, robustCV_ratio_of_means=cr, robustCV_marker_mean=cm, robustCV_FPI_ratio=rcv(FPI))).to_csv(f'{RES}/fpik_definition_decision_trainHealthy.tsv', sep='\t', index=False)

# ---- FPI vs FPI_K 비교
rows = []; pair = []
for sp in ['Train', 'Test']:
    a = T[(T.split == sp) & (T.metric == 'FPI')].set_index('stratum'); b = T[(T.split == sp) & (T.metric == 'FPI_K')].set_index('stratum'); hh = T[(T.split == sp) & (T.metric == 'h_floor')].set_index('stratum')
    j = a[['delta', 'q_BH', 'mean_healthy', 'mean_HCC']].join(b[['delta', 'q_BH', 'ci_lo', 'ci_hi', 'mean_healthy', 'mean_HCC']], lsuffix='_FPI', rsuffix='_FPIK').join(hh[['delta', 'mean_healthy']].rename(columns=dict(delta='delta_h', mean_healthy='h_healthy'))).dropna(subset=['delta_FPI', 'delta_FPIK'])
    j['split'] = sp; j['ratio'] = j.delta_FPIK / j.delta_FPI; pair.append(j.reset_index())
    rows.append(dict(split=sp, n_strata=len(j), FPI_pos=int((j.delta_FPI > 0).sum()), FPIK_pos=int((j.delta_FPIK > 0).sum()), sign_agree=int((np.sign(j.delta_FPI) == np.sign(j.delta_FPIK)).sum()),
                     FPI_q05_pos=int(((j.delta_FPI > 0) & (j.q_BH_FPI < .05)).sum()), FPIK_q05_pos=int(((j.delta_FPIK > 0) & (j.q_BH_FPIK < .05)).sum()), FPIK_CI_pos=int((j.ci_lo > 0).sum()),
                     ratio_median=j.ratio.median(), ratio_IQR_lo=j.ratio.quantile(.25), ratio_IQR_hi=j.ratio.quantile(.75), pooled_ratio_medians=j.delta_FPIK.median() / j.delta_FPI.median(),
                     spearman=stats.spearmanr(j.delta_FPI, j.delta_FPIK)[0], median_delta_FPI=j.delta_FPI.median(), median_delta_FPIK=j.delta_FPIK.median(),
                     median_delta_h=j.delta_h.median(), median_h_healthy=j.h_healthy.median(), median_theory_ratio_1_over_1mh=(1 / (1 - j.h_healthy)).median()))
CMP = pd.DataFrame(rows); CMP.to_csv(f'{RES}/fpik_FPI_vs_FPIK_comparison.tsv', sep='\t', index=False); print(CMP.T.to_string())
pd.concat(pair).to_csv(f'{RES}/fpik_FPI_vs_FPIK_by_stratum.tsv', sep='\t', index=False)

# ---- 회귀 (A2와 같은 설계): 모델 A group+depth+age+sex, 모델 B + mean_K + mean_fraglen, 같은 표본에서
def ols(y, X, g_idx=1):
    o = np.isfinite(y) & np.isfinite(X).all(1); y = y[o]; X = X[o]; n, k = X.shape
    if n < k + 10 or X[:, g_idx].std() == 0: return np.nan, np.nan, np.nan, np.nan, n
    XtXi = np.linalg.pinv(X.T @ X); beta = XtXi @ X.T @ y; e = y - X @ beta; h = np.einsum('ij,jk,ik->i', X, XtXi, X); w = (e / (1 - h)) ** 2
    cov = XtXi @ (X.T * w) @ X @ XtXi; se = np.sqrt(cov[g_idx, g_idx]); df = n - k; t = beta[g_idx] / se; p = 2 * stats.t.sf(abs(t), df); tc = stats.t.ppf(.975, df)
    return beta[g_idx], beta[g_idx] - tc * se, beta[g_idx] + tc * se, p, n
rows = []
REGM = {k: v for k, v in ARR.items() if k in ('FPI', 'FPI_K', 'h_floor', 'FPI_K_cglobal')}
for sp in ['Train', 'Test']:
    m = (smp.split == sp).to_numpy(); s_ = smp[m]; g = s_.hcc.to_numpy(float); dep = s_.log_depth.to_numpy(); age = s_.age.to_numpy(); sex = s_.sex.to_numpy(float)
    for nm, Y in REGM.items():
        for j, s in enumerate(strata):
            y = Y[m, j]
            if np.isfinite(y).sum() < 30: continue
            Xa = np.column_stack([np.ones(len(g)), g, dep, age, sex]); Xb = np.column_stack([Xa, K[m, j], FL[m, j]])
            okk = np.isfinite(K[m, j]) & np.isfinite(FL[m, j]) & np.isfinite(y); ca = ols(y[okk], Xa[okk]); cb = ols(y[okk], Xb[okk])
            rows.append(dict(split=sp, metric=nm, stratum=s, n=ca[4], raw_diff=np.nanmean(y[g == 1]) - np.nanmean(y[g == 0]), coef_A=ca[0], A_lo=ca[1], A_hi=ca[2], A_p=ca[3], coef_B=cb[0], B_lo=cb[1], B_hi=cb[2], B_p=cb[3],
                             retained_pct=100 * cb[0] / ca[0] if np.isfinite(ca[0]) and ca[0] != 0 else np.nan))
REG = pd.DataFrame(rows)
for (sp, mt), g in REG.groupby(['split', 'metric']):
    REG.loc[g.index, 'B_q'] = bh(g.B_p.to_numpy()); REG.loc[g.index, 'A_q'] = bh(g.A_p.to_numpy())
REG.to_csv(f'{RES}/fpik_regression_by_stratum.tsv', sep='\t', index=False)
rs = []
for (sp, mt), g in REG.groupby(['split', 'metric']):
    g = g[np.isfinite(g.coef_A) & np.isfinite(g.coef_B)]; sgn = np.sign(g.coef_A.median())
    aa = g[g.stratum == 'All|All']
    rs.append(dict(split=sp, metric=mt, n_strata=len(g), median_coef_A=g.coef_A.median(), median_coef_B=g.coef_B.median(), median_retained_pct=g.retained_pct.median(), pooled_retained_pct=100 * g.coef_B.median() / g.coef_A.median(),
                   n_A_q05_same_sign=int(((g.A_q < .05) & (np.sign(g.coef_A) == sgn)).sum()), n_B_q05_same_sign=int(((g.B_q < .05) & (np.sign(g.coef_B) == sgn)).sum()), n_sign_kept=int((np.sign(g.coef_A) == np.sign(g.coef_B)).sum()),
                   n_B_pos=int((g.coef_B > 0).sum()), n_B_neg=int((g.coef_B < 0).sum()), n_B_CI_excl0_same=int((((g.B_lo > 0) | (g.B_hi < 0)) & (np.sign(g.coef_B) == sgn)).sum()),
                   allall_coef_A=aa.coef_A.squeeze() if len(aa) else np.nan, allall_coef_B=aa.coef_B.squeeze() if len(aa) else np.nan, allall_B_lo=aa.B_lo.squeeze() if len(aa) else np.nan, allall_B_hi=aa.B_hi.squeeze() if len(aa) else np.nan, allall_B_p=aa.B_p.squeeze() if len(aa) else np.nan))
RS = pd.DataFrame(rs); RS.to_csv(f'{RES}/fpik_regression_summary.tsv', sep='\t', index=False); print(RS.to_string())

# ---- healthy Train FPI vs h 관계 (표본 수준 All|All, 층 수준 중앙값)
j0 = strata.index('All|All'); hm = trh & np.isfinite(FPI[:, j0])
rel = []
x = HW[hm, j0]; y = FPI[hm, j0]; sl = stats.linregress(x, y)
rel.append(dict(level='sample (healthy Train, All|All)', n=int(hm.sum()), spearman=stats.spearmanr(x, y)[0], slope=sl.slope, intercept=sl.intercept, mean_h=x.mean(), mean_FPI=y.mean(), mean_FPI_K=FPIK[hm, j0].mean(), sd_FPI=y.std(), sd_FPI_K=FPIK[hm, j0].std(), cv_FPI=y.std() / y.mean(), cv_FPI_K=FPIK[hm, j0].std() / FPIK[hm, j0].mean()))
xs = np.nanmean(HW[trh], 0); ys = np.nanmean(FPI[trh], 0); okk = np.isfinite(xs) & np.isfinite(ys); sl = stats.linregress(xs[okk], ys[okk])
rel.append(dict(level='stratum (healthy Train mean, 59층)', n=int(okk.sum()), spearman=stats.spearmanr(xs[okk], ys[okk])[0], slope=sl.slope, intercept=sl.intercept, mean_h=xs[okk].mean(), mean_FPI=ys[okk].mean(), mean_FPI_K=np.nanmean(FPIK[trh], 0)[okk].mean()))
pd.DataFrame(rel).to_csv(f'{RES}/fpik_FPI_vs_h_healthyTrain.tsv', sep='\t', index=False); print(pd.DataFrame(rel).to_string())
pd.concat([smp[['AnalysisID', 'group', 'split', 'age', 'sex', 'log_depth']], pd.DataFrame({k: v[:, j0] for k, v in ARR.items()}), pd.DataFrame(dict(K=K[:, j0], fraglen=FL[:, j0]))], axis=1).to_csv(f'{RES}/fpik_sample_AllAll_values.tsv', sep='\t', index=False)
print('FPI_K analysis done')
