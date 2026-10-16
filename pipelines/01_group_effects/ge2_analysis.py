"""군 비교 단계 2단계: agg_sums.npy(표본 x 59층 x 열 합계)로 group-effect tables 분석표를 만든다.

규칙: 결정(R_B 정의, bin 가중, margin 목록)은 Train에서만 정하고 Test에는 같은 설정을 한 번 적용한다.
강건 신뢰구간: 군별 표본 부트스트랩(2000회, 다항 가중, 고정 seed).
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import json, sys, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings('ignore')
D = CFG.DATA; A = CFG.work('group_effects'); RES = A + '/results'
NB = 2000; MINM = 20

sums = np.load(f'{RES}/agg_sums.npy'); meta = json.load(open(f'{RES}/agg_meta.json'))
names = meta['names']; strata = meta['strata']; col = {n: i for i, n in enumerate(names)}
smp = pd.read_csv(f'{RES}/agg_samples.tsv', sep='\t')
N, S, _ = sums.shape
smp['has_strict'] = sums[:, 0, col['n5s|cnt']] > 0   # All|All 층에서 strict 통계 유효 marker가 있는 표본
smp['hcc'] = (smp.group == 'liver').astype(int)
assert len(smp) == N


def cnt(tag):
    return sums[:, :, col[f'{tag}|cnt']]


def mean(tag, nm, minm=MINM):
    c = cnt(tag); x = sums[:, :, col[f'{tag}|{nm}']] / np.where(c > 0, c, np.nan)
    return np.where(c >= minm, x, np.nan)


def bh(p):
    p = np.asarray(p, float); q = np.full(len(p), np.nan); ok = np.isfinite(p); pp = p[ok]
    if not len(pp): return q
    o = np.argsort(pp); r = pp[o] * len(pp) / (np.arange(len(pp)) + 1); r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(len(pp)); out[o] = np.minimum(r, 1); q[ok] = out; return q


rng = np.random.default_rng(20261014)
WCACHE = {}


def boot_w(n, key):
    if (n, key) not in WCACHE:
        r = np.random.default_rng(1000 + n * 7 + key)
        WCACHE[(n, key)] = r.multinomial(n, np.ones(n) / n, size=NB).astype(float)
    return WCACHE[(n, key)]


def bmean(X, W):
    ok = np.isfinite(X); num = W @ np.where(ok, X, 0); den = W @ ok
    return num / np.where(den > 0, den, np.nan)


def contrast(Xc, Xh, key=0):
    """열(층)별 HCC(Xc) - healthy(Xh) 차이: 평균차, 중앙값차, bootstrap 95% CI, rank-sum p."""
    P = Xc.shape[1]; Wc = boot_w(Xc.shape[0], key); Wh = boot_w(Xh.shape[0], key + 1)
    d = np.nanmean(Xc, 0) - np.nanmean(Xh, 0); md = np.nanmedian(Xc, 0) - np.nanmedian(Xh, 0)
    bd = bmean(Xc, Wc) - bmean(Xh, Wh)
    lo, hi = np.nanpercentile(bd, [2.5, 97.5], axis=0)
    p = np.full(P, np.nan)
    for j in range(P):
        a = Xc[:, j]; b = Xh[:, j]; a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
        if len(a) >= 5 and len(b) >= 5: p[j] = stats.mannwhitneyu(a, b, alternative='two-sided').pvalue
    return d, md, lo, hi, p, np.isfinite(Xc).sum(0), np.isfinite(Xh).sum(0), np.nanmean(Xc, 0), np.nanmean(Xh, 0)


def subsets(strict):
    """(label, 표본 마스크 healthy, HCC). Test: 전원. Train: strict 필요 지표는 strict 보유 표본만."""
    out = []
    for sp in ['Train', 'Test']:
        m = (smp.split == sp).to_numpy()
        if strict: m = m & smp.has_strict.to_numpy().astype(bool)
        out.append((sp + ('(strict보유)' if strict and sp == 'Train' else ''), m))
    return out


def effect_table(metric_arrays, strict_flags, extra=None):
    """metric_arrays: {name: N x S 배열}. 반환 long table (split x metric x stratum)."""
    rows = []
    for nm, X in metric_arrays.items():
        for lab, m in subsets(strict_flags[nm]):
            h = X[m & (smp.hcc == 0).to_numpy()]; c = X[m & (smp.hcc == 1).to_numpy()]
            d, md, lo, hi, p, nc, nh, mc_, mh_ = contrast(c, h, key={'Train': 0, 'Test': 10, 'Train(strict보유)': 20}[lab])
            q = bh(p)
            for j, s in enumerate(strata):
                rows.append(dict(split=lab, metric=nm, stratum=s, n_healthy=nh[j], n_HCC=nc[j], mean_healthy=mh_[j], mean_HCC=mc_[j], delta=d[j], median_delta=md[j], ci_lo=lo[j], ci_hi=hi[j], p=p[j], q_BH=q[j]))
    return pd.DataFrame(rows)


def summ(t, up=True):
    """metric x split 별 증가(또는 감소) 층 수, 유의 층 수, 층 중앙 효과."""
    out = []
    for (sp, mt), g in t.groupby(['split', 'metric']):
        g = g[np.isfinite(g.delta)]
        pos = (g.delta > 0).sum(); sig_pos = ((g.delta > 0) & (g.q_BH < .05)).sum(); sig_neg = ((g.delta < 0) & (g.q_BH < .05)).sum()
        ci_pos = (g.ci_lo > 0).sum(); ci_neg = (g.ci_hi < 0).sum()
        out.append(dict(split=sp, metric=mt, n_strata=len(g), n_delta_pos=pos, n_delta_neg=len(g) - pos, n_q05_pos=sig_pos, n_q05_neg=sig_neg, n_CI_pos=ci_pos, n_CI_neg=ci_neg,
                        median_delta=g.delta.median(), n_median_delta_pos=(g.median_delta > 0).sum(), allall_delta=g.loc[g.stratum == 'All|All', 'delta'].squeeze() if (g.stratum == 'All|All').any() else np.nan))
    return pd.DataFrame(out)


# ---------------------------------------------------------------- 기본 지표 (n>=5, 층당 marker>=20)
B = mean('n5', 'B'); U = mean('n5', 'U'); V = mean('n5', 'V'); MU = mean('n5', 'mu'); K = mean('n5', 'K'); FL = mean('n5', 'fl'); SP = mean('span', 'span')
Rmk = mean('rmk', 'Rmk'); Rr = B / V
SM = mean('n5s', 'sM'); SU = mean('n5s', 'sU'); MIX = mean('n5s', 'mixed')
OPP = mean('full_m0', 'opp'); CONC = mean('full_m0', 'conc')
OPPT = mean('trH_m0', 'opp')
B10 = mean('n10', 'B'); U10 = mean('n10', 'U'); V10 = mean('n10', 'V'); R10 = B10 / V10
strict_of = dict(mu=False, U=False, B=False, V=False, R_ratio=False, R_markermean=False, mixed=True, strictM=True, strictU=True, opp=True, conc=True, K=False, fraglen=False, span=False, B_n10=False, R_ratio_n10=False)
ARR = dict(mu=MU, U=U, B=B, V=V, R_ratio=Rr, R_markermean=Rmk, mixed=MIX, strictM=SM, strictU=SU, opp=OPP, conc=CONC, K=K, fraglen=FL, span=SP, B_n10=B10, R_ratio_n10=R10)
T0 = effect_table(ARR, strict_of)
T0.to_csv(f'{RES}/ge_strata_effects_all_metrics.tsv', sep='\t', index=False)
S0 = summ(T0); S0.to_csv(f'{RES}/ge_effect_summary.tsv', sep='\t', index=False)
print(S0.to_string())

# ---------------------------------------------------------------- group-effect step : earlier analysis(기존)과 대조
prior = pd.read_csv(CFG.inp('prior_strata_effects'), sep='\t').set_index('stratum')
rows = []
cmap = dict(B='dB', U='dU', mixed='dMixedStrict', opp='dHoppStrict')
for sp in ['Train', 'Test', 'Train(strict보유)']:
    for mt, pc in cmap.items():
        g = T0[(T0.split == sp) & (T0.metric == mt)].set_index('stratum')
        if not len(g): continue
        j = g.join(prior[pc].rename('prior'), how='inner'); j = j[np.isfinite(j.delta) & np.isfinite(j.prior)]
        rows.append(dict(split=sp, metric=mt, n_strata=len(j), sign_agree=(np.sign(j.delta) == np.sign(j.prior)).sum(), spearman=stats.spearmanr(j.delta, j.prior)[0],
                         ratio_median=(j.delta / j.prior).median(), n_delta_expected_sign=((j.delta < 0) if mt in ('U', 'mixed') else (j.delta > 0)).sum(),
                         n_prior_expected_sign=((j.prior < 0) if mt in ('U', 'mixed') else (j.prior > 0)).sum()))
A0cmp = pd.DataFrame(rows); A0cmp.to_csv(f'{RES}/ge_prior_comparison.tsv', sep='\t', index=False); print(A0cmp.to_string())

# 짝지은 Test 73쌍 (기존 방식) 재현
man = pd.read_csv(CFG.inp('prior_engine_manifest'), sep='\t')
hm = man[man.role.eq('healthy_matched')]
pairs = [(lv, h) for h, lv in zip(hm.sid, hm.matched_to) if isinstance(lv, str)]
idx = {a: i for i, a in enumerate(smp.AnalysisID)}
pp = [(idx[l], idx[h]) for l, h in pairs if l in idx and h in idx and smp.group[idx[l]] == 'liver' and smp.split[idx[l]] == 'Test']
GE = pd.read_csv(CFG.inp('prior_group_effects'), sep='\t'); GE = GE[GE.contrast == 'TestLiver_vs_matchedHealthy']
rows = []
for mt, arr, pm in [('B', B, 'between_read_B_emp'), ('U', U, 'within_read_U'), ('mixed', MIX, 'mixed_strict'), ('opp', OPP, 'hopp_strict')]:
    li = [a for a, b in pp]; hi_ = [b for a, b in pp]; d = arr[li] - arr[hi_]
    g = GE[GE.metric == pm].set_index('stratum')
    for j, s in enumerate(strata):
        x = d[:, j]; x = x[np.isfinite(x)]
        rows.append(dict(metric=mt, stratum=s, n_pairs=len(x), delta_paired_core=x.mean() if len(x) else np.nan, prior_paired=g.mean_diff.get(s, np.nan)))
PC = pd.DataFrame(rows); PC['abs_diff'] = (PC.delta_paired_core - PC.prior_paired).abs(); PC.to_csv(f'{RES}/ge_Test_paired73_vs_prior.tsv', sep='\t', index=False)
print('paired pairs', len(pp)); print(PC.groupby('metric').agg(n=('abs_diff', 'count'), max_abs_diff=('abs_diff', 'max'), med_abs_diff=('abs_diff', 'median'), mean_pairs=('n_pairs', 'mean')))

# ---------------------------------------------------------------- group-effect step : R_B 정의 결정 (healthy Train만)
trh = ((smp.split == 'Train') & (smp.hcc == 0)).to_numpy()
def robust_cv(X):
    x = X[trh]; med = np.nanmedian(x, 0); mad = np.nanmedian(np.abs(x - med), 0) * 1.4826
    return mad / np.abs(med)
cv_r = robust_cv(Rr); cv_m = robust_cv(Rmk)
dec = pd.DataFrame(dict(stratum=strata, robustCV_ratio_of_means=cv_r, robustCV_marker_mean=cv_m, valid_ratio=np.isfinite(Rr[trh]).mean(0), valid_markermean=np.isfinite(Rmk[trh]).mean(0)))
n_r_better = (cv_r < cv_m).sum()
primary_R = 'R_ratio' if n_r_better >= 30 else 'R_markermean'
dec.attrs['primary'] = primary_R
dec.to_csv(f'{RES}/ge_R_definition_decision_trainHealthy.tsv', sep='\t', index=False)
json.dump(dict(n_strata_ratio_lower_CV=int(n_r_better), n_strata=int(len(strata)), median_CV_ratio=float(np.nanmedian(cv_r)), median_CV_markermean=float(np.nanmedian(cv_m)), primary=primary_R,
               rule='healthy Train 표본 간 robust CV(MAD*1.4826/|중앙값|)가 더 작은 정의를 59층 중 다수(>=30)에서 보이면 ratio-of-means를 1차로 채택; 아니면 marker평균'),
          open(f'{RES}/ge_R_definition_decision.json', 'w'), ensure_ascii=False, indent=1)
print('R primary:', primary_R, n_r_better, np.nanmedian(cv_r), np.nanmedian(cv_m))

# ---------------------------------------------------------------- group-effect step : ΔB 분해 (평균 이동 vs R_B)
rows = []
for sp in ['Train', 'Test']:
    m = (smp.split == sp).to_numpy(); h = m & (smp.hcc == 0).to_numpy(); c = m & (smp.hcc == 1).to_numpy()
    for j, s in enumerate(strata):
        b = B[:, j]; v = V[:, j]; r = Rr[:, j]
        ok = np.isfinite(b) & np.isfinite(v) & np.isfinite(r)
        hh = h & ok; cc = c & ok
        if hh.sum() < 5 or cc.sum() < 5: continue
        dB = b[cc].mean() - b[hh].mean(); dV = v[cc].mean() - v[hh].mean(); dR = r[cc].mean() - r[hh].mean()
        Rm = (r[cc].mean() + r[hh].mean()) / 2; Vm = (v[cc].mean() + v[hh].mean()) / 2
        tV = dV * Rm; tR = dR * Vm; res = dB - tV - tR
        rows.append(dict(split=sp, stratum=s, delta_B=dB, delta_V_mu1mu=dV, delta_R=dR, mean_R_healthy=r[hh].mean(), mean_R_HCC=r[cc].mean(), term_V_shift=tV, term_R_change=tR, term_interaction_resid=res,
                         frac_V_shift=tV / dB if dB else np.nan, frac_R_change=tR / dB if dB else np.nan))
DEC = pd.DataFrame(rows); DEC.to_csv(f'{RES}/ge_deltaB_decomposition.tsv', sep='\t', index=False)
dsum = []
for sp, g in DEC.groupby('split'):
    gp = g[g.delta_B > 0]
    dsum.append(dict(split=sp, n_strata=len(g), n_deltaB_pos=len(gp), n_deltaR_pos=(g.delta_R > 0).sum(), median_frac_V_shift=gp.frac_V_shift.median(), median_frac_R_change=gp.frac_R_change.median(),
                     n_R_term_gt_V_term=(gp.term_R_change.abs() > gp.term_V_shift.abs()).sum(), n_V_shift_pos=(g.delta_V_mu1mu > 0).sum(), n_V_shift_neg=(g.delta_V_mu1mu < 0).sum()))
pd.DataFrame(dsum).to_csv(f'{RES}/ge_decomposition_summary.tsv', sep='\t', index=False); print(pd.DataFrame(dsum).to_string())

# ---------------------------------------------------------------- group-effect step : 표본 수준 공변량 (K, 길이, span) 군 차이 + 회귀
cov_rows = []
for nm in ['K', 'fraglen', 'span']:
    pass
sub = T0[T0.metric.isin(['K', 'fraglen', 'span']) & T0.split.isin(['Train', 'Test'])]
sub.to_csv(f'{RES}/ge_sample_covariate_group_diffs.tsv', sep='\t', index=False)
print(sub[sub.stratum == 'All|All'][['split', 'metric', 'mean_healthy', 'mean_HCC', 'delta', 'ci_lo', 'ci_hi', 'p']].to_string())


def ols(y, X, g_idx=1):
    """HC3 로버스트 OLS. 반환 group 계수, CI, p."""
    ok = np.isfinite(y) & np.isfinite(X).all(1); y = y[ok]; X = X[ok]; n, k = X.shape
    if n < k + 10 or X[:, g_idx].std() == 0: return np.nan, np.nan, np.nan, np.nan, n
    XtXi = np.linalg.pinv(X.T @ X); beta = XtXi @ X.T @ y; e = y - X @ beta
    h = np.einsum('ij,jk,ik->i', X, XtXi, X); w = (e / (1 - h)) ** 2
    cov = XtXi @ (X.T * w) @ X @ XtXi; se = np.sqrt(cov[g_idx, g_idx]); df = n - k
    t = beta[g_idx] / se; p = 2 * stats.t.sf(abs(t), df); tc = stats.t.ppf(.975, df)
    return beta[g_idx], beta[g_idx] - tc * se, beta[g_idx] + tc * se, p, n


rows = []
outs = dict(B=B, U=U, R_ratio=Rr, mixed=MIX, opp=OPP)
for sp in ['Train', 'Test']:
    for nm, Y in outs.items():
        st = strict_of[nm]; m = (smp.split == sp).to_numpy() & ((smp.has_strict.to_numpy().astype(bool)) if st else True)
        lab = sp + ('(strict보유)' if st and sp == 'Train' else '')
        s_ = smp[m]; g = s_.hcc.to_numpy(float); dep = s_.log_depth.to_numpy(); age = s_.age.to_numpy(); sex = s_.sex.to_numpy(float)
        for j, s in enumerate(strata):
            y = Y[m, j]
            if np.isfinite(y).sum() < 30: continue
            Xa = np.column_stack([np.ones(len(g)), g, dep, age, sex])                         # 모델 A: 군 + depth + age + sex
            Xb = np.column_stack([Xa, K[m, j], FL[m, j]])                                      # 모델 B: + mean_K + mean_fraglen
            cr = ols(y, Xa)
            # 같은 표본(K,fraglen 유한)에서 비교하도록 A도 동일 ok 집합 위에서 재적합
            okk = np.isfinite(K[m, j]) & np.isfinite(FL[m, j]) & np.isfinite(y)
            ca = ols(y[okk], Xa[okk]); cb = ols(y[okk], Xb[okk])
            # 군 평균차(원 비교)
            d0 = np.nanmean(y[g == 1]) - np.nanmean(y[g == 0])
            rows.append(dict(split=lab, metric=nm, stratum=s, n=ca[4], raw_diff=d0, coef_A=ca[0], A_lo=ca[1], A_hi=ca[2], A_p=ca[3], coef_B=cb[0], B_lo=cb[1], B_hi=cb[2], B_p=cb[3],
                             retained_pct=100 * cb[0] / ca[0] if ca[0] not in (0, np.nan) and np.isfinite(ca[0]) else np.nan))
REG = pd.DataFrame(rows)
for (sp, mt), g in REG.groupby(['split', 'metric']):
    REG.loc[g.index, 'B_q'] = bh(g.B_p.to_numpy()); REG.loc[g.index, 'A_q'] = bh(g.A_p.to_numpy())
REG.to_csv(f'{RES}/ge_regression_adjusted_by_K_fraglen.tsv', sep='\t', index=False)
rs = []
for (sp, mt), g in REG.groupby(['split', 'metric']):
    g = g[np.isfinite(g.coef_A) & np.isfinite(g.coef_B)]
    same = np.sign(g.coef_A) == np.sign(g.coef_B)
    rs.append(dict(split=sp, metric=mt, n_strata=len(g), median_coef_A=g.coef_A.median(), median_coef_B=g.coef_B.median(), median_retained_pct=g.retained_pct.median(),
                   pooled_retained_pct=100 * g.coef_B.median() / g.coef_A.median(), n_A_q05_same_sign_as_median=((g.A_q < .05) & (np.sign(g.coef_A) == np.sign(g.coef_A.median()))).sum(),
                   n_B_q05_same_sign_as_median=((g.B_q < .05) & (np.sign(g.coef_B) == np.sign(g.coef_A.median()))).sum(), n_sign_kept=same.sum(),
                   n_B_pos=(g.coef_B > 0).sum(), n_B_neg=(g.coef_B < 0).sum(), n_B_CI_excl0_same=(((g.B_lo > 0) | (g.B_hi < 0)) & (np.sign(g.coef_B) == np.sign(g.coef_A.median()))).sum()))
RS = pd.DataFrame(rs); RS.to_csv(f'{RES}/ge_regression_summary.tsv', sep='\t', index=False); print(RS.to_string())

# ---- group-effect step(a) marker 수준 bin 표준화 (K x fraglen 4x4; 경계와 가중은 Train에서만)
kb = json.load(open(f'{RES}/bin_edges.json'))
def binarr(prefix, nm):
    c = np.stack([cnt(f'{prefix}{k}') for k in range(16)], -1); x = np.stack([sums[:, :, col[f'{prefix}{k}|{nm}']] for k in range(16)], -1)
    return c, x
def standardize(prefix, nm, strict):
    c, x = binarr(prefix, nm); trm = (smp.split == 'Train').to_numpy() & ((smp.has_strict.to_numpy().astype(bool)) if strict else True)
    w = c[trm].sum(0); w = w / w.sum(-1, keepdims=True)                              # Train(두 군 합산) 층 x bin 가중
    out = np.full((N, S), np.nan); cover = np.zeros((N, S))
    ok = c >= 3; xm = np.where(ok, x / np.where(c > 0, c, 1), 0.0)
    num = (xm * w[None] * ok).sum(-1); wsum = (w[None] * ok).sum(-1)
    out = np.where(wsum >= .8, num / np.where(wsum > 0, wsum, np.nan), np.nan)
    tot = c.sum(-1); out = np.where(tot >= MINM, out, np.nan)
    return out
Bs = standardize('bin', 'B', False); Vs = standardize('bin', 'V', False); Us = standardize('bin', 'U', False)
Rs = Bs / Vs; Ms = standardize('bins', 'mixed', True); Os = standardize('bins', 'opp', True)
ARRs = {'B_binstd': Bs, 'V_binstd': Vs, 'U_binstd': Us, 'R_ratio_binstd': Rs, 'mixed_binstd': Ms, 'opp_binstd': Os}
TS = effect_table(ARRs, {k: k.startswith(('mixed', 'opp')) for k in ARRs}); TS.to_csv(f'{RES}/ge_binstandardized_strata_effects.tsv', sep='\t', index=False)
SS = summ(TS); SS.to_csv(f'{RES}/ge_binstandardized_summary.tsv', sep='\t', index=False); print(SS.to_string())
# 유지율: bin 표준화 delta / 비표준화 delta (같은 split, 같은 층)
ret = []
for mt, raw in [('B_binstd', 'B'), ('V_binstd', 'V'), ('U_binstd', 'U'), ('R_ratio_binstd', 'R_ratio'), ('mixed_binstd', 'mixed'), ('opp_binstd', 'opp')]:
    for sp in T0.split.unique():
        a = TS[(TS.split == sp) & (TS.metric == mt)].set_index('stratum'); b = T0[(T0.split == sp) & (T0.metric == raw)].set_index('stratum')
        if not len(a) or not len(b): continue
        j = a[['delta']].join(b[['delta']], lsuffix='_std', rsuffix='_raw', how='inner').dropna()
        ret.append(dict(split=sp, metric=mt, n_strata=len(j), median_retained_pct=(100 * j.delta_std / j.delta_raw).median(), pooled_retained_pct=100 * j.delta_std.median() / j.delta_raw.median(),
                        n_sign_kept=(np.sign(j.delta_std) == np.sign(j.delta_raw)).sum()))
pd.DataFrame(ret).to_csv(f'{RES}/ge_binstandardized_retention.tsv', sep='\t', index=False); print(pd.DataFrame(ret).to_string())

# ---------------------------------------------------------------- group-effect step : 방향 margin
rows = []
for bn in ['full', 'trH']:
    for mg in [0, 1, 2, 3]:
        tag = f'{bn}_m{mg}'
        arrs = {f'opp@{tag}': mean(tag, 'opp'), f'conc@{tag}': mean(tag, 'conc'), f'B@{tag}': mean(tag, 'B'), f'mixed@{tag}': mean(tag, 'mixed'),
                f'Bsub@{tag}d': mean(tag + 'd', 'Bsub'), f'Bdel@{tag}d': mean(tag + 'd', 'Bdel')}
        TT = effect_table(arrs, {k: True for k in arrs}); TT['baseline'] = bn; TT['margin'] = mg / 10; TT['metric'] = TT.metric.str.split('@').str[0]
        rows.append(TT)
M3 = pd.concat(rows, ignore_index=True); M3.to_csv(f'{RES}/ge_margin_strata_effects.tsv', sep='\t', index=False)
# 층별 감소 비율
red = []
for (bn, mg, sp), g in M3.groupby(['baseline', 'margin', 'split']):
    e = g.pivot_table(index='stratum', columns='metric', values='delta'); q = g.pivot_table(index='stratum', columns='metric', values='q_BH')
    e['red'] = 1 - e.Bdel / e.Bsub; e.loc[~(e.Bsub > 0), 'red'] = np.nan
    ev = e.dropna(subset=['opp']);
    red.append(dict(baseline=bn, margin=mg, split=sp, n_strata_eval=int(ev.opp.notna().sum()), n_dOpp_pos=int((ev.opp > 0).sum()), n_dOpp_q05_pos=int(((e.opp > 0) & (q.opp < .05)).sum()),
                    median_dOpp=ev.opp.median(), n_dB_pos=int((e.B > 0).sum()), n_dB_q05_pos=int(((e.B > 0) & (q.B < .05)).sum()), median_dB=e.B.median(),
                    n_red_evaluable=int(e.red.notna().sum()), median_reduction_pct=100 * e.red.median(), n_reduction_ge50=int((e.red >= .5).sum()),
                    n_dMixed_neg=int((e.mixed < 0).sum()), n_dConc_pos=int((e.conc > 0).sum()), n_dConc_neg=int((e.conc < 0).sum())))
RD = pd.DataFrame(red); RD.to_csv(f'{RES}/ge_margin_summary.tsv', sep='\t', index=False); print(RD.to_string())
# 층별 reduction
rr = []
for (bn, mg, sp), g in M3.groupby(['baseline', 'margin', 'split']):
    e = g.pivot_table(index='stratum', columns='metric', values='delta'); e['reduction'] = 1 - e.Bdel / e.Bsub; e.loc[~(e.Bsub > 0), 'reduction'] = np.nan
    e['baseline'] = bn; e['margin'] = mg; e['split'] = sp; rr.append(e.reset_index())
pd.concat(rr).to_csv(f'{RES}/ge_margin_reduction_by_stratum.tsv', sep='\t', index=False)

# baseline 진단
mb = pd.read_parquet(f'{RES}/marker_baselines.parquet'); ok = mb.baseline_trainH.notna()
dg = dict(n_markers=int(len(mb)), n_trainH_defined=int(ok.sum()),
          corr_full_vs_trainH=float(np.corrcoef(mb.baseline_fullTrain[ok], mb.baseline_trainH[ok])[0, 1]),
          corr_full_vs_TestHealthy_diag=float(np.corrcoef(mb.baseline_fullTrain[mb.mean_mu_TestHealthy_diag.notna()], mb.mean_mu_TestHealthy_diag.dropna())[0, 1]),
          corr_trainH_vs_TestHealthy_diag=float(mb[ok & mb.mean_mu_TestHealthy_diag.notna()][['baseline_trainH', 'mean_mu_TestHealthy_diag']].corr().iloc[0, 1]),
          direction_agree_full_trainH=float(((mb.baseline_fullTrain[ok] < .5) == (mb.baseline_trainH[ok] < .5)).mean()),
          direction_agree_full_TestHealthy=float(((mb.baseline_fullTrain < .5)[mb.mean_mu_TestHealthy_diag.notna()] == (mb.mean_mu_TestHealthy_diag.dropna() < .5)).mean()),
          direction_agree_trainH_TestHealthy=float(((mb.baseline_trainH < .5)[ok & mb.mean_mu_TestHealthy_diag.notna()] == (mb.mean_mu_TestHealthy_diag[ok & mb.mean_mu_TestHealthy_diag.notna()] < .5)).mean()),
          mean_abs_diff_full_trainH=float((mb.baseline_fullTrain[ok] - mb.baseline_trainH[ok]).abs().mean()),
          mean_abs_diff_full_TestHealthy=float((mb.baseline_fullTrain - mb.mean_mu_TestHealthy_diag).abs().mean()),
          mean_abs_diff_trainH_TestHealthy=float((mb.baseline_trainH - mb.mean_mu_TestHealthy_diag).abs().mean()))
for m in [0, .1, .2, .3]:
    dg[f'n_markers_margin{m}_full'] = int((np.abs(mb.baseline_fullTrain - .5) >= m).sum()); dg[f'n_markers_margin{m}_trainH'] = int((np.abs(mb.baseline_trainH - .5) >= m).sum())
json.dump(dg, open(f'{RES}/ge_baseline_diagnostics.json', 'w'), indent=1); print(dg)

# 표본 단위 파일 (하류 그림용)
cols = {}
for nm, X in dict(B=B, U=U, V=V, R_ratio=Rr, opp=OPP, mixed=MIX, K=K, fraglen=FL).items():
    cols[nm] = X[:, strata.index('All|All')]
pd.concat([smp[['AnalysisID', 'group', 'split', 'age', 'sex', 'log_depth', 'has_strict']], pd.DataFrame(cols)], axis=1).to_csv(f'{RES}/ge_sample_AllAll_values.tsv', sep='\t', index=False)
print('analysis done')
