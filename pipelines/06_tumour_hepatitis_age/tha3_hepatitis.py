"""Hepatitis: healthy / HCC / hepatitis(sample_lists.hepatitis_included) 세 군의 직접 상태 비교 (같은 marker 층, 59층). 세 대비 HCC-H, Hep-H, HCC-Hep.
표본 집합: full(건강인 1009 / HCC 177 / 간염 목록, marker 행렬 + 간염 suff 같은 공식), Test(336/73/간염 목록), priorset(이전 분석 표본 249개 중 건강인·HCC, 본 계산값),
priorown(이전 분석 표(inputs.prior_strata_samples)의 원값 hopp_strict/strictM/strictU/mixed_strict/between_read_B_emp/within_read_U).
bootstrap: 군별 표본 복원추출 2000 (seed 고정). 보조: log_depth 보정 OLS.
"""
import numpy as np, pandas as pd
from scipy import stats
from tha_common import *
W, nm, S = load(); NB = 2000
rng = np.random.default_rng(SEED)
STR59 = pd.read_csv(CFG.inp('strata_definition'), sep='\t').query('used').stratum.tolist()
P = pd.read_parquet(CFG.inp('prior_strata_samples')); P = P[P['mode'] == 'raw']
mp = {'opp': 'hopp_strict', 'strictM': 'strictM', 'strictU': 'strictU', 'mixed': 'mixed_strict', 'B': 'between_read_B_emp', 'U': 'within_read_U', 'mu': 'mu'}
Pw = P[P.metric.isin(mp.values())].pivot_table(index='sid', columns=['stratum', 'metric'], values='value')
inv = {v: k for k, v in mp.items()}; Pw.columns = pd.MultiIndex.from_tuples([(a, inv[b]) for a, b in Pw.columns])
for s in Pw.columns.get_level_values(0).unique(): Pw[(s, 'BV')] = Pw[(s, 'B')] / (Pw[(s, 'B')] + Pw[(s, 'U')])
Pw = Pw.sort_index(axis=1)
p103 = set(P.sid.unique())
g = S.group
HEP = set(CFG.ids('hepatitis_included'))   # 간염은 이 목록만 사용
def ids(grp, split=None, within=None):
    m = (g == grp) if split is None else ((g == grp) & (S.split == split))
    if grp == 'hepatitis': m = m & S.index.isin(HEP)
    i = S.index[m]
    return [x for x in i if within is None or x in within]
SETS = {
 'full': (W, ids('healthy'), ids('liver'), ids('hepatitis')),
 'Test': (W, ids('healthy', 'Test'), ids('liver', 'Test'), ids('hepatitis')),
 'priorset': (W, ids('healthy', within=p103), ids('liver', within=p103), ids('hepatitis')),
 'priorown': (Pw, ids('healthy', within=p103), ids('liver', within=p103), ids('hepatitis')),
}
STRP = [s for s in STR59 if s in set(W.columns.get_level_values(0))]
def cols(X): return [(s, m) for s in (STR59 + ['U|nonIsland', 'M|nonIsland', 'mid|nonIsland', 'All|nonIsland']) for m in METRICS if (s, m) in X.columns]
def bmeans(X):
    v = X.to_numpy(float); val = np.isfinite(v); f = np.where(val, v, 0.0); n = len(v)
    cnt = rng.multinomial(n, np.full(n, 1 / n), size=NB).astype(float)
    return (f.mean(0) * 0 + (f.sum(0) / np.maximum(val.sum(0), 1))), (cnt @ f) / np.maximum(cnt @ val.astype(float), 1), val.sum(0)
rows = []; dirrows = []; means = []
for sname, (X, h, c, p) in SETS.items():
    cc = cols(X); Xc = X[cc]
    est = {}
    for lab, ii in [('healthy', h), ('HCC', c), ('hepatitis', p)]:
        m, b, n = bmeans(Xc.loc[ii]); est[lab] = (m, b, n)
        sd = Xc.loc[ii].std(ddof=1).to_numpy()
        for k, col in enumerate(cc): means.append(dict(subset=sname, group=lab, stratum=col[0], metric=col[1], n=int(n[k]), mean=m[k], sd=sd[k]))
    shsd = Xc.loc[h].std(ddof=1).to_numpy()
    for cn, (a, b) in {'HCC_minus_healthy': ('HCC', 'healthy'), 'hepatitis_minus_healthy': ('hepatitis', 'healthy'), 'HCC_minus_hepatitis': ('HCC', 'hepatitis')}.items():
        d = est[a][0] - est[b][0]; db = est[a][1] - est[b][1]; lo, hi = np.nanpercentile(db, [2.5, 97.5], axis=0)
        for k, col in enumerate(cc):
            rows.append(dict(subset=sname, contrast=cn, stratum=col[0], metric=col[1], diff=d[k], lo=lo[k], hi=hi[k], n_a=int(est[a][2][k]), n_b=int(est[b][2][k]), healthy_SD=shsd[k], diff_SD=d[k] / shsd[k] if shsd[k] > 0 else np.nan))
E_ = pd.DataFrame(rows); E_.to_csv(RES + '/hep_contrasts.tsv', sep='\t', index=False)
pd.DataFrame(means).to_csv(RES + '/hep_group_means.tsv', sep='\t', index=False)
# 방향 일치 집계 (59층)
for sname in SETS:
    for m in METRICS:
        e = E_[(E_.subset == sname) & (E_.metric == m) & E_.stratum.isin(STR59)]
        if m not in set(E_[E_.subset == sname].metric): continue
        w = e.pivot(index='stratum', columns='contrast', values=['diff', 'lo', 'hi', 'diff_SD'])
        H, Hp, C = 'HCC_minus_healthy', 'hepatitis_minus_healthy', 'HCC_minus_hepatitis'
        ok = w['diff'][[H, Hp]].notna().all(1); sig = lambda c: (w['lo'][c] > 0) | (w['hi'][c] < 0)
        sH, sP = np.sign(w['diff'][H]), np.sign(w['diff'][Hp]); both = ok & sig(H) & sig(Hp)
        rho = stats.spearmanr(w['diff_SD'][H][ok], w['diff_SD'][Hp][ok])[0] if ok.sum() > 5 else np.nan
        dirrows.append(dict(subset=sname, metric=m, n_strata=int(ok.sum()), same_sign=int((ok & (sH == sP)).sum()), opposite_sign=int((ok & (sH != sP)).sum()),
                            HCC_CI_excl0=int((ok & sig(H)).sum()), HCC_up=int((ok & sig(H) & (w['diff'][H] > 0)).sum()), HCC_down=int((ok & sig(H) & (w['diff'][H] < 0)).sum()),
                            hep_CI_excl0=int((ok & sig(Hp)).sum()), hep_up=int((ok & sig(Hp) & (w['diff'][Hp] > 0)).sum()), hep_down=int((ok & sig(Hp) & (w['diff'][Hp] < 0)).sum()),
                            both_CI_excl0=int(both.sum()), both_same_sign=int((both & (sH == sP)).sum()), both_opposite_sign=int((both & (sH != sP)).sum()),
                            HCCvsHep_CI_excl0=int((ok & sig(C)).sum()), HCCvsHep_up=int((ok & sig(C) & (w['diff'][C] > 0)).sum()), HCCvsHep_down=int((ok & sig(C) & (w['diff'][C] < 0)).sum()),
                            spearman_effect_vectors_HCC_vs_hep=rho))
pd.DataFrame(dirrows).to_csv(RES + '/hep_direction_counts.tsv', sep='\t', index=False)
# log_depth 보정 OLS (full 집합, 세 군 동시): y ~ HCC + hep + log_depth ; 군별 bootstrap 1000
h, c, p = SETS['full'][1:]
allids = h + c + p; gdum = np.column_stack([np.r_[np.zeros(len(h)), np.ones(len(c)), np.zeros(len(p))], np.r_[np.zeros(len(h)), np.zeros(len(c)), np.ones(len(p))]])
dep = S.loc[allids, 'log_depth'].to_numpy(float); grp = np.r_[np.zeros(len(h)), np.ones(len(c)), 2 * np.ones(len(p))]
Xd = np.column_stack([np.ones(len(allids)), gdum, dep])
adj = []; gi = [np.flatnonzero(grp == k) for k in range(3)]
for col in cols(W):
    if col[0] not in STR59: continue
    y = W.loc[allids, col].to_numpy(float); ok = np.isfinite(y)
    def fit(ix):
        ix = ix[ok[ix]]; return np.linalg.lstsq(Xd[ix], y[ix], rcond=None)[0]
    b0 = fit(np.arange(len(y))); bb = np.array([fit(np.concatenate([rng.choice(gi[k], len(gi[k])) for k in range(3)])) for _ in range(300)])
    for nmc, v, vb in [('HCC_minus_healthy', b0[1], bb[:, 1]), ('hepatitis_minus_healthy', b0[2], bb[:, 2]), ('HCC_minus_hepatitis', b0[1] - b0[2], bb[:, 1] - bb[:, 2])]:
        adj.append(dict(stratum=col[0], metric=col[1], contrast=nmc, diff_depthadj=v, lo=np.percentile(vb, 2.5), hi=np.percentile(vb, 97.5)))
pd.DataFrame(adj).to_csv(RES + '/hep_contrasts_depthadj.tsv', sep='\t', index=False)
# 표본 구성 이력
comp = S.loc[allids].assign(in_prior=[i in p103 for i in allids]).groupby(['group', 'split']).agg(n=('age', 'size'), age_mean=('age', 'mean'), male_frac=('sex', lambda x: (x == 1).mean()), log_depth_mean=('log_depth', 'mean'), n_prior=('in_prior', 'sum'))
comp.to_csv(RES + '/hep_sample_composition.tsv', sep='\t')
print(comp)
print(pd.DataFrame(dirrows).round(2).to_string())
