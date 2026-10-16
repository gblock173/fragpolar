"""component-origin step 분석: 두 성분(P 양극화, H 균질화)의 59층 효과 벡터 vs 기원 후보 시뮬레이션/조직 차이.
단위: 모든 벡터는 건강인 Train 표본 SD(층, 지표별) 단위(z)로 맞춘다. 유사도 = 코사인, 층별 Pearson, 최소제곱 배율 alpha_hat.
usage: co5_analysis.py [full|t167]   (PAT 변형; 기본 full)
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, os, json, glob, re, warnings; warnings.filterwarnings('ignore')
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.stats import spearmanr

VAR = sys.argv[1] if len(sys.argv) > 1 else 'full'
SUF = '' if VAR == 'full' else '_t'
D = CFG.DATA; AN = CFG.WORK; E = AN + '/component_origin'; THA_DIR = AN + '/tumour_hepatitis_age'; PB_DIR = AN + '/pbmc'
RES = E + '/results'; NM = 190218; SEED = 20261014
rng = np.random.default_rng(SEED)
st = pd.read_parquet(CFG.inp('marker_annotation'))
sd = pd.read_csv(CFG.inp('strata_definition'), sep='\t'); sd = sd[sd.used]
masks = []
for r in sd.itertuples():
    bm = np.ones(NM, bool) if r.baseline == 'All' else (st.band == r.baseline).to_numpy()
    am = np.ones(NM, bool) if r.axis == 'all' else (st.annot == r.annotation).to_numpy() if r.axis == 'annot' else \
        (st.genomic_context == r.annotation).to_numpy() if r.axis == 'context' else st['atlas_' + r.annotation].to_numpy()
    masks.append(bm & am)
NAMES = sd.stratum.tolist(); assert len(NAMES) == 59
S = sp.csr_matrix(np.stack(masks).astype(np.float64))
b0 = st.baseline_fullTrain.to_numpy(np.float64); bsel = np.where(b0 < .5, 1, np.where(b0 > .5, 0, -1))
MET = ['opp', 'conc', 'strictM', 'strictU', 'mixed', 'U', 'B', 'BV', 'mu']
PSET = ['opp', 'B', 'BV']; HSET = ['mixed', 'U', 'conc', 'strictM']
mi = {m: i for i, m in enumerate(MET)}


def strata_from_moments(mu, p2, sM, sU, valid):
    """marker 수준 모멘트 -> (59, len(MET)) 층 값. 규칙: 층 내 유효 marker >= 20, marker 동일 가중 평균 (표본 수준 규칙과 동일)."""
    v = valid & np.isfinite(mu) & np.isfinite(p2) & np.isfinite(sM) & np.isfinite(sU)
    U = np.maximum(mu - p2, 0); B = np.maximum(p2 - mu * mu, 0)
    opp = np.where(bsel == 1, sM, np.where(bsel == 0, sU, np.nan)); conc = np.where(bsel == 1, sU, np.where(bsel == 0, sM, np.nan))
    vo = v & np.isfinite(opp)
    out = np.full((59, len(MET)), np.nan)
    cnt = S @ v.astype(float); cnto = S @ vo.astype(float)
    def mean(x, vm, c): return np.where(c >= 20, (S @ np.where(vm, np.nan_to_num(x), 0)) / np.maximum(c, 1), np.nan)
    d = dict(opp=mean(opp, vo, cnto), conc=mean(conc, vo, cnto), strictM=mean(sM, v, cnt), strictU=mean(sU, v, cnt), mixed=mean(1 - sM - sU, v, cnt),
             U=mean(U, v, cnt), B=mean(B, v, cnt), mu=mean(mu, v, cnt))
    d['BV'] = d['B'] / (d['U'] + d['B'])
    for m in MET: out[:, mi[m]] = d[m]
    return out


# ---------- 관측: 표본 x 층 x 지표 ----------
L = pd.read_parquet(f'{THA_DIR}/results/sample_strata_long.parquet'); L = L[L.stratum.isin(NAMES)]
W = L.pivot_table(index='sid', columns=['stratum', 'metric'], values='value', dropna=False)
X = np.full((len(W), 59, len(MET)), np.nan)
for j, s in enumerate(NAMES):
    for m in MET:
        X[:, j, mi[m]] = (W[(s, 'B')] / W[(s, 'V')]).to_numpy() if m == 'BV' else W[(s, m)].to_numpy()
sid = W.index.to_numpy()
ms = pd.read_csv(CFG.inp('samples_table'), sep='\t').set_index('AnalysisID').loc[sid]
fs = pd.read_csv(f'{THA_DIR}/results/tf_f_samples.tsv', sep='\t').set_index('AnalysisID')
isH = ((ms.group == 'healthy') & (ms.split == 'Train')).to_numpy(); isHCC = (ms.group == 'liver').to_numpy()
isF = isHCC & np.isin(sid, fs.index); isNF = isHCC & ~isF
isHep = ((ms.group == 'hepatitis') & ms.index.isin(CFG.ids('hepatitis_included'))).to_numpy()   # 간염: 시료 품질 기준을 충족한 표본만
assert isH.sum() == 673 and isHCC.sum() == 177 and isF.sum() == 46 and isNF.sum() == 131 and isHep.sum() == 10, (isH.sum(), isHCC.sum(), isF.sum(), isNF.sum(), isHep.sum())
SDH = np.nanstd(X[isH], axis=0, ddof=1); MUH = np.nanmean(X[isH], axis=0)
fvals = fs.loc[sid[isF], 'f'].to_numpy(); Xf = X[isF]


def slope_f(Xs, f):
    fc = f - f.mean(); w = np.isfinite(Xs); y = np.where(w, Xs, 0); yb = y.sum(0) / np.maximum(w.sum(0), 1)
    return (fc[:, None, None] * np.where(w, Xs - yb, 0)).sum(0) / (fc ** 2).sum()


def obs_vectors(idx_nf, idx_f, idx_all, idx_hep, Xh, Xf_, Xnf, Xall, Xhep, f):
    mh = np.nanmean(Xh, 0)
    return {'P_slope': slope_f(Xf_, f) * 0.1 / SDH,                    # f 0.1 증가당 건강인 SD
            'P_f46delta': (np.nanmean(Xf_, 0) - mh) / SDH,
            'H_nof131': (np.nanmean(Xnf, 0) - mh) / SDH,
            'H_HCCall': (np.nanmean(Xall, 0) - mh) / SDH,
            'Hepatitis_10': (np.nanmean(Xhep, 0) - mh) / SDH}


OBS = obs_vectors(None, None, None, None, X[isH], Xf, X[isNF], X[isHCC], X[isHep], fvals)
# PBMC - plasma (PBMC step)
Lc = pd.read_parquet(f'{PB_DIR}/results/sample_strata_long.parquet'); Lc = Lc[Lc.stratum.isin(NAMES)]
Wc = Lc[Lc.kind == 'PBMC'].pivot_table(index='sid', columns=['stratum', 'metric'], values='value', dropna=False)
XP = np.full((len(Wc), 59, len(MET)), np.nan)
for j, s in enumerate(NAMES):
    for m in MET: XP[:, j, mi[m]] = (Wc[(s, 'B')] / Wc[(s, 'V')]).to_numpy() if m == 'BV' else Wc[(s, m)].to_numpy()
OBS['PBMC_minus_plasma'] = (np.nanmean(XP, 0) - np.nanmean(X[isH], 0)) / SDH


# ---------- 원천 분포 (marker 수준 모멘트) ----------
def moments_from_df(d, suf='', nmin=3):
    n = d['n_reads' + suf].to_numpy(float)
    with np.errstate(invalid='ignore', divide='ignore'):
        mu = d['sum_mean_meth' + suf].to_numpy() / n; p2 = d['sum_mean_meth_sq' + suf].to_numpy() / n
        sM = d['n_strictM' + suf].to_numpy() / n; sU = d['n_strictU' + suf].to_numpy() / n
        K = d['sum_ncpg' + suf].to_numpy() / n
    ok = n >= nmin
    return n, np.where(ok, mu, np.nan), np.where(ok, p2, np.nan), np.where(ok, sM, np.nan), np.where(ok, sU, np.nan), np.where(ok, K, np.nan)


def pool(frames, suf, nmin_each=3, nmin_pool=10):
    """표본(공여자/library) 동일 가중 평균. 반환: dict(mu,p2,sM,sU,K,valid,n)."""
    A = [moments_from_df(f, suf, nmin_each) for f in frames]
    n = np.sum([np.nan_to_num(a[0]) * np.isfinite(a[1]) for a in A], 0)
    out = {k: np.nanmean(np.stack([a[i] for a in A]), 0) for k, i in [('mu', 1), ('p2', 2), ('sM', 3), ('sU', 4), ('K', 5)]}
    out['valid'] = n >= nmin_pool; out['n'] = n
    return out


def load_pat(ids):
    return [pd.read_parquet(f'{E}/data/pat/{i}.parquet') for i in ids]


tasks = pd.read_csv(f'{E}/data/pat_tasks.tsv', sep='\t')
tasks = tasks[[os.path.exists(f'{E}/data/pat/{i}.parquet') for i in tasks.id]]
SRC = {}
for ct, g in tasks.groupby('celltype'):
    SRC[ct] = pool(load_pat(g.id), SUF if True else '')
# 조직
tt = pd.read_csv(f'{E}/data/tissue_tasks.tsv', sep='\t', dtype=str, keep_default_na=False)
tt = tt[[os.path.exists(f'{E}/data/tissue_suff/{s}.parquet') for s in tt.sid]]
TIS = {s: pd.read_parquet(f'{E}/data/tissue_suff/{s}.parquet').set_index('marker_id').reindex(st.marker_id).reset_index() for s in tt.sid}
for role, nm in [('tumour', 'TISSUE_tumour'), ('normal', 'TISSUE_nontumour')]:
    ids = tt[tt.role == role].sid.tolist()
    if ids: SRC[nm] = pool([TIS[i] for i in ids], '', 5, 20)

# 기준: 건강인 Train marker 평균
G = pd.read_parquet(f'{E}/data/marker_obs_groups.parquet')
BASE = dict(mu=G['H_Train__mu'].to_numpy(), p2=(G['H_Train__mu'] - G['H_Train__U']).to_numpy(), sM=G['H_Train__sM'].to_numpy(), sU=G['H_Train__sU'].to_numpy(), valid=(G['H_Train__nsamp'].to_numpy() >= 300))


def mix(src, alpha, ref=None):
    """분포 수준 혼합: M' = M_base + alpha (M_src - M_ref). ref=None이면 순수 혼합 (1-alpha) base + alpha src. 선형 모멘트(mu, p2, sM, sU), U, B는 모멘트에서 재계산 -> B에 성분 간 분산항 포함."""
    s = SRC[src]; r = BASE if ref is None else SRC[ref]
    valid = BASE['valid'] & s['valid'] & (r['valid'] if ref is not None else True)
    out = {}
    for k in ['mu', 'p2', 'sM', 'sU']:
        rk = BASE[k] if ref is None else r[k]
        out[k] = BASE[k] + alpha * (s[k] - rk)
    return out, valid


def sim_delta(src, alpha, ref=None):
    m, v = mix(src, alpha, ref)
    a = strata_from_moments(m['mu'], m['p2'], m['sM'], m['sU'], v); b = strata_from_moments(BASE['mu'], BASE['p2'], BASE['sM'], BASE['sU'], v)
    return (a - b) / SDH, v


# ---------- 유사도 ----------
def flat(Z, mets): return Z[:, [mi[m] for m in mets]].reshape(-1)


def sim_stats(zo, zs, mets):
    a, b = flat(zo, mets), flat(zs, mets); ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    cos = float(a @ b / np.sqrt((a @ a) * (b @ b))) if (a @ a) > 0 and (b @ b) > 0 else np.nan
    k = float(a @ b / (b @ b)) if (b @ b) > 0 else np.nan
    r2 = float(1 - ((a - k * b) ** 2).sum() / (a @ a)) if np.isfinite(k) else np.nan
    return cos, k, r2


def per_metric(zo, zs, mets):
    o = {}
    for m in mets:
        a, b = zo[:, mi[m]], zs[:, mi[m]]; ok = np.isfinite(a) & np.isfinite(b)
        o[m + '_pearson'] = float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() > 5 and np.std(b[ok]) > 0 and np.std(a[ok]) > 0 else np.nan
        o[m + '_signagree'] = float((np.sign(a[ok]) == np.sign(b[ok])).mean()) if ok.sum() else np.nan
    return o


# 부트스트랩 (표본 재추출)로 관측 벡터 분포
def boot_obs(B=1000):
    out = {k: [] for k in OBS if k != 'PBMC_minus_plasma'}
    Xh, Xnf, Xall, Xhep = X[isH], X[isNF], X[isHCC], X[isHep]
    for _ in range(B):
        ih = rng.integers(0, len(Xh), len(Xh)); inf = rng.integers(0, len(Xnf), len(Xnf)); ia = rng.integers(0, len(Xall), len(Xall))
        iF = rng.integers(0, len(Xf), len(Xf)); ihp = rng.integers(0, len(Xhep), len(Xhep))
        v = obs_vectors(None, None, None, None, Xh[ih], Xf[iF], Xnf[inf], Xall[ia], Xhep[ihp], fvals[iF])
        for k in out: out[k].append(v[k])
    return {k: np.stack(v) for k, v in out.items()}


BO = boot_obs(1000)
COMP = {'P': ('P_slope', PSET), 'P_alt(f46 Delta)': ('P_f46delta', PSET), 'H': ('H_nof131', HSET), 'H_alt(HCC all)': ('H_HCCall', HSET)}

# ---------- 후보 목록 ----------
cands = {}   # name -> (zvec, group, params)
def add(name, z, grp, **kw): cands[name] = (z, grp, kw)
ALPHA = 0.02
for a in [0.005, 0.01, 0.02, 0.03, 0.05]:
    if 'Liver-Hepatocytes' in SRC: add(f'hepatocyte_a{a}', sim_delta('Liver-Hepatocytes', a)[0], '1 hepatocyte', alpha=a)
blood_types = [c for c in ['Blood-Granulocytes', 'Blood-Monocytes', 'Blood-B', 'Blood-NK', 'Blood-T-CD4', 'Blood-T-CD8'] if c in SRC]
for c in blood_types: add(f'{c}_a{ALPHA}', sim_delta(c, ALPHA)[0], '2 blood', alpha=ALPHA)
if 'CNVS_WBC' in SRC: add(f'WBC_bulk(CNVS)_a{ALPHA}', sim_delta('CNVS_WBC', ALPHA)[0], '2 blood', alpha=ALPHA)
# 혈구 구성: 호중구 <-> 림프구
LYM = [c for c in ['Blood-T-CD4', 'Blood-T-CD8', 'Blood-B', 'Blood-NK'] if c in SRC]
if len(LYM) == 4 and 'Blood-Granulocytes' in SRC:
    w = dict(zip(LYM, [.13, .08, .06, .05])); tot = sum(w.values())
    SRC['LYM_pool'] = {k: sum(w[c] * SRC[c][k] for c in LYM) / tot for k in ['mu', 'p2', 'sM', 'sU', 'K']}
    SRC['LYM_pool']['valid'] = np.all([SRC[c]['valid'] for c in LYM], 0); SRC['LYM_pool']['n'] = 0
    for d_ in [0.05, 0.10]:
        add(f'neutrophil_up_lymphocyte_down_d{d_}', sim_delta('Blood-Granulocytes', d_, ref='LYM_pool')[0], '2 blood composition', alpha=d_)
        add(f'lymphocyte_up_neutrophil_down_d{d_}', sim_delta('LYM_pool', d_, ref='Blood-Granulocytes')[0], '2 blood composition', alpha=d_)
if 'Blood-Monocytes' in SRC and 'Blood-Granulocytes' in SRC:
    add('monocyte_up_neutrophil_down_d0.05', sim_delta('Blood-Monocytes', .05, ref='Blood-Granulocytes')[0], '2 blood composition', alpha=.05)
for a in [0.02, 0.10]:
    if 'TISSUE_tumour' in SRC: add(f'tumour_tissue_a{a}', sim_delta('TISSUE_tumour', a)[0], '3 tumour tissue', alpha=a)
if 'TISSUE_tumour' in SRC and 'TISSUE_nontumour' in SRC:
    add('tumour_minus_nontumour_a0.10', sim_delta('TISSUE_tumour', .10, ref='TISSUE_nontumour')[0], '3 tumour-specific', alpha=.10)
if 'TISSUE_nontumour' in SRC: add(f'nontumour_liver_tissue_a{ALPHA}', sim_delta('TISSUE_nontumour', ALPHA)[0], '4 nontumour liver tissue', alpha=ALPHA)
OTHER = [c for c in tasks.celltype.unique() if not c.startswith(('Blood', 'CNVS', 'Liver-Hep'))]
for c in OTHER:
    if c in SRC: add(f'{c}_a{ALPHA}', sim_delta(c, ALPHA)[0], '5 other cell type (control)', alpha=ALPHA)
# 관측 비교 벡터
for k in ['Hepatitis_10', 'PBMC_minus_plasma', 'P_slope', 'H_nof131']: add('OBS_' + k, OBS[k], '6 observed', alpha=np.nan)

# 조직 직접 차이 (3),(4): 층 값 (조직 library별), 쌍 차이
def tissue_strata(sid_):
    d = TIS[sid_]; n, mu, p2, sM, sU, K = moments_from_df(d, '', 5)
    return strata_from_moments(mu, p2, sM, sU, np.isfinite(mu)), np.nanmean(K)
TS = {s: tissue_strata(s) for s in TIS}
tt_ = tt.set_index('sid'); mkd = lambda role: {s: v[0] for s, v in TS.items() if tt_.loc[s, 'role'] == role}
TUM, NOR = mkd('tumour'), mkd('normal')
donor = tt_.donor
pairs = []
for dn in sorted(set(donor[donor != ''])):
    ts = [s for s in TUM if donor[s] == dn]; ns = [s for s in NOR if donor[s] == dn]
    if ts and ns: pairs.append((dn, ts, ns))
print('tissue pairs', len(pairs), 'tumour libs', len(TUM), 'normal libs', len(NOR), flush=True)
dT = np.stack([np.nanmean([TUM[s] for s in ts], 0) - np.nanmean([NOR[s] for s in ns], 0) for dn, ts, ns in pairs])   # (pairs,59,met)
TN = np.nanmean(dT, 0) / SDH
add('TISSUE_tumour_minus_adjacent(paired)', TN, '3 tissue direct', alpha=np.nan)
NP = (np.nanmean(list(NOR.values()), 0) - MUH) / SDH
add('TISSUE_nontumour_minus_plasma', NP, '4 tissue direct', alpha=np.nan)
if 'Liver-Hepatocytes' in SRC:
    m, v = mix('Liver-Hepatocytes', 1.0); hp = strata_from_moments(m['mu'], m['p2'], m['sM'], m['sU'], v)
    add('hepatocyte_minus_plasma(100%)', (hp - strata_from_moments(BASE['mu'], BASE['p2'], BASE['sM'], BASE['sU'], v)) / SDH, '4 tissue direct', alpha=1.0)
# 조직 쌍 차이 부트스트랩
dTb = []
for _ in range(1000):
    ix = rng.integers(0, len(dT), len(dT)); dTb.append(np.nanmean(dT[ix], 0) / SDH)
dTb = np.stack(dTb)
NORb = []
nl = np.stack(list(NOR.values()))
for _ in range(1000):
    ix = rng.integers(0, len(nl), len(nl)); NORb.append((np.nanmean(nl[ix], 0) - MUH) / SDH)
NORb = np.stack(NORb)

rows = []
for cname, (z, grp, kw) in cands.items():
    for comp, (ok_, mets) in COMP.items():
        zo = OBS[ok_]; cos, k, r2 = sim_stats(zo, z, mets)
        # 부트스트랩 CI: 관측 벡터 재추출 (후보 벡터 고정; 조직 직접 차이는 같이 재추출)
        bo = BO[ok_]; cb = []
        for b in range(len(bo)):
            zs_b = z
            if cname.startswith('TISSUE_tumour_minus_adjacent'): zs_b = dTb[b]
            elif cname == 'TISSUE_nontumour_minus_plasma': zs_b = NORb[b]
            cb.append(sim_stats(bo[b], zs_b, mets)[0])
        lo, hi = np.nanpercentile(cb, [2.5, 97.5])
        rows.append(dict(candidate=cname, group=grp, component=comp, cosine=cos, cos_lo=lo, cos_hi=hi, alpha_hat_per_unit=k * kw.get('alpha', np.nan) if np.isfinite(kw.get('alpha', np.nan)) else k,
                         R2_through_origin=r2, **per_metric(zo, z, mets)))
Rdf = pd.DataFrame(rows); Rdf.insert(0, 'pat_variant', VAR)
Rdf.to_csv(f'{RES}/similarity_{VAR}.tsv', sep='\t', index=False)
# 벡터 저장
vec = []
for nm, z in list(OBS.items()) + [(c, v[0]) for c, v in cands.items()]:
    for j, s in enumerate(NAMES):
        vec.append(dict(vector=nm, stratum=s, **{m: z[j, mi[m]] for m in MET}))
pd.DataFrame(vec).drop_duplicates(['vector', 'stratum']).to_csv(f'{RES}/effect_vectors_z_{VAR}.tsv', sep='\t', index=False)

# 관측 벡터 사이 비교
oo = []
for a_, ma in [('P_slope', PSET), ('H_nof131', HSET)]:
    for b_ in ['Hepatitis_10', 'PBMC_minus_plasma', 'P_slope', 'H_nof131', 'H_HCCall']:
        if a_ == b_: continue
        cos, k, r2 = sim_stats(OBS[a_], OBS[b_], ma)
        bo = BO.get(b_); cb = [sim_stats(BO[a_][i], BO[b_][i] if b_ in BO else OBS[b_], ma)[0] for i in range(1000)]
        oo.append(dict(component=a_, vs=b_, metrics=','.join(ma), cosine=cos, lo=np.nanpercentile(cb, 2.5), hi=np.nanpercentile(cb, 97.5), **per_metric(OBS[a_], OBS[b_], ma)))
pd.DataFrame(oo).to_csv(f'{RES}/observed_vs_observed_{VAR}.tsv', sep='\t', index=False)

# 그룹 평균 층 합의: 성분 안의 부호 방향 (층 수)
dirs = []
for k_, mets in [('P_slope', PSET), ('P_f46delta', PSET), ('H_nof131', HSET), ('H_HCCall', HSET)]:
    for m in mets:
        v = OBS[k_][:, mi[m]]; dirs.append(dict(obs=k_, metric=m, n_pos=int((v > 0).sum()), n_neg=int((v < 0).sum()), median_z=float(np.nanmedian(v))))
pd.DataFrame(dirs).to_csv(f'{RES}/observed_component_directions.tsv', sep='\t', index=False)

# ---------- 대조군 특이성: 단일 세포 유형 2% 혼합 중 간세포 코사인의 순위 ----------
single = Rdf[Rdf.candidate.str.endswith('_a0.02') & Rdf.group.isin(['1 hepatocyte', '2 blood', '5 other cell type (control)']) & ~Rdf.candidate.str.contains('CNVS')].copy()
single['rank_among_single_celltypes'] = single.groupby('component').cosine.rank(ascending=False)
single['n_single'] = single.groupby('component').cosine.transform('size')
single.to_csv(f'{RES}/specificity_single_celltype_rank_{VAR}.tsv', sep='\t', index=False)
# ---------- 대조군 특이성(옛 표) ----------
ctrl = Rdf[(Rdf.group.str.startswith('5')) | (Rdf.candidate == f'hepatocyte_a{ALPHA}') | (Rdf.group == '2 blood')]
ctrl = ctrl[ctrl.candidate.str.endswith(f'a{ALPHA}') | ctrl.candidate.str.contains('CNVS')]
ctrl.to_csv(f'{RES}/specificity_controls_{VAR}.tsv', sep='\t', index=False)

# ---------- marker 수준 ----------
if VAR == 'full':
    g = lambda k: G[k].to_numpy()
    def mk_metrics(prefix):
        sM, sU = g(prefix + '__sM'), g(prefix + '__sU'); mu = g(prefix + '__mu'); U = g(prefix + '__U'); Bv = g(prefix + '__B')
        opp = np.where(bsel == 1, sM, np.where(bsel == 0, sU, np.nan)); conc = np.where(bsel == 1, sU, np.where(bsel == 0, sM, np.nan))
        return dict(opp=opp, conc=conc, strictM=sM, mixed=1 - sM - sU, U=U, B=Bv, mu=mu)
    mH = mk_metrics('H_Train')
    def mdelta(grp): m = mk_metrics(grp); return {k: m[k] - mH[k] for k in m}
    mo = {'HCC_nof131': mdelta('HCC_nof131'), 'HCC_all': mdelta('HCC_all'), 'HCC_Train': mdelta('HCC_Train'), 'HCC_Test': mdelta('HCC_Test'), 'PBMC': mdelta('PBMC'), 'Hepatitis': mdelta('Hepatitis')}
    FS = pd.read_parquet(f'{E}/data/marker_obs_fslope.parquet')
    sM_s, sU_s = FS['slope_sM'].to_numpy(), FS['slope_sU'].to_numpy()
    mo['P_slope'] = dict(opp=np.where(bsel == 1, sM_s, np.where(bsel == 0, sU_s, np.nan)), B=FS['slope_B'].to_numpy(), mixed=-(sM_s + sU_s), U=FS['slope_U'].to_numpy(), conc=np.where(bsel == 1, sU_s, np.where(bsel == 0, sM_s, np.nan)), strictM=sM_s, mu=FS['slope_mu'].to_numpy())
    def sim_marker(src, alpha, ref=None):
        m, v = mix(src, alpha, ref); sM, sU = m['sM'], m['sU']; mu = m['mu']; p2 = m['p2']
        U = np.maximum(mu - p2, 0); Bv = np.maximum(p2 - mu * mu, 0)
        bU = np.maximum(BASE['mu'] - BASE['p2'], 0); bB = np.maximum(BASE['p2'] - BASE['mu'] ** 2, 0)
        opp = lambda a, b: np.where(bsel == 1, a, np.where(bsel == 0, b, np.nan))
        r = dict(opp=opp(sM, sU) - opp(BASE['sM'], BASE['sU']), conc=opp(sU, sM) - opp(BASE['sU'], BASE['sM']), strictM=sM - BASE['sM'], mixed=-(sM - BASE['sM']) - (sU - BASE['sU']), U=U - bU, B=Bv - bB, mu=mu - BASE['mu'])
        return {k: np.where(v, x, np.nan) for k, x in r.items()}
    ms_rows = []
    sims = {'hepatocyte': ('Liver-Hepatocytes', None), 'WBC_bulk': ('CNVS_WBC', None), 'granulocyte': ('Blood-Granulocytes', None), 'tumour_tissue': ('TISSUE_tumour', None),
            'nontumour_tissue': ('TISSUE_nontumour', None), 'neutrophil_vs_lymph': ('Blood-Granulocytes', 'LYM_pool'), 'tumour_minus_nontumour': ('TISSUE_tumour', 'TISSUE_nontumour'),
            'monocyte_vs_neutrophil': ('Blood-Monocytes', 'Blood-Granulocytes')}
    for c_ in OTHER + ['Blood-Monocytes', 'Blood-B', 'Blood-NK', 'Blood-T-CD4', 'Blood-T-CD8']: sims['ctrl:' + c_] = (c_, None)
    for sname, (src, ref) in sims.items():
        if src not in SRC or (ref and ref not in SRC): continue
        sm_ = sim_marker(src, 0.02, ref)
        for ok_, mets in [('P_slope', PSET[:2]), ('HCC_nof131', ['mixed', 'U', 'conc', 'strictM']), ('HCC_all', ['opp', 'B', 'mixed', 'U', 'conc', 'strictM']), ('PBMC', ['mixed', 'U', 'conc', 'strictM', 'opp', 'B']), ('Hepatitis', ['mixed', 'U', 'conc', 'strictM', 'opp', 'B'])]:
            for m in mets:
                a, b = mo[ok_][m], sm_[m]; okm = np.isfinite(a) & np.isfinite(b)
                if okm.sum() > 1000: ms_rows.append(dict(sim=sname, obs=ok_, metric=m, n_markers=int(okm.sum()), spearman=float(spearmanr(a[okm], b[okm])[0])))
    # 신뢰도 상한: HCC Train vs Test 효과의 marker 상관
    for m in ['opp', 'B', 'mixed', 'U', 'conc', 'strictM']:
        a, b = mo['HCC_Train'][m], mo['HCC_Test'][m]; okm = np.isfinite(a) & np.isfinite(b)
        ms_rows.append(dict(sim='(reliability) HCC Train vs Test', obs='HCC', metric=m, n_markers=int(okm.sum()), spearman=float(spearmanr(a[okm], b[okm])[0])))
    # 조직 직접 marker 차이: 종양-정상 (library 동일 가중 평균), 비종양-혈장
    pd.DataFrame(ms_rows).to_csv(f'{RES}/marker_level_spearman.tsv', sep='\t', index=False)

# ---------- 파이프라인 검증 ----------
if VAR == 'full':
    val = []
    def mk_src(name):
        s = SRC[name]; U = np.maximum(s['mu'] - s['p2'], 0); Bv = np.maximum(s['p2'] - s['mu'] ** 2, 0)
        return dict(mu=s['mu'], U=U, B=Bv, strictM=s['sM'], mixed=1 - s['sM'] - s['sU'], K=s['K'])
    def corr(a, b):
        ok = np.isfinite(a) & np.isfinite(b); return float(spearmanr(a[ok], b[ok])[0]), int(ok.sum())
    if 'CNVS_cfDNA' in SRC:
        c = mk_src('CNVS_cfDNA'); h = mk_metrics('H_Train')
        for m in ['mu', 'B', 'U', 'mixed', 'strictM']:
            r, n = corr(c[m], h[m] if m in h else mH[m]); val.append(dict(test='CNVS cfDNA(PAT) vs 우리 건강인 Train 혈장', metric=m, spearman=r, n_markers=n))
    if 'CNVS_WBC' in SRC:
        c = mk_src('CNVS_WBC'); p = mk_metrics('PBMC')
        for m in ['mu', 'B', 'U', 'mixed', 'strictM']:
            r, n = corr(c[m], p[m]); val.append(dict(test='CNVS WBC(PAT) vs 우리 PBMC 30', metric=m, spearman=r, n_markers=n))
    if 'CNVS_cfDNA' in SRC and 'CNVS_WBC' in SRC:   # Loyfer 자체 WBC-cfDNA 층 벡터 vs 우리 PBMC-혈장 층 벡터
        z = sim_delta('CNVS_WBC', 1.0)[0]; cos, k, r2 = sim_stats(OBS['PBMC_minus_plasma'], z, HSET + PSET); val.append(dict(test='CNVS WBC 100%혼합 vs 우리 PBMC-혈장 (층 z 벡터, H+P 지표)', metric='cosine', spearman=cos, n_markers=np.nan))
    # K/길이
    kr = []
    for ct in SRC:
        s = SRC[ct]; kr.append(dict(source=ct, mean_K_marker_mean=float(np.nanmean(s['K'])), valid_markers=int(s['valid'].sum())))
    pd.DataFrame(kr).to_csv(f'{RES}/source_meanK_validmarkers.tsv', sep='\t', index=False)
    pd.DataFrame(val).to_csv(f'{RES}/pipeline_validation.tsv', sep='\t', index=False)
print('done', VAR)
