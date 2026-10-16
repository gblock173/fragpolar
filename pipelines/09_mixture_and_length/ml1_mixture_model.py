"""mixture step: "반대 균일 fragment가 BFV를 만든다"가 정의상 자명한가?  혼합 분산식으로 기대 ΔBFV를 계산하고 관측 ΔBFV와 비교한다.

입력(읽기 전용): marker_matrices(mu,U,B,strictM,strictU,n), marker annotation, samples table, strict_recompute(Train 422).
처리
 1) 표본 x marker 행렬에서 (군 x split) = healthy/HCC x Train/Test 별 marker 평균을 만든다(유효: n>=5 & 전 항 유한).
    marker 포함 조건: 두 군 모두 유효 표본이 군 크기의 50% 이상, baseline_fullTrain != 0.5 (opposite 방향 정의 가능).
 2) 건강인 기준 marker 분포(mu_h, WFV_h=U, BFV_h=B, strictM/U)에 반대 방향 완전 균일 fragment(p_o=1 if baseline<.5 else 0)를
    비율 δ 만큼 섞는다. δ = 그 층의 관측 ΔsOUF(층 평균, HCC-건강). 모형(전부 표본 단위 평균 기준):
    M1   : 건강 fragment 를 무작위로 δ 만큼 p_o 로 대체. ΔB = δ(1-δ)E_i[(p_o-mu_i)^2] - δ*B_h  (정확식, 표본별 mu_i 평균 반영)
    M1p  : M1 을 군 평균 분포(pooled)로 계산 (M1과의 차이 = 표본 간 mu 변동의 영향)
    M1m  : M1 이지만 δ 를 marker 별 관측 ΔsOUF_m(부호 있음)로 사용
    M1c  : 대조. 같은 δ 를 건강 방향 균일 fragment(p_c=1-p_o)로 넣음
    M2   : 건강 strict-mixed fragment 를 δ 만큼 opposite 균일로 바꿈 (혼합 class 의 평균/2차 모멘트는 건강 분포에서 역산)
    M2c  : M2 의 대조(건강 방향 균일로 바꿈)
    M3   : 관측된 strict 조성 변화 전부(ΔsOUF, Δ건강방향 균일 → strict-mixed 감소분)를 반영한 모형. 관측 BFV 와 비교
 3) marker -> 층(59층, 군 비교 단계와 동일 층 정의) 평균, ΔFPI = 층 평균 B / 층 평균 (U+B) 의 군 차이.
 4) 표본 bootstrap(군x split 내 복원추출) B회.
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, json, time
import numpy as np, pandas as pd, scipy.sparse as sp

D = CFG.DATA; V = CFG.WORK; OUT = V + '/mixture_and_length/results'
NS, NM = 2044, 190218
NBOOT = int(sys.argv[1]) if len(sys.argv) > 1 else 300
TEST = len(sys.argv) > 2 and sys.argv[2] == 'test'
t0 = time.time()
MM = {f: np.memmap(f'{CFG.inp("marker_matrices")}/{f}.f32', dtype=np.float32, mode='r', shape=(NS, NM)) for f in ['mu', 'U', 'B', 'strictM', 'strictU', 'n']}
cs = pd.read_csv(CFG.inp('marker_matrix_samples'), sep='\t', low_memory=False, usecols=['AnalysisID', 'group', 'split']); cs['col'] = np.arange(len(cs))
ms = pd.read_csv(CFG.inp('samples_table'), sep='\t')
ms = ms[ms.group.isin(['healthy', 'liver']) & ms.split.isin(['Train', 'Test'])]
sel = cs.merge(ms[['AnalysisID']], on='AnalysisID', how='inner').sort_values('col').reset_index(drop=True)
assert len(sel) == 1186
st = pd.read_parquet(CFG.inp('marker_annotation'))
sd = pd.read_csv(CFG.inp('strata_definition'), sep='\t'); sd = sd[sd.used]
masks = []
for r in sd.itertuples():
    bm = np.ones(NM, bool) if r.baseline == 'All' else (st.band == r.baseline).to_numpy()
    am = np.ones(NM, bool) if r.axis == 'all' else (st.annot == r.annotation).to_numpy() if r.axis == 'annot' else \
        (st.genomic_context == r.annotation).to_numpy() if r.axis == 'context' else st['atlas_' + r.annotation].to_numpy()
    masks.append(bm & am)
masks = np.stack(masks); STRATA = sd.stratum.tolist(); assert len(STRATA) == 59 and (masks.sum(1) == sd.n_markers.to_numpy()).all()
SM = sp.csr_matrix(masks.astype(np.float64))
b = st.baseline_fullTrain.to_numpy(np.float64)
dir_ok = np.isfinite(b) & (b != 0.5)
po = np.where(b < 0.5, 1.0, 0.0)          # 건강 기준 반대 방향 완전 균일 fragment 의 p
pc = 1.0 - po
rec_dir = V + '/tumour_hepatitis_age/results/strict_recompute/'
GROUPS = [('healthy', 'Train'), ('healthy', 'Test'), ('liver', 'Train'), ('liver', 'Test')]
F = 7   # mu, U, B, mu2, sM, sU, valid
X = {}; nsamp = {}; nfill = 0
for g in GROUPS:
    ss = sel[(sel.group == g[0]) & (sel.split == g[1])]
    if TEST: ss = ss.head(20)
    A = np.zeros((len(ss), F * NM), np.float32)
    for k, r in enumerate(ss.itertuples()):
        i = int(r.col); gf = lambda f: np.asarray(MM[f][i], dtype=np.float64)
        n, mu, U, B, sM, sU = gf('n'), gf('mu'), gf('U'), gf('B'), gf('strictM'), gf('strictU')
        if not np.isfinite(sM).any():     # Train strict 결측 표본: raw 재집계(strict_recompute)로 대체. n, mu 패리티 확인
            s = pd.read_parquet(rec_dir + r.AnalysisID + '.parquet').set_index('marker_id').reindex(st.marker_id)
            n2 = s.n_reads.to_numpy(float)
            with np.errstate(invalid='ignore', divide='ignore'):
                mu2 = s.sum_mean_meth.to_numpy() / n2; sM = s.n_strictM.to_numpy(float) / n2; sU = s.n_strictU.to_numpy(float) / n2
            ok = np.isfinite(n) & (n2 > 0); assert np.nanmax(np.abs(n2[ok] - n[ok])) == 0 and np.nanmax(np.abs(mu2[ok] - mu[ok])) < 1e-5
            nfill += 1
        v = (n >= 5) & np.isfinite(mu) & np.isfinite(U) & np.isfinite(B) & np.isfinite(sM) & np.isfinite(sU)
        for f, x in enumerate([mu, U, B, mu * mu, sM, sU, np.ones(NM)]):
            A[k, f * NM:(f + 1) * NM] = np.where(v, x, 0.0)
        A[k, 6 * NM:7 * NM] = v
    X[g] = A; nsamp[g] = len(ss)
    print('loaded', g, A.shape, round(time.time() - t0), flush=True)
print('strict filled', nfill, flush=True)


def gsum(g, w):
    s = (w.astype(np.float32) @ X[g]).astype(np.float64)
    return s.reshape(F, NM)


def means(s):
    cnt = s[6]
    with np.errstate(invalid='ignore', divide='ignore'):
        return cnt, s[0] / cnt, s[1] / cnt, s[2] / cnt, s[3] / cnt, s[4] / cnt, s[5] / cnt


def smean(x, elig):
    """층 평균 (유효 marker >=20)."""
    w = elig.astype(np.float64); cnt = SM @ w
    v = SM @ np.where(elig, x, 0.0)
    return np.where(cnt >= 20, v / np.maximum(cnt, 1), np.nan), cnt


def model(sh, sc, nh, nc, elig):
    """sh, sc: 건강/HCC 군 합 (F x NM). 반환: dict(name -> 59 vector). elig 는 고정."""
    _, mu_h, U_h, B_h, mu2_h, sM_h, sU_h = means(sh)
    _, mu_c, U_c, B_c, mu2_c, sM_c, sU_c = means(sc)
    opp_h = np.where(po == 1, sM_h, sU_h); con_h = np.where(po == 1, sU_h, sM_h)
    opp_c = np.where(po == 1, sM_c, sU_c); con_c = np.where(po == 1, sU_c, sM_c)
    el = elig & np.isfinite(mu_h + U_h + B_h + mu2_h + sM_h + sU_h + mu_c + U_c + B_c + sM_c + sU_c)
    z = lambda x: np.where(el, x, 0.0)
    out = {}
    out['n_markers'] = smean(np.ones(NM), el)[1]
    dopp_m = opp_c - opp_h; dcon_m = con_c - con_h
    d_o = smean(dopp_m, el)[0]; d_c = smean(dcon_m, el)[0]            # 층 δ (ΔsOUF, Δ건강방향 균일)
    out['delta_sOUF'] = d_o; out['delta_conc'] = d_c
    out['obs_dB'] = smean(B_c - B_h, el)[0]; out['obs_dU'] = smean(U_c - U_h, el)[0]
    out['obs_dmixed'] = smean(-(sM_c + sU_c) + (sM_h + sU_h), el)[0]
    out['obs_dmu'] = smean(mu_c - mu_h, el)[0]
    out['B_h'] = smean(B_h, el)[0]; out['U_h'] = smean(U_h, el)[0]; out['B_c'] = smean(B_c, el)[0]; out['U_c'] = smean(U_c, el)[0]
    fpi = lambda B, U: smean(B, el)[0] / (smean(B, el)[0] + smean(U, el)[0])
    out['FPI_h'] = fpi(B_h, U_h); out['obs_dFPI'] = fpi(B_c, U_c) - out['FPI_h']
    # 정확식 계수 (표본 단위 평균)
    Qo = po - 2 * po * mu_h + mu2_h; Qc = pc - 2 * pc * mu_h + mu2_h          # E_i[(p-mu_i)^2], p in {0,1} -> p^2=p
    Qo_p = (po - mu_h) ** 2; Qc_p = (pc - mu_h) ** 2
    Qo_s, Qc_s, Bh_s, Uh_s = smean(Qo, el)[0], smean(Qc, el)[0], smean(B_h, el)[0], smean(U_h, el)[0]
    Qop_s, Qcp_s = smean(Qo_p, el)[0], smean(Qc_p, el)[0]
    Bp_h = B_h + (mu2_h - mu_h ** 2)                                           # pooled BFV (= B_h + Var_i(mu_i))
    d = d_o
    out['exp_M1'] = d * (1 - d) * Qo_s - d * Bh_s
    out['exp_M1p'] = d * (1 - d) * Qop_s - d * smean(Bp_h, el)[0]
    out['exp_M1c'] = d * (1 - d) * Qc_s - d * Bh_s                              # 대조: 같은 δ 를 건강 방향 균일 fragment 로
    # M1m: marker 별 δ_m (부호 있음)
    out['exp_M1m'] = smean(dopp_m * (1 - dopp_m) * Qo - dopp_m * B_h, el)[0]
    out['exp_M1_dU'] = -d * Uh_s                                                # WFV 기대 변화
    Bexp = Bh_s + out['exp_M1']; Uexp = Uh_s * (1 - d)
    out['exp_M1_dFPI'] = Bexp / (Bexp + Uexp) - out['FPI_h']
    out['exp_M1c_dFPI'] = (Bh_s + out['exp_M1c']) / ((Bh_s + out['exp_M1c']) + Uh_s * (1 - d)) - out['FPI_h']
    # M2 / M3 : 건강 strict-mixed class 에서 역산
    a_o, a_c = opp_h, con_h; a_x = 1 - a_o - a_c
    E2 = mu_h - U_h                                                            # E[p^2] (pooled) = mu - U
    okx = el & (a_x > 0.02)
    with np.errstate(invalid='ignore', divide='ignore'):
        mux = (mu_h - a_o * po - a_c * pc) / a_x; s2x = (E2 - a_o * po - a_c * pc) / a_x
    Bp0 = E2 - mu_h ** 2
    def comp(do_s, dc_s, tag):
        # 층마다 δ 가 다르므로 층 루프 (59 x 190k, 가벼움)
        res = np.full(59, np.nan); resU = np.full(59, np.nan); cap = np.full(59, np.nan); fb = np.full(59, np.nan)
        for s in range(59):
            m = masks[s] & okx
            if m.sum() < 20 or not np.isfinite(do_s[s]) or not np.isfinite(dc_s[s]): continue
            do_, dc_ = do_s[s], dc_s[s]
            # 층 δ 가 음수이면 (건강 쪽으로 감소) 해당 class 에서 되돌리는 것으로 형식 외삽한다
            tot = max(do_, 0) + max(dc_, 0)
            sc_ = np.minimum(1.0, a_x[m] / np.maximum(tot, 1e-12))
            do_m = do_ * np.where(do_ > 0, sc_, 1.0); dc_m = dc_ * np.where(dc_ > 0, sc_, 1.0)
            mu_n = mu_h[m] + do_m * (po[m] - mux[m]) + dc_m * (pc[m] - mux[m])
            E2_n = E2[m] + do_m * (po[m] - s2x[m]) + dc_m * (pc[m] - s2x[m])
            dB = (E2_n - mu_n ** 2) - Bp0[m]
            Ux = U_h[m] / a_x[m]; U_n = U_h[m] - (do_m + dc_m) * Ux
            res[s] = dB.mean(); resU[s] = (U_n - U_h[m]).mean(); cap[s] = (sc_ < 1).mean(); fb[s] = m.sum()
        return res, resU, cap, fb
    zero = np.zeros(59)
    r2 = comp(d_o, zero, 'M2'); out['exp_M2'] = r2[0]; out['exp_M2_dU'] = r2[1]; out['M2_capfrac'] = r2[2]; out['M2_n_markers'] = r2[3]
    r2c = comp(zero, d_o, 'M2c'); out['exp_M2c'] = r2c[0]
    r3 = comp(d_o, d_c, 'M3'); out['exp_M3'] = r3[0]; out['exp_M3_dU'] = r3[1]; out['M3_capfrac'] = r3[2]
    # M2 와 대응하는 관측 ΔB 는 동일 marker 집합에서 계산(okx)이 아님 -> 참고: obs_dB(전체 elig)로 비교
    return out


def elig_fix(g_h, g_c):
    s1 = gsum(g_h, np.ones(nsamp[g_h])); s2 = gsum(g_c, np.ones(nsamp[g_c]))
    return dir_ok & (s1[6] >= 0.5 * nsamp[g_h]) & (s2[6] >= 0.5 * nsamp[g_c]), s1, s2


def summarize(o):
    obs = o['obs_dB']; res = {}
    for k in ['exp_M1', 'exp_M1p', 'exp_M1m', 'exp_M1c', 'exp_M2', 'exp_M2c', 'exp_M3']:
        e = o[k]; ok = np.isfinite(obs) & np.isfinite(e)
        res[f'AllAll_ratio_{k}'] = obs[0] / e[0] if e[0] != 0 else np.nan
        pos = ok & (e > 0)
        res[f'median_ratio_{k}'] = np.median(obs[pos] / e[pos]) if pos.sum() else np.nan
        res[f'slope0_{k}'] = (obs[ok] * e[ok]).sum() / (e[ok] ** 2).sum() if ok.sum() else np.nan
        res[f'pearson_{k}'] = np.corrcoef(obs[ok], e[ok])[0, 1] if ok.sum() > 2 else np.nan
        res[f'sum_ratio_{k}'] = obs[ok].sum() / e[ok].sum() if ok.sum() else np.nan
        res[f'median_resid_{k}'] = np.median(obs[ok] - e[ok]) if ok.sum() else np.nan
    for k in ['exp_M1_dU']:
        ok = np.isfinite(o['obs_dU']) & np.isfinite(o[k]); res['AllAll_ratio_dU_M1'] = o['obs_dU'][0] / o[k][0]; res['sum_ratio_dU_M1'] = o['obs_dU'][ok].sum() / o[k][ok].sum()
    ok = np.isfinite(o['obs_dFPI']) & np.isfinite(o['exp_M1_dFPI']); res['sum_ratio_dFPI_M1'] = o['obs_dFPI'][ok].sum() / o['exp_M1_dFPI'][ok].sum(); res['AllAll_ratio_dFPI_M1'] = o['obs_dFPI'][0] / o['exp_M1_dFPI'][0]
    return res


rng = np.random.default_rng(20261014)
points = {}; boots = {}
for split in ['Train', 'Test']:
    gh, gc = ('healthy', split), ('liver', split)
    elig, sh, sc = elig_fix(gh, gc)
    print(split, 'eligible markers', int(elig.sum()), 'of', NM, flush=True)
    o = model(sh, sc, nsamp[gh], nsamp[gc], elig)
    points[split] = o
    pd.DataFrame(dict(split=split, stratum=STRATA, **o)).to_csv(f'{OUT}/mix_strata_point{"_TEST" if TEST else ""}.tsv' if split == 'Train' else f'{OUT}/mix_strata_point_{split}{"_TEST" if TEST else ""}.tsv', sep='\t', index=False)
    bs = []
    for it in range(NBOOT):
        wh = np.bincount(rng.integers(0, nsamp[gh], nsamp[gh]), minlength=nsamp[gh]); wc = np.bincount(rng.integers(0, nsamp[gc], nsamp[gc]), minlength=nsamp[gc])
        ob = model(gsum(gh, wh), gsum(gc, wc), nsamp[gh], nsamp[gc], elig)
        bs.append(dict(summarize(ob), **{f'AA_{k}': ob[k][0] for k in ['delta_sOUF', 'delta_conc', 'obs_dB', 'exp_M1', 'exp_M1c', 'exp_M2', 'exp_M3', 'obs_dU', 'exp_M1_dU', 'obs_dFPI', 'exp_M1_dFPI', 'obs_dmixed']}))
        if it % 25 == 0: print(split, 'boot', it, round(time.time() - t0), flush=True)
    boots[split] = pd.DataFrame(bs)
rows = []
for split in ['Train', 'Test']:
    ps = summarize(points[split]); pa = {f'AA_{k}': points[split][k][0] for k in ['delta_sOUF', 'delta_conc', 'obs_dB', 'exp_M1', 'exp_M1c', 'exp_M2', 'exp_M3', 'obs_dU', 'exp_M1_dU', 'obs_dFPI', 'exp_M1_dFPI', 'obs_dmixed']}
    for k, v in {**ps, **pa}.items():
        bb = boots[split][k].to_numpy(); rows.append(dict(split=split, quantity=k, point=v, boot_lo=np.nanpercentile(bb, 2.5), boot_hi=np.nanpercentile(bb, 97.5), boot_n=len(bb)))
    boots[split].to_csv(f'{OUT}/mix_boot_{split}{"_TEST" if TEST else ""}.tsv', sep='\t', index=False)
pd.DataFrame(rows).to_csv(f'{OUT}/mix_summary_boot{"_TEST" if TEST else ""}.tsv', sep='\t', index=False)
json.dump(dict(strict_filled_samples=nfill, nsamp={f'{a}|{b_}': v for (a, b_), v in nsamp.items()}, nboot=NBOOT, seconds=time.time() - t0), open(f'{OUT}/mix_meta{"_TEST" if TEST else ""}.json', 'w'))
print('done', time.time() - t0, flush=True)
