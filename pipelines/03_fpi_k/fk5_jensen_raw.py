"""FPI_K step: Train healthy 3표본의 raw cpg_calls(a6_2 산출)로 marker 수준 정확한 E[1/K] vs 1/mean_K 비교.
(1) raw 재계산 mu/B/mean_K가 marker_matrices와 같은지 확인, (2) marker 수준 h 오차, (3) 59층 FPI_K 값이 h=1/mean_K vs 정확한 E[1/K]에서 얼마나 달라지는지,
(4) c = ΣV·E[1/K]/ΣV·(1/mean_K) (All|All, V 가중) -> A6_raw_jensen_summary.json (분석의 c_global 민감도에 사용; Train healthy에서만 정한다)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, json
import numpy as np, pandas as pd
A = CFG.WORK
sys.path.insert(0, CFG.pipeline_dir('group_effects')); import ge1_aggregate as a1
RES = A + '/fpi_k/results'; sids = sys.argv[1:]
a1.setup(); masks = a1.G['masks']; strata = a1.G['strata']; sel = a1.G['sel'].set_index('AnalysisID')
rows = []; srows = []; scat = []; cs = []
for sid in sids:
    z = np.load(f'{RES}/raw_invK_{sid}.npz'); acc = z['acc']; n = acc[:, 0]
    with np.errstate(all='ignore'):
        K = acc[:, 1] / n; hex_ = acc[:, 2] / n; mu = acc[:, 3] / n; B = acc[:, 4] / n - mu ** 2; V = mu * (1 - mu); h1 = 1 / K
    i = int(sel.loc[sid, 'col']); cm = {f: np.asarray(a1.MM[f][i], float) for f in ['n', 'mu', 'B', 'mean_K']}
    v = (n >= 5) & (cm['n'] >= 5) & np.isfinite(cm['mu'])
    chk = dict(sid=sid, n_valid=int(v.sum()), max_abs_n=np.nanmax(np.abs(n[v] - cm['n'][v])), max_abs_mu=np.nanmax(np.abs(mu[v] - cm['mu'][v])), max_abs_B=np.nanmax(np.abs(B[v] - cm['B'][v])), max_abs_K=np.nanmax(np.abs(K[v] - cm['mean_K'][v])),
               frac_markers_n_equal=float((n[v] == cm['n'][v]).mean()))
    v = (n >= 5) & np.isfinite(B) & np.isfinite(h1)
    ratio = hex_[v] / h1[v]; err = hex_[v] - h1[v]
    chk.update(marker_c_median=np.median(ratio), marker_c_mean=ratio.mean(), marker_c_Vweighted=(V[v] * hex_[v]).sum() / (V[v] * h1[v]).sum(), marker_c_IQR_lo=np.quantile(ratio, .25), marker_c_IQR_hi=np.quantile(ratio, .75),
               mean_h_1overK=h1[v].mean(), mean_h_exact=hex_[v].mean(), mean_abs_err_h=np.abs(err).mean(), pooled_c_fragment=(acc[:, 2].sum() / acc[:, 0].sum()) * (acc[:, 1].sum() / acc[:, 0].sum()), n_markers=int(v.sum()))
    rows.append(chk); cs.append(chk['marker_c_Vweighted'])
    f = lambda m_, h_: (np.where(m_, B, 0).sum() - 0) # placeholder
    for j, s in enumerate(strata):
        mk = masks[j] & v
        if mk.sum() < 20: continue
        sb, sv, sh1, she = B[mk].sum(), V[mk].sum(), (h1 * V)[mk].sum(), (hex_ * V)[mk].sum()
        fpi = sb / sv; fk1 = (sb - sh1) / (sv - sh1); fke = (sb - she) / (sv - she)
        srows.append(dict(sid=sid, stratum=s, n_markers=int(mk.sum()), FPI=fpi, FPI_K_1overK=fk1, FPI_K_exact=fke, diff_exact_minus_1overK=fke - fk1, h_w_1overK=sh1 / sv, h_w_exact=she / sv, c_V=she / sh1))
    # 산점도 데이터 (n>=20 marker 무작위 3만개)
    vv = np.flatnonzero(v & (n >= 20)); rng = np.random.default_rng(1); pick = rng.choice(vv, min(30000, len(vv)), replace=False)
    scat.append(pd.DataFrame(dict(sid=sid, marker=pick, n=n[pick], h_1overK=h1[pick], h_exact=hex_[pick], FPI_m=B[pick] / V[pick], mu=mu[pick])))
C = pd.DataFrame(rows); C.to_csv(f'{RES}/fpik_raw_jensen_by_sample.tsv', sep='\t', index=False); print(C.T.to_string())
S = pd.DataFrame(srows); S.to_csv(f'{RES}/fpik_raw_jensen_by_stratum.tsv', sep='\t', index=False)
g = S.groupby('stratum').agg(diff_mean=('diff_exact_minus_1overK', 'mean'), diff_sd_between_samples=('diff_exact_minus_1overK', 'std'), FPI=('FPI', 'mean'), FPIK=('FPI_K_1overK', 'mean'), FPIKexact=('FPI_K_exact', 'mean')).reset_index()
g.to_csv(f'{RES}/fpik_raw_jensen_by_stratum_mean.tsv', sep='\t', index=False)
aa = S[S.stratum == 'All|All']
summ = dict(sids=sids, c_marker_Vweighted_AllAll_mean=float(np.mean(cs)), c_marker_Vweighted_AllAll_each=cs, n_strata=int(g.shape[0]), AllAll_FPI=float(aa.FPI.mean()), AllAll_FPIK_1overK=float(aa.FPI_K_1overK.mean()), AllAll_FPIK_exact=float(aa.FPI_K_exact.mean()),
            AllAll_diff=float(aa.diff_exact_minus_1overK.mean()), strata_diff_median=float(g['diff_mean'].median()), strata_diff_min=float(g['diff_mean'].min()), strata_diff_max=float(g['diff_mean'].max()),
            strata_diff_sd_between_strata=float(g['diff_mean'].std()), median_between_sample_sd=float(g.diff_sd_between_samples.median()), spearman_FPIK_vs_FPIKexact_strata=float(g[['FPIK', 'FPIKexact']].corr('spearman').iloc[0, 1]))
json.dump(summ, open(f'{RES}/fpik_raw_jensen_summary.json', 'w'), indent=1); print(json.dumps(summ, indent=1))
pd.concat(scat).to_parquet(f'{RES}/fpik_marker_scatter_data.parquet', index=False)
