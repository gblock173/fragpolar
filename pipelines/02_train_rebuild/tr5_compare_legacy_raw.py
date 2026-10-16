"""legacy(marker_matrices, legacy 경로) vs raw 재생성 값: 422개 Train 표본의 marker 수준/표본 수준 차이 (B, mu, U, mean_K, n)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, numpy as np, pandas as pd
from scipy import stats
TR_DIR = CFG.work('train_rebuild')
import importlib.util
spec = importlib.util.spec_from_file_location('a5agg', str(_pl.Path(__file__).resolve().parent / 'tr2_aggregate.py')); ag = importlib.util.module_from_spec(spec); spec.loader.exec_module(ag)
ag.setup(); sel = ag.G['sel']; ids = open(TR_DIR + '/data/legacy422_ids.txt').read().split()
sel = sel[sel.AnalysisID.isin(ids)]; rows = []
F = ['mu', 'U', 'B', 'mean_K', 'n']
for r in sel.itertuples():
    rb = ag.rebuilt_rows(r.AnalysisID); lg = {f: np.asarray(ag.MM[f][r.col], np.float64) for f in F}; rw = {f: np.asarray(rb[f], np.float64) for f in F}
    d = dict(AnalysisID=r.AnalysisID, group=r.group)
    d['markers_n5_legacy'] = int(((lg['n'] >= 5) & np.isfinite(lg['B'])).sum()); d['markers_n5_raw'] = int(((rw['n'] >= 5) & np.isfinite(rw['B'])).sum())
    ok = (lg['n'] >= 5) & (rw['n'] >= 5) & np.isfinite(lg['B']) & np.isfinite(rw['B']); d['markers_n5_both'] = int(ok.sum())
    for f in F:
        x = lg[f][ok]; y = rw[f][ok]
        d[f'{f}_mean_legacy'] = x.mean(); d[f'{f}_mean_raw'] = y.mean(); d[f'{f}_mean_diff_raw_minus_legacy'] = (y - x).mean()
        d[f'{f}_median_abs_marker_diff'] = np.median(np.abs(y - x)); d[f'{f}_corr_marker'] = np.corrcoef(x, y)[0, 1]
    d['n_ratio_raw_over_legacy_median'] = np.median(rw['n'][ok] / lg['n'][ok]); d['frac_markers_n_equal'] = float((rw['n'][ok] == lg['n'][ok]).mean())
    rows.append(d); 
T = pd.DataFrame(rows); T.to_csv(TR_DIR + '/results/legacy_vs_raw_per_sample.tsv', sep='\t', index=False)
out = []
for f in F:
    dd = T[f'{f}_mean_diff_raw_minus_legacy']; out.append(dict(metric=f, n_samples=len(T), mean_legacy=T[f'{f}_mean_legacy'].mean(), mean_raw=T[f'{f}_mean_raw'].mean(), mean_diff_raw_minus_legacy=dd.mean(),
        rel_diff_pct=100 * dd.mean() / T[f'{f}_mean_legacy'].mean(), sd_sample_diff=dd.std(), sample_corr_legacy_raw=np.corrcoef(T[f'{f}_mean_legacy'], T[f'{f}_mean_raw'])[0, 1],
        median_marker_abs_diff=T[f'{f}_median_abs_marker_diff'].median(), median_marker_corr=T[f'{f}_corr_marker'].median()))
pd.DataFrame(out).to_csv(TR_DIR + '/results/legacy_vs_raw_summary.tsv', sep='\t', index=False)
# 군 비교(HCC vs healthy) in 422 : legacy vs raw
g = []
for f in ['B', 'mu', 'U', 'mean_K']:
    for nm in ['legacy', 'raw']:
        a = T[T.group == 'liver'][f'{f}_mean_{nm}']; b = T[T.group == 'healthy'][f'{f}_mean_{nm}']
        g.append(dict(metric=f, path=nm, n_HCC=len(a), n_healthy=len(b), mean_HCC=a.mean(), mean_healthy=b.mean(), delta=a.mean() - b.mean(), p=stats.mannwhitneyu(a, b).pvalue))
pd.DataFrame(g).to_csv(TR_DIR + '/results/legacy_vs_raw_group_delta_422.tsv', sep='\t', index=False)
print(pd.DataFrame(out).round(5).to_string()); print(pd.DataFrame(g).round(5).to_string())
