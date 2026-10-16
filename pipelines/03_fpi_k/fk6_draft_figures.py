"""FPI_K step: 그림 2개. (1) healthy Train에서 FPI vs h (표본 All|All, marker 수준), (2) 층별 HCC-healthy ΔFPI vs ΔFPI_K."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
A = CFG.work('fpi_k')
for f in ['NanumGothic', 'Noto Sans CJK KR', 'Malgun Gothic', 'AppleGothic']:
    if any(f in x.name for x in fm.fontManager.ttflist): plt.rcParams['font.family'] = f; break
plt.rcParams['axes.unicode_minus'] = False
v = pd.read_csv(A + '/results/fpik_sample_AllAll_values.tsv', sep='\t'); th = v[(v.split == 'Train') & (v.group == 'healthy')]
sc = pd.read_parquet(A + '/results/fpik_marker_scatter_data.parquet')
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
a = ax[0]; a.scatter(th.h_floor, th.FPI, s=9, alpha=.5, color='#3b6ea5', label=f'healthy Train 표본 (n={len(th)})')
xs = np.linspace(th.h_floor.min() * .97, th.h_floor.max() * 1.03, 10); a.plot(xs, xs, 'k--', lw=1, label='FPI = h (CpG 독립, r=0)')
a.set_xlabel('h = V 가중 평균 1/mean_K (All|All)'); a.set_ylabel('FPI = BFV/μ(1-μ) (All|All)'); a.set_title('표본 수준: FPI는 h 바닥 위에 있다'); a.legend(fontsize=8)
a = ax[1]; s = sc.sample(min(len(sc), 40000), random_state=0); hb = a.hexbin(s.h_1overK, s.FPI_m.clip(0, 1.2), gridsize=45, bins='log', cmap='viridis', mincnt=1)
bins = np.quantile(sc.h_1overK, np.linspace(0, 1, 21)); idx = np.clip(np.digitize(sc.h_1overK, bins) - 1, 0, 19)
g = sc.assign(b=idx).groupby('b').agg(h=('h_1overK', 'mean'), f=('FPI_m', 'mean')); a.plot(g.h, g.f, 'r-o', ms=3, lw=1.3, label='h 분위 구간 평균')
xs = np.linspace(.1, .25, 5); a.plot(xs, xs, 'w--', lw=1.2, label='FPI = h'); a.set_xlim(.09, .26); a.set_ylim(0, 1.0)
a.set_xlabel('h = 1/mean_K (marker)'); a.set_ylabel('marker FPI = B/μ(1-μ)'); a.set_title('marker 수준 (Train healthy 3표본, n≥20)'); a.legend(fontsize=8, loc='upper left'); fig.colorbar(hb, ax=a, label='log10 marker 수')
fig.tight_layout(); fig.savefig(A + '/figures/Fig_fpik_1_FPI_vs_h_healthyTrain.png', dpi=200); plt.close(fig)
p = pd.read_csv(A + '/results/fpik_FPI_vs_FPIK_by_stratum.tsv', sep='\t')
fig, ax = plt.subplots(1, 2, figsize=(10, 4.6), sharex=True, sharey=True)
for a, sp, c in zip(ax, ['Train', 'Test'], ['#3b6ea5', '#c0504d']):
    q = p[p.split == sp]; a.scatter(q.delta_FPI, q.delta_FPIK, s=16, color=c, alpha=.75)
    lim = [min(0, q.delta_FPI.min(), q.delta_FPIK.min()) - .003, max(q.delta_FPI.max(), q.delta_FPIK.max()) + .004]; a.plot(lim, lim, 'k--', lw=1, label='y = x')
    sl = np.sum(q.delta_FPI * q.delta_FPIK) / np.sum(q.delta_FPI ** 2); a.plot(lim, [sl * l for l in lim], '-', color='gray', lw=1, label=f'원점 회귀 기울기 {sl:.2f}')
    a.axhline(0, color='k', lw=.5); a.axvline(0, color='k', lw=.5); a.set_title(f'{sp}: ΔFPI>0 {int((q.delta_FPI>0).sum())}층, ΔFPI_K>0 {int((q.delta_FPIK>0).sum())}층 (/{len(q)})', fontsize=10)
    a.set_xlabel('ΔFPI (HCC - healthy)'); a.legend(fontsize=8, loc='upper left')
ax[0].set_ylabel('ΔFPI_K (HCC - healthy)'); fig.tight_layout(); fig.savefig(A + '/figures/Fig_fpik_2_dFPI_vs_dFPIK_by_stratum.png', dpi=200)
print('figs done')
