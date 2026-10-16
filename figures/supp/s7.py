from sc import *
R = AN / 'mixture_and_length/results/'
pt = pd.concat([pd.read_csv(R / 'mix_strata_point.tsv', sep='\t'), pd.read_csv(R / 'mix_strata_point_Test.tsv', sep='\t')])
pt = pt.rename(columns={'obs_dB': 'obs_dBFV', 'obs_dU': 'obs_dWFV', 'exp_M1_dU': 'exp_M1_dWFV', 'exp_M3_dU': 'exp_M3_dWFV'})
wt(pt[['split', 'stratum', 'delta_sOUF', 'obs_dBFV', 'exp_M1', 'exp_M1p', 'exp_M1m', 'exp_M2', 'exp_M3', 'exp_M1c', 'exp_M2c', 'obs_dWFV', 'exp_M1_dWFV', 'exp_M3_dWFV', 'obs_dFPI', 'exp_M1_dFPI']], 'S7a.tsv')
sb = pd.concat([pd.read_csv(R / 'mix_summary_boot.tsv', sep='\t')]); sb = sb[sb.quantity.str.startswith(('sum_ratio_exp_', 'pearson_exp_'))]
wt(sb[['split', 'quantity', 'point', 'boot_lo', 'boot_hi']], 'S7b.tsv')
A, B = rt('S7a.tsv'), rt('S7b.tsv'); rng = np.random.default_rng(1)
MOD = [('exp_M1', 'M1 exact'), ('exp_M1p', 'M1p pooled'), ('exp_M1m', 'M1m signed δ per marker'), ('exp_M2', 'M2 replace mixed'), ('exp_M3', 'M3 full composition')]
cv = SC(116, 7)
LIM = max(12, float(max(A[c].max() for c, _ in MOD) * 1e3), float(A.obs_dBFV.max() * 1e3)) + .6
cv.panel('a', 1, 1); cv.title(8, 2.2, 'Observed ΔBFV against expectation from alternative mixing models')
for k, (c, nm) in enumerate(MOD):
    ax = cv.ax(12 + k * 34, 13, 26, 26)
    for sp, mk in (('Train', 'o'), ('Test', 's')):
        g = A[A.split == sp]; ax.scatter(g[c] * 1e3, g.obs_dBFV * 1e3, s=5, marker=mk, color=SPLIT_C[sp], alpha=.7, lw=0)
    ax.plot([-1, LIM], [-1, LIM], color=PALE, ls='--', lw=1.0, zorder=0); ax.set_xlim(-1, LIM); ax.set_ylim(-1, LIM)
    ax.set_title(nm, fontsize=5.5, pad=2, loc='left'); ax.set_xlabel('Expected ΔBFV (×10$^{-3}$)', labelpad=1, fontsize=5.5)
    if k == 0: ax.set_ylabel('Observed ΔBFV (×10$^{-3}$)', labelpad=1)
    if k == 0: ax.legend(frameon=False, fontsize=5.5, loc='upper left', handletextpad=.1, borderaxespad=0, labels=['Train', 'Test'])
cv.panel('b', 1, 47); cv.title(8, 48.2, 'Observed / expected ΔBFV (sum over 59 strata)')
ax = cv.ax(26, 58, 36, 38); lab = ['M1', 'M1p', 'M1m', 'M2', 'M3']; names = ['exp_M1', 'exp_M1p', 'exp_M1m', 'exp_M2', 'exp_M3']
for sp, off in (('Train', -.17), ('Test', .17)):
    for i, n in enumerate(names):
        r = B[(B.split == sp) & (B.quantity == 'sum_ratio_' + n)].iloc[0]; ax.plot([r.boot_lo, r.boot_hi], [i + off] * 2, color=SPLIT_C[sp], lw=1.0); ax.scatter([r.point], [i + off], s=9, color=SPLIT_C[sp], lw=0, zorder=3)
ax.axvline(1, color=PALE, ls='--', lw=1.0, zorder=0); ax.set_yticks(range(5)); ax.set_yticklabels(lab); ax.invert_yaxis(); ax.set_xlabel('Ratio (bootstrap 95% CI)', labelpad=1); ax.set_xlim(.6, 1.5)
cv.panel('c', 72, 47); cv.title(79, 48.2, 'Same δ as same-direction uniform')
ax = cv.ax(88, 58, 36, 38)
for sp, mk in (('Train', 'o'), ('Test', 's')):
    g = A[A.split == sp]; ax.scatter(g.exp_M1 * 1e3, g.exp_M1c * 1e3, s=5, marker=mk, color=SPLIT_C[sp], alpha=.7, lw=0)
ax.axhline(0, color=PALE, lw=1.0); ax.set_xlabel('Expected ΔBFV, opposite uniform (×10$^{-3}$)', labelpad=1, fontsize=5.5); ax.set_ylabel('Same-direction (×10$^{-3}$)', labelpad=1, fontsize=5.5)
cv.panel('d', 134, 47); cv.title(141, 48.2, 'ΔWFV and ΔFPI vs sOUF-only')
ax = cv.ax(150, 58, 26, 15); ax2 = cv.ax(150, 82, 26, 15)
for sp, mk in (('Train', 'o'), ('Test', 's')):
    g = A[A.split == sp]; ax.scatter(g.exp_M1_dWFV * 1e3, g.obs_dWFV * 1e3, s=4, marker=mk, color=SPLIT_C[sp], alpha=.7, lw=0); ax2.scatter(g.exp_M1_dFPI, g.obs_dFPI, s=4, marker=mk, color=SPLIT_C[sp], alpha=.7, lw=0)
for a_, lo, hi in ((ax, -20, 5), (ax2, -.01, .09)): a_.plot([lo, hi], [lo, hi], color=PALE, ls='--', lw=1.0, zorder=0)
ax.set_xlabel('Exp. ΔWFV (×10$^{-3}$)', labelpad=1, fontsize=5.5); ax.set_ylabel('Obs. ΔWFV (×10$^{-3}$)', labelpad=1, fontsize=5.5)
ax2.set_xlabel('Exp. ΔFPI', labelpad=1, fontsize=5.5); ax2.set_ylabel('Obs. ΔFPI', labelpad=1, fontsize=5.5)
cv.text(134, 106, 'Expected from sOUF alone (M1).', 5.5, color=MUTED)
cv.save()
