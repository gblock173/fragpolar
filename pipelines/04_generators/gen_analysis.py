"""[rebuilt variant: Train strict states from train-rebuild step (all 777 Train samples), polarization also at every f]
generator step analysis: generator response (Train seed) vs observed Train / Test HCC-minus-healthy vectors; per-stratum cosine; strict polarization rule.
Fixed f = 0.05 is the primary result; 'best f' (the original rule, maximum over f per generator and stratum) is reported for comparison only."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import glob, json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
R = Path(CFG.DATA); GEN_DIR = Path(CFG.work('generators'))
DAT = GEN_DIR / 'data'; RES = GEN_DIR / 'results'
PF = ['amf', 'pdr', 'pm', 'me', 'mhl', 'uhl', 'fdrp', 'qfdrp', 'within_read_U', 'between_read_B_emp', 'hof', 'arrangement_excess_per_link']
GENS = ['G1', 'G2', 'G3']; FG = [0.01, 0.02, 0.05, 0.1, 0.2]; FIX = 0.05
SEEDS = [l.strip() for l in open(CFG.ids_file('generator_base')) if l.strip()]


def bh(p):
    p = np.asarray(p, float); o = np.argsort(p); q = np.empty_like(p); n = len(p); r = p[o] * n / (np.arange(n) + 1); r = np.minimum.accumulate(r[::-1])[::-1]; q[o] = np.minimum(r, 1); return q


def cos(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 4 or not np.any(a[ok]) or not np.any(b[ok]): return np.nan
    return float(a[ok] @ b[ok] / np.linalg.norm(a[ok]) / np.linalg.norm(b[ok]))


# ---------------------------------------------------------------- generator responses
def response(files, sids=None):
    L = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True); Z = L[L.scale.eq('z')]
    raw = Z[Z.condition.eq('raw')].set_index(['sid', 'stratum', 'feature']).value; rows = []
    for c, g in Z[Z.condition.ne('raw')].groupby('condition'):
        gi = g.set_index(['sid', 'stratum', 'feature']); d = gi.value - raw.reindex(gi.index).to_numpy()
        a = d.groupby(level=['stratum', 'feature']).agg(['mean', 'std', 'count']).reset_index(); a.insert(0, 'condition', c); rows.append(a)
    return pd.concat(rows, ignore_index=True)


tr_files = sorted(glob.glob(str(DAT / 'panel/*.parquet'))); assert len(tr_files) == 20, len(tr_files)
assert sorted(Path(f).stem for f in tr_files) == sorted(SEEDS)
RT = response(tr_files); RT.to_csv(RES / 'gen_response_TrainSeed20_rebuilt.tsv', sep='\t', index=False)
RO = CFG.unmap_features(pd.read_csv(CFG.inp('prior_perturbation_response'), sep='\t')); RO = RO[RO.condition.str.match(r'^G[123]_')]       # original Test seeds (read-only)
SN = list(dict.fromkeys(RT.stratum))
def vec(Rr, cond, s):
    return Rr[(Rr.condition == cond) & (Rr.stratum == s)].set_index('feature')['mean'].reindex(PF).to_numpy(float)

# ---------------------------------------------------------------- observed vectors
ms = pd.read_csv(DAT / 'obs_samples.tsv', sep='\t').set_index('AnalysisID'); PZ = pd.read_parquet(DAT / 'obs_panelz.parquet'); ST = pd.read_parquet(DAT / 'obs_strict_rebuilt.parquet')
def wide(X, cols): return {c: X.pivot(index='AnalysisID', columns='stratum', values=c).reindex(columns=SN) for c in cols}
WZ = wide(PZ, PF); WS = wide(ST, ['within_read_U', 'between_read_B_emp', 'Ushare', 'mixed_strict'])
liv = lambda sp: ms.index[(ms.group == 'liver') & (ms.split == sp)]; hea = lambda sp: ms.index[(ms.group == 'healthy') & (ms.split == sp)]
man = pd.read_csv(CFG.inp('prior_engine_manifest'), sep='\t'); hm = man[man.role.eq('healthy_matched')]
pairs_test = [(lv, h) for h, lv in zip(hm.sid, hm.matched_to) if lv in ms.index and h in ms.index]
# Train 1:1 depth/era matching, same algorithm as one-pass engine (nearest log-depth within era, without replacement, liver sorted by id); the 20 seed ids are not matching candidates
H = ms.loc[[h for h in hea('Train') if h not in set(SEEDS)]]; used = set(); pairs_train = []
for a in sorted(liv('Train')):
    c = H[~H.index.isin(used) & H.era.eq(ms.era[a])]
    if not len(c): c = H[~H.index.isin(used)]
    j = (c.log_depth - ms.log_depth[a]).abs().idxmin(); used.add(j); pairs_train.append((a, j))
CONTR = {'Test_paired': pairs_test, 'Train_matched': pairs_train, 'Test_unmatched': (list(liv('Test')), list(hea('Test'))), 'Train_unmatched': (list(liv('Train')), list(hea('Train')))}
pd.DataFrame([dict(contrast=k, liver=a, healthy=b) for k in ('Test_paired', 'Train_matched') for a, b in CONTR[k]]).to_csv(RES / 'gen_matched_pairs_rebuilt.tsv', sep='\t', index=False)


def effect(W, col, spec):
    """mean difference (HCC - healthy) and p per stratum for one metric table W[col]; paired -> Wilcoxon, unmatched -> Mann-Whitney."""
    T = W[col]; out = []
    for s in SN:
        if isinstance(spec, list):
            d = T.loc[[p[0] for p in spec], s].to_numpy(float) - T.loc[[p[1] for p in spec], s].to_numpy(float); d = d[np.isfinite(d)]
            out.append((float(d.mean()) if len(d) else np.nan, float(stats.wilcoxon(d).pvalue) if len(d) >= 10 and np.any(d) else np.nan, len(d)))
        else:
            a = T.loc[spec[0], s].to_numpy(float); b = T.loc[spec[1], s].to_numpy(float); a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
            out.append((float(a.mean() - b.mean()) if len(a) and len(b) else np.nan, float(stats.mannwhitneyu(a, b).pvalue) if len(a) >= 10 and len(b) >= 10 else np.nan, len(a)))
    return pd.DataFrame(out, columns=['delta', 'p', 'n'], index=SN)

OBS = {}; CELLS = []; EFF = {}
for k, spec in CONTR.items():
    v = {}; rows = []
    for f in PF:
        e = effect(WZ, f, spec); v[f] = e.delta; rows += [dict(contrast=k, stratum=s, feature=f, delta_z=r.delta, p=r.p, n=r.n) for s, r in e.iterrows()]
    C = pd.DataFrame(rows); ok = C.p.notna(); C['q'] = np.nan; C.loc[ok, 'q'] = bh(C.p[ok]); CELLS.append(C); OBS[k] = pd.DataFrame(v)[PF]
    EFF[k] = {m: effect(WS, m, spec) for m in ['within_read_U', 'between_read_B_emp', 'Ushare', 'mixed_strict']}
    for m in EFF[k]: e = EFF[k][m]; ok = e.p.notna(); e['q'] = np.nan; e.loc[ok, 'q'] = bh(e.p[ok])
pd.concat(CELLS).to_csv(RES / 'gen_observed_cells_deltaZ_rebuilt.tsv', sep='\t', index=False)
pd.concat([e.assign(contrast=k, metric=m).reset_index(names='stratum') for k, d in EFF.items() for m, e in d.items()]).to_csv(RES / 'gen_observed_strict_effects_rebuilt.tsv', sep='\t', index=False)

# parity of the observed side: Test paired (v1 matrices) vs original earlier analysis panel cells / group effects
pc = CFG.unmap_features(pd.read_csv(CFG.inp('prior_panel_cells'), sep='\t')); pcw = pc.pivot(index='stratum', columns='feature', values='delta_z').reindex(SN)[PF]
par = dict(panelz_cos_per_stratum_median=float(np.nanmedian([cos(OBS['Test_paired'].loc[s].to_numpy(), pcw.loc[s].to_numpy()) for s in SN])),
           panelz_pearson_all_cells=float(np.corrcoef(*[np.c_[OBS['Test_paired'].to_numpy().ravel(), pcw.to_numpy().ravel()][np.isfinite(OBS['Test_paired'].to_numpy().ravel()) & np.isfinite(pcw.to_numpy().ravel())].T])[0, 1]),
           panelz_max_abs_diff=float(np.nanmax(np.abs(OBS['Test_paired'].to_numpy() - pcw.to_numpy()))), panelz_n_pairs=len(pairs_test))
ge = pd.read_csv(CFG.inp('prior_group_effects'), sep='\t'); ge = ge[ge.contrast.eq('TestLiver_vs_matchedHealthy')]
for m in ['within_read_U', 'between_read_B_emp', 'Ushare', 'mixed_strict']:
    g = ge[ge.metric == m].set_index('stratum').reindex(SN); x = EFF['Test_paired'][m]
    par[f'strict_{m}_max_abs_diff_vs_prior'] = float(np.nanmax(np.abs(x.delta.to_numpy() - g.mean_diff.to_numpy())))
json.dump(par, open(RES / 'gen_observed_parity_vs_prior_rebuilt.json', 'w'), indent=1)
print('observed parity', par)

# ---------------------------------------------------------------- cosines
rows = []
for seedset, Rr in [('TrainSeed20', RT), ('TestSeed20_orig', RO)]:
    for k, O in OBS.items():
        for s in SN:
            o = O.loc[s].to_numpy(float); r = dict(seedset=seedset, contrast=k, stratum=s)
            for g in GENS:
                for f in FG: r[f'{g}_{f:g}'] = cos(o, vec(Rr, f'{g}_{f:g}', s))
                r[f'{g}_fixed'] = r[f'{g}_{FIX:g}']; cs = [r[f'{g}_{f:g}'] for f in FG]; r[f'{g}_best'] = np.nanmax(cs) if np.isfinite(cs).any() else np.nan
                r[f'{g}_bestf'] = FG[int(np.nanargmax(cs))] if np.isfinite(cs).any() else np.nan
            rows.append(r)
CS = pd.DataFrame(rows); CS.to_csv(RES / 'gen_cosine_by_stratum_rebuilt.tsv', sep='\t', index=False)
def q3(x): x = x.dropna(); return dict(median=x.median(), q25=x.quantile(.25), q75=x.quantile(.75), min=x.min(), max=x.max(), n=len(x), frac_pos=(x > 0).mean())
summ = []
for (sd, k), g in CS.groupby(['seedset', 'contrast']):
    for g_ in GENS:
        for how in ['fixed', 'best']: summ.append(dict(seedset=sd, contrast=k, gen=g_, f_rule=how, **q3(g[f'{g_}_{how}'])))
        for f in FG: summ.append(dict(seedset=sd, contrast=k, gen=g_, f_rule=f'f={f:g}', **q3(g[f'{g_}_{f:g}'])))
    summ.append(dict(seedset=sd, contrast=k, gen='G1>G2', f_rule='fixed', median=np.nan, n=len(g), frac_pos=(g.G1_fixed > g.G2_fixed).mean()))
    summ.append(dict(seedset=sd, contrast=k, gen='G1>G2', f_rule='best', median=np.nan, n=len(g), frac_pos=(g.G1_best > g.G2_best).mean()))
pd.DataFrame(summ).to_csv(RES / 'gen_cosine_summary_rebuilt.tsv', sep='\t', index=False)

# ---------------------------------------------------------------- strict polarization rule (registered C2_3 rule from s78c_mechanism.py)
sig = lambda e, s, sgn: bool(np.isfinite(e.delta[s]) and np.isfinite(e.q[s]) and e.q[s] < .05 and np.sign(e.delta[s]) == sgn)
pol = []
for (sd, k), g in CS.groupby(['seedset', 'contrast']):
    E = EFF[k]
    for how in ['fixed', 'best'] + [f'{f:g}' for f in FG]:
        for s in g.stratum:
            r = g[g.stratum == s].iloc[0]; c1 = np.nan_to_num(r[f'G1_{how}'], nan=-1); c2 = np.nan_to_num(r[f'G2_{how}'], nan=-1)
            sB = sig(E['between_read_B_emp'], s, 1); sUs = sig(E['Ushare'], s, -1); sUp = sig(E['within_read_U'], s, 1); sUsp = sig(E['Ushare'], s, 1); sMx = sig(E['mixed_strict'], s, 1)
            mx_na = not np.isfinite(E['mixed_strict'].delta[s])
            pol.append(dict(seedset=sd, contrast=k, f_rule=how, stratum=s, cosG1=c1, cosG2=c2, sigB_pos=sB, sigUshare_neg=sUs, notsig_mixed_pos=not sMx, cos1_ge_08=c1 >= .8, cos1_gt_cos2_02=c1 >= c2 + .2,
                            polarization=bool(sB and sUs and not sMx and c1 >= .8 and c1 >= c2 + .2),
                            within_fragment=bool(sUp and sUsp and sMx and c2 >= .8 and c2 >= c1 + .2)))
PL = pd.DataFrame(pol); PL.to_csv(RES / 'gen_polarization_by_stratum_rebuilt.tsv', sep='\t', index=False)
PS = PL.groupby(['seedset', 'contrast', 'f_rule']).agg(n_strata=('stratum', 'size'), polarization=('polarization', 'sum'), within_fragment=('within_fragment', 'sum'), sigB_pos=('sigB_pos', 'sum'),
                                                       sigUshare_neg=('sigUshare_neg', 'sum'), notsig_mixed_pos=('notsig_mixed_pos', 'sum'), cos1_ge_08=('cos1_ge_08', 'sum'), cos1_gt_cos2_02=('cos1_gt_cos2_02', 'sum')).reset_index()
PS.to_csv(RES / 'gen_polarization_summary_rebuilt.tsv', sep='\t', index=False); print(PS.to_string()); print(pd.DataFrame(summ).query("f_rule in ['fixed','best'] and gen!='G1>G2'").round(3).to_string())
