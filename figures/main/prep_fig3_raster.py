"""Fig 3a source data. Pre-specified rule (pipelines/05_direct_fragments/df6_raster.py), pairs from gen_matched_pairs_rebuilt.tsv:
Train pair of the first HCC by sorted AnalysisID; layer All|All; first two markers (index order) in the HCC sample with n>=40,
>=3 opposite, >=1 concordant, >=1 mixed fragment; <=45 K>=4 fragments per sample (seed 20261014); rows ordered opposite/mixed/concordant."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys; sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from common import *
import pyarrow.parquet as pq
R = Path(CFG.DATA); RAW = Path(CFG.CPG_CALLS)
P = pd.read_csv(AN / 'generators/results/gen_matched_pairs_rebuilt.tsv', sep='\t'); P = P[P.contrast == 'Train_matched'].sort_values('liver')
pr = P.iloc[0]; hcc, hl = pr.liver, pr.healthy; print('pair', hcc, hl)
st = pd.read_parquet(CFG.inp('marker_annotation'), columns=['marker_id', CFG.col('healthy_direction')])
z = np.load(AN / 'direct_fragments/work/ml' / (hcc + '.npz'))
cand = np.flatnonzero((z['n'] >= 40) & (z['r'] >= 3) & (z['c_con'] >= 1) & (z['c_mix'] >= 1))[:2]; print('markers', cand, st.marker_id.iloc[cand].tolist())
rows, summ = [], []
for mi in cand:
    mid = st.marker_id.iloc[mi]; dirn = int(st[CFG.col('healthy_direction')].iloc[mi])
    for sid, lab in ((hcc, 'HCC'), (hl, 'Healthy')):
        d = pq.read_table(RAW / (sid + '.cpg_calls.parquet'), columns=['marker_id', 'fragment_id', 'pos', 'call'], filters=[('marker_id', '=', mid)]).to_pandas()
        d = d[d.call.isin([0, 1])].drop_duplicates(['fragment_id', 'pos'])
        g = d.groupby('fragment_id'); k = g.call.size(); m = g.call.sum(); keep = k[k >= 4].index; d = d[d.fragment_id.isin(keep)]
        kk, mm = k[keep], m[keep]; opp = (mm == kk) if dirn == 1 else (mm == 0); con = (mm == 0) if dirn == 1 else (mm == kk)
        cls = np.where(opp, 0, np.where(con, 2, 1)); p = (mm / kk); nfr = len(kk)
        summ.append(dict(sample=sid, group=lab, marker=mid, direction=dirn, n_fragments=nfr, n_opposite=int(opp.sum()), n_mixed=int((cls == 1).sum()), n_concordant=int(con.sum()), BFV=float(p.var(ddof=0))))
        rng = np.random.default_rng(20261014); sel = rng.choice(nfr, min(45, nfr), replace=False)
        ids = np.array(keep)[sel]; cl = cls[sel]; pp = p.to_numpy()[sel]; o = np.lexsort((pp, cl)); ids, cl = ids[o], cl[o]
        pos = np.sort(d.pos.unique()); pidx = {q: j for j, q in enumerate(pos)}
        sub = d[d.fragment_id.isin(ids)].copy(); order = {f: i for i, f in enumerate(ids)}
        sub['row'] = sub.fragment_id.map(order); sub['col'] = sub.pos.map(pidx); sub['class'] = pd.Series(cl, index=range(len(cl))).reindex(sub.row).to_numpy()
        sub['class'] = sub['class'].map({0: 'opposite', 1: 'mixed', 2: 'concordant'}); sub['sample'] = sid; sub['group'] = lab; sub['marker'] = mid
        rows.append(sub[['sample', 'group', 'marker', 'row', 'class', 'col', 'pos', 'call']])
out = pd.concat(rows); out.to_csv(SD / 'Fig3a.tsv', sep='\t', index=False)
pd.DataFrame(summ).to_csv(SD / 'Fig3a_summary.tsv', sep='\t', index=False); print(pd.DataFrame(summ))
