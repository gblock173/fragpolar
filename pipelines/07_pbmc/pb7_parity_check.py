"""재추출 parity: 새 cpg_calls vs 기존 cpg_calls(cpg_calls_dir) (190,218 primary marker 한정) + 새 suff vs 기존 core suff."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, json, numpy as np, pandas as pd, pyarrow.parquet as pq, pyarrow.compute as pc, pyarrow as pa
sid = sys.argv[1]
D = CFG.DATA
PB_DIR = CFG.work('pbmc')
old = CFG.cpg_calls(sid)
new = f'{PB_DIR}/data/parity/cpg_calls/{sid}.cpg_calls.parquet'
prim = pd.read_parquet(CFG.inp('marker_regions'))
prim = prim[prim.primary_eligible]; allowed = pa.array(prim.marker_id)
po, pn = pq.ParquetFile(old), pq.ParquetFile(new)
print('old schema', po.schema_arrow.names == pn.schema_arrow.names, po.metadata.num_rows, pn.metadata.num_rows)
cols = ['marker_id', 'fragment_id', 'pos', 'call', 'fragment_start', 'fragment_end', 'conversion_context', 'is_paired', 'is_complete_pair', 'n_alignments']
res = dict(sample=sid, chroms={}); tot = dict(old=0, new=0, only_old=0, only_new=0, call_diff=0, other_diff=0)
for ch in ['chr%d' % i for i in range(1, 23)]:
    to = pq.read_table(old, columns=cols, filters=[('chrom', '=', ch)]); to = to.filter(pc.is_in(to['marker_id'], value_set=allowed)).to_pandas()
    tn = pq.read_table(new, columns=cols, filters=[('chrom', '=', ch)]).to_pandas()
    k = ['marker_id', 'fragment_id', 'pos']
    m = to.merge(tn, on=k, how='outer', suffixes=('_o', '_n'), indicator=True, validate='one_to_one')
    oo = int((m._merge == 'left_only').sum()); on = int((m._merge == 'right_only').sum()); b = m[m._merge == 'both']
    cd = int((b.call_o != b.call_n).sum())
    od = int(sum((b[c + '_o'] != b[c + '_n']).sum() for c in cols[4:] if c != 'call'))
    res['chroms'][ch] = dict(old=len(to), new=len(tn), only_old=oo, only_new=on, call_diff=cd, other_col_diff=od)
    for key, v in [('old', len(to)), ('new', len(tn)), ('only_old', oo), ('only_new', on), ('call_diff', cd), ('other_diff', od)]: tot[key] += v
    print(ch, res['chroms'][ch], flush=True)
res['total'] = tot
# suff 비교
cs = pd.read_csv(CFG.inp('marker_matrix_samples'), sep='\t', low_memory=False, usecols=['AnalysisID', 'suff']).set_index('AnalysisID')
so = pd.read_parquet(cs.loc[sid, 'suff']).set_index('marker_id'); sn = pd.read_parquet(f'{PB_DIR}/data/parity/suff/{sid}.parquet').set_index('marker_id')
res['suff'] = dict(rows_old=len(so), rows_new=len(sn), index_equal=bool(so.index.sort_values().equals(sn.index.sort_values())))
sn = sn.reindex(so.index)
res['suff']['max_abs_diff'] = {c: float(np.nanmax(np.abs(so[c].to_numpy(float) - sn[c].to_numpy(float)))) for c in so.columns if c in sn.columns}
res['suff']['n_markers_with_any_diff'] = int((np.abs(so[[c for c in so.columns if c in sn.columns]].to_numpy(float) - sn[[c for c in so.columns if c in sn.columns]].to_numpy(float)) > 1e-9).any(1).sum())
json.dump(res, open(f'{PB_DIR}/results/parity_{sid}.json', 'w'), indent=1)
print(json.dumps(res['total']), json.dumps(res['suff']))
