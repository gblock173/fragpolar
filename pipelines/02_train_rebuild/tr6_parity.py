"""parity: 재생성 sufficient parquet vs 기존(inputs.prior_sufficient_dir) parquet, 그리고 marker_matrices 행(float32). ids = argv[1:] (raw 경로 표본) 또는 --all-new 로 legacy 비교."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, json, numpy as np, pandas as pd
sys.argv_ = sys.argv
TR_DIR = CFG.work('train_rebuild')
OLD = CFG.inp('prior_sufficient_dir')
sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
import importlib.util
spec = importlib.util.spec_from_file_location('a5agg', str(_pl.Path(__file__).resolve().parent / 'tr2_aggregate.py')); ag = importlib.util.module_from_spec(spec); spec.loader.exec_module(ag)
ag.setup()
col = dict(zip(ag.G['sel'].AnalysisID, ag.G['sel'].col))
rows = []
for sid in sys.argv[1:]:
    new = pd.read_parquet(f'{TR_DIR}/data/suff/{sid}.parquet').sort_values('marker_id').reset_index(drop=True)
    old = pd.read_parquet(f'{OLD}/{sid}.parquet').sort_values('marker_id').reset_index(drop=True)
    r = dict(sample=sid, rows_new=len(new), rows_old=len(old), same_markers=bool((new.marker_id.values == old.marker_id.values).all()) if len(new) == len(old) else False)
    if r['same_markers']:
        mx = {}
        for c in new.columns:
            if c == 'marker_id': continue
            a = new[c].to_numpy(float); b = old[c].to_numpy(float)
            mx[c] = 0.0 if np.array_equal(a, b, equal_nan=True) else float(np.nanmax(np.abs(a - b)))
        r['max_abs_diff_all_cols'] = max(mx.values()); r['exact_parquet'] = r['max_abs_diff_all_cols'] == 0
        r['cols_differing'] = ','.join(k for k, v in mx.items() if v > 0)
    ag.LEG.add(sid); rb = ag.rebuilt_rows(sid); i = col[sid]; md = {}
    for f in ag.MATS:
        a = np.asarray(rb[f], np.float64); b = np.asarray(ag.MM[f][i], np.float64)
        md[f] = 0.0 if np.array_equal(a, b, equal_nan=True) else float(np.nanmax(np.abs(a - b)))
    r['matrix_max_abs_diff'] = max(md.values()); r['matrix_exact'] = r['matrix_max_abs_diff'] == 0
    rows.append(r); print(r, flush=True)
pd.DataFrame(rows).to_csv(TR_DIR + '/results/parity_raw_samples.tsv', sep='\t', index=False)
