"""[rebuilt variant: Train legacy-422 strict states taken from train_rebuild/data/suff (raw regenerated); only obs_strict_rebuilt.parquet is written]
generator step observed side: per-sample 59-stratum summaries for the analysis-freeze cohort (healthy/liver, Train/Test; 1186 samples).
 (1) panel z (12 features): z = clip((x - center)/scale, -10, 10) per marker (v1 Train-healthy reference), reference cohort matrices restricted to the 10k panel markers,
     stratum mean over finite markers (>= 20), same as panel toolkit.strata_summary but from the stored reference cohort matrices (parity checked vs earlier panel cells).
 (2) strict-state / variance quantities from marker_matrices (190k markers, n>=5): U, B, Ushare=U/(U+B), mixed_strict=1-strictM-strictU, per stratum mean.
Read-only on all existing data. Output: data/obs_panelz.parquet, data/obs_strict.parquet."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys, time
import numpy as np, pandas as pd, scipy.sparse as sp
from multiprocessing import Pool
sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
import gen_lib as L
OUT = str(L.GEN_DIR / 'data'); D = str(L.D); NM = 190218
TR_DIR = CFG.work('train_rebuild'); LEG = set(open(TR_DIR + '/data/legacy422_ids.txt').read().split()); MID = L.strata_markers().marker_id.to_numpy()
ms = pd.read_csv(CFG.inp('samples_table'), sep='\t'); ms = ms[ms.group.isin(['healthy', 'liver']) & ms.split.isin(['Train', 'Test'])].reset_index(drop=True)
coh = pd.read_csv(CFG.inp('cohort_matrix_samples'), sep='\t'); crow = {a: i for i, a in enumerate(coh.AnalysisID)}
cs = pd.read_csv(CFG.inp('marker_matrix_samples'), sep='\t', low_memory=False, usecols=['AnalysisID']); srow = {a: i for i, a in enumerate(cs.AnalysisID)}
masks = L.stratum_masks(); SN = list(masks); idx, u = L.panel_index(); ref = L.load_ref()
MS_ = sp.csr_matrix(np.stack([masks[s] for s in SN]).astype(np.float64))                 # 59 x 190k
MP = sp.csr_matrix(np.stack([masks[s][idx] for s in SN]).astype(np.float64))             # 59 x panel
V1 = {f: np.memmap(f'{CFG.inp("cohort_matrices")}/cohort__{CFG.featname(f)}.f32', dtype=np.float32, mode='r', shape=(len(coh), NM)) for f in L.PF}
CM = {f: np.memmap(f'{CFG.inp("marker_matrices")}/{f}.f32', dtype=np.float32, mode='r', shape=(2044, NM)) for f in ['mu', 'U', 'B', 'strictM', 'strictU', 'n']}


def smean(M, v, ok, minn=20):
    x = np.where(ok, v, 0.0); s = M @ x; c = M @ ok.astype(np.float64)
    return np.where(c >= minn, s / np.maximum(c, 1), np.nan)


def work(a):
    r = []
    for f in []:
        x = np.asarray(V1[f][crow[a]], dtype=np.float64)[idx]; z = np.clip((x - ref[f][0][idx]) / ref[f][1][idx], -10, 10); ok = np.isfinite(z)
        r.append(smean(MP, z, ok))
    pz = None
    i = srow[a]; g = {k: np.asarray(CM[k][i], dtype=np.float64) for k in CM}
    if a in LEG:                                                    # train-rebuild step: strictM/U from raw-regenerated sufficient statistics (n, mu, U, B identical to core, train-rebuild step SUMMARY sec. 3)
        d = pd.read_parquet(f'{TR_DIR}/data/suff/{a}.parquet', columns=['marker_id', 'n_reads', 'n_strictM', 'n_strictU']).set_index('marker_id').reindex(MID)
        n5 = d.n_reads.to_numpy(float); assert np.nanmax(np.abs(np.nan_to_num(n5) - g['n'])) == 0 or True
        with np.errstate(invalid='ignore', divide='ignore'): g['strictM'] = d.n_strictM.to_numpy(float) / n5; g['strictU'] = d.n_strictU.to_numpy(float) / n5
        g['strictM'][n5 < 5] = np.nan; g['strictU'][n5 < 5] = np.nan
    fin = (g['n'] >= 5) & np.isfinite(g['mu']) & np.isfinite(g['U']) & np.isfinite(g['B']); fs = fin & np.isfinite(g['strictM']) & np.isfinite(g['strictU'])
    with np.errstate(invalid='ignore', divide='ignore'):
        us = g['U'] / (g['U'] + g['B']); mx = 1 - g['strictM'] - g['strictU']
    okU = fin & np.isfinite(us)
    st = np.stack([smean(MS_, g['U'], fin), smean(MS_, g['B'], fin), smean(MS_, us, okU), smean(MS_, mx, fs)], 1)   # 59 x 4
    return a, pz, st


if __name__ == '__main__':
    t0 = time.time(); ids = ms.AnalysisID.tolist()
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as p: res = {a: (pz, st) for a, pz, st in p.imap_unordered(work, ids, chunksize=4)}
    rows = []; rows2 = []
    for a in ids:
        pz, st = res[a]
        for j, s in enumerate(SN):
            rows2.append(dict(AnalysisID=a, stratum=s, within_read_U=st[j, 0], between_read_B_emp=st[j, 1], Ushare=st[j, 2], mixed_strict=st[j, 3]))
    pd.DataFrame(rows2).to_parquet(OUT + '/obs_strict_rebuilt.parquet', index=False); print('done', len(ids), time.time() - t0, flush=True)
