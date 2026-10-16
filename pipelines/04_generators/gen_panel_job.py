"""generator step per-sample generator job: copy of the generator panel step per-sample cfDNA base-pool branch (raw + G1/G2/G3 at f in {1,2,5,10,20}%, one draw per (g, f),
rng = default_rng(SEED + 7800 + b) shared across g and f in the same order as the original). Tissue spike-ins (N/T) are not run.
--engdir selects where panel_calls.parquet lives; --b is the base index (0-19) used in the rng seed; --out the result parquet."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import argparse, time, json
from pathlib import Path
import gen_lib as L
from gen_lib import *

FG = [.01, .02, .05, .10, .20]


def run(mats, idx, dirA, masks, ref, cond):
    F = features(mats, idx, dirA); S_ = strata_summary(F, idx, masks, ref); S_.insert(0, 'condition', cond); return F, S_


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--sid', required=True); ap.add_argument('--engdir', required=True); ap.add_argument('--b', type=int, required=True)
    ap.add_argument('--out', required=True); a = ap.parse_args(); t0 = time.time()
    L.ENG = Path(a.engdir)
    idx, u = panel_index(); dirA = healthy_dir(); ref = load_ref(); masks = stratum_masks(); out = []
    mats = load_mats(a.sid); F, S_ = run(mats, idx, dirA, masks, ref, 'raw'); out.append(S_)
    rng = np.random.default_rng(SEED + 7800 + a.b)
    for g in ('G1', 'G2', 'G3'):
        for f in FG:
            pm = {mk: ((g_flip(M, f, rng) if g == 'G3' else g_add(M, f, int(dirA[mk]), rng, g)), p) for mk, (M, p) in mats.items()}
            out.append(run(pm, idx, dirA, masks, ref, f'{g}_{f:g}')[1])
    R = pd.concat(out, ignore_index=True); R.insert(0, 'sid', a.sid); R.to_parquet(a.out, index=False)
    print(json.dumps(dict(sid=a.sid, b=a.b, seconds=time.time() - t0, conditions=sorted(R.condition.unique().tolist()))), flush=True)


if __name__ == '__main__':
    main()
