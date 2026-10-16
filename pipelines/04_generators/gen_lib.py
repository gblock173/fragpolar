"""generator step isolated copy of the panel toolkit toolkit (verbatim function bodies from the panel toolkit toolkit + the few helpers it takes from panel toolkit).
Only difference: no import of the p9/pa/lib_v2 chain; ENG (panel_calls root) is set by the caller (generator step data/engine)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from onepass_engine import gold_summary, trim_keep, marker_table

D = Path(CFG.DATA)
REFD = Path(CFG.inp('z_reference_dir'))
GEN_DIR = Path(CFG.work('generators')); ENG = GEN_DIR / 'data' / 'engine'
SEED = 20261008                                                     # pa_common.SEED
PF = ['amf', 'pdr', 'pm', 'me', 'mhl', 'uhl', 'fdrp', 'qfdrp', 'within_read_U', 'between_read_B_emp', 'hof', 'arrangement_excess_per_link']
AUXP = ['mu', 'fracM', 'fracU', 'strictM', 'strictU', 'mixed_strict', 'hopp_strict', 'n_k4']


def healthy_dir():
    return np.load(CFG.inp('healthy_direction_npy'))


def strata_markers():
    return pd.read_parquet(CFG.inp('strata_markers'))


def stratum_masks():
    st = strata_markers(); sd = pd.read_csv(CFG.inp('reference_strata_definition'), sep='\t'); sd = sd[sd.used]; M = {}
    for r in sd.itertuples():
        bm = np.ones(len(st), bool) if r.baseline == 'All' else (st.band == r.baseline).to_numpy()
        am = np.ones(len(st), bool) if r.axis == 'all' else (st.annot == r.annotation).to_numpy() if r.axis == 'annot' else \
            (st.genomic_context == r.annotation).to_numpy() if r.axis == 'context' else st['atlas_' + r.annotation].to_numpy()
        M[r.stratum] = bm & am
    return M


def panel_index():
    u = marker_table(); p = set(pd.read_csv(CFG.inp('generator_panel'), sep='\t', usecols=['marker_id']).marker_id)
    return np.flatnonzero(u.marker_id.isin(p).to_numpy()), u


def load_mats(sid, trim=False):
    """{global marker index: M} from panel_calls (optionally trimmed with the engine's deterministic window)."""
    d = pd.read_parquet(ENG / sid / 'panel_calls.parquet')
    if trim:
        q = np.load(CFG.inp('span_quantiles')); h = d.frag.to_numpy(np.uint64)
        u1 = (h & np.uint64(0xFFFFFFFF)).astype(np.float64) / 2 ** 32; u2 = (h >> np.uint64(32)).astype(np.float64) / 2 ** 32
        L = np.interp(u1, np.linspace(0, 1, len(q)), q); span = d.span.to_numpy(np.float64); fs = d.fs.to_numpy(np.float64)
        off = fs + np.floor(u2 * np.maximum(span - L, 0)); pos = d.pos.to_numpy(np.float64); d = d[(span <= L) | ((pos >= off) & (pos < off + L))]
    u = marker_table(); st = u.start.to_numpy(); en = u.end.to_numpy(); out = {}
    for mk, g in d.groupby('marker', sort=False):
        p = g.pos.to_numpy(np.int64); inr = (p >= st[mk]) & (p < en[mk]); g = g[inr]; p = p[inr]
        if not len(g): continue
        fi = pd.factorize(g.frag)[0]; pos = np.unique(p); M = np.full((fi.max() + 1, len(pos)), -1, np.int8); M[fi, np.searchsorted(pos, p)] = g.call.to_numpy(np.int8)
        out[int(mk)] = (M, pos)
    return out


def suff_from_M(M, dirA):
    cov = M >= 0; k = cov.sum(1); keep = k >= 4; n = int(keep.sum()); out = dict.fromkeys(['within_read_U', 'between_read_B_emp', 'hof', 'arrangement_excess_per_link'] + AUXP[:-1], np.nan)
    out['n_k4'] = n
    if n < 5: return out
    Mk = M[keep]; ck = cov[keep]; k = k[keep]; m = (Mk == 1).sum(1); p = m / k
    r, c = np.nonzero(ck); v = Mk[r, c]; sw = np.bincount(r[1:][(v[1:] != v[:-1]) & (r[1:] == r[:-1])], minlength=n)
    E = 2 * m * (k - m) / k; mu = p.mean(); vv = (p * p).mean(); fM = np.mean(p >= .7); fU = np.mean(p <= .3); sM = np.mean(m == k); sU = np.mean(m == 0)
    out.update(within_read_U=max(mu - vv, 0), between_read_B_emp=max(vv - mu * mu, 0), arrangement_excess_per_link=(E.sum() - sw.sum()) / (k - 1).sum(),
               hof=fM if dirA == 1 else (fU if dirA == 0 else np.nan), mu=mu, fracM=fM, fracU=fU, strictM=sM, strictU=sU, mixed_strict=1 - sM - sU,
               hopp_strict=sM if dirA == 1 else (sU if dirA == 0 else np.nan))
    return out


def features(mats, idx_all, dirA_all):
    """per-marker DataFrame (index = global marker index) of PF + AUXP."""
    rows = {}
    for mk, (M, _) in mats.items():
        g = gold_summary(M, 4, 4, 10, 10, 20); s = suff_from_M(M, dirA_all[mk])
        rows[mk] = dict(zip(PF[:8], g), **s)
    return pd.DataFrame.from_dict(rows, orient='index').reindex(idx_all)


def strata_summary(F, idx, masks, ref):
    """mean z per (stratum, feature) for PF and raw mean for AUXP; masks on the 190k axis restricted to the panel index."""
    out = []
    for f in PF:
        z = (F[f].to_numpy(float) - ref[f][0][idx]) / ref[f][1][idx]; z = np.clip(z, -10, 10)
        for sname, m in masks.items():
            x = z[m[idx]]; ok = np.isfinite(x); out.append((sname, f, 'z', float(x[ok].mean()) if ok.sum() >= 20 else np.nan, int(ok.sum())))
    for f in PF + AUXP[:-1]:
        x0 = F[f].to_numpy(float)
        for sname, m in masks.items():
            x = x0[m[idx]]; ok = np.isfinite(x); out.append((sname, f, 'raw', float(x[ok].mean()) if ok.sum() >= 20 else np.nan, int(ok.sum())))
    return pd.DataFrame(out, columns=['stratum', 'feature', 'scale', 'value', 'n'])


def load_ref():
    return {f: (np.load(REFD / f'{CFG.featname(f)}.npz')['center'][0], np.load(REFD / f'{CFG.featname(f)}.npz')['scale'][0]) for f in PF}


# ---------------------------------------------------------------- perturbations
def g_add(M, f, dirA, rng, kind):
    """G1 healthy-opposite uniform / G2 within-fragment discordant fragments, n_add = round(f/(1-f) * n) with positions borrowed from real fragments."""
    if dirA not in (0, 1) and kind == 'G1': return M
    cov = M >= 0; src = np.flatnonzero(cov.sum(1) >= 4); src = src if len(src) else np.flatnonzero(cov.any(1))
    n_add = int(round(f / (1 - f) * len(M)))
    if n_add == 0 or not len(src): return M
    pick = rng.choice(src, n_add, replace=True); A = np.full((n_add, M.shape[1]), -1, np.int8); C = cov[pick]
    if kind == 'G1':
        A[C] = dirA
    else:
        r, c = np.nonzero(C); first = np.r_[True, r[1:] != r[:-1]]; rank = np.arange(len(r)) - np.maximum.accumulate(np.where(first, np.arange(len(r)), 0))
        A[r, c] = ((rng.integers(0, 2, n_add)[r] + rank) % 2).astype(np.int8)
    return np.vstack([M, A])


def g_flip(M, f, rng):
    X = M.copy(); cov = X >= 0; fl = cov & (rng.random(X.shape) < f); X[fl] = 1 - X[fl]; return X


def mix_pos(base, srcs, fracs, rng):
    """base, srcs: (M, positions) tuples; returns mixed M on the union of positions with total fragments = base n."""
    pos = np.unique(np.concatenate([base[1]] + [s[1] for s in srcs if s is not None]))
    def lift(M, p):
        X = np.full((len(M), len(pos)), -1, np.int8); X[:, np.searchsorted(pos, p)] = M; return X
    B = lift(*base); n = len(B); parts = []; n_rep = 0
    for s, f in zip(srcs, fracs):
        if s is None or f <= 0 or not len(s[0]): continue
        k = int(rng.binomial(n, f)); n_rep += k
        if k: parts.append(lift(*s)[rng.integers(0, len(s[0]), k)])
    keep = rng.permutation(n)[:max(n - n_rep, 0)]
    return np.vstack([B[keep]] + parts) if parts else B
