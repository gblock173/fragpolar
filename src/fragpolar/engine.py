"""One streaming pass over a sample's CpG calls (row groups + carry-over of the open marker).

Copy of the analysis engine (one-pass engine; generator step variant with explicit --path/--outroot instead of a manifest).
Only the hard-coded directories were replaced by configuration lookups; calculations are unchanged.

Outputs in <outroot>/<sample>/:
  suff_{mode}.parquet  190,218 primary markers: extract.aggregate-equivalent sums over K>=4 fragments (same valid/duplicate/conflict rule)
                       + strict states (M=K, M=0) + AMF numerator/denominator (region calls of all fragments, production AMF definition)
  panel_calls.parquet  10k-panel cleaned call rows (marker index, fragment hash, pos, call, span) for simulations (untrimmed)
  gold_{mode}.parquet  10k-panel production features (k=4, min reads 10/10, min pairs 20) on uncleaned rows (gold_summary)
  bins.parquet         tissue only: unique fragments per 1 Mb bin (fragment_start) for the CNV identity check
  meta.json            QA counts, K and span histograms, timing
mode 'raw' for every sample; tissue additionally 'trim': each fragment is cut to a window whose length is drawn from the pooled healthy cfDNA span
distribution at a uniform offset, both from a deterministic hash of fragment_id (same window for the same fragment in every marker).
"""
import argparse, json, sys, time, resource, math
from pathlib import Path
import numpy as np, pandas as pd, pyarrow.parquet as pq
from scipy import sparse
from fragpolar import config as CFG


def _SPANQ():          # pooled healthy span quantiles (tissue 'trim' mode only)
    return CFG.inp('span_quantiles')


def _PANEL():          # 10,000-marker panel (marker_id column)
    return CFG.inp('generator_panel')


FIELDS = ['marker_id', 'fragment_id', 'chrom', 'pos', 'call', 'fragment_start', 'fragment_end', 'is_complete_pair']
GOLD = ['amf', 'pdr', 'pm', 'me', 'mhl', 'uhl', 'fdrp', 'qfdrp']
SUMS = ['n_reads', 'sum_ncpg', 'sum_switches', 'sum_expswitch', 'sum_mean_meth', 'sum_mean_meth_sq', 'n_full', 'n_fullU', 'n_strictM', 'n_strictU',
        'sum_fraglen', 'n_complete', 'amf_m', 'amf_c']


def marker_table():
    u = pd.read_parquet(CFG.inp('marker_regions'), columns=['marker_id', 'start', 'end', 'primary_eligible'])
    return u[u.primary_eligible].reset_index(drop=True)


def trim_keep(fid, fs, fe, pos, q):
    """Deterministic cfDNA-length window per fragment: length = healthy span quantile at u1, offset uniform at u2 (u1, u2 from hash(fragment_id))."""
    h = pd.util.hash_pandas_object(fid, index=False).to_numpy(np.uint64)
    u1 = (h & np.uint64(0xFFFFFFFF)).astype(np.float64) / 2 ** 32; u2 = (h >> np.uint64(32)).astype(np.float64) / 2 ** 32
    L = np.interp(u1, np.linspace(0, 1, len(q)), q); span = (fe - fs).astype(np.float64)
    off = fs + np.floor(u2 * np.maximum(span - L, 0)); return (span <= L) | ((pos >= off) & (pos < off + L))


def clean_sort(d, mkcode):
    """extract.aggregate cleaning: valid calls, drop (marker, fragment, pos) with conflicting calls, keep one of identical duplicates; sort."""
    d = d[d.call.isin([0, 1]) & d.pos.notna()]
    mk = mkcode(d.marker_id); fr = pd.factorize(d.fragment_id)[0]; pos = d.pos.to_numpy(np.int64)
    o = np.lexsort((pos, fr, mk)); d = d.iloc[o]; mk = mk[o]; fr = fr[o]; pos = pos[o]; call = d.call.to_numpy(np.int8)
    same = np.r_[False, (mk[1:] == mk[:-1]) & (fr[1:] == fr[:-1]) & (pos[1:] == pos[:-1])]
    n_dup = n_conf = 0
    if same.any():
        grp = np.cumsum(~same) - 1; g_n = np.bincount(grp); g_m = np.bincount(grp, weights=call)
        bad_g = (g_m > 0) & (g_m < g_n); dupg = g_n > 1; n_dup = int(g_n[dupg].sum()); n_conf = int(bad_g.sum())
        keep = ~same & ~bad_g[grp]; d = d[keep]; mk = mk[keep]; fr = fr[keep]; pos = pos[keep]; call = call[keep]
    return d, mk, fr, pos, call, n_dup, n_conf


def aggregate(mk, fr, pos, call, span, comp, inreg, nm, khist, shist):
    """per-marker sums (rows sorted by marker, fragment, pos and cleaned)."""
    out = np.zeros((nm, len(SUMS)))
    if not len(mk): return out
    new = np.r_[True, (mk[1:] != mk[:-1]) | (fr[1:] != fr[:-1])]; idx = np.flatnonzero(new); rid = np.cumsum(new) - 1
    c = call.astype(np.float64); k = np.bincount(rid); m = np.bincount(rid, weights=c)
    sw = np.bincount(rid, weights=(np.r_[False, call[1:] != call[:-1]] & ~new).astype(np.float64))
    out[:, SUMS.index('amf_m')] = np.bincount(mk, weights=c * inreg, minlength=nm); out[:, SUMS.index('amf_c')] = np.bincount(mk, weights=inreg.astype(np.float64), minlength=nm)
    keep = k >= 4; g = mk[idx][keep]; k = k[keep]; m = m[keep]; sw = sw[keep]; p = m / k; sp = span[idx][keep]; cp = comp[idx][keep]
    khist += np.bincount(np.minimum(k, 63), minlength=64); shist += np.bincount(np.clip(sp, 0, 999) // 10, minlength=100)
    vals = [np.ones(len(k)), k, sw, 2 * m * (k - m) / k, p, p * p, p >= .7, p <= .3, m == k, m == 0, sp, cp]
    for j, v in enumerate(vals): out[:, j] = np.bincount(g, weights=np.asarray(v, np.float64), minlength=nm)
    return out


def gold_summary(M, k=4, min_pdr=4, min_pm=10, min_mhl=10, min_pairs=20, min_common=2):
    """production compute_methylation_features_fast summary values (amf, pdr, pm, me, mhl, uhl, fdrp, qfdrp) without the per-window FDRP,
    pattern strings and pair lists that the summary does not use; same formulas and thresholds (vectorised). M: reads x region CpGs, -1 uncovered."""
    cov = M >= 0; N, C = M.shape; out = [np.nan] * 8
    tot = int(cov.sum()); out[0] = (M == 1).sum() / tot if tot else np.nan
    v = cov.sum(1) >= min_pdr; disc = v & (M == 1).any(1) & (M == 0).any(1)
    nc_ = (cov & v[:, None]).sum(0); nd_ = (cov & disc[:, None]).sum(0); ok = nc_ >= 1
    out[1] = float(np.mean(nd_[ok] / nc_[ok])) if ok.any() else np.nan
    W = C - k + 1
    if W > 0:
        Ms = np.maximum(M, 0).astype(np.int64); cs = np.zeros((N, C + 1), np.int32); cs[:, 1:] = np.cumsum(cov, 1)
        valid = (cs[:, k:] - cs[:, :W]) == k; pw = 1 << np.arange(k - 1, -1, -1, dtype=np.int64)
        pat = sum(Ms[:, b:b + W] * pw[b] for b in range(k)); nvw = valid.sum(0)
        cnt = np.bincount((np.arange(W)[None, :] * (1 << k) + pat)[valid], minlength=W * (1 << k)).reshape(W, 1 << k)
        use = nvw >= min_pm; den = int(nvw[use].sum()); pm_s = me_s = 0.0
        if use.any():
            f = cnt[use] / nvw[use, None]; lg = np.log2(np.where(f > 0, f, 1.0))
            pm_s = float(np.sum((1.0 - (f ** 2).sum(1)) * nvw[use])); me_s = float(np.sum((-1.0 / k) * (f * lg).sum(1) * nvw[use]))
        poss = np.zeros(k + 1); fm = np.zeros(k + 1); fu = np.zeros(k + 1); nvt = int(valid.sum())
        for L in range(1, k + 1):
            a1 = (1 << L) - 1
            for o in range(k - L + 1):
                sb = (pat >> (k - o - L)) & a1; poss[L] += nvt; fm[L] += int(((sb == a1) & valid).sum()); fu[L] += int(((sb == 0) & valid).sum())
        if den: out[2] = pm_s / den; out[3] = me_s / den
        Ls = [L for L in range(1, k + 1) if poss[L] >= min_mhl]
        if Ls:
            out[4] = sum(L * fm[L] / poss[L] for L in Ls) / sum(Ls); out[5] = sum(L * fu[L] / poss[L] for L in Ls) / sum(Ls)
    sub = cov.any(1); Mb = M[sub]; n = len(Mb)
    if n >= 2:
        if n * C > 40000:                                                 # large marker: only overlapping fragment pairs share CpGs (sparse)
            cf = sparse.csr_matrix((Mb >= 0).astype(np.float64)); mf = sparse.csr_matrix((Mb == 1).astype(np.float64)); iv = cf - mf
            S = sparse.triu(cf @ cf.T, 1).tocoo(); D = (cf @ cf.T - mf @ mf.T - iv @ iv.T).tocsr()
            r, c, nc = S.row, S.col, S.data; nd = np.rint(np.asarray(D[r, c]).ravel()); pv = nc >= min_common
            r, c, nc, nd = r[pv], c[pv], nc[pv], nd[pv]
            mk = lambda w: sparse.csr_matrix((w, (r, c)), shape=(n, n))
            npair = np.asarray(cf.multiply(mk(np.ones(len(r))) @ cf).sum(0)).ravel(); ndisc = np.asarray(cf.multiply(mk((nd > 0).astype(float)) @ cf).sum(0)).ravel()
            nham = np.asarray(cf.multiply(mk(nd / nc) @ cf).sum(0)).ravel()
        else:
            cf = (Mb >= 0).astype(np.float64); mf = np.where(Mb >= 0, Mb, 0).astype(np.float64); iv = cf - mf
            nc = cf @ cf.T; nd = np.rint(nc - mf @ mf.T - iv @ iv.T); up = np.triu(np.ones((n, n), bool), 1)
            pv = up & (nc >= min_common); Wp = pv.astype(np.float64); Wd = (pv & (nd > 0)).astype(np.float64)
            Wh = np.where(pv, nd / np.where(nc > 0, nc, 1), 0.0)
            npair = (cf * (Wp @ cf)).sum(0); ndisc = (cf * (Wd @ cf)).sum(0); nham = (cf * (Wh @ cf)).sum(0)   # pairs i<j covering c
        okp = npair >= min_pairs
        if okp.any(): out[6] = float(np.mean(ndisc[okp] / npair[okp])); out[7] = float(np.mean(nham[okp] / npair[okp]))
    return out


def gold(rows, lims):
    """per panel marker on uncleaned rows in file order (last write wins within a fragment, as production)."""
    out = []
    for mid, g in rows.groupby('marker_id', sort=False):
        ch, se = mid.split(':'); s, e = map(int, se.split('-')); p = g.pos.to_numpy(np.int64); r = g.call.to_numpy(np.int8)
        fi = pd.factorize(g.fragment_id)[0]; inr = (p >= s) & (p < e); pos = np.unique(p[inr]); M = np.full((fi.max() + 1, len(pos)), -1, np.int8)
        if len(pos): M[fi[inr], np.searchsorted(pos, p[inr])] = r[inr]
        nwc = int((M >= 0).any(1).sum())
        out.append([mid] + (gold_summary(M, 4, 4, *lims) if len(pos) else [np.nan] * 8) + [nwc])
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--sid', required=True); ap.add_argument('--limit', type=int); ap.add_argument('--outroot', required=True); ap.add_argument('--path', required=True); a = ap.parse_args(); t0 = time.time()   # generator step: no manifest, explicit cfDNA path/outroot
    from types import SimpleNamespace; row = SimpleNamespace(kind='cfDNA', role='healthy_base', path=a.path); tissue = False
    od = Path(a.outroot) / a.sid; od.mkdir(parents=True, exist_ok=True)
    u = marker_table(); ids = u.marker_id.to_numpy(); nm = len(u); pos_of = pd.Index(ids)
    mkcode = lambda s: pos_of.get_indexer(s.to_numpy())
    mstart = u.start.to_numpy(np.int64); mend = u.end.to_numpy(np.int64)
    panel = set(pd.read_csv(_PANEL(), sep='\t', usecols=['marker_id']).marker_id); pmask = pos_of.isin(panel)
    modes = ['raw', 'trim'] if tissue else ['raw']
    q = np.load(_SPANQ()) if tissue else None
    acc = {md: np.zeros((nm, len(SUMS))) for md in modes}; kh = {md: np.zeros(64) for md in modes}; sh = {md: np.zeros(100) for md in modes}
    gold_rows = {md: [] for md in modes}; pcalls = []; bins = {}; qa = dict(rows=0, rows_primary=0, dup_rows=0, conflict_sites=0, trimmed_rows_dropped=0)
    pf = pq.ParquetFile(row.path); nrg = pf.metadata.num_row_groups if a.limit is None else min(a.limit, pf.metadata.num_row_groups)
    carry = None

    def process(d):
        qa['rows'] += len(d)
        if tissue:
            fr1 = d.drop_duplicates('fragment_id'); b = fr1.chrom.astype(str) + ':' + (fr1.fragment_start.to_numpy() // 1_000_000).astype(str)
            for k_, v_ in b.value_counts().items(): bins[k_] = bins.get(k_, 0) + int(v_)
        d = d[d.marker_id.isin(pos_of)]; qa['rows_primary'] += len(d)
        if not len(d): return
        for md in modes:
            dd = d
            if md == 'trim':
                keep = trim_keep(d.fragment_id, d.fragment_start.to_numpy(np.int64), d.fragment_end.to_numpy(np.int64), d.pos.to_numpy(np.int64), q)
                qa['trimmed_rows_dropped'] += int((~keep).sum()); dd = d[keep]
            pr = dd[pmask[mkcode(dd.marker_id)]]
            if len(pr): gold_rows[md] += gold(pr, (10, 10, 20))
            c, mk, fr, pos, call, ndup, nconf = clean_sort(dd, mkcode)
            if md == 'raw':
                qa['dup_rows'] += ndup; qa['conflict_sites'] += nconf
                pm_ = pmask[mk]
                if pm_.any():
                    pcalls.append(pd.DataFrame(dict(marker=mk[pm_].astype(np.int32), frag=pd.util.hash_pandas_object(c.fragment_id[pm_], index=False).to_numpy(np.uint64),
                                                    pos=pos[pm_].astype(np.int32), call=call[pm_], span=(c.fragment_end.to_numpy() - c.fragment_start.to_numpy())[pm_].astype(np.int32), fs=c.fragment_start.to_numpy()[pm_].astype(np.int64))))
            span = (c.fragment_end.to_numpy() - c.fragment_start.to_numpy()).astype(np.int64); comp = c.is_complete_pair.to_numpy(bool)
            inreg = (pos >= mstart[mk]) & (pos < mend[mk])
            acc[md] += aggregate(mk, fr, pos, call, span, comp, inreg, nm, kh[md], sh[md])

    for j in range(nrg):
        d = pf.read_row_group(j, columns=FIELDS, use_threads=False).to_pandas()
        if carry is not None: d = pd.concat([carry, d], ignore_index=True)
        if not len(d): carry = None; continue
        last = d.marker_id.iloc[-1]; end = d.marker_id.eq(last).to_numpy(); carry = d[end]; d = d[~end]
        if len(d):
            process(d)
        if j % 25 == 0: print(json.dumps(dict(rg=j, of=nrg, sec=round(time.time() - t0))), flush=True)
    if carry is not None and len(carry): process(carry)
    for md in modes:
        s = pd.DataFrame(acc[md], columns=SUMS); s.insert(0, 'marker_id', ids); s.to_parquet(od / f'suff_{md}.parquet', index=False, compression='zstd')
        gdf = pd.DataFrame(gold_rows[md], columns=['marker_id'] + GOLD + ['n_reads_with_calls'])
        assert not gdf.marker_id.duplicated().any(), 'panel marker split across chunks'
        gdf.to_parquet(od / f'gold_{md}.parquet', index=False)
    (pd.concat(pcalls, ignore_index=True) if pcalls else pd.DataFrame()).to_parquet(od / 'panel_calls.parquet', index=False, compression='zstd')
    if tissue: pd.Series(bins, name='n_fragments').rename_axis('bin').reset_index().to_parquet(od / 'bins.parquet', index=False)
    n = acc['raw'][:, 0]
    meta = dict(sid=a.sid, kind=row.kind, role=row.role, path=row.path, row_groups=nrg, seconds=time.time() - t0,
                peak_RSS_GiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20, QA=qa, markers_with_K4=int((n > 0).sum()),
                median_n_K4=float(np.median(n[n > 0])) if (n > 0).any() else 0.0,
                K_hist={md: kh[md].tolist() for md in modes}, span_hist_10bp={md: sh[md].tolist() for md in modes},
                gold_markers={md: len(gold_rows[md]) for md in modes})
    (od / 'meta.json').write_text(json.dumps(meta)); print(json.dumps({k: v for k, v in meta.items() if k not in ('K_hist', 'span_hist_10bp')}), flush=True)


if __name__ == '__main__':
    main()
