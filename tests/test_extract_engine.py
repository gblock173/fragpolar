"""Tests of the CpG-call cleaning / aggregation code (fragpolar.extract, fragpolar.engine) on synthetic CpG-call tables."""
import numpy as np, pandas as pd, pytest
from fragpolar import extract, engine, feature_definitions as FD

QA0 = lambda: {x: 0 for x in ['invalid_calls', 'duplicate_rows', 'conflict_sites', 'eligible_fragment_marker_pairs']}


def synthetic_calls(rng, n_markers=3, n_frag=60, max_k=8):
    rows = []
    for mk in range(n_markers):
        mid = f'chr1:{mk * 1000}-{mk * 1000 + 500}'
        for f in range(n_frag):
            k = int(rng.integers(1, max_k + 1)); state = rng.random()
            calls = (rng.random(k) < (1. if state < .2 else 0. if state < .4 else .5)).astype(int)
            fs = mk * 1000 + int(rng.integers(0, 100)); fe = fs + int(rng.integers(120, 300))
            for j, c in enumerate(calls):
                rows.append(dict(marker_id=mid, fragment_id=f'f{mk}_{f}', chrom='chr1', pos=mk * 1000 + 10 + 2 * j, call=c,
                                 fragment_start=fs, fragment_end=fe, is_complete_pair=bool(f % 3)))
    d = pd.DataFrame(rows)
    # inject an identical duplicate (kept once), a conflicting duplicate (site removed) and an invalid call
    dup = d.iloc[[0]]; conf = d.iloc[[5]].assign(call=1 - d.call.iloc[5]); bad = d.iloc[[7]].assign(call=2)
    return pd.concat([d, dup, conf, bad], ignore_index=True)


def direct_stats(d):
    """Reference: clean as documented, then per marker over K >= 4 fragments."""
    d = d[d.call.isin([0, 1]) & d.pos.notna()]
    g = d.groupby(['marker_id', 'fragment_id', 'pos']).call
    ok = g.transform('nunique') == 1
    d = d[ok].drop_duplicates(['marker_id', 'fragment_id', 'pos'])
    fr = d.groupby(['marker_id', 'fragment_id']).call.agg(['size', 'sum']).rename(columns={'size': 'k', 'sum': 'm'})
    fr = fr[fr.k >= 4]; fr['p'] = fr.m / fr.k
    out = fr.groupby(level=0).agg(n=('p', 'size'), S1=('p', 'sum'), S2=('p', lambda x: (x * x).sum()),
                                  sM=('p', lambda x: (x == 1).sum()), sU=('p', lambda x: (x == 0).sum()), K=('k', 'sum'))
    return out


def test_upstream_synthetic_selfcheck():
    cases = extract.synthetic_checks()
    assert cases['mixture']['B'] == .25 and cases['mixture']['U'] == 0


def test_extract_aggregate_matches_direct():
    rng = np.random.default_rng(10); d = synthetic_calls(rng); qa = QA0()
    a = extract.aggregate(d, qa).set_index('marker_id'); r = direct_stats(d)
    assert qa['conflict_sites'] == 1 and qa['invalid_calls'] == 1
    a = a.loc[r.index]
    assert np.allclose(a.n_reads, r.n) and np.allclose(a.sum_mean_meth, r.S1) and np.allclose(a.sum_mean_meth_sq, r.S2)
    assert np.allclose(a.n_strictM, r.sM) and np.allclose(a.n_strictU, r.sU) and np.allclose(a.sum_ncpg, r.K)
    info = extract.check_metrics(a.reset_index())                     # WFV + BFV = mu(1 - mu) per marker
    assert info['identity_max_error'] < 1e-10


def test_engine_aggregate_matches_extract():
    rng = np.random.default_rng(11); d = synthetic_calls(rng)
    ids = np.array(sorted(d.marker_id.unique())); pos_of = pd.Index(ids); mkcode = lambda s: pos_of.get_indexer(s.to_numpy())
    c, mk, fr, pos, call, ndup, nconf = engine.clean_sort(d, mkcode)
    span = (c.fragment_end.to_numpy() - c.fragment_start.to_numpy()).astype(np.int64); comp = c.is_complete_pair.to_numpy(bool)
    out = engine.aggregate(mk, fr, pos, call, span, comp, np.ones(len(mk), bool), len(ids), np.zeros(64), np.zeros(100))
    e = pd.DataFrame(out, columns=engine.SUMS, index=ids)
    a = extract.aggregate(d, QA0()).set_index('marker_id').loc[ids]
    for ce, ca in [('n_reads', 'n_reads'), ('sum_ncpg', 'sum_ncpg'), ('sum_mean_meth', 'sum_mean_meth'), ('sum_mean_meth_sq', 'sum_mean_meth_sq'),
                   ('n_strictM', 'n_strictM'), ('n_strictU', 'n_strictU'), ('sum_switches', 'sum_switches'), ('sum_expswitch', 'sum_expswitch')]:
        assert np.allclose(e[ce].to_numpy(), a[ca].to_numpy()), ce
    assert nconf == 1


@pytest.mark.parametrize('feat', ['AMF', 'PDR', 'PM', 'ME', 'MHL', 'UHL', 'FDRP', 'qFDRP'])
def test_gold_summary_matches_reference_definitions(feat):
    """engine.gold_summary (vectorised production summary) vs feature_definitions.compute_features (gate-free reference)
    on the 24 toy fragment sets, with all gates set to 1."""
    idx = ['AMF', 'PDR', 'PM', 'ME', 'MHL', 'UHL', 'FDRP', 'qFDRP'].index(feat)
    for key, (lo, hi, _) in FD.EXAMPLES.items():
        for reads in (lo, hi):
            M = FD._matrix(reads).astype(np.int8)
            g = engine.gold_summary(M, k=4, min_pdr=4, min_pm=1, min_mhl=1, min_pairs=1)
            v, _ = FD.compute_features(reads)
            assert (np.isnan(g[idx]) and np.isnan(v[feat])) or abs(g[idx] - v[feat]) < 1e-12, (key, reads, g[idx], v[feat])
