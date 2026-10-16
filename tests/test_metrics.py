"""Unit tests of the core fragment-state formulas on synthetic fragments (no data needed; runs in seconds)."""
import itertools, math
from pathlib import Path
import numpy as np, pandas as pd, pytest
from fragpolar import metrics as M

REPO = Path(__file__).resolve().parents[1]


def random_fragments(rng, n=400, K=6, p_meth=.3, frac_uniform=.2):
    """fragments x CpG 0/1 matrix: a share of fully methylated / unmethylated fragments plus independent CpG calls."""
    X = (rng.random((n, K)) < p_meth).astype(int)
    u = rng.random(n) < frac_uniform; X[u] = (rng.random(u.sum()) < .5)[:, None].astype(int)
    return X


# ---------------------------------------------------------------- variance decomposition
def test_decomposition_identity_random():
    rng = np.random.default_rng(1)
    for _ in range(20):
        p = rng.random(rng.integers(5, 500))
        d = M.decompose(p)
        assert abs(d['WFV'] + d['BFV'] - d['mu'] * (1 - d['mu'])) < 1e-12
        assert d['WFV'] >= -1e-15 and d['BFV'] >= -1e-15
        assert abs(d['BFV'] - np.var(p)) < 1e-12                       # population variance (denominator n)


def test_decomposition_extremes():
    d = M.decompose([0, 0, 1, 1]); assert d['WFV'] == 0 and d['BFV'] == .25          # opposite uniform mixture: all between
    d = M.decompose([.5] * 4); assert d['BFV'] == 0 and d['WFV'] == .25              # within-fragment mixing only


@pytest.mark.parametrize('K', [4, 5, 8])
def test_fpi_independent_cpgs_equals_one_over_K(K):
    # exact binomial enumeration, p = 1/2 (same identity as the extraction self-test)
    ps = []; w = []
    for m in range(K + 1):
        ps.append(m / K); w.append(math.comb(K, m) * .5 ** K)
    ps, w = np.array(ps), np.array(w); mu = (w * ps).sum(); bfv = (w * ps * ps).sum() - mu ** 2
    assert abs(M.fpi_ratio([bfv], [mu * (1 - mu)]) - 1 / K) < 1e-12
    assert abs(M.fpi_k([bfv], [mu * (1 - mu)], [K])) < 1e-12                       # FPI_K = 0 without within-fragment correlation


def test_fpi_k_one_for_uniform_fragments():
    p = np.array([0, 1, 1, 0, 1]); d = M.decompose(p)
    assert abs(M.fpi_ratio([d['BFV']], [d['V']]) - 1) < 1e-12 and abs(M.fpi_k([d['BFV']], [d['V']], [6]) - 1) < 1e-12


def test_fpi_ratio_of_means():
    B = np.array([.01, .05, .02]); V = np.array([.1, .2, .15])
    assert abs(M.fpi_ratio(B, V) - B.sum() / V.sum()) < 1e-15


def test_delta_bfv_decomposition_sums_to_total():
    Vc, Fc, Vh, Fh = .21, .42, .23, .38; Bc, Bh = Vc * Fc, Vh * Fh
    d = M.delta_bfv_decomposition(Vc, Fc, Vh, Fh, Bc, Bh)
    assert abs(d['term_mean_shift'] + d['term_fpi_change'] + d['residual'] - (Bc - Bh)) < 1e-15
    assert abs(d['residual']) < 1e-15                                            # with B = V * F the midpoint split is exact


# ---------------------------------------------------------------- fragment states / sOUF
def test_strict_states_and_souf_direction():
    X = np.array([[1, 1, 1, 1], [0, 0, 0, 0], [0, 0, 0, 0], [0, 1, 0, 0], [1, 1, 1, 1]])
    K, m, p = M.fragment_p(X)
    s = M.strict_states(m, K, healthy_mu=.2)          # healthy unmethylated -> sOUF = fully methylated
    assert s['sOUF'] == .4 and s['concordant'] == .4 and abs(s['mixed'] - .2) < 1e-15
    s = M.strict_states(m, K, healthy_mu=.8)          # healthy methylated -> sOUF = fully unmethylated
    assert s['sOUF'] == .4 and s['concordant'] == .4


def test_fragment_p_ignores_unobserved():
    X = np.array([[1, -1, 0, 1, 1], [-1, -1, 0, 0, 0]])
    K, m, p = M.fragment_p(X)
    assert list(K) == [4, 3] and list(m) == [3, 0] and abs(p[0] - .75) < 1e-15


# ---------------------------------------------------------------- sOUF removal (analytic vs direct)
@pytest.mark.parametrize('p0', [0., 1.])
def test_analytic_removal_equals_direct(p0):
    rng = np.random.default_rng(2)
    X = random_fragments(rng); K, m, p = M.fragment_p(X)
    opp = p == p0; f = opp.mean(); d = M.decompose(p)
    direct = np.var(p[~opp])
    analytic = M.remove_opposite_analytic(d['mu'], d['BFV'], f, p0)
    assert abs(direct - analytic) < 1e-12
    # same quantity from sums (matched-dose code path, r removed fragments with p = p0)
    n = len(p); r = opp.sum()
    assert abs(M.bfv_after_removal(p.sum(), (p * p).sum(), n, r, p0) - direct) < 1e-12


def test_random_removal_from_sums():
    rng = np.random.default_rng(3); p = rng.random(300); rem = rng.choice(300, 40, replace=False); keep = np.setdiff1d(np.arange(300), rem)
    assert abs(M.var_from_sums(p.sum() - p[rem].sum(), (p * p).sum() - (p[rem] ** 2).sum(), 260) - np.var(p[keep])) < 1e-12


def test_reduction():
    assert M.reduction(0.01, 0.003) == pytest.approx(.7)


# ---------------------------------------------------------------- read-split hash
def test_readsplit_deterministic_and_balanced():
    ids = [f'frag{i}' for i in range(20000)]
    a = M.readsplit_half('S1', ids); b = M.readsplit_half('S1', ids)
    assert (a == b).all() and set(np.unique(a)) <= {0, 1}
    assert abs(a.mean() - .5) < .02                                             # ~50/50
    c = M.readsplit_half('S2', ids)
    assert 0.4 < (a == c).mean() < 0.6                                          # sample identifier enters the hash
    # the same fragment gets the same half whatever marker it appears in (marker is not hashed)
    sub = M.readsplit_half('S1', ids[100:200]); assert (sub == a[100:200]).all()


# ---------------------------------------------------------------- mixture-expected BFV
@pytest.mark.parametrize('p_o', [0., 1.])
def test_mixture_expected_bfv_matches_weighted_mixture(p_o):
    rng = np.random.default_rng(4)
    p_h = rng.beta(2, 5, 1000); w_h = np.full(1000, 1 / 1000)
    mu_h = (w_h * p_h).sum(); B_h = (w_h * p_h ** 2).sum() - mu_h ** 2
    for delta in [.005, .02, .1]:
        p = np.r_[p_h, p_o]; w = np.r_[(1 - delta) * w_h, delta]
        mu = (w * p).sum(); B_mix = (w * p * p).sum() - mu ** 2
        assert abs((B_mix - B_h) - M.expected_dbfv_mixture(delta, mu_h, B_h, p_o)) < 1e-12


# ---------------------------------------------------------------- the formulas above are the ones written in the pipelines
@pytest.mark.parametrize('rel, text', [
    ('pipelines/01_group_effects/ge1_aggregate.py', 'mu2 = (mu - f * p0) / (1 - f); m2 = (B + mu * mu - f * p0) / (1 - f)'),
    ('pipelines/01_group_effects/ge1_aggregate.py', 'Bdel = np.maximum(m2 - mu2 ** 2, 0)'),
    ('pipelines/05_direct_fragments/df2_sample.py', "Bopp=np.where(okd,var(W['S1']-r*p0f,W['S2']-r*p0f,n2),np.nan)"),
    ('pipelines/05_direct_fragments/df2_sample.py', 'mu=S1/n;return np.maximum(S2/n-mu*mu,0)'),
    ('pipelines/05_direct_fragments/df1_extract.py', "hb=((pd.util.hash_pandas_object(pd.Series(uq).radd(sid+'|'+SALT+'|'),index=False,hash_key='P03readsplitv1__').to_numpy(np.uint64)>>np.uint64(63))&np.uint64(1)).astype(np.int8)"),
    ('pipelines/05_direct_fragments/df1_extract.py', "SALT='P03_readsplit_v1'"),
    ('pipelines/05_direct_fragments/df4_stats.py', 'red=1-po/pf'),
    ('pipelines/09_mixture_and_length/ml1_mixture_model.py', "out['exp_M1'] = d * (1 - d) * Qo_s - d * Bh_s"),
    ('pipelines/09_mixture_and_length/ml1_mixture_model.py', 'Qo = po - 2 * po * mu_h + mu2_h'),
    ('pipelines/01_group_effects/ge2_analysis.py', 'tV = dV * Rm; tR = dR * Vm; res = dB - tV - tR'),
    ('pipelines/03_fpi_k/fk5_jensen_raw.py', 'fk1 = (sb - sh1) / (sv - sh1)'),
])
def test_formula_text_in_pipeline(rel, text):
    assert text in (REPO / rel).read_text(), f'{rel}: expression changed; update fragpolar/metrics.py and this test'
