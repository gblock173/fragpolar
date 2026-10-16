"""Small pure functions that restate, expression for expression, the core formulas written inline in the pipeline scripts.

The pipeline scripts are kept as close as possible to the code that produced the manuscript, so they do not import this
module; it exists so that the formulas can be unit-tested on synthetic input (tests/). Each function names the script
and the expression it mirrors.
"""
import numpy as np
import pandas as pd


# ---------------------------------------------------------------- fragment states and variance decomposition
def fragment_p(M):
    """M: fragments x CpG matrix (1 methylated, 0 unmethylated, -1 not observed) -> (K, M_count, p) for every row."""
    obs = M >= 0
    K = obs.sum(1); m = (M == 1).sum(1)
    with np.errstate(invalid='ignore', divide='ignore'):
        p = m / K
    return K, m, p


def decompose(p):
    """mu = E[p], WFV = E[p(1-p)], BFV = E[p^2] - mu^2 (population variance, denominator n); WFV + BFV = mu(1-mu).
    Same expressions as extract.check_metrics (u = mu - q, b = q - mu^2) and gen_lib.suff_from_M."""
    p = np.asarray(p, float); mu = p.mean(); q = (p * p).mean()
    return dict(mu=mu, WFV=mu - q, BFV=q - mu * mu, V=mu * (1 - mu))


def strict_states(m, K, healthy_mu):
    """strict-mixed (0 < M < K), sOUF (uniform opposite to healthy) and healthy-direction uniform fractions.
    healthy_mu < 0.5 -> healthy unmethylated -> sOUF = fraction with M = K (ge1_aggregate: opp = strictM if baseline < .5 else strictU)."""
    m = np.asarray(m); K = np.asarray(K)
    sM = np.mean(m == K); sU = np.mean(m == 0)
    opp, conc = (sM, sU) if healthy_mu < .5 else (sU, sM)
    return dict(strictM=sM, strictU=sU, mixed=1 - sM - sU, sOUF=opp, concordant=conc)


# ---------------------------------------------------------------- FPI and FPI_K (stratum level, ratio of marker means)
def fpi_ratio(B, V):
    """FPI of a stratum = sum(BFV) / sum(mu(1-mu)) over valid markers (ge2_analysis: Rr = B / V on marker means)."""
    B = np.asarray(B, float); V = np.asarray(V, float); return B.sum() / V.sum()


def fpi_k(B, V, meanK):
    """FPI_K = sum(B - h V) / sum((1 - h) V), h = 1/mean_K per marker (fk4_analysis / fk5_jensen_raw: (sb - sh1)/(sv - sh1))."""
    B = np.asarray(B, float); V = np.asarray(V, float); h = 1.0 / np.asarray(meanK, float)
    sb, sv, sh = B.sum(), V.sum(), (h * V).sum()
    return (sb - sh) / (sv - sh)


def delta_bfv_decomposition(V_c, F_c, V_h, F_h, B_c, B_h):
    """dB = dV * F_bar + dF * V_bar + residual (ge2_analysis group-effect step decomposition: tV = dV * Rm, tR = dR * Vm, res = dB - tV - tR)."""
    dV, dF, dB = V_c - V_h, F_c - F_h, B_c - B_h; Fm, Vm = (F_c + F_h) / 2, (V_c + V_h) / 2
    tV, tR = dV * Fm, dF * Vm
    return dict(term_mean_shift=tV, term_fpi_change=tR, residual=dB - tV - tR)


# ---------------------------------------------------------------- analytic sOUF removal (ge1_aggregate margin analysis)
def remove_opposite_analytic(mu, B, f, p0):
    """Remove a fraction f of fragments with p = p0 (opposite strict-uniform) from a marker with mean mu and BFV B.
    ge1_aggregate: mu2 = (mu - f p0)/(1 - f); m2 = (B + mu^2 - f p0)/(1 - f); Bdel = max(m2 - mu2^2, 0)  (p0^2 = p0)."""
    with np.errstate(invalid='ignore', divide='ignore'):
        mu2 = (mu - f * p0) / (1 - f); m2 = (B + mu * mu - f * p0) / (1 - f)
        return np.maximum(m2 - mu2 ** 2, 0)


# ---------------------------------------------------------------- matched-dose removal from sums (df2_sample)
def var_from_sums(S1, S2, n):
    """df2_sample.var / df3_stratum.var: max(S2/n - (S1/n)^2, 0)."""
    with np.errstate(invalid='ignore', divide='ignore'):
        mu = S1 / n; return np.maximum(S2 / n - mu * mu, 0)


def bfv_after_removal(S1, S2, n, r, p_removed):
    """BFV after removing r fragments that all have p = p_removed in {0, 1} (df2_sample: var(S1 - r p0, S2 - r p0, n - r))."""
    return var_from_sums(S1 - r * p_removed, S2 - r * p_removed, n - r)


def reduction(d_before, d_after):
    """1 - dBFV(after removal) / dBFV(before) (df4_stats: red = 1 - po/pf)."""
    return 1 - d_after / d_before


# ---------------------------------------------------------------- read-split hash (df1_extract)
SALT = 'P03_readsplit_v1'
HASH_KEY = 'P03readsplitv1__'


def readsplit_half(sid, fragment_ids):
    """Half (0/1) of each fragment = top bit of pandas siphash(hash_key) of 'sid|SALT|fragment_id' (marker and state not hashed).
    Same expression as df1_extract.process."""
    uq = pd.Series(list(fragment_ids))
    h = pd.util.hash_pandas_object(uq.radd(sid + '|' + SALT + '|'), index=False, hash_key=HASH_KEY).to_numpy(np.uint64)
    return ((h >> np.uint64(63)) & np.uint64(1)).astype(np.int8)


# ---------------------------------------------------------------- mixture-expected BFV (ml1_mixture_model, model M1)
def expected_dbfv_mixture(delta, mu_h, B_h, p_o):
    """Mix healthy fragments with a fraction delta of opposite strict-uniform fragments (p = p_o):
    BFV_mix = (1 - delta) B_h + delta (1 - delta)(p_o - mu_h)^2, so dBFV = delta (1 - delta)(p_o - mu_h)^2 - delta B_h
    (ml1_mixture_model exp_M1 with a single healthy distribution; Qo = p_o - 2 p_o mu_h + E[mu^2] reduces to (p_o - mu_h)^2)."""
    return delta * (1 - delta) * (p_o - mu_h) ** 2 - delta * B_h
