"""Reference implementation of the 12 fragment summaries on small toy fragment sets (Fig. 1b definition cards).

Verbatim copy of sections 1-2 of the figure-gallery module (compute_features, example fragment sets, value formatting,
derivation strings). Plotting, layout and cross-check code of that module are not included.

compute_features() is an independent implementation of the formula definitions (AMF, PDR, PM, ME, MHL, UHL, FDRP,
qFDRP on the call-cache definitions; U (= WFV), B (= BFV), A (= HOF), J on K >= 4 fragments). Production depth gates
(PM/ME n_w >= 10, MHL/UHL N_l >= 10, FDRP/qFDRP >= 20 pairs per CpG) default to 1 here (formula value of the reads given).
J = (sum E - sum S)/sum(K - 1), S = observed switches, E = 2M(K - M)/K.
"""
import itertools, math
from collections import Counter
import numpy as np

F12 = ['AMF', 'PDR', 'PM', 'ME', 'MHL', 'UHL', 'FDRP', 'qFDRP', 'U', 'B', 'A', 'J']


# ════════════════════════════════════════════════════════════════════════════════════════════════════════════════
# 1. independent implementation of the 12 formulas (FEATURE_FORMULAS_KO.md §1-2)
# ════════════════════════════════════════════════════════════════════════════════════════════════════════════════
def _matrix(reads):
    """reads: equal-length strings over {'0','1','.'} in genomic CpG order -> int matrix (1 meth, 0 unmeth, -1 no call)."""
    return np.array([[1 if c == '1' else 0 if c == '0' else -1 for c in r] for r in reads], dtype=int)


def compute_features(reads, healthy_unmethylated=True, k=4, min_calls_pdr=4, min_window_reads=1, min_substrings=1,
                     min_pairs=1, min_shared=2):
    """Return (values, aux). Gates default to 1 (= formula value of the reads given); production: 10 / 10 / 20."""
    X = _matrix(reads); obs = X >= 0; N, C = X.shape; v = {}; aux = {}
    # AMF — methylated calls / all calls (every fragment, every CpG)
    ncall, nmeth = int(obs.sum()), int((X == 1).sum())
    v['AMF'] = nmeth / ncall if ncall else math.nan; aux.update(n_calls=ncall, n_meth=nmeth)
    # PDR — per CpG: among reads with >= 4 calls covering it, share holding both states; CpG-equal mean
    K = obs.sum(1); elig = K >= min_calls_pdr
    mixed = elig & (X == 1).any(1) & (X == 0).any(1)
    per = [mixed[elig & obs[:, j]].mean() for j in range(C) if (elig & obs[:, j]).any()]
    v['PDR'] = float(np.mean(per)) if per else math.nan; aux.update(mixed_rows=np.flatnonzero(mixed).tolist(), n_pdr_reads=int(elig.sum()))
    # PM / ME / MHL / UHL — consecutive k-windows over the observed CpGs; reads covering the whole window
    cols = [j for j in range(C) if obs[:, j].any()]
    pm_s = me_s = 0.0; den = 0; NL = np.zeros(k + 1); ML = np.zeros(k + 1); UL = np.zeros(k + 1); npat = []
    for w in range(len(cols) - k + 1):
        cw = cols[w:w + k]; full = obs[:, cw].all(1); P = X[np.ix_(full, cw)]; n = len(P)
        if n == 0:
            continue
        f = np.array(list(Counter(map(tuple, P)).values()), float) / n; npat.append(len(f))
        if n >= min_window_reads:
            pm_s += n * (1.0 - float(np.sum(f * f))); me_s += n * (-float(np.sum(f * np.log2(f))) / k); den += n
        for row in P:
            for l in range(1, k + 1):
                for o in range(k - l + 1):
                    s = row[o:o + l]; NL[l] += 1; ML[l] += bool((s == 1).all()); UL[l] += bool((s == 0).all())
    v['PM'] = pm_s / den if den else math.nan; v['ME'] = me_s / den if den else math.nan
    ls = [l for l in range(1, k + 1) if NL[l] > 0 and NL[l] >= min_substrings]
    v['MHL'] = sum(l * ML[l] / NL[l] for l in ls) / sum(ls) if ls else math.nan
    v['UHL'] = sum(l * UL[l] / NL[l] for l in ls) / sum(ls) if ls else math.nan
    aux.update(n_patterns=npat)
    # FDRP / qFDRP — read pairs sharing >= 2 CpGs; per CpG covered by both: any-difference share / mean differing fraction
    rows = np.flatnonzero(obs.any(1)); pc = np.zeros(C); pdis = np.zeros(C); pham = np.zeros(C); npair = ndis = 0; hsum = 0.0
    for a, b in itertools.combinations(rows, 2):
        sh = obs[a] & obs[b]; ns = int(sh.sum())
        if ns < min_shared:
            continue
        nd = int((X[a, sh] != X[b, sh]).sum())
        pc[sh] += 1; pdis[sh] += nd > 0; pham[sh] += nd / ns; npair += 1; ndis += nd > 0; hsum += nd / ns
    ok = (pc > 0) & (pc >= min_pairs)
    v['FDRP'] = float(np.mean(pdis[ok] / pc[ok])) if ok.any() else math.nan
    v['qFDRP'] = float(np.mean(pham[ok] / pc[ok])) if ok.any() else math.nan
    aux.update(n_pairs=npair, n_disc_pairs=ndis, ham_sum=hsum)
    # U / B / A / J — K >= 4 fragments, fragment-equal weights
    r4 = np.flatnonzero(K >= 4)
    if len(r4):
        Kr = K[r4].astype(float); Mr = np.array([(X[r] == 1).sum() for r in r4], float); p = Mr / Kr
        S = np.array([int(np.sum(X[r][obs[r]][1:] != X[r][obs[r]][:-1])) for r in r4], float)
        E = 2 * Mr * (Kr - Mr) / Kr
        mu = float(p.mean()); p2 = float(np.mean(p * p))
        v['U'] = max(mu - p2, 0.0); v['B'] = max(p2 - mu * mu, 0.0)
        v['A'] = float(np.mean(p >= .7)) if healthy_unmethylated else float(np.mean(p <= .3))
        L = float(np.sum(Kr - 1)); v['J'] = (float(E.sum()) - float(S.sum())) / L if L > 0 else math.nan
        aux.update(p=p.tolist(), S=S.tolist(), E=E.tolist(), mu=mu, Q=float(S.sum()) / L if L else math.nan,
                   Q0=float(E.sum()) / L if L else math.nan, n_K4=len(r4))
    else:
        v.update(U=math.nan, B=math.nan, A=math.nan, J=math.nan); aux.update(n_K4=0)
    return v, aux


PROD_GATES = dict(min_window_reads=10, min_substrings=10, min_pairs=20)

# ════════════════════════════════════════════════════════════════════════════════════════════════════════════════
# 2. gallery content
# ════════════════════════════════════════════════════════════════════════════════════════════════════════════════
# (Low reads, High reads, highlight). Slide-3 patterns kept where correct; changed: AMF High (14/16 instead of 16/16, shows
# a fraction), A (healthy-like Low; High mixes a fully and a p = .75 methylated read: p >= .7, not 'fully'), qFDRP Low
# (pairs differing at 1 of 4 CpGs, shows the quantitative weighting; slide Low = identical reads).
EXAMPLES = {
    'AMF': (['0000', '0010', '0100', '0000'], ['1111', '1011', '1101', '1111'], None),
    'PDR': (['0000', '0000', '1111', '1111'], ['0001', '1110', '0001', '1110'], 'mixed'),
    'U': (['0000', '0000', '1111', '1111'], ['0011', '0011', '0011', '0011'], 'mixed'),
    'MHL': (['0101'] * 4, ['0111'] * 4, 'run1'),
    'UHL': (['0101'] * 4, ['0001'] * 4, 'run0'),
    'J': (['0101'] * 4, ['0011'] * 4, 'switch'),
    'A': (['0000', '0001', '0000', '0100'], ['0000', '1111', '0111', '0000'], 'tail'),
    'PM': (['0011'] * 4, ['0111', '0011', '1010', '1100'], None),
    'ME': (['0011'] * 4, ['0011', '1011', '0101', '1100'], None),
    'B': (['0011'] * 4, ['0000', '0000', '1111', '1111'], None),
    'FDRP': (['0011'] * 4, ['0011', '0101', '0100', '1100'], (0, 1)),
    'qFDRP': (['0011', '0011', '0111', '0011'], ['0011', '0011', '1100', '1100'], (1, 2)),
}
HEALTHY_UNMETH = True     # stated baseline for A (U-marker convention)

DEFINITION = {
    'AMF': 'Fraction of observed CpG calls that are methylated',
    'PDR': 'Fraction of reads (≥4 CpGs) that carry both methylation states',
    'U': 'Mean within-read mixing, mean of p(1−p)',
    'MHL': 'Length-weighted load of fully methylated runs',
    'UHL': 'Length-weighted load of fully unmethylated runs',
    'J': 'Switch deficit vs random CpG order, (E−S)/(K−1); S observed, E expected switches',
    'A': 'Fraction of reads in the tail opposite to healthy (healthy unmethylated: p ≥ 0.7)',
    'PM': 'Diversity of complete 4-CpG patterns, 1 − Σf²',
    'ME': 'Diversity and evenness of 4-CpG patterns, −¼ Σ f·log₂f',
    'B': 'Spread of read methylation fractions p across reads (variance)',
    'FDRP': 'Fraction of read pairs that differ at any shared CpG',
    'qFDRP': 'Mean fraction of shared CpGs that differ within read pairs',
}

ROWS = {6: [['AMF', 'PDR', 'U', 'MHL', 'UHL', 'J'], ['A', 'PM', 'ME', 'B', 'FDRP', 'qFDRP']],
        4: [['AMF', 'MHL', 'UHL', 'J'], ['PDR', 'U', 'PM', 'ME'], ['B', 'FDRP', 'qFDRP', 'A']],
        3: [['AMF', 'PDR', 'U'], ['MHL', 'UHL', 'J'], ['A', 'PM', 'ME'], ['B', 'FDRP', 'qFDRP']]}


def fmt_value(x, signed=False):
    if not np.isfinite(x):
        return 'NA'
    r = round(float(x), 3); r = 0.0 if r == 0 else r
    s = f'{abs(r):.3f}'.rstrip('0').rstrip('.')
    sign = ('+' if r > 0 else '−' if r < 0 else '') if signed else ('−' if r < 0 else '')
    return sign + s


def _fmt_p(p):
    return '0' if p == 0 else '1' if p == 1 else f'{p:.2f}'.rstrip('0').lstrip('0')


def _max_run(reads, state):
    best = 0
    for r in reads:
        cur = 0
        for c in r:
            cur = cur + 1 if c == state else 0; best = max(best, cur)
    return best


def derivation(feat, reads, aux):
    if feat == 'AMF':
        return f"{aux['n_meth']}/{aux['n_calls']} calls"
    if feat in ('PDR', 'U') and feat == 'PDR':
        return f"{len(aux['mixed_rows'])}/{aux['n_pdr_reads']} reads mixed"
    if feat in ('U', 'B'):
        return 'p ' + ' '.join(_fmt_p(p) for p in aux['p'])
    if feat == 'A':
        return f"{int(round(sum(p >= .7 for p in aux['p'])))}/{aux['n_K4']} p ≥ 0.7"
    if feat == 'MHL':
        return f'longest 1-run {_max_run(reads, "1")}'
    if feat == 'UHL':
        return f'longest 0-run {_max_run(reads, "0")}'
    if feat == 'J':
        S, E = aux['S'], aux['E']
        if len(set(S)) == 1 and len(set(E)) == 1:
            return f'S {S[0]:g} vs E {E[0]:g}'
        return f'ΣS {sum(S):g} vs ΣE {sum(E):g}'
    if feat in ('PM', 'ME'):
        n = aux['n_patterns'][0]
        return f'{n} pattern' + ('s' if n > 1 else '')
    if feat == 'FDRP':
        return f"{aux['n_disc_pairs']}/{aux['n_pairs']} pairs differ"
    if feat == 'qFDRP':
        return f"{aux['ham_sum']:g}/{aux['n_pairs']} pairs"
    raise KeyError(feat)


