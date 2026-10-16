"""Fig. 1b definition cards: Low/High example values recomputed from the definitions.

12 summaries: compute_features() of fragpolar/feature_definitions.py (independent implementation of the formula definitions, gate-free)
on the example fragments, cross-checked (optional) with the stored gallery value table (config paths.gallery_values_tsv).
Interpretive indices (00_TERMS_KO.md), computed here on K >= 4 fragments with equal weights:
  FPI   = BFV / mu(1 - mu)
  FPI_K = (FPI - E[1/K]) / (1 - E[1/K])
  sOUF  = fraction of fragments fully opposite to healthy (healthy unmethylated: p = 1)
HOF / sOUF use the stated healthy reference 'unmethylated' (healthy mu < 0.5), as in earlier.
Writes source_data/Fig1b.tsv. Runs in < 1 s.
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys
from pathlib import Path
import numpy as np, pandas as pd
from fragpolar import feature_definitions as GAL      # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parent))
from style_v8 import SD, METRIC_CAT, INTERP            # noqa: E402

INTERNAL = {'WFV': 'U', 'BFV': 'B', 'HOF': 'A'}         # display -> internal symbol (gallery / production tables)
CARDS = [  # display, full name, definition, category
    ('AMF', 'Average methylation fraction', 'Fraction of observed CpG calls that are methylated'),
    ('PDR', 'Proportion of discordant fragments', 'Fraction of fragments (≥4 CpGs) carrying both methylation states'),
    ('WFV', 'Within-fragment variance', 'Mean within-fragment mixing, E[p(1\u00a0−\u00a0p)]'),
    ('PM', 'Epipolymorphism', 'Diversity of complete 4-CpG patterns, 1 − Σf²'),
    ('ME', 'Methylation entropy', 'Diversity and evenness of 4-CpG patterns, −¼ Σ f·log₂f'),
    ('MHL', 'Methylation haplotype load', 'Length-weighted load of fully methylated runs'),
    ('UHL', 'Unmethylated haplotype load', 'Length-weighted load of fully unmethylated runs'),
    ('J', 'Arrangement coherence', 'Switch deficit vs random CpG order, (ΣE − ΣS)/Σ(K − 1)'),
    ('HOF', 'Healthy-opposite fragment fraction', 'Near-uniform fragments opposite to healthy (here p\u00a0≥\u00a00.7)'),
    ('sOUF', 'Strict opposite-uniform fraction', 'Fragments fully opposite to healthy (here p\u00a0=\u00a01)'),
    ('BFV', 'Between-fragment variance', 'Spread of fragment methylation p across fragments, Var(p)'),
    ('FDRP', 'Fraction of discordant fragment pairs', 'Fraction of fragment pairs that differ at any shared CpG'),
    ('qFDRP', 'Quantitative FDRP', 'Mean fraction of shared CpGs that differ within fragment pairs'),
    ('FPI', 'Fragment polarisation index', 'Share of total variance between fragments, BFV/μ(1 − μ)'),
    ('FPI_K', 'K-corrected FPI', '(FPI − E[1/K])/(1 − E[1/K]); K = 4 here, so E[1/K] = 0.25'),
]
EXTRA = {  # interpretive cards: (Low, High, highlight)
    'FPI': (['0011', '0101', '1010', '1100'], ['0000', '1111', '0000', '1111'], None),
    'FPI_K': (['0100', '0010', '1101', '1011'], ['0000', '0000', '1111', '1111'], None),
    'sOUF': (['0000', '0001', '0000', '0100'], ['0000', '1111', '0111', '0000'], 'p1'),
}


def interp(reads):
    X = GAL._matrix(reads); obs = X >= 0; K = obs.sum(1); r4 = K >= 4
    p = np.array([(X[i] == 1).sum() / K[i] for i in np.flatnonzero(r4)], float); Kr = K[r4].astype(float)
    mu = p.mean(); bfv = p.var(); tot = mu * (1 - mu)
    fpi = bfv / tot if tot > 0 else np.nan; e1k = float(np.mean(1 / Kr))
    return dict(FPI=fpi, FPI_K=(fpi - e1k) / (1 - e1k), sOUF=float(np.mean(p == 1)), mu=mu, BFV=bfv, tot=tot, e1k=e1k, p=p)


def fmt(x, signed=False):
    return GAL.fmt_value(x, signed=signed)


rows = []
_gv = CFG.path('gallery_values_tsv', required=False)   # optional cross-check table
GV = pd.read_csv(_gv, sep='\t') if _gv else None
for disp, name, dfn in CARDS:
    cat = METRIC_CAT.get(disp, INTERP)
    if disp in EXTRA:
        lo, hi, hl = EXTRA[disp]
    else:
        lo, hi, hl = GAL.EXAMPLES[INTERNAL.get(disp, disp)]
    for ex, reads in (('Low', lo), ('High', hi)):
        if disp in EXTRA:
            d = interp(reads); v = d[disp]
            if disp == 'FPI':
                der = f"BFV {fmt(d['BFV'])} / {fmt(d['tot'])}"
            elif disp == 'FPI_K':
                der = f"FPI {fmt(d['FPI'])}"
            else:
                der = f"{int(round((d['p'] == 1).sum()))}/{len(d['p'])} p = 1"
            vt = fmt(v)
        else:
            key = INTERNAL.get(disp, disp)
            vals, aux = GAL.compute_features(reads, healthy_unmethylated=True); v = vals[key]
            vt = fmt(v, signed=key == 'J'); der = GAL.derivation(key, reads, aux)
            for a_, b_ in ((' reads mixed', ' mixed'), (' pairs differ', ' pairs differ'), ('longest 1-run', 'max 1-run'), ('longest 0-run', 'max 0-run')):
                der = der.replace(a_, b_)
            if GV is not None:
                t = GV[GV.feature.eq(key) & GV.example.eq(ex)]
                assert len(t) == 1 and abs(float(t.value.iloc[0]) - v) < 1e-12 and t.reads.iloc[0] == ' '.join(reads), (disp, ex)
        rows.append(dict(card=disp, full_name=name, category=cat, definition=dfn, example=ex, fragments=' '.join(reads),
                         value=float(v), value_text=vt, derivation=der,
                         highlight=str(hl) if hl is not None else ''))
out = pd.DataFrame(rows)
# sanity: HOF and sOUF on the same High fragments differ (0.5 vs 0.25); FPI_K(Low) = 0 for the binomial-like set
assert abs(out.query("card=='HOF' & example=='High'").value.iloc[0] - .5) < 1e-12
assert abs(out.query("card=='sOUF' & example=='High'").value.iloc[0] - .25) < 1e-12
assert abs(out.query("card=='FPI_K' & example=='Low'").value.iloc[0]) < 1e-12
out.to_csv(SD / 'Fig1b.tsv', sep='\t', index=False, float_format='%.6g')
print(out[['card', 'example', 'fragments', 'value_text', 'derivation']].to_string())
print('12 summaries == stored gallery values (24 examples): OK' if GV is not None else 'gallery value cross-check skipped (paths.gallery_values_tsv not set)')
