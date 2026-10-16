import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import sys; sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from style_v8 import *
V7 = Path(CFG.path('concept_examples_dir'))
# a: constructed sets
pt = pd.read_csv(V7 / 'S0_Figure1_concepts_20261014/constructed_example_patterns.tsv', sep='\t')
names = {'Uniform-state mixture': 'Opposite uniform', 'Within-fragment mixing': 'Within-fragment mixed', 'Combined distribution': 'Combined'}
pt['set'] = pt.example.map(names); pt['binary_pattern'] = pt.binary_pattern.astype(str).str.zfill(4)
rows = []
for s, g in pt.groupby('set', sort=False):
    p = g.binary_pattern.map(lambda x: sum(int(c) for c in x) / len(x)).to_numpy()
    mu = p.mean(); wfv = (p * (1 - p)).mean(); bfv = p.var(); rows.append((s, mu, wfv, bfv, bfv / (mu * (1 - mu))))
summ = pd.DataFrame(rows, columns=['set', 'mu', 'WFV', 'BFV', 'FPI'])
out = pt[['set', 'constructed_fragment_index', 'binary_pattern']].merge(summ, on='set')
out.to_csv(SD / 'Fig1a.tsv', sep='\t', index=False); print(summ)
# c: reversal
fw = {'0001': 10, '0011': 20, '0111': 30, '1111': 40}
pos = []
for d, cnt in (('Original', fw), ('Reversed', {k[::-1]: v for k, v in fw.items()})):
    n = sum(cnt.values())
    for i in range(4): pos.append((d, i + 1, sum(v for k, v in cnt.items() if k[i] == '1') / n))
pd.DataFrame(pos, columns=['order', 'cpg_position', 'methylated_fraction']).to_csv(SD / 'Fig1c_positions.tsv', sep='\t', index=False)
rf = pd.read_csv(V7 / 'S0_Figure1_identifiability_20261014/constructed_reversal_features.tsv', sep='\t')
rf = rf.rename(columns={'direction': 'healthy_reference_direction', 'forward': 'original', 'reversed': 'reversed_order'})
rf.to_csv(SD / 'Fig1c.tsv', sep='\t', index=False); print(rf.absolute_difference.max())
# d: generator response, f = 5 %, Train seed
r = pd.read_csv(AN / 'generators/results/gen_response_TrainSeed20_rebuilt.tsv', sep='\t')
strata = ['All|All', 'U|All', 'mid|All', 'M|All']
r = r[r.condition.isin(['G1_0.05', 'G2_0.05', 'G3_0.05']) & r.stratum.isin(strata)].copy()
r['generator'] = r.condition.str[:2]
disp = {m[0]: m[1] for m in METRICS}; cat = {m[0]: m[2] for m in METRICS}
r['metric'] = r.feature.map(disp); r['category'] = r.feature.map(cat)
r['scaled'] = r['mean'] / r.groupby('feature')['mean'].transform(lambda x: x.abs().max())
r.rename(columns={'mean': 'mean_delta_z', 'count': 'n_seed_samples'})[['stratum', 'generator', 'metric', 'category', 'mean_delta_z', 'std', 'n_seed_samples', 'scaled']].to_csv(SD / 'Fig1d.tsv', sep='\t', index=False)
# b: definitions (text table)
D = [('AMF', 'Average level', 'Mean methylation, μ = E[p]'),
 ('PDR', 'Within-fragment mixing', 'Fraction of mixed fragments (0 < p < 1)'),
 ('WFV', 'Within-fragment mixing', 'E[p(1 − p)]'),
 ('MHL', 'Continuity and arrangement', 'Methylated-run fraction, length-weighted'),
 ('UHL', 'Continuity and arrangement', 'Unmethylated-run fraction, length-weighted'),
 ('J', 'Continuity and arrangement', 'Expected minus observed switches per CpG link'),
 ('HOF', 'Reference deviation', 'Fraction with p ≥ 0.7 or p ≤ 0.3, opposite to healthy'),
 ('PM', 'Pattern diversity', 'Pattern-frequency diversity'),
 ('ME', 'Pattern diversity', 'Pattern-frequency diversity'),
 ('BFV', 'Between-fragment difference', 'Var(p); WFV + BFV = μ(1 − μ)'),
 ('FDRP', 'Between-fragment difference', 'Pairwise pattern disagreement'),
 ('qFDRP', 'Between-fragment difference', 'Sitewise pair disagreement')]
I = [('FPI', 'Interpretive index', 'BFV / μ(1 − μ), range 0–1'),
     ('FPI_K', 'Interpretive index', '(FPI − E[1/K]) / (1 − E[1/K])'),
     ('sOUF', 'Interpretive index', 'Fraction of fragments with p = 1 (p = 0) when healthy μ < 0.5 (> 0.5)')]
# design pass 2026-10-15: Fig1b.tsv is now written by prep_fig1b_cards.py (definition cards with Low/High values)
# e: study design
E_ = [('Healthy cfDNA', 1009, 673, 336), ('HCC cfDNA', 177, 104, 73), ('HCC with tumour-fraction estimate', 46, 28, 18),
      ('Hepatitis cfDNA', 10, None, None), ('Tumour tissue libraries', 26, None, None), ('Adjacent non-tumour tissue libraries', 27, None, None),
      ('Tissue donor pairs', 22, None, None), ('PBMC', 30, None, None), ('Matched HCC-healthy pairs (year, depth)', None, 104, 73)]
# design pass 2026-10-15: study design moved to Supplementary Fig. 1c -> supp/source_data/S1c.tsv
pd.DataFrame(E_, columns=['group', 'n_total', 'n_train', 'n_test']).to_csv(FIG / 'supp/source_data/S1c.tsv', sep='\t', index=False)
