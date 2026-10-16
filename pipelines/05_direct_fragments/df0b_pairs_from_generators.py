"""direct-fragment step: direct-fragment analysis on the generator step matched pairs (pair set shared with the generator analysis; manuscript main result).

RECONSTRUCTED SCRIPT (the original conversion was not kept as code). Writes, under <work_dir>/direct_fragments/results_A4pairs/:
  B0_pairs.tsv          split, pair_id, hcc, healthy, era, era_match, ld_hcc, ld_healthy  (pairs in the order of
                        generators/results/gen_matched_pairs_rebuilt.tsv; Train_matched -> Train, Test_paired -> Test)
  B0_extra_samples.tsv  sid, task: pair samples that were not extracted for the original direct-fragment step pairs (sorted), for b1/b2 jobs
  B0_union_samples.tsv  sid, task: original B0_samples.tsv order followed by the extra samples (input of df3_stratum.py)
Only split/hcc/healthy are used downstream (df4_stats.py). Then run, with SAMPLES/PKL/PAIRS/OUTD/FIGD set as in
README, "Run order" (df1_extract.py and df2_sample.py on the extra samples; df3_stratum.py, df4_stats.py, df5_draft_figures.py on the union).
"""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import pandas as pd, numpy as np
from pathlib import Path
O = Path(CFG.work('direct_fragments')); OUT = O / 'results_genpairs'; OUT.mkdir(parents=True, exist_ok=True)
a = pd.read_csv(Path(CFG.work('generators', 'results')) / 'gen_matched_pairs_rebuilt.tsv', sep='\t')
s = pd.read_csv(CFG.inp('samples_table'), sep='\t').set_index('AnalysisID')   # era = collection-year column (as in generator step matching)
sp = a.contrast.map({'Train_matched': 'Train', 'Test_paired': 'Test'})
P = pd.DataFrame(dict(split=sp, pair_id=sp + '_' + a.liver, hcc=a.liver, healthy=a.healthy))
P['era'] = s.era.reindex(P.hcc).to_numpy(); P['era_match'] = s.era.reindex(P.healthy).to_numpy() == P.era.to_numpy()
P['ld_hcc'] = s.log_depth.reindex(P.hcc).to_numpy(); P['ld_healthy'] = s.log_depth.reindex(P.healthy).to_numpy()
P.to_csv(OUT / 'df_pairs.tsv', sep='\t', index=False)
old = pd.read_csv(O / 'results' / 'df_samples.tsv', sep='\t').sid.tolist()
extra = sorted((set(P.hcc) | set(P.healthy)) - set(old))
pd.DataFrame(dict(sid=extra, task=np.arange(len(extra)) + 1)).to_csv(OUT / 'df_extra_samples.tsv', sep='\t', index=False)
u = old + extra
pd.DataFrame(dict(sid=u, task=np.arange(len(u)) + 1)).to_csv(OUT / 'df_union_samples.tsv', sep='\t', index=False)
print(P.split.value_counts().to_dict(), 'era match', P.groupby('split').era_match.mean().round(3).to_dict(), 'extra', len(extra), 'union', len(u))
