import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import pandas as pd
c=pd.read_csv(CFG.path('sample_provenance_tsv'),sep='\t')
ms=pd.read_csv(CFG.inp('samples_table'),sep='\t');ms=ms[ms.group.isin(['healthy','liver'])&(ms.split=='Train')]
c=c[c.sample_id.isin(ms.AnalysisID)]
P=CFG.work('train_rebuild')
leg=c[c.source==CFG.col('legacy_source_label')].sample_id.tolist(); raw=c[c.source=='raw_CpG_calls'].sample_id.tolist()
open(P+'/data/legacy422_ids.txt','w').write('\n'.join(leg)+'\n'); open(P+'/data/raw355_ids.txt','w').write('\n'.join(raw)+'\n')
print(len(leg),len(raw))
