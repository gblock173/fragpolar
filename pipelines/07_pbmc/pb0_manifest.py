"""PBMC step 입력 목록: PBMC 30 library BAM 경로 + 패리티용 혈장 1개. input audit 감사 산출물 읽기 전용."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import pandas as pd
A=CFG.path('pbmc_input_audit_dir')+'/'
O=CFG.work('pbmc')+'/'
m=pd.read_csv(A+'pbmc_library_donor_manifest.tsv',sep='\t')
r=pd.read_csv(A+'raw_input_availability.tsv',sep='\t')
b=r[(r.specimen=='PBMC')&(r.input_type=='raw_BAM')][['sample_id','path','bytes_now']].rename(columns={'sample_id':'pbmc_lab_id','path':'bam','bytes_now':'bam_bytes'})
m=m.merge(b,on='pbmc_lab_id',validate='one_to_one')
assert len(m)==30
m['task']=range(30)
m[['task','pbmc_lab_id','pbmc_library_id','bam','bam_bytes','allowed_healthy_pair','healthy_cfdna_lab_id','original_split','legacy_11_subset']].to_csv(O+'data/pbmc_tasks.tsv',sep='\t',index=False)
print(m[['task','pbmc_lab_id','bam_bytes']].head())
