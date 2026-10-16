"""Train 중 strict 통계 보유/결측 표본 구성과 B 효과 비교 (가벼운 표 계산; All|All 층)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import pandas as pd, numpy as np
from scipy import stats
A=CFG.work('group_effects','results')
s=pd.read_csv(A+'/ge_sample_AllAll_values.tsv',sep='\t')
c=pd.read_csv(CFG.inp('marker_matrix_samples'),sep='\t',low_memory=False,usecols=['AnalysisID','source'])
s=s.merge(c,on='AnalysisID')
s['hcc']=(s.group=='liver').astype(int)
rows=[]
for sp,sub in [('Train',s[s.split=='Train']),('Train_strict보유',s[(s.split=='Train')&s.has_strict]),('Train_strict결측',s[(s.split=='Train')&~s.has_strict]),('Test',s[s.split=='Test'])]:
    for m in ['B','R_ratio','U','V','fraglen']:
        a=sub[sub.hcc==1][m].dropna(); b=sub[sub.hcc==0][m].dropna()
        rows.append(dict(subset=sp,metric=m,n_HCC=len(a),n_healthy=len(b),mean_HCC=a.mean(),mean_healthy=b.mean(),delta=a.mean()-b.mean(),p=stats.mannwhitneyu(a,b).pvalue if len(a)>4 and len(b)>4 else np.nan))
pd.DataFrame(rows).to_csv(A+'/ge_Train_strict_subset_check_AllAll.tsv',sep='\t',index=False)
comp=s.groupby(['split','has_strict','group','source']).size().rename('n').reset_index()
comp.to_csv(A+'/ge_sample_composition_by_source_strict.tsv',sep='\t',index=False)
pd.set_option('display.width',200); print(pd.DataFrame(rows).round(5).to_string()); print(comp[comp.split=='Train'].to_string())
