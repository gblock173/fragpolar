"""Fig 4c rho and 95% CI taken from analysis/tumour_hepatitis_age/results/C1_dose_response.tsv (5,000-sample patient bootstrap), same source as the text."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import pandas as pd
R=CFG.WORK+'/'; FIGD=CFG.FIG+'/'
c=pd.read_csv(R+'tumour_hepatitis_age/results/tf_dose_response.tsv',sep='\t')
c=c[(c.stratum=='All|All')&(c.ref=='median')&c.endpoint.isin(['opp','B'])&c.subset.isin(['All46','Train28','Test18'])]
o=pd.DataFrame({'subset':c.subset.map({'All46':'all','Train28':'Train','Test18':'Test'}),'n':c.n,'rho':c.rho,'ci_lo':c.rho_lo,'ci_hi':c.rho_hi,'y':c.endpoint.map({'opp':'sOUF_diff_pp','B':'BFV_diff_x1e2'}),'source':'tf_dose_response.tsv (bootstrap 5000)'})
o.to_csv(FIGD+'source_data/Fig4c_stats.tsv',sep='\t',index=False,float_format='%.6g'); print(o)
