"""PAT 처리 목록. 간세포 6, 혈구(과립구, 단핵구, B, NK, T-CD4, T-CD8, T-CD3), 검증용 CNVS-NORM cfDNA/WBC 짝 8쌍."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import glob,os,re,collections
L=CFG.path('pat_dir')
E=CFG.work('component_origin')
want=['Liver-Hepatocytes','Blood-Granulocytes','Blood-Monocytes','Blood-B','Blood-NK','Blood-T-CD4','Blood-T-CD8','Blood-T-CD3']
rows=[]
for f in sorted(glob.glob(L+'/GSM*.hg38.pat.gz')):
    b=os.path.basename(f)
    if 'CNVS-NORM' in b:
        m=re.match(r'GSM\d+_CNVS-NORM-(cfDNA|WBC)-(\d+)-',b); rows.append((b.replace('.hg38.pat.gz',''),'CNVS_'+m.group(1),f)); continue
    m=re.match(r'(GSM\d+)_(.+)-(Z[0-9A-Za-z]+)\.hg38\.pat\.gz',b); gsm,ct,z=m.groups()
    if ct in want: rows.append((gsm+'_'+ct+'-'+z,ct,f))
don=lambda r:re.search(r'-(\d{9})-',r[0]).group(1)
cf=[r for r in rows if r[1]=='CNVS_cfDNA']; wb=[r for r in rows if r[1]=='CNVS_WBC']
common=sorted(set(map(don,cf))&set(map(don,wb)))[:8]
rows=[r for r in rows if not r[1].startswith('CNVS') or don(r) in common]
with open(E+'/data/pat_tasks.tsv','w') as o:
    o.write('task\tid\tcelltype\tpath\n')
    for i,r in enumerate(rows,1): o.write(f'{i}\t{r[0]}\t{r[1]}\t{r[2]}\n')
print(len(rows), collections.Counter(r[1] for r in rows))
