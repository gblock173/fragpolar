import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
from pathlib import Path
import os,json,time,resource,argparse,traceback,hashlib
import pandas as pd,numpy as np,pyarrow as pa,pyarrow.parquet as pq,pyarrow.compute as pc
WOUT=Path(CFG.work('train_rebuild','data'))  # train-rebuild step: only write path changed
COLS=['n_reads','sum_ncpg','sum_switches','sum_expswitch','sum_mean_meth','sum_mean_meth_sq','n_full','n_fullU','sum_fraglen','n_complete','sum_switch_rate']
EXTRA=['n_strictM','n_strictU','n_K5','n_short','n_long','sum_complete_span']

def aggregate(tb,qa):
 if len(tb)==0:return pd.DataFrame(columns=['marker_id']+COLS+EXTRA)
 valid=tb.call.isin([0,1])&tb.pos.notna();qa['invalid_calls']+=int((~valid).sum());tb=tb[valid].copy()
 tb=tb.sort_values(['marker_id','fragment_id','pos'],kind='mergesort')
 du=tb.duplicated(['marker_id','fragment_id','pos'],keep=False)
 if du.any():
  z=tb[du].groupby(['marker_id','fragment_id','pos']).call.nunique();bad=z[z>1].index
  qa['duplicate_rows']+=int(du.sum());qa['conflict_sites']+=len(bad)
  ix=pd.MultiIndex.from_frame(tb[['marker_id','fragment_id','pos']]);tb=tb[~ix.isin(bad)].drop_duplicates(['marker_id','fragment_id','pos'])
 if not len(tb):return pd.DataFrame(columns=['marker_id']+COLS+EXTRA)
 mk,_=pd.factorize(tb.marker_id,sort=False);fr,_=pd.factorize(tb.fragment_id,sort=False);pos=tb.pos.to_numpy();call=tb.call.to_numpy(float)
 new=np.r_[True,(mk[1:]!=mk[:-1])|(fr[1:]!=fr[:-1])];idx=np.flatnonzero(new);rid=np.cumsum(new)-1
 k=np.bincount(rid);m=np.bincount(rid,weights=call);sw=np.bincount(rid,weights=np.r_[False,call[1:]!=call[:-1]]&~new)
 keep=k>=4;k=k[keep];m=m[keep];sw=sw[keep];ix=idx[keep];p=m/k
 marker=tb.marker_id.to_numpy()[ix];g,names=pd.factorize(marker,sort=False);n=len(names)
 span=(tb.fragment_end.to_numpy()-tb.fragment_start.to_numpy())[ix];complete=tb.is_complete_pair.to_numpy()[ix].astype(bool)
 vals=[np.ones(len(k)),k,sw,2*m*(k-m)/k,p,p*p,p>=.7,p<=.3,span,complete,sw/(k-1),m==k,m==0,k>=5,complete&(span>=80)&(span<=150),complete&(span>=151)&(span<=220),span*complete]
 d={'marker_id':names}
 for c,v in zip(COLS+EXTRA,vals):d[c]=np.bincount(g,weights=v,minlength=n)
 qa['eligible_fragment_marker_pairs']+=len(k)
 return pd.DataFrame(d)

def raw_extract(path,regions,limit=None):
 pf=pq.ParquetFile(path);fields=['marker_id','fragment_id','pos','call','fragment_start','fragment_end','is_complete_pair']
 assert set(fields)<=set(pf.schema.names)
 allowed=pa.array(regions.marker_id);parts=[];carry=None;seen=set();qa={x:0 for x in ['invalid_calls','duplicate_rows','conflict_sites','eligible_fragment_marker_pairs']};qa['raw_rows']=pf.metadata.num_rows;qa['row_groups']=pf.metadata.num_row_groups
 for j in range(pf.metadata.num_row_groups if limit is None else min(limit,pf.metadata.num_row_groups)):
  t=pf.read_row_group(j,columns=fields,use_threads=False);t=t.filter(pc.is_in(t['marker_id'],value_set=allowed))
  d=t.to_pandas();del t
  if carry is not None:d=pd.concat([carry,d],ignore_index=True)
  if len(d)==0:carry=None;continue
  last=d.marker_id.iloc[-1];end=d.marker_id.eq(last);carry=d[end].copy();d=d[~end]
  if len(d):
   cur=set(d.marker_id);assert not(cur&seen),'Noncontiguous marker encountered; streaming incomplete';seen.update(cur)
   parts.append(aggregate(d,qa))
  if j%25==0:print(json.dumps({'rowgroup':j,'of':pf.metadata.num_row_groups,'path':Path(path).name}),flush=True)
 if carry is not None and len(carry):
  assert not(set(carry.marker_id)&seen);parts.append(aggregate(carry,qa))
 d=pd.concat(parts,ignore_index=True) if parts else pd.DataFrame(columns=['marker_id']+COLS+EXTRA)
 assert not d.marker_id.duplicated().any();return d,qa

def check_metrics(d):
 n=d.n_reads.to_numpy(float);mu=d.sum_mean_meth.to_numpy()/n;q=d.sum_mean_meth_sq.to_numpy()/n;u=mu-q;b=q-mu*mu
 assert np.nanmin(n)>=1 and np.nanmin(u)>-1e-7 and np.nanmin(b)>-1e-7
 assert np.nanmax(abs(u+b-mu*(1-mu)))<1e-10
 assert np.all(d.n_full<=n)&np.all(d.n_fullU<=n)
 return {'identity_max_error':float(np.nanmax(abs(u+b-mu*(1-mu)))),'zero_arrangement_denominator_markers':int((d.sum_expswitch==0).sum()),'markers':len(d)}

def process(row,force_raw=False,limit=None):
 sid=row.AnalysisID;folder='pilot_raw_checks' if force_raw else 'sufficient';out=WOUT/'suff';out.mkdir(parents=True,exist_ok=True);p=out/(sid+'.parquet');jp=p.with_suffix('.json')
 if p.exists() and jp.exists():return {'AnalysisID':sid,'status':'already_complete'}
 start=time.time();regions=pd.read_parquet(CFG.inp('extraction_regions'))
 if str(row.legacy_path) not in ('','nan') and not force_raw:
  d=pd.read_parquet(row.legacy_path,columns=['marker_id','lenbin']+COLS,filters=[('lenbin','==','ALL')]).drop(columns='lenbin');source=CFG.col('legacy_source_label');qa={}
  for c in EXTRA:d[c]=np.nan
 else:d,qa=raw_extract(row.raw_path,regions,limit);source='raw_CpG_calls'
 qa.update(check_metrics(d));assert d.marker_id.isin(regions.marker_id).all()
 # Store double sums: cancellation-sensitive variance must not lose precision.
 tmp=p.with_suffix('.tmp');d.to_parquet(tmp,index=False,compression='zstd');os.replace(tmp,p)
 rec={'AnalysisID':sid,'status':'complete','source':source,'source_path':row.legacy_path if source.startswith('reuse') else row.raw_path,'source_stat':{'size':Path(row.legacy_path if source.startswith('reuse') else row.raw_path).stat().st_size},'seconds':time.time()-start,'peak_RSS_GiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20,'output_bytes':p.stat().st_size,'QA':qa}
 jp.write_text(json.dumps(rec,indent=2));print(json.dumps(rec),flush=True);return rec

def synthetic():
 qa={x:0 for x in ['invalid_calls','duplicate_rows','conflict_sites','eligible_fragment_marker_pairs']}
 def table(patterns):
  return pd.DataFrame([{'marker_id':'chr1:0-100','fragment_id':str(i),'pos':j*2,'call':int(c),'fragment_start':0,'fragment_end':150,'is_complete_pair':True} for i,s in enumerate(patterns) for j,c in enumerate(s)])
 cases={}
 for name,patterns in {'allM':['11111'],'allU':['00000'],'mixture':['11111','00000'],'internal':['11000','10001'],'alternating':['10101'],'K4':['1100'],'empty':[]}.items():
  d=aggregate(table(patterns),qa) if patterns else pd.DataFrame()
  if len(d):
   check_metrics(d);r=d.iloc[0];n=r.n_reads;mu=r.sum_mean_meth/n;u=(r.sum_mean_meth-r.sum_mean_meth_sq)/n;b=r.sum_mean_meth_sq/n-mu*mu
   score=1-r.sum_switches/r.sum_expswitch if r.sum_expswitch>0 else None
   cases[name]={'patterns':patterns,'mu':mu,'U':u,'B':b,'A_hyper':r.n_full/n,'A_hypo':r.n_fullU/n,'arrangement_score':score}
  else:cases[name]={'patterns':patterns,'status':'NA_no_reads'}
 assert cases['allM']['arrangement_score'] is None and cases['mixture']['B']==.25 and cases['mixture']['U']==0
 assert abs(cases['internal']['arrangement_score']-.375)<1e-12
 t=table(['11000']);t=pd.concat([t,t.iloc[:1]]);d=aggregate(t,qa);assert d.iloc[0].sum_ncpg==5
 t=table(['11000']);t=pd.concat([t,t.iloc[:1].assign(call=0)]);d=aggregate(t,qa);assert d.iloc[0].sum_ncpg==4
 import math
 for k in [4,5,8]:
  b=sum(math.comb(k,m)*.5**k*(m/k-.5)**2 for m in range(k+1));assert abs(b/.25-1/k)<1e-12
 (WOUT/'metric_definition_examples.json').write_text(json.dumps(cases,indent=2));(WOUT/'synthetic_QA.json').write_text(json.dumps({'PASS':True,'cases':len(cases),'duplicate_and_conflict':qa,'zero_denominator':'NA','independent_CpG_R':'1/K'},indent=2))
 print('SYNTHETIC PASS',flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--sid');ap.add_argument('--scope',choices=['Train','Other','all'],default='all');ap.add_argument('--pilot',action='store_true');ap.add_argument('--synthetic',action='store_true');ap.add_argument('--force-raw',action='store_true');ap.add_argument('--limit-rowgroups',type=int);ap.add_argument('--shard',type=int,default=0);ap.add_argument('--shards',type=int,default=1);args=ap.parse_args()
 if args.synthetic:synthetic()
 else:
  m=pd.read_csv(CFG.inp('extraction_pilot_manifest' if args.pilot else 'extraction_manifest'),sep='\t',keep_default_na=False);m=m[m.raw_path.ne('')|m.legacy_path.ne('')]
  if args.scope!='all':m=m[m.Split_AgeSexStage.eq('Train') if args.scope=='Train' else ~m.Split_AgeSexStage.eq('Train')]
  if args.sid:m=m[m.AnalysisID.eq(args.sid)]
  m=m.sort_values(['figure_train','CancerType','AnalysisID'],ascending=[False,True,True]).iloc[args.shard::args.shards]
  failures=0
  for r in m.itertuples():
   try:process(r,args.force_raw,args.limit_rowgroups)
   except Exception as e:
    failures+=1;err=WOUT/'errors';err.mkdir(parents=True,exist_ok=True);(err/(r.AnalysisID+'.txt')).write_text(traceback.format_exc());print(traceback.format_exc(),flush=True)
  if failures:raise SystemExit(f'{failures} samples failed; explicit error files retained')
