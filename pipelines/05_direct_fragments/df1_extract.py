"""direct-fragment step: raw cpg_calls -> per-sample histogram of K>=4 fragments by (marker, read-split half, K, M).
Cleaning follows fragpolar/engine.py (one-pass engine) (engine rule that produced the existing 59-stratum numbers): valid call 0/1, pos not NA, marker in 190,218 primary set,
(marker,fragment,pos) conflicting calls dropped, identical duplicates kept once, fragment = (marker, fragment_id), K>=4.
Half = top bit of pandas siphash(hash_key='P03readsplitv1__') of 'sid|P03_readsplit_v1|fragment_id' (marker/state not hashed)."""
import sys as _sys, pathlib as _pl; _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[2] / 'src'))  # code_release/src
from fragpolar import config as CFG  # noqa: E402
import os,sys,time,json
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np,pandas as pd,pyarrow.parquet as pq,resource
from pathlib import Path
R=Path(CFG.DATA);O=Path(CFG.work('direct_fragments'))
RAW=Path(CFG.CPG_CALLS)
SALT='P03_readsplit_v1'
def run(sid,limit=None):
 t=time.time();out=O/'work/hist';out.mkdir(parents=True,exist_ok=True)
 st=pd.read_parquet(CFG.inp('marker_annotation'),columns=['marker_id']);lookup=pd.Index(st.marker_id)
 pf=pq.ParquetFile(RAW/(sid+'.cpg_calls.parquet'));cols=['marker_id','fragment_id','pos','call']
 qa=dict(rows=0,rows_primary=0,dup_rows=0,conflict_sites=0,fragments_k4=0);parts=[];carry=None;seen=set()
 def process(d):
  qa['rows']+=len(d);d=d[d.call.isin([0,1])&d.pos.notna()]
  mk=lookup.get_indexer(d.marker_id);ok=mk>=0;d=d[ok];mk=mk[ok];qa['rows_primary']+=len(d)
  if not len(d):return
  fr,uq=pd.factorize(d.fragment_id);hb=((pd.util.hash_pandas_object(pd.Series(uq).radd(sid+'|'+SALT+'|'),index=False,hash_key='P03readsplitv1__').to_numpy(np.uint64)>>np.uint64(63))&np.uint64(1)).astype(np.int8)
  pos=d.pos.to_numpy(np.int64);o=np.lexsort((pos,fr,mk));mk=mk[o];fr=fr[o];pos=pos[o];call=d.call.to_numpy(np.int8)[o]
  same=np.r_[False,(mk[1:]==mk[:-1])&(fr[1:]==fr[:-1])&(pos[1:]==pos[:-1])]
  if same.any():
   grp=np.cumsum(~same)-1;gn=np.bincount(grp);gm=np.bincount(grp,weights=call);bad=(gm>0)&(gm<gn)
   qa['dup_rows']+=int(gn[gn>1].sum());qa['conflict_sites']+=int(bad.sum());keep=~same&~bad[grp];mk=mk[keep];fr=fr[keep];call=call[keep]
  cur=set(np.unique(mk).tolist());assert not(seen&cur),'noncontiguous marker';seen.update(cur)
  new=np.r_[True,(mk[1:]!=mk[:-1])|(fr[1:]!=fr[:-1])];rid=np.cumsum(new)-1;starts=np.flatnonzero(new)
  k=np.bincount(rid);m=np.bincount(rid,weights=call).astype(np.int64);pick=k>=4
  k=k[pick];m=m[pick];mm=mk[starts][pick];h=hb[fr[starts][pick]];qa['fragments_k4']+=len(k)
  assert k.max()<256
  key=((mm.astype(np.int64)*2+h)*256+k)*256+m
  u,c=np.unique(key,return_counts=True);parts.append(pd.DataFrame(dict(key=u,cnt=c.astype(np.int32))))
 n=pf.num_row_groups if limit is None else min(limit,pf.num_row_groups)
 for j in range(n):
  d=pf.read_row_group(j,columns=cols,use_threads=False).to_pandas()
  if carry is not None:d=pd.concat([carry,d],ignore_index=True)
  last=d.marker_id.iloc[-1];sel=d.marker_id.eq(last);carry=d[sel];d=d[~sel]
  if len(d):process(d)
  if len(parts)>40:parts=[pd.concat(parts).groupby('key',as_index=False).cnt.sum()]
  if j%20==0:print(sid,'RG',j,'/',n,round(time.time()-t),'s rss',round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20,2),'GB',flush=True)
 if carry is not None and len(carry):process(carry)
 H=pd.concat(parts).groupby('key',as_index=False).cnt.sum();key=H.key.to_numpy()
 H=pd.DataFrame(dict(marker=(key>>17).astype(np.int32),half=((key>>16)&1).astype(np.int8),k=((key>>8)&255).astype(np.int16),m=(key&255).astype(np.int16),cnt=H.cnt.to_numpy()))
 H.to_parquet(out/(sid+('' if limit is None else f'.limit{limit}')+'.parquet'),index=False,compression='zstd')
 qa.update(sid=sid,hist_rows=len(H),seconds=time.time()-t,peak_rss_gb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20,limit=limit)
 (out/(sid+'.QA.json')).write_text(json.dumps(qa,indent=1));print(json.dumps(qa),flush=True)
if __name__=='__main__':
 sid=sys.argv[1] if not sys.argv[1].isdigit() else pd.read_csv(Path(os.environ.get('SAMPLES',O/'results/df_samples.tsv')),sep='\t').set_index('task').sid[int(sys.argv[1])]
 run(sid,int(sys.argv[2]) if len(sys.argv)>2 else None)
