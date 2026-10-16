"""CpG-call cleaning and per-marker sufficient statistics (K >= 4 fragments).

Verbatim copy of aggregate(), raw_extract() and check_metrics() of the upstream extraction module (extract.py; see
docs/CODE_STRUCTURE_KO.md, internal section, for the source path). Only the module-level output path was removed;
the per-sample driver (process / __main__) lives in pipelines/02_train_rebuild/tr1_extract.py.

Rules: valid calls 0/1 with a position; (marker, fragment, position) duplicates with identical calls kept once,
conflicting calls removed; fragment = (marker, fragment_id); only fragments with K >= 4 CpG calls are aggregated.
"""
from pathlib import Path
import os,json,time,resource,argparse,traceback,hashlib
import pandas as pd,numpy as np,pyarrow as pa,pyarrow.parquet as pq,pyarrow.compute as pc
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



def synthetic_checks():
    """Assertions of the upstream synthetic() routine without writing files; returns the case table."""
    qa = {x: 0 for x in ['invalid_calls', 'duplicate_rows', 'conflict_sites', 'eligible_fragment_marker_pairs']}
    def table(patterns):
        return pd.DataFrame([{'marker_id': 'chr1:0-100', 'fragment_id': str(i), 'pos': j * 2, 'call': int(c), 'fragment_start': 0, 'fragment_end': 150,
                              'is_complete_pair': True} for i, s in enumerate(patterns) for j, c in enumerate(s)])
    cases = {}
    for name, patterns in {'allM': ['11111'], 'allU': ['00000'], 'mixture': ['11111', '00000'], 'internal': ['11000', '10001'], 'alternating': ['10101'],
                           'K4': ['1100'], 'empty': []}.items():
        d = aggregate(table(patterns), qa) if patterns else pd.DataFrame()
        if len(d):
            check_metrics(d); r = d.iloc[0]; n = r.n_reads; mu = r.sum_mean_meth / n; u = (r.sum_mean_meth - r.sum_mean_meth_sq) / n; b = r.sum_mean_meth_sq / n - mu * mu
            score = 1 - r.sum_switches / r.sum_expswitch if r.sum_expswitch > 0 else None
            cases[name] = {'patterns': patterns, 'mu': mu, 'U': u, 'B': b, 'A_hyper': r.n_full / n, 'A_hypo': r.n_fullU / n, 'arrangement_score': score}
        else:
            cases[name] = {'patterns': patterns, 'status': 'NA_no_reads'}
    assert cases['allM']['arrangement_score'] is None and cases['mixture']['B'] == .25 and cases['mixture']['U'] == 0
    assert abs(cases['internal']['arrangement_score'] - .375) < 1e-12
    t = table(['11000']); t = pd.concat([t, t.iloc[:1]]); d = aggregate(t, qa); assert d.iloc[0].sum_ncpg == 5
    t = table(['11000']); t = pd.concat([t, t.iloc[:1].assign(call=0)]); d = aggregate(t, qa); assert d.iloc[0].sum_ncpg == 4
    import math
    for k in [4, 5, 8]:
        b = sum(math.comb(k, m) * .5 ** k * (m / k - .5) ** 2 for m in range(k + 1)); assert abs(b / .25 - 1 / k) < 1e-12
    return cases
