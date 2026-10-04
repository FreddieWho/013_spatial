#!/usr/bin/env python3
"""Complete native-field orientation references; significance is separately gated."""
from pathlib import Path
import argparse,importlib.util,json,time
import numpy as np
from scipy import sparse,spatial
from scipy.sparse import csgraph

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'infra/tls_multiscale_gobp_20261002'
def mod(name):
 s=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
M=mod('tls_multiscale_gobp');F=mod('tls_multiscale_forward');ANGLES=np.arange(4.5,180,4.5)


def residual_parts(weight,Z,visible,centre):
    training=visible&~centre
    if training.sum()<max(50,2*Z.shape[1]):return None,None
    X=Z[training];coef=np.linalg.solve(X.T@X+np.eye(X.shape[1])*1e-6,Z.T@weight)
    local=np.zeros(len(weight));local[centre]=Z[centre]@coef
    return coef,local


def composition(W,labels):
    # Metadata sensitivity only; no expression signal enters reference qualification.
    positive=np.maximum(W,0);negative=-np.minimum(W,0)
    return np.c_[positive@(labels=='TUM'),negative@(labels=='TUM'),
                 positive@np.isin(labels,['','UNASSIGNED']),negative@np.isin(labels,['','UNASSIGNED'])]


def reference_geometry(xy,labels,points,length,visible,width,true_counts,true_meta,span,angle):
    mid=(points[0]+points[-1])/2;t=np.deg2rad(angle);R=np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]])
    rotated=(points-mid)@R.T+mid
    along,side,_=F.path_coordinates(xy,rotated[0],rotated[-1],polyline=rotated)
    result,status,counts=F.corridor_weights(along,side,length,visible,width,span)
    if status!='ELIGIBLE':return None,1,False
    counts=np.array(counts);ratio=counts/np.asarray(true_counts)
    if np.any((ratio<.5)|(ratio>2)):return None,2,False
    W,V,centre=result
    matched=bool(np.all(np.abs(composition(W,labels)-true_meta)<=.2))
    return (W,centre),0,matched


def build(sid):
    directory=OUT/'controls'/sid;directory.mkdir(parents=True,exist_ok=True)
    if (directory/'geometry_complete.json').exists():return
    source=OUT/'forward'/sid;complete=json.loads((source/'operator_complete.json').read_text());n=complete['evaluable_geometries']
    if not n:
        M.js(directory/'geometry_complete.json',dict(status='NO_EVALUABLE_GEOMETRY',rows=0,matched_references=0));return
    rows=M.tab(source/'effect_rows.tsv');pairs=M.tab(source/'pairs.tsv');pairmap={(r['pair_id'],r['mode']):r for r in pairs}
    stored_segments=sparse.load_npz(source/'segments.npz')
    endpoints=M.tab(OUT/'geometry'/f'{sid}_endpoints.tsv')
    ends={int(e['endpoint_id']):e for e in endpoints}
    g=np.load(OUT/'geometry'/f'{sid}.npz');xy=g['xy'];labels=g['labels'];graph=sparse.load_npz(OUT/'geometry'/f'{sid}_graph.npz')
    dtls=spatial.cKDTree(xy[labels=='TLS']).query(xy)[0]
    starts=sorted({int(ends[int(p['source_endpoint'])]['spot_index']) for p in pairs})
    _,pred=csgraph.dijkstra(graph,directed=False,indices=starts,return_predecessors=True);predmap={s:pred[i] for i,s in enumerate(starts)}
    # Bounded chunks avoid keeping hundreds of thousands of rows in RAM.
    group_ids={};group_meta=[];pending_w=[];pending_local=[];pending_coef=[];pending_meta=[];chunk=0;all_meta=[]
    flags=np.full((n,len(ANGLES)),255,np.uint8);tissue=np.zeros_like(flags,bool)
    last_key=None
    def flush():
        nonlocal pending_w,pending_local,pending_coef,pending_meta,chunk
        if not pending_meta:return
        sparse.save_npz(directory/f'chunk_{chunk:05}_segments.npz',sparse.vstack(pending_w,format='csr'))
        sparse.save_npz(directory/f'chunk_{chunk:05}_local.npz',sparse.vstack(pending_local,format='csr'))
        np.save(directory/f'chunk_{chunk:05}_coef.npy',np.array(pending_coef))
        M.tsv(directory/f'chunk_{chunk:05}_rows.tsv',pending_meta)
        all_meta.extend(pending_meta);chunk+=1;pending_w=[];pending_local=[];pending_coef=[];pending_meta=[]
    for row in rows:
        i=int(row['effect_row']);rad=int(row['mask_um']);width=float(row['width_um']);pair=pairmap[(row['pair_id'],row['mode'])]
        a=ends[int(pair['source_endpoint'])];b=ends[int(pair['target_endpoint'])];start=int(a['spot_index']);end=int(b['spot_index'])
        key=(row['pair_id'],rad)
        if key!=last_key:
            srcids=g['endpoint_members'][g['endpoint_indptr'][int(a['endpoint_id'])]:g['endpoint_indptr'][int(a['endpoint_id'])+1]]
            d1=spatial.cKDTree(xy[srcids]).query(xy)[0];d2=np.linalg.norm(xy-xy[end],axis=1);visible=(dtls>rad)&(d2>rad)
            Z=np.column_stack([np.ones(len(xy))]+[np.exp(-d1/l) for l in [100,250,500,1000,2000,4000]]+[np.exp(-d2/l) for l in [100,250,500,1000,2000,4000]])
            last_key=key
        if key not in group_ids:
            group_ids[key]=len(group_ids);group_meta.append(dict(group_id=group_ids[key],pair_id=row['pair_id'],mask_um=rad,source_endpoint=a['endpoint_id'],target_endpoint=b['endpoint_id']))
        group=group_ids[key]
        if row['mode']=='straight':points=xy[[start,end]]
        else:points=xy[F.trace_path(predmap[start],start,end)]
        along,side,length=F.path_coordinates(xy,xy[start],xy[end],polyline=points if row['mode']!='straight' else None)
        true,status,_=F.corridor_weights(along,side,length,visible,width)
        if status!='ELIGIBLE':raise RuntimeError('stored geometry no longer reproduces: '+sid+' '+str(i))
        W,_,centre=true
        np.testing.assert_allclose(W,stored_segments[i*5:(i+1)*5].toarray(),rtol=0,atol=1e-14,
                                   err_msg='reference geometry differs from scored forward operator')
        span=(float(along[centre].min()),float(along[centre].max()));meta=composition(W,labels);counts=json.loads(row['segment_counts'])
        for j,angle in enumerate(ANGLES):
            rr,status,matched=reference_geometry(xy,labels,points,length,visible,width,counts,meta,span,angle)
            flags[i,j]=status;tissue[i,j]=matched
            if status:continue
            WW,cc=rr;coef,local=residual_parts(WW.mean(axis=0),Z,visible,cc)
            if coef is None:flags[i,j]=3;continue
            rid=len(all_meta)+len(pending_meta)
            pending_meta.append(dict(reference_id=rid,effect_row=i,angle_index=j,angle_degrees=float(angle),group_id=group,tissue_matched=matched))
            pending_w.extend([sparse.csr_matrix(w) for w in WW]);pending_local.append(sparse.csr_matrix(local));pending_coef.append(coef)
            if len(pending_meta)>=256:flush()
        if i%500==0:print('references',sid,i,'/',n,'accepted',len(all_meta)+len(pending_meta),flush=True)
    flush();M.tsv(directory/'groups.tsv',group_meta)
    np.savez_compressed(directory/'qualification.npz',status=flags,tissue_matched=tissue,angles=ANGLES)
    M.js(directory/'geometry_complete.json',dict(status='COMPLETE',rows=n,angles=len(ANGLES),attempted_references=n*len(ANGLES),matched_references=len(all_meta),chunks=chunk,groups=len(group_ids),codes={'0':'ELIGIBLE','1':'INSUFFICIENT_SUPPORT','2':'COUNT_MISMATCH','3':'ENDPOINT_BASELINE_UNAVAILABLE','255':'UNPROCESSED'}))
    print('references complete',sid,n,len(all_meta),flush=True)


def source_groups(sid,cache):
    directory=OUT/'controls'/sid;g=np.load(OUT/'geometry'/f'{sid}.npz');xy=g['xy'];labels=g['labels']
    endpoints={int(r['endpoint_id']):r for r in M.tab(OUT/'geometry'/f'{sid}_endpoints.tsv')}
    dtls=spatial.cKDTree(xy[labels=='TLS']).query(xy)[0];S=cache
    T=[];SD=[]
    for row in M.tab(directory/'groups.tsv'):
        a=endpoints[int(row['source_endpoint'])];b=endpoints[int(row['target_endpoint'])];rad=int(row['mask_um']);target=int(b['spot_index'])
        ids=g['endpoint_members'][g['endpoint_indptr'][int(a['endpoint_id'])]:g['endpoint_indptr'][int(a['endpoint_id'])+1]]
        d1=spatial.cKDTree(xy[ids]).query(xy)[0];d2=np.linalg.norm(xy-xy[target],axis=1);vis=(dtls>rad)&(d2>rad)
        Z=np.column_stack([np.ones(len(xy))]+[np.exp(-d1/l) for l in [100,250,500,1000,2000,4000]]+[np.exp(-d2/l) for l in [100,250,500,1000,2000,4000]])
        T.append(Z[vis].T@S[vis]);SD.append(S[vis].std(axis=0,ddof=1))
    return np.array(T),np.array(SD)


def score(sid):
    directory=OUT/'controls'/sid;meta=json.loads((directory/'geometry_complete.json').read_text())
    if (directory/'scoring_complete.json').exists():return
    n=meta['matched_references']
    if not n:M.js(directory/'scoring_complete.json',dict(status='NO_REFERENCES',n_references=0));return
    cache=np.load(OUT/'scores'/f'{sid}.npz');S=np.nan_to_num(cache['scores']).astype(float);valid=cache['valid'];P=S.shape[1]
    T,SD=source_groups(sid,S)
    shape=(n,P)
    arrays={k:np.lib.format.open_memmap(directory/f'{k}.npy',mode='w+',dtype='float32' if k in ['raw','adjusted'] else 'uint8',shape=shape) for k in ['raw','adjusted','positive_segments','negative_segments']}
    for ch in range(meta['chunks']):
        rows=M.tab(directory/f'chunk_{ch:05}_rows.tsv');idx=np.array([int(r['reference_id']) for r in rows]);groups=np.array([int(r['group_id']) for r in rows])
        W=sparse.load_npz(directory/f'chunk_{ch:05}_segments.npz');L=sparse.load_npz(directory/f'chunk_{ch:05}_local.npz');C=np.load(directory/f'chunk_{ch:05}_coef.npy')
        R=np.asarray(W@S).reshape(len(rows),5,P);mean=R.mean(axis=1);adjust=mean+L@S-np.einsum('nd,ndp->np',C,T[groups],optimize=True)
        ok=valid[None,:]&(SD[groups]>1e-12);scale=np.where(ok,SD[groups],np.nan)
        arrays['raw'][idx]=mean/scale;arrays['adjusted'][idx]=adjust/scale
        pos=(R>0).sum(axis=1);neg=(R<0).sum(axis=1);pos[~ok]=0;neg[~ok]=0
        arrays['positive_segments'][idx]=pos;arrays['negative_segments'][idx]=neg
        if ch%20==0:print('reference scoring',sid,ch,'/',meta['chunks'],flush=True)
    for a in arrays.values():a.flush()
    M.js(directory/'scoring_complete.json',dict(status='COMPLETE_NATIVE_REFERENCE_NOT_YET_CALIBRATED',n_references=n,n_gobp=1691,n_controls=51))
    print('reference scoring complete',sid,n,flush=True)

def run_source(task):
    stage,sid=task
    if stage in ['build','all']:build(sid)
    if stage in ['score','all']:score(sid)
    return sid


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['build','score','all']);ap.add_argument('--section',default='');ap.add_argument('--workers',type=int,default=1);a=ap.parse_args()
    tasks=[(a.stage,info['section_id']) for info in json.loads((OUT/'contract.json').read_text())['sources'] if not a.section or a.section==info['section_id']]
    if not tasks:raise ValueError('No matching section')
    if a.workers==1:
        for task in tasks:run_source(task)
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            for sid in pool.map(run_source,tasks):print('source finished',a.stage,sid,flush=True)
