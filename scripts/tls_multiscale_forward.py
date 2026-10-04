#!/usr/bin/env python3
"""Full 1691-panel forward operators; no formal p-values or geometry-selected claims."""
from pathlib import Path
import argparse,csv,importlib.util,json,time
import numpy as np
from scipy import spatial,sparse
from scipy.sparse import csgraph

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'infra/tls_multiscale_gobp_20261002'
def mod(name):
 spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
M=mod('tls_multiscale_gobp')


def path_coordinates(xy,start,end,path=None,polyline=None):
    if path is None and polyline is None:
        delta=end-start;length=float(np.linalg.norm(delta))
        if length==0:return np.full(len(xy),np.nan),np.full(len(xy),np.nan),0.
        unit=delta/length;rel=xy-start
        return rel@unit,rel[:,0]*unit[1]-rel[:,1]*unit[0],length
    points=xy[path] if polyline is None else polyline;steps=np.linalg.norm(np.diff(points,axis=0),axis=1);arc=np.r_[0,np.cumsum(steps)]
    _,nearest=spatial.cKDTree(points).query(xy)
    # Project to either segment adjacent to the closest vertex, retaining exact segment distances.
    best=np.full(len(xy),np.inf);along=np.zeros(len(xy));side=np.zeros(len(xy))
    for offset in [-1,0]:
        k=np.clip(nearest+offset,0,len(points)-2);a=points[k];b=points[k+1];v=b-a;length=steps[k]
        u=np.clip(np.sum((xy-a)*v,axis=1)/(length*length),0,1)
        foot=a+u[:,None]*v;dist=np.linalg.norm(xy-foot,axis=1);cross=(xy[:,0]-a[:,0])*v[:,1]-(xy[:,1]-a[:,1])*v[:,0]
        use=dist<best;along[use]=(arc[k]+u*length)[use];side[use]=(dist*np.where(cross>=0,1,-1))[use];best[use]=dist[use]
    return along,side,float(arc[-1])


def corridor_weights(along,side,length,visible,width,span=None):
    central=visible&(along>0)&(along<length)&(np.abs(side)<=width)
    if central.sum()<25:return None,'CENTER_LT25',None
    lo=float(along[central].min());hi=float(along[central].max())
    if span is not None:
        lo,hi=span
        central &= (along>=lo)&(along<=hi)
    if hi-lo<100:return None,'RESIDUAL_SPAN_LT_PITCH',None
    seg=np.minimum(np.floor((along-lo)/(hi-lo)*5).astype(int),4)
    within=(along>=lo)&(along<=hi)&visible
    rows=[];counts=[];side_rows=[]
    for k in range(5):
        a=within&(seg==k)&(np.abs(side)<=width)
        b=within&(seg==k)&(side>width)&(side<=3*width)
        c=within&(seg==k)&(side<-width)&(side>=-3*width)
        ns=(int(a.sum()),int(b.sum()),int(c.sum()));counts.append(ns)
        if min(ns)<5:return None,'SEGMENT_OR_SIDEBAND_LT5',counts
        rows.append(a/ns[0]-.5*b/ns[1]-.5*c/ns[2]);side_rows.append(b/ns[1]-c/ns[2])
    return (np.array(rows),np.array(side_rows),central),'ELIGIBLE',counts


def endpoint_adjust(weight,design,training):
    if training.sum()<max(50,2*design.shape[1]):return None
    Z=design[training];rhs=design.T@weight
    coeff=np.linalg.lstsq(Z.T@Z+np.eye(Z.shape[1])*1e-6,rhs,rcond=None)[0]
    adjusted=weight.copy();adjusted[training]-=Z@coeff
    return adjusted


def trace_path(pred,start,end):
    ids=[end];node=end
    while node!=start:
        node=int(pred[node])
        if node<0:return None
        ids.append(node)
        if len(ids)>len(pred):raise RuntimeError('path predecessor cycle')
    return np.array(ids[::-1])


def scale_id(length):
    return int(np.searchsorted([500,1000,2000,4000,8000],length,side='right'))


def pairs_for(endpoints):
    tls=[e for e in endpoints if e['family']=='TLS'];pairs=[]
    for i,a in enumerate(tls):
        for b in tls[i+1:]:pairs.append(('TLS_TLS',a,b))
        for b in endpoints:
            if b['family'].startswith('TUMOR'):
                pairs.append(('TLS_EDGE' if b['family']=='TUMOR_EDGE' else 'TLS_REGION',a,b))
    return pairs


def operators(sid):
    directory=OUT/'forward'/sid;directory.mkdir(parents=True,exist_ok=True)
    if (directory/'operator_complete.json').exists():return
    g=np.load(OUT/'geometry'/f'{sid}.npz');xy=g['xy'];labels=g['labels'];graph=sparse.load_npz(OUT/'geometry'/f'{sid}_graph.npz')
    info=next(r for r in M.tab(OUT/'sections.tsv') if r['section_id']==sid)
    endpoints=M.tab(OUT/'geometry'/f'{sid}_endpoints.tsv')
    for e in endpoints:
        for k in ['endpoint_id','spot_index','n_spots']:e[k]=int(e[k])
        for k in ['x_um','y_um','diameter_um']:e[k]=float(e[k])
    tlsmask=labels=='TLS';dtls=spatial.cKDTree(xy[tlsmask]).query(xy)[0]
    op_rows=[];metadata=[];side_rows=[];adjusted_rows=[];row_ids=[]
    eligibility=[];pairmeta=[];pindex=0
    sources=sorted({a['spot_index'] for _,a,b in pairs_for(endpoints)})
    distances,preds=csgraph.dijkstra(graph,directed=False,indices=sources,return_predecessors=True) if sources else ([],[])
    source_lookup={s:i for i,s in enumerate(sources)}
    for family,a,b in pairs_for(endpoints):
        source=a['spot_index'];target=b['spot_index'];start=xy[source];end=xy[target]
        srcids=g['endpoint_members'][g['endpoint_indptr'][a['endpoint_id']]:g['endpoint_indptr'][a['endpoint_id']+1]]
        d1=spatial.cKDTree(xy[srcids]).query(xy)[0];d2=np.linalg.norm(xy-end,axis=1)
        Z=np.column_stack([np.ones(len(xy))]+[np.exp(-d1/l) for l in [100,250,500,1000,2000,4000]]+[np.exp(-d2/l) for l in [100,250,500,1000,2000,4000]])
        for mode in ['straight','tissue_graph_shortest']:
            path=None
            if mode=='tissue_graph_shortest':
                path=trace_path(preds[source_lookup[source]],source,target)
                if path is None:
                    for radius in [130,230]:
                        for width in [100,250,500,1000]:eligibility.append(dict(pair_id=pindex,family=family,mode=mode,mask_um=radius,width_um=width,status='DISCONNECTED_TISSUE',length_um='',length_stratum='',effect_row=''))
                    continue
            along,side,length=path_coordinates(xy,start,end,path)
            pairmeta.append(dict(pair_id=pindex,family=family,source_endpoint=a['endpoint_id'],target_endpoint=b['endpoint_id'],target_family=b['family'],mode=mode,length_um=length,length_over_source_diameter=length/max(a['diameter_um'],100),source_size_spots=a['n_spots'],target_size_spots=b['n_spots']))
            for radius in [130,230]:
                visible=(dtls>radius)&(d2>radius)
                for width in [100,250,500,1000]:
                    rec=dict(pair_id=pindex,family=family,mode=mode,mask_um=radius,width_um=width,status='',length_um=length,length_stratum=scale_id(length),effect_row='')
                    if info['geometry_status']!='SOURCE_GEOMETRY_AVAILABLE':rec['status']='KNOWN_TLS_MISSING';eligibility.append(rec);continue
                    result,status,counts=corridor_weights(along,side,length,visible,width)
                    rec['status']=status
                    if status!='ELIGIBLE':eligibility.append(rec);continue
                    W,V,center=result;mean=W.mean(axis=0);adj=endpoint_adjust(mean,Z,visible&~center)
                    effect_id=len(metadata);rec['effect_row']=effect_id;eligibility.append(rec)
                    metadata.append(dict(**{k:rec[k] for k in ['pair_id','family','mode','mask_um','width_um','length_um','length_stratum']},effect_row=effect_id,n_visible=int(visible.sum()),endpoint_baseline_status='AVAILABLE' if adj is not None else 'NOT_TESTABLE',segment_counts=json.dumps(counts),target_spot=target))
                    # Sparse linear operators ensure every GO set uses exactly the same geometric support.
                    for row in W:op_rows.append(sparse.csr_matrix(row))
                    side_rows.append(sparse.csr_matrix(V.mean(axis=0)))
                    adjusted_rows.append(sparse.csr_matrix(adj if adj is not None else np.zeros(len(xy))))
        pindex+=1
    M.tsv(directory/'eligibility.tsv',eligibility);M.tsv(directory/'pairs.tsv',pairmeta);M.tsv(directory/'effect_rows.tsv',metadata)
    if op_rows:
        sparse.save_npz(directory/'segments.npz',sparse.vstack(op_rows,format='csr'))
        sparse.save_npz(directory/'side_asymmetry.npz',sparse.vstack(side_rows,format='csr'))
        sparse.save_npz(directory/'endpoint_adjusted.npz',sparse.vstack(adjusted_rows,format='csr'))
    M.js(directory/'operator_complete.json',dict(status='COMPLETE',pairs=pindex,attempted_geometries=len(eligibility),evaluable_geometries=len(metadata),go_scored=False))
    print('operators',sid,pindex,'pairs',len(metadata),'evaluable /',len(eligibility),flush=True)


def score_forward(sid):
    directory=OUT/'forward'/sid
    if (directory/'scoring_complete.json').exists():return
    done=json.loads((directory/'operator_complete.json').read_text());g=np.load(OUT/'geometry'/f'{sid}.npz');xy=g['xy'];labels=g['labels']
    cache=np.load(OUT/'scores'/f'{sid}.npz');S=np.nan_to_num(cache['scores']).astype(np.float64);names=cache['set_ids'];valid=cache['valid']
    dtls=spatial.cKDTree(xy[labels=='TLS']).query(xy)[0];info=next(r for r in M.tab(OUT/'sections.tsv') if r['section_id']==sid)
    # All radial and edge bins are retained as dense program arrays with explicit counts.
    radial=[];radmeta=[];edge=[];edgemeta=[]
    if info['geometry_status']=='SOURCE_GEOMETRY_AVAILABLE':
        for radius in [130,230]:
            vis=dtls>radius
            edges=[float(x) for x in [radius,230,400,800,1600,3200,6400] if x>=radius]
            edges=sorted(set(edges));edges=[x for x in edges if x<float(dtls.max())]+[float(dtls.max())+1e-6]
            for lo,hi in zip(edges[:-1],edges[1:]):
                ix=vis&(dtls>=lo)&(dtls<hi);n=int(ix.sum());radmeta.append(dict(mask_um=radius,lo_um=lo,hi_um=hi,n_spots=n,status='EVALUABLE' if n>=30 else 'NOT_TESTABLE_LT30'))
                radial.append(S[ix].mean(axis=0) if n>=30 else np.full(len(names),np.nan))
            if g['edge'].any():
                de=spatial.cKDTree(xy[g['edge']]).query(xy)[0]
                for comp in ['TUM','EXPLICIT_NON_TUMOR']:
                    labmask=labels=='TUM' if comp=='TUM' else np.isin(labels,['NOR','INFL','LN'])
                    for lo,hi in [(0,230),(230,500),(500,1000),(1000,2000),(2000,4000),(4000,float(de.max())+1e-6)]:
                        if hi<=lo:continue
                        ix=vis&labmask&(de>=lo)&(de<hi);n=int(ix.sum());edgemeta.append(dict(mask_um=radius,compartment=comp,lo_um=lo,hi_um=hi,n_spots=n,status='EVALUABLE' if n>=30 else 'NOT_TESTABLE_LT30'))
                        edge.append(S[ix].mean(axis=0) if n>=30 else np.full(len(names),np.nan))
    if radial:
        A=np.array(radial);A[:,~valid]=np.nan;np.savez_compressed(directory/'radial.npz',means=A.astype(np.float32),set_ids=names);M.tsv(directory/'radial_rows.tsv',radmeta)
    if edge:
        A=np.array(edge);A[:,~valid]=np.nan;np.savez_compressed(directory/'edge.npz',means=A.astype(np.float32),set_ids=names);M.tsv(directory/'edge_rows.tsv',edgemeta)
    n=done['evaluable_geometries']
    if n:
        rows=M.tab(directory/'effect_rows.tsv');W=sparse.load_npz(directory/'segments.npz');V=sparse.load_npz(directory/'side_asymmetry.npz');A=sparse.load_npz(directory/'endpoint_adjusted.npz')
        shape=(n,len(names));raw=np.lib.format.open_memmap(directory/'raw_effect.npy',mode='w+',dtype='float32',shape=shape);adjusted=np.lib.format.open_memmap(directory/'endpoint_residual.npy',mode='w+',dtype='float32',shape=shape);asym=np.lib.format.open_memmap(directory/'side_asymmetry.npy',mode='w+',dtype='float32',shape=shape);cont=np.lib.format.open_memmap(directory/'same_sign_segments.npy',mode='w+',dtype='uint8',shape=shape)
        sums={r:(S[dtls>r].sum(axis=0),(S[dtls>r]**2).sum(axis=0),int(sum(dtls>r))) for r in [130,230]}
        # Per-target SD excludes the endpoint disk as well as all TLS cores.
        sd_cache={}
        for start in range(0,n,128):
            end=min(start+128,n);R=np.asarray(W[start*5:end*5]@S).reshape(end-start,5,-1);av=np.asarray(A[start:end]@S);vv=np.asarray(V[start:end]@S)
            for k,row in enumerate(rows[start:end]):
                rad=int(row['mask_um']);target=int(row['target_spot']);key=(rad,target)
                if key not in sd_cache:
                    omit=(dtls>rad)&(np.linalg.norm(xy-xy[target],axis=1)<=rad);ss,ss2,nn=sums[rad];nn-=int(omit.sum());aa=ss-S[omit].sum(axis=0);bb=ss2-(S[omit]**2).sum(axis=0);sd_cache[key]=np.sqrt(np.maximum((bb-aa*aa/nn)/max(nn-1,1),0))
                sd=sd_cache[key];ok=valid&(sd>1e-12);scale=np.where(ok,sd,np.nan)
                raw[start+k]=(R[k].mean(axis=0)/scale).astype(np.float32);adjusted[start+k]=(av[k]/scale).astype(np.float32) if row['endpoint_baseline_status']=='AVAILABLE' else np.nan;asym[start+k]=(vv[k]/scale).astype(np.float32)
                cnt=np.maximum(np.sum(R[k]>0,axis=0),np.sum(R[k]<0,axis=0));cnt[~ok]=0;cont[start+k]=cnt
        for arr in [raw,adjusted,asym,cont]:arr.flush()
        M.js(directory/'program_axis.json',names.tolist())
    M.js(directory/'scoring_complete.json',dict(status='COMPLETE_EXPLORATORY',go_terms=1691,controls=51,attempted_bridge_geometries=done['attempted_geometries'],evaluable_bridge_geometries=n,radial_rows=len(radmeta),edge_rows=len(edgemeta),formal_p='NOT_CALIBRATED',tumor_status=info['tumor_status'],geometry_status=info['geometry_status']))
    print('forward scored',sid,n,'geometries x',len(names),'programs',flush=True)

def continuity(sid):
    directory=OUT/'forward'/sid
    done=json.loads((directory/'operator_complete.json').read_text());n=done['evaluable_geometries']
    if not n:
        M.js(directory/'continuity_complete.json',dict(status='NO_EVALUABLE_GEOMETRY',rows=0));return
    if (directory/'continuity_complete.json').exists():return
    cache=np.load(OUT/'scores'/f'{sid}.npz');S=np.nan_to_num(cache['scores']).astype(float)
    W=sparse.load_npz(directory/'segments.npz')
    raw=np.load(directory/'raw_effect.npy',mmap_mode='r')
    pos=np.lib.format.open_memmap(directory/'positive_segments.npy',mode='w+',dtype='uint8',shape=raw.shape)
    neg=np.lib.format.open_memmap(directory/'negative_segments.npy',mode='w+',dtype='uint8',shape=raw.shape)
    for a in range(0,n,128):
        b=min(n,a+128);R=np.asarray(W[a*5:b*5]@S).reshape(b-a,5,-1)
        pp=np.sum(R>0,axis=1);nn=np.sum(R<0,axis=1);invalid=~np.isfinite(raw[a:b])
        pp[invalid]=0;nn[invalid]=0;pos[a:b]=pp;neg[a:b]=nn
    pos.flush();neg.flush();M.js(directory/'continuity_complete.json',dict(status='COMPLETE',rows=n))
    print('directional continuity',sid,n,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['operators','score','continuity']);ap.add_argument('--section',default='');args=ap.parse_args()
    contract=json.loads((OUT/'contract.json').read_text())
    for info in contract['sources']:
        if args.section and args.section!=info['section_id']:continue
        {'operators':operators,'score':score_forward,'continuity':continuity}[args.stage](info['section_id'])
