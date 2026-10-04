#!/usr/bin/env python3
"""Uniform-query masked field-to-TLS benchmark; no GT geometry in expression features."""
from pathlib import Path
import argparse,importlib.util,json
import numpy as np
from scipy import sparse,spatial
from scipy.sparse import csgraph

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'infra/tls_multiscale_gobp_20261002'
def mod(name):
 spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
M=mod('tls_multiscale_gobp');RADII=[230,500,1000]


def query_labels(xy,labels,graph,radius):
    row=np.rint(xy[:,1]/(100*np.sqrt(3)/2)).astype(int);col=np.rint(xy[:,0]/50).astype(int)
    grid=np.flatnonzero((row%2==0)&(col%4==0));tls=np.flatnonzero(labels=='TLS')
    _,c=csgraph.connected_components(graph[tls][:,tls],directed=False)
    member={tls[j]:tls[c==c[j]] for j in range(len(tls))};rows=[]
    for i in grid:
        lab=labels[i];y=-1;status='UNKNOWN_CENTRE_LABEL'
        if lab=='TLS':
            ids=member[i]
            if len(ids)<2:status='SINGLETON_TLS_COMPONENT'
            elif np.linalg.norm(xy[ids]-xy[i],axis=1).max()+130>radius+1e-6:status='TARGET_CORE_NOT_FULLY_HIDDEN'
            else:y=1;status='LABELLED'
        elif lab in ['NOR','TUM','INFL','LN','NO_TLS']:y=0;status='LABELLED'
        dist=np.linalg.norm(xy-xy[i],axis=1);bins=np.searchsorted([radius,radius+250,radius+750,radius+1750,radius+3750],dist,side='right')-1
        visible=dist>radius;n0=sum(visible&(bins==0));n1=sum(visible&(bins==1))
        if status=='LABELLED' and (n0<20 or n1<20 or sum(visible)<100):status='INSUFFICIENT_VISIBLE_SUPPORT';y=-1
        rows.append(dict(spot_index=int(i),x_um=float(xy[i,0]),y_um=float(xy[i,1]),mask_um=radius,label=y,centre_annotation=lab or 'UNKNOWN',status=status,n_first_ring=int(n0),n_second_ring=int(n1)))
    return rows


def pool_query_features(xy,S,lib,queries,radius,chunk=32):
    """S has independent per-spot scores. All score moments exclude each query disk."""
    P=S.shape[1];features=np.empty((len(queries),P,15),dtype=np.float32);baseline=np.zeros((len(queries),17),np.float32)
    for a in range(0,len(queries),chunk):
        sub=queries[a:a+chunk];ri=[];ci=[];vv=[];counts=np.zeros((len(sub),5,8),int);means=[];stds=[]
        for qi,q in enumerate(sub):
            center=xy[q['spot_index']];rel=xy-center;dist=np.linalg.norm(rel,axis=1);vis=dist>radius
            b=np.searchsorted([radius,radius+250,radius+750,radius+1750,radius+3750],dist,side='right')-1
            ang=np.mod(np.arctan2(rel[:,1],rel[:,0]),2*np.pi);sector=np.minimum((ang/(2*np.pi)*8).astype(int),7)
            use=np.flatnonzero(vis);key=b[use]*8+sector[use];count=np.bincount(key,minlength=40);counts[qi]=count.reshape(5,8)
            ri.extend((qi*40+key).tolist());ci.extend(use.tolist());vv.extend((1/np.maximum(count[key],1)).tolist())
            mu=S[vis].mean(axis=0);sd=S[vis].std(axis=0,ddof=1);means.append(mu);stds.append(sd)
            baseline[a+qi,:5]=np.log1p(count.reshape(5,8).sum(axis=1));baseline[a+qi,5:10]=(count.reshape(5,8)>=5).sum(axis=1)/8
            for bi in range(5):
                ids=vis&(b==bi);baseline[a+qi,10+bi]=np.mean(np.log1p(lib[ids])) if ids.any() else 0
            baseline[a+qi,15:]=(center-xy.min(axis=0))/np.maximum(np.ptp(xy,axis=0),1)
        op=sparse.csr_matrix((np.array(vv,np.float32),(ri,ci)),shape=(len(sub)*40,len(xy)))
        values=np.asarray(op@S).reshape(len(sub),5,8,P)
        mu=np.array(means);sd=np.array(stds);safe=np.where(sd>1e-12,sd,1)
        theta=(np.arange(8)+.5)*2*np.pi/8
        ring=(values*counts[:,:,:,None]).sum(axis=2)/np.maximum(counts.sum(axis=2)[:,:,None],1)
        radial=(ring-mu[:,None,:])/safe[:,None,:];radial[counts.sum(axis=2)<20]=0
        centered=values-values.mean(axis=2,keepdims=True)
        amp=[]
        for order in [1,2]:
            re=np.mean(centered*np.cos(order*theta)[None,None,:,None],axis=2);im=np.mean(centered*np.sin(order*theta)[None,None,:,None],axis=2)
            v=2*np.sqrt(re*re+im*im)/safe[:,None,:];v[np.any(counts<5,axis=2)]=0;amp.append(v)
        F=np.concatenate([radial,*amp],axis=1).transpose(0,2,1)
        F[sd<=1e-12]=0
        features[a:a+len(sub)]=F
    return features,baseline


def prepare():
    directory=OUT/'reverse';directory.mkdir(exist_ok=True)
    if (directory/'contract.json').exists():raise FileExistsError('reverse contract exists')
    contract=json.loads((OUT/'contract.json').read_text());allrows=[]
    for info in contract['sources']:
        sid=info['section_id'];g=np.load(OUT/'geometry'/f'{sid}.npz');graph=sparse.load_npz(OUT/'geometry'/f'{sid}_graph.npz');meta=next(r for r in M.tab(OUT/'sections.tsv') if r['section_id']==sid)
        for rad in RADII:
            rows=query_labels(g['xy'],g['labels'],graph,rad)
            if meta['geometry_status']!='SOURCE_GEOMETRY_AVAILABLE':
                for r in rows:r['label']=-1;r['status']='KNOWN_TLS_MISSING_IN_SOURCE'
            for r in rows:r['section_id']=sid;r['cohort']=info['cohort']
            M.tsv(directory/f'{sid}_{rad}_queries.tsv',rows)
            allrows.append(dict(section_id=sid,cohort=info['cohort'],mask_um=rad,n_grid=len(rows),n_positive=sum(r['label']==1 for r in rows),n_negative=sum(r['label']==0 for r in rows),n_excluded=sum(r['label']<0 for r in rows)))
    M.tsv(directory/'eligibility.tsv',allrows)
    M.js(directory/'contract.json',dict(status='FROZEN_BEFORE_REVERSE_FEATURES',decision='D-173',grid='row%2==0 and doubled-column%4==0; independent of labels',radii_um=RADII,features='5 radial normalized means + 5 angular order1 amplitudes + 5 angular order2 amplitudes',
        mask='same circular mask for every query; all moments from outside circle; positive full annotated component plus130um hidden',target='central spot TLS vs explicit central non-TLS; not presence/absence anywhere in the window',
        models={'single_program':'covariance standardized ridge lambda=1,15 features','combined':'training-fold top10 GO by training AUC + 17 geometry/depth features; ridge lambda=1','baseline':'17 geometry/depth features','fixed_controls':51},
        split='leave-one-section-out and leave-cohort-out, train-only selection; not guaranteed distinct patients',missing='inadequate bins encode0 plus baseline support features; UNKNOWN target excluded, not relabelled negative',n_gobp=1691,n_controls=51,formal_patient_inference='NOT_TESTABLE_IDENTITY'))


def features():
    directory=OUT/'reverse';assert (directory/'contract.json').exists()
    main=json.loads((OUT/'contract.json').read_text())
    for info in main['sources']:
        sid=info['section_id'];g=np.load(OUT/'geometry'/f'{sid}.npz');cache=np.load(OUT/'scores'/f'{sid}.npz');S=np.nan_to_num(cache['scores']);lib=cache['libsize']
        for rad in RADII:
            out=directory/f'{sid}_{rad}_features.npz'
            if out.exists():continue
            rows=M.tab(directory/f'{sid}_{rad}_queries.tsv');qs=[]
            for r in rows:
                if int(r['label'])>=0:qs.append({'spot_index':int(r['spot_index']),'label':int(r['label'])})
            if qs:F,B=pool_query_features(g['xy'],S,lib,qs,rad)
            else:F=np.zeros((0,S.shape[1],15),np.float32);B=np.zeros((0,17),np.float32)
            F[:,~cache['valid']]=0
            with out.with_suffix('.partial').open('wb') as f:np.savez_compressed(f,features=F,baseline=B,y=np.array([q['label'] for q in qs]),spot_index=np.array([q['spot_index'] for q in qs]),valid_program=cache['valid'],set_ids=cache['set_ids'])
            out.with_suffix('.partial').replace(out)
            print('reverse features',sid,rad,len(qs),F.shape,flush=True)
    M.js(directory/'features_complete.json',dict(status='COMPLETE',source_mask_combinations=26*3,n_gobp=1691,n_controls=51))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare','features']);a=ap.parse_args();{'prepare':prepare,'features':features}[a.stage]()
