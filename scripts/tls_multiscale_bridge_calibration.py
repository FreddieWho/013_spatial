#!/usr/bin/env python3
"""Necessary exchangeability gate with identical geometry searches at 40 angles.

Passing this univariate gate does not validate the full joint real GO null.
"""
from pathlib import Path
import argparse,hashlib,importlib.util,json
import numpy as np
from scipy import sparse,spatial
from scipy.stats import beta

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'infra/tls_multiscale_gobp_20261002'
def mod(name):
    s=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
M=mod('tls_multiscale_gobp')
MODELS=['gaussian_150','gaussian_500','gaussian_1500','gaussian_4000','mixture','skew_heterogeneous_gradient']


def selected_statistics(E,positive,negative,target=None):
    # E: geometry x angle x field; identical candidates at every angle.
    mu=E.mean(axis=1,keepdims=True);sd=E.std(axis=1,ddof=1,keepdims=True)
    z=np.divide(E-mu,sd,out=np.zeros_like(E),where=np.isfinite(sd)&(sd>1e-12))
    use=((positive>=4)&(E>0)&(z>0))|((negative>=4)&(E<0)&(z<0))
    statistic=np.where(use,np.abs(z),0)
    maximum=statistic.max(axis=0)
    return maximum,statistic[target,0] if target is not None else None


def diagnostic_rank(maximum,target=None):
    if target is None:return np.sum(maximum>=maximum[0],axis=0)/len(maximum)
    return (1+np.sum(maximum[1:]>=target,axis=0))/len(maximum)


def bound(k,n,upper):
    if not n:return None
    if upper:return 1. if k==n else float(beta.ppf(.95,k+1,n-k))
    return 0. if k==0 else float(beta.ppf(.05,k,n-k+1))


def prepare(sid):
    dest=OUT/'bridge_calibration'/sid;dest.mkdir(parents=True,exist_ok=True)
    if (dest/'operators_complete.json').exists():return
    ctl=OUT/'controls'/sid;info=json.loads((ctl/'geometry_complete.json').read_text())
    if not info['rows']:
        M.js(dest/'operators_complete.json',dict(status='NO_EVALUABLE_GEOMETRY',n_geometry=0));return
    forward=OUT/'forward'/sid;rows=M.tab(forward/'effect_rows.tsv')
    q=np.load(ctl/'qualification.npz')['status'];complete=np.all(q==0,axis=1)
    ids=np.array([int(r['effect_row']) for r in rows if complete[int(r['effect_row'])] and float(r['length_um'])>=2*float(r['width_um'])],dtype=int)
    if not len(ids):
        M.js(dest/'operators_complete.json',dict(status='NO_COMMON_39_DIRECTION_GEOMETRY',n_geometry=0,attempted_geometry=len(rows)));return
    lookup={i:g for g,i in enumerate(ids)};selected=[rows[i] for i in ids];n=len(ids)
    W=sparse.load_npz(forward/'segments.npz');A=sparse.load_npz(forward/'endpoint_adjusted.npz')
    Wtrue=W[(ids[:,None]*5+np.arange(5)).ravel()];Atrue=A[ids]
    sparse.save_npz(dest/'true_segments.npz',Wtrue);sparse.save_npz(dest/'true_adjusted.npz',Atrue)
    Ws=[];Ls=[];Cs=[];slots=[];groups=[];reference_ids=[]
    for path in sorted(ctl.glob('chunk_*_rows.tsv')):
        meta=M.tab(path);keep=[i for i,r in enumerate(meta) if int(r['effect_row']) in lookup]
        if not keep:continue
        stem=path.name.replace('_rows.tsv','');ix=np.asarray(keep)
        Ws.append(sparse.load_npz(ctl/f'{stem}_segments.npz')[(ix[:,None]*5+np.arange(5)).ravel()])
        Ls.append(sparse.load_npz(ctl/f'{stem}_local.npz')[ix]);Cs.append(np.load(ctl/f'{stem}_coef.npy')[ix])
        for i in keep:
            r=meta[i];slots.append(lookup[int(r['effect_row'])]*40+int(r['angle_index'])+1)
            groups.append(int(r['group_id']));reference_ids.append(int(r['reference_id']))
    assert set(slots)=={g*40+a for g in range(n) for a in range(1,40)} and len(slots)==n*39
    sparse.save_npz(dest/'reference_segments.npz',sparse.vstack(Ws,format='csr'))
    sparse.save_npz(dest/'reference_local.npz',sparse.vstack(Ls,format='csr'))
    original_groups=sorted(set(groups));remap={g:i for i,g in enumerate(original_groups)}
    geometry=np.load(OUT/'geometry'/f'{sid}.npz');xy=geometry['xy'];labels=geometry['labels']
    endpoints={int(r['endpoint_id']):r for r in M.tab(OUT/'geometry'/f'{sid}_endpoints.tsv')}
    dtls=spatial.cKDTree(xy[labels=='TLS']).query(xy)[0];Zs=[]
    group_rows=M.tab(ctl/'groups.tsv')
    for gid in original_groups:
        row=group_rows[gid];a=endpoints[int(row['source_endpoint'])];b=endpoints[int(row['target_endpoint'])];rad=int(row['mask_um'])
        members=geometry['endpoint_members'][geometry['endpoint_indptr'][int(a['endpoint_id'])]:geometry['endpoint_indptr'][int(a['endpoint_id'])+1]]
        d1=spatial.cKDTree(xy[members]).query(xy)[0];d2=np.linalg.norm(xy-xy[int(b['spot_index'])],axis=1)
        vis=(dtls>rad)&(d2>rad)
        Z=np.column_stack([np.ones(len(xy))]+[np.exp(-d1/l) for l in [100,250,500,1000,2000,4000]]+[np.exp(-d2/l) for l in [100,250,500,1000,2000,4000]])
        Zs.append((Z*vis[:,None]).T)
    np.save(dest/'moment_operator.npy',np.concatenate(Zs,axis=0))
    ref_slots=np.array(slots);ref_groups=np.array([remap[g] for g in groups]);true_groups=np.zeros(n,int)
    for g in range(n):true_groups[g]=ref_groups[np.flatnonzero(ref_slots==g*40+1)[0]]
    np.savez_compressed(dest/'metadata.npz',effect_ids=ids,ref_slots=ref_slots,ref_groups=ref_groups,true_groups=true_groups,
                        reference_ids=np.array(reference_ids),coefficients=np.concatenate(Cs))
    M.tsv(dest/'common_geometries.tsv',selected)
    representatives={}
    for i,r in enumerate(selected):
        key=tuple(r[k] for k in ['family','mode','mask_um','length_stratum','width_um']);representatives.setdefault(key,[]).append(i)
    reps=[indices[len(indices)//2] for _,indices in sorted(representatives.items())]
    np.save(dest/'power_representatives.npy',np.asarray(reps))
    M.js(dest/'operators_complete.json',dict(status='COMMON_CANDIDATES_FROZEN',n_geometry=n,attempted_geometry=len(rows),n_power_groups=len(reps),n_angles=40))
    print('calibration operators',sid,n,'geometries',len(reps),'power groups',flush=True)


def load_operators(sid):
    d=OUT/'bridge_calibration'/sid;meta=dict(np.load(d/'metadata.npz'))
    for key,filename in [('W','true_segments'),('A','true_adjusted'),('RW','reference_segments'),('L','reference_local')]:meta[key]=sparse.load_npz(d/f'{filename}.npz')
    meta['Z']=np.load(d/'moment_operator.npy');meta['n']=len(meta['effect_ids'])
    return meta


def evaluate(S,op):
    n=op['n'];P=S.shape[1];R=np.empty((n*40,5,P));E=np.empty((n*40,P));true_slots=np.arange(n)*40
    R[true_slots]=np.asarray(op['W']@S).reshape(n,5,P);E[true_slots]=op['A']@S
    moments=(op['Z']@S).reshape(-1,13,P)
    # Chunk sparse reference products to bound transient arrays for full simulations.
    for lo in range(0,len(op['ref_slots']),256):
        hi=min(lo+256,len(op['ref_slots']));slots=op['ref_slots'][lo:hi]
        r=np.asarray(op['RW'][lo*5:hi*5]@S).reshape(hi-lo,5,P)
        correction=np.einsum('nd,ndp->np',op['coefficients'][lo:hi],moments[op['ref_groups'][lo:hi]],optimize=True)
        R[slots]=r;E[slots]=r.mean(axis=1)+op['L'][lo:hi]@S-correction
    return R.reshape(n,40,5,P),E.reshape(n,40,P)


def simulate(xy,model,seed,n=400):
    rng=np.random.default_rng(seed);x=xy-xy.mean(axis=0)
    def basis(length,count):
        frequencies=rng.normal(size=(2,count))/length;phase=rng.uniform(0,2*np.pi,count)
        return np.sqrt(2/count)*np.cos(x@frequencies+phase)
    if model.startswith('gaussian_'):B=basis(float(model.split('_')[1]),512)
    else:B=np.concatenate([np.sqrt(w)*basis(l,128) for l,w in zip([150,500,1500,4000],[.4,.3,.2,.1])],axis=1)
    S=B@rng.normal(size=(B.shape[1],n))
    if model=='skew_heterogeneous_gradient':
        loc=x/np.maximum(np.ptp(x,axis=0),1)
        S=(.5+np.abs(loc[:,0,None]))*(np.exp(.8*S)-np.exp(.32))+.75*loc[:,0,None]+.5*loc[:,1,None]
    return S


def real_statistics(sid):
    d=OUT/'bridge_calibration'/sid
    if (d/'real_complete.json').exists():return
    meta=json.loads((d/'operators_complete.json').read_text());n=meta['n_geometry']
    if not n:M.js(d/'real_complete.json',dict(status=meta['status'],n_geometry=0));return
    assert (OUT/'controls'/sid/'scoring_complete.json').exists()
    ix=np.load(d/'metadata.npz');ids=ix['effect_ids'];slots=ix['ref_slots'];rids=ix['reference_ids'];P=1742
    E=np.empty((n*40,P));pos=np.empty((n*40,P),np.uint8);neg=np.empty_like(pos);actual=np.arange(n)*40
    for arr,f,g in [(E,'endpoint_residual','adjusted'),(pos,'positive_segments','positive_segments'),(neg,'negative_segments','negative_segments')]:
        arr[actual]=np.load(OUT/'forward'/sid/f'{f}.npy',mmap_mode='r')[ids]
        arr[slots]=np.load(OUT/'controls'/sid/f'{g}.npy',mmap_mode='r')[rids]
    E=E.reshape(n,40,P);T,_=selected_statistics(E,pos.reshape(n,40,P),neg.reshape(n,40,P));ranks=diagnostic_rank(T)
    valid=np.any(np.all(np.isfinite(E),axis=1),axis=0);ranks[~valid]=np.nan
    names=np.load(OUT/'scores'/f'{sid}.npz')['set_ids']
    np.savez_compressed(d/'real_statistics.npz',max_by_direction=T,diagnostic_rank=ranks,set_ids=names,valid=valid)
    rows=[dict(program=str(name),diagnostic_rank=float(ranks[j]) if valid[j] else '',statistic=float(T[0,j]) if valid[j] else '',
               n_geometry=n,status='DIAGNOSTIC_NOT_CALIBRATED' if valid[j] else 'NOT_TESTABLE',formal_p='') for j,name in enumerate(names)]
    M.tsv(d/'all_program_diagnostics.tsv',rows)
    joint=T[:,:1691].max(axis=1)[:,None]
    M.js(d/'real_complete.json',dict(status='FULL_PANEL_DIAGNOSTIC_NOT_FORMAL_P',n_geometry=n,n_gobp=1691,n_controls=51,
        joint_1691_diagnostic_rank=float(diagnostic_rank(joint)[0]),patient_confirmation=False))


def calibrate(sid):
    d=OUT/'bridge_calibration'/sid
    if (d/'calibration_complete.json').exists():return
    meta=json.loads((d/'operators_complete.json').read_text());n=meta['n_geometry']
    if not n:M.js(d/'calibration_complete.json',dict(status=meta['status'],n_geometry=0));return
    op=load_operators(sid);xy=np.load(OUT/'geometry'/f'{sid}.npz')['xy'];reps=np.load(d/'power_representatives.npy')
    rows=M.tab(d/'common_geometries.tsv')
    # Every injected pattern is frozen by geometry; no real GO expression is inspected.
    H=np.column_stack([(np.asarray(op['W'][g*5:(g+1)*5].sum(axis=0)).ravel()>0).astype(float) for g in reps])
    for j,g in enumerate(reps):H[:,j]/=float(np.mean(op['W'][g*5:(g+1)*5]@H[:,j]))
    dR,dE=evaluate(H,op)
    allnull=[];allpower=[]
    for model in MODELS:
        receipt=d/f'{model}_complete.json'
        if receipt.exists():
            allnull.extend(M.tab(d/f'{model}_null.tsv'));allpower.extend(M.tab(d/f'{model}_power.tsv'));continue
        seed=int.from_bytes(hashlib.sha256((sid+'|'+model+'|D177').encode()).digest()[:8],'little')
        S=simulate(xy,model,seed);R,E=evaluate(S,op)
        T,_=selected_statistics(E,(R>0).sum(axis=2),(R<0).sum(axis=2));p=diagnostic_rank(T);k=int((p<=.05).sum())
        nr=dict(section_id=sid,model=model,n_replicates=400,n_false_positive=k,rate=k/400,upper95=bound(k,400,True),passed=bound(k,400,True)<=.08,n_common_geometry=n)
        allnull.append(nr);M.tsv(d/f'{model}_null.tsv',[nr]);np.savez_compressed(d/f'{model}_null.npz',max_by_direction=T,diagnostic_rank=p)
        V=op['Z'][::13];counts=V.sum(axis=1)[:,None];sums=V@S[:,:200]
        variance=np.maximum((V@(S[:,:200]**2)-sums*sums/counts)/(counts-1),0);SD=np.sqrt(variance)
        power=[]
        for j,g in enumerate(reps):
            for amplitude in [.5,1.]:
                gain=amplitude*SD[op['true_groups'][g]]
                r=R[:,:,:,:200]+dR[:,:,:,j,None]*gain;e=E[:,:,:200]+dE[:,:,j,None]*gain
                selected,own=selected_statistics(e,(r>0).sum(axis=2),(r<0).sum(axis=2),target=int(g))
                pp=diagnostic_rank(selected,own);hits=int(((pp<=.05)&(own>0)).sum());rr=rows[g]
                power.append(dict(section_id=sid,model=model,effect_row=rr['effect_row'],family=rr['family'],mode=rr['mode'],mask_um=rr['mask_um'],
                    length_stratum=rr['length_stratum'],halfwidth_um=rr['width_um'],amplitude_sd=amplitude,n_replicates=200,n_detected=hits,
                    power=hits/200,lower95=bound(hits,200,False),passed=bound(hits,200,False)>=.8))
        M.tsv(d/f'{model}_power.tsv',power);allpower.extend(power)
        M.js(receipt,dict(status='COMPLETE',seed=seed,n_null=400,n_power_cases=len(power),joint_GO_null_validated=False))
        print('calibration',sid,model,'FP',k,'/400','power cases',len(power),flush=True)
    def truth(value):return value is True or value=='True'
    null_pass=all(truth(r['passed']) for r in allnull)
    power_half=[r for r in allpower if float(r['amplitude_sd'])==.5]
    M.js(d/'calibration_complete.json',dict(status='NECESSARY_UNIVARIATE_GATE_PASSED_ONLY' if null_pass else 'ROTATION_EXCHANGEABILITY_GATE_FAILED',
        n_geometry=n,null_cells=len(allnull),failed_null_cells=sum(not truth(r['passed']) for r in allnull),power_cells=len(allpower),
        failed_power_half_sd=sum(not truth(r['passed']) for r in power_half),joint_GO_null_validated=False,formal_p_allowed=False))


def run_source(task):
    sid,stage=task
    prepare(sid)
    if stage in ['real','all']:real_statistics(sid)
    if stage in ['calibrate','all']:calibrate(sid)
    return sid


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare','real','calibrate','all']);ap.add_argument('--section',default='');ap.add_argument('--workers',type=int,default=1);ap.add_argument('--ready-only',action='store_true');a=ap.parse_args()
    tasks=[(r['section_id'],a.stage) for r in json.loads((OUT/'contract.json').read_text())['sources'] if not a.section or a.section==r['section_id']]
    if a.ready_only:
        tasks=[t for t in tasks if (OUT/'controls'/t[0]/('scoring_complete.json' if a.stage in ['real','all'] else 'geometry_complete.json')).exists()]
    if a.workers==1:
        for task in tasks:run_source(task)
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            for sid in pool.map(run_source,tasks):print('finished',a.stage,sid,flush=True)
