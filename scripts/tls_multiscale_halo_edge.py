#!/usr/bin/env python3
"""Full-panel halo/interface translation references and necessary null diagnostics."""
from pathlib import Path
import argparse,csv,hashlib,importlib.util,json
import numpy as np
from scipy import sparse,spatial

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'infra/tls_multiscale_gobp_20261002'
def mod(name):
    s=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
M=mod('tls_multiscale_gobp');B=mod('tls_multiscale_bridge_calibration')
FAMILIES=['halo','edge'];RADII=[130,230];LENGTHS=[100,250,500,1000,2000,4000]


def tsv(path,rows):
    if not rows:return
    with Path(path).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)),delimiter='\t');w.writeheader();w.writerows(rows)


def translations():
    y=100*np.sqrt(3)/2
    primitive=np.array([[100,0],[50,y],[-50,y],[-100,0],[-50,-y],[50,-y],
                        [150,y],[0,2*y],[-150,y],[-150,-y],[0,-2*y],[150,-y]],float)
    return np.vstack([np.zeros((1,2))]+[primitive*k for k in [1,2,3,5,8,13,20]])


def masks(xy,labels,dtls,distance,rows,family):
    vis=[(dtls>r)&((distance>r) if family=='halo' else True) for r in RADII]
    bins=[]
    for row in rows:
        mask=vis[RADII.index(int(row['mask_um']))].copy()
        mask&=(distance>=float(row['lo_um']))&(distance<float(row['hi_um']))
        if family=='edge':mask&=(labels=='TUM') if row['compartment']=='TUM' else np.isin(labels,['NOR','INFL','LN'])
        bins.append(mask)
    return np.array(bins),np.array(vis)


def contrasts(rows,family):
    cases=[];groups={}
    for i,r in enumerate(rows):groups.setdefault((r['mask_um'],r.get('compartment','ALL')),[]).append(i)
    for (radius,comp),indices in groups.items():
        indices.sort(key=lambda i:float(rows[i]['lo_um']))
        for a,b in zip(indices[:-1],indices[1:]):
            assert abs(float(rows[a]['hi_um'])-float(rows[b]['lo_um']))<1e-5
            cases.append(dict(kind='adjacent_distance',mask_um=radius,compartment=comp,a=a,b=b))
    if family=='edge':
        lookup={(r['mask_um'],r['lo_um'],r['hi_um'],r['compartment']):i for i,r in enumerate(rows)}
        for a,r in enumerate(rows):
            if r['compartment']!='TUM':continue
            key=(r['mask_um'],r['lo_um'],r['hi_um'],'EXPLICIT_NON_TUMOR')
            if key in lookup:cases.append(dict(kind='across_interface',mask_um=r['mask_um'],compartment='TUM_MINUS_EXPLICIT_NON_TUMOR',a=a,b=lookup[key]))
    for i,c in enumerate(cases):
        a,b=rows[c['a']],rows[c['b']];c.update(case_id=i,a_lo_um=a['lo_um'],a_hi_um=a['hi_um'],b_lo_um=b['lo_um'],b_hi_um=b['hi_um'],
            status='EVALUABLE' if a['status']==b['status']=='EVALUABLE' else 'NOT_TESTABLE_BIN_SUPPORT')
    return cases


def prepare(sid,family):
    dest=OUT/'halo_edge'/sid/family;dest.mkdir(parents=True,exist_ok=True)
    if (dest/'operators_complete.json').exists():return
    g=np.load(OUT/'geometry'/f'{sid}.npz');xy=g['xy'];labels=g['labels']
    source=next(r for r in M.tab(OUT/'sections.tsv') if r['section_id']==sid)
    filename='radial_rows.tsv' if family=='halo' else 'edge_rows.tsv';path=OUT/'forward'/sid/filename
    if not path.exists():
        reason='KNOWN_TLS_MISSING' if source['geometry_status']!='SOURCE_GEOMETRY_AVAILABLE' else ('NO_DISTANCE_BINS' if family=='halo' else ('NO_TUMOR_LABELS' if not np.any(labels=='TUM') else 'NO_EXPLICIT_ADJACENT_INTERFACE'))
        M.js(dest/'operators_complete.json',dict(status=reason,n_bins=0,n_cases=0,n_common_references=0));return
    rows=M.tab(path);cases=contrasts(rows,family);n=len(rows);offsets=translations();anchors=xy[labels=='TLS' if family=='halo' else g['edge']]
    tree=spatial.cKDTree(xy);dtls=spatial.cKDTree(xy[labels=='TLS']).query(xy)[0]
    original=np.array([int(r['n_spots']) for r in rows]);original_ok=np.array([r['status']=='EVALUABLE' for r in rows])
    Ws=[];Vs=[];counts=[];flags=[];coverage=[]
    for j,offset in enumerate(offsets):
        moved=anchors+offset;fraction=float(np.mean(tree.query(moved)[0]<=51));coverage.append(fraction)
        distance=spatial.cKDTree(moved).query(xy)[0];bins,visible=masks(xy,labels,dtls,distance,rows,family)
        number=bins.sum(axis=1);counts.append(number);status=np.zeros(n,np.uint8)
        status[~original_ok]=4;status[original_ok&(number<30)]=2
        ratio=number/np.maximum(original,1);status[original_ok&(number>=30)&((ratio<.5)|(ratio>2))]=3
        if fraction<.9:status[original_ok]=1
        if j==0:np.testing.assert_array_equal(number,original)
        flags.append(status);weights=bins/np.maximum(number[:,None],1)
        Ws.append(sparse.csr_matrix(weights));Vs.append(sparse.csr_matrix(visible.astype(float)))
    flags=np.array(flags);counts=np.array(counts)
    valid_cases=np.array([i for i,c in enumerate(cases) if c['status']=='EVALUABLE'],dtype=int)
    needed=sorted({c[key] for c in cases if c['status']=='EVALUABLE' for key in ['a','b']})
    common=np.flatnonzero(np.all(flags[:,needed]==0,axis=1)) if needed else np.array([],int)
    if len(common):assert common[0]==0
    sparse.save_npz(dest/'bin_operators.npz',sparse.vstack(Ws,format='csr'));sparse.save_npz(dest/'visibility.npz',sparse.vstack(Vs,format='csr'))
    np.savez_compressed(dest/'geometry.npz',status=flags,counts=counts,offsets=offsets,anchor_coverage=coverage,valid_cases=valid_cases,common_positions=common,
                        bin_masks=np.array([RADII.index(int(r['mask_um'])) for r in rows]))
    tsv(dest/'bins.tsv',rows);tsv(dest/'cases.tsv',cases)
    pos=[dict(position_id=i,dx_um=float(v[0]),dy_um=float(v[1]),anchor_coverage=coverage[i],is_common=i in common,n_matched_bins=int(sum(flags[i]==0))) for i,v in enumerate(offsets)]
    tsv(dest/'positions.tsv',pos)
    M.js(dest/'operators_complete.json',dict(status='COMPLETE_GEOMETRY',n_bins=n,n_cases=len(cases),n_evaluable_cases=len(valid_cases),n_common_references=max(0,len(common)-1),n_anchors=len(anchors),n_positions=len(offsets),
        flags={'0':'EVALUABLE','1':'ANCHOR_COVERAGE_LT90_PERCENT','2':'BIN_LT30','3':'COUNT_MISMATCH','4':'ORIGINAL_BIN_NOT_TESTABLE'}))


def load(sid,family):
    d=OUT/'halo_edge'/sid/family;op=dict(np.load(d/'geometry.npz'));op.update(W=sparse.load_npz(d/'bin_operators.npz'),V=sparse.load_npz(d/'visibility.npz'),rows=M.tab(d/'bins.tsv'),cases=M.tab(d/'cases.tsv') if (d/'cases.tsv').exists() else [])
    return op


def evaluate(S,op):
    J,n=op['status'].shape;P=S.shape[1];means=(op['W']@S).reshape(J,n,P)
    counts=np.asarray(op['V'].sum(axis=1));sums=op['V']@S;square=op['V']@(S*S)
    variance=np.maximum((square-sums*sums/np.maximum(counts,1))/np.maximum(counts-1,1),0)
    sd=np.sqrt(variance).reshape(J,2,P)
    cases=op['cases'];E=np.empty((len(cases),J,P))
    for i,c in enumerate(cases):
        a,b=int(c['a']),int(c['b']);scale=sd[:,RADII.index(int(c['mask_um']))]
        ok=(op['status'][:,a]==0)&(op['status'][:,b]==0)
        E[i]=np.divide(means[:,a]-means[:,b],scale,out=np.full((J,P),np.nan),where=ok[:,None]&(scale>1e-12))
    means=np.where(op['status'][:,:,None]==0,means,np.nan)
    return means,E,sd


def search(E,target=None):
    mu=E.mean(axis=1,keepdims=True);sd=E.std(axis=1,ddof=1,keepdims=True)
    Z=np.divide(E-mu,sd,out=np.zeros_like(E),where=np.isfinite(sd)&(sd>1e-12))
    stat=np.abs(Z);return stat.max(axis=0),stat[target,0] if target is not None else None


def score(sid,family):
    d=OUT/'halo_edge'/sid/family
    if (d/'scoring_complete.json').exists():return
    meta=json.loads((d/'operators_complete.json').read_text());cache=np.load(OUT/'scores'/f'{sid}.npz');names=cache['set_ids']
    if not meta['n_bins']:
        tsv(d/'all_program_diagnostics.tsv',[dict(program=str(name),diagnostic_rank='',status=meta['status'],formal_p='') for name in names])
        M.js(d/'scoring_complete.json',dict(status=meta['status'],n_gobp=1691,n_controls=51));return
    op=load(sid,family);S=np.nan_to_num(cache['scores']).astype(float);means,E,sd=evaluate(S,op)
    means[:,:,~cache['valid']]=np.nan;E[:,:,~cache['valid']]=np.nan
    # Original identity-position bins must exactly reproduce the existing full scan.
    original=np.load(OUT/'forward'/sid/('radial.npz' if family=='halo' else 'edge.npz'))['means']
    np.testing.assert_allclose(means[0],original,rtol=1e-6,atol=1e-7,equal_nan=True)
    np.savez_compressed(d/'full_panel.npz',bin_means=means.astype(np.float32),effect_sd=E.astype(np.float32),set_ids=names,visible_sd=sd.astype(np.float32))
    diagnostic=[];common=op['common_positions'];ids=op['valid_cases'];T=None
    if len(common)>=20 and len(ids):
        selected=E[ids][:,common];T,_=search(selected);rank=B.diagnostic_rank(T);valid=np.any(np.all(np.isfinite(selected),axis=1),axis=0)
        for j,name in enumerate(names):diagnostic.append(dict(program=str(name),diagnostic_rank=float(rank[j]) if valid[j] else '',status='DIAGNOSTIC_NOT_CALIBRATED' if valid[j] else 'NOT_TESTABLE_PROGRAM',formal_p=''))
        np.savez_compressed(d/'diagnostics.npz',max_by_position=T,diagnostic_rank=np.where(valid,rank,np.nan),set_ids=names)
    else:
        for name in names:diagnostic.append(dict(program=str(name),diagnostic_rank='',status='NOT_TESTABLE_COMMON_REFERENCE_RESOLUTION',formal_p=''))
    tsv(d/'all_program_diagnostics.tsv',diagnostic)
    M.js(d/'scoring_complete.json',dict(status='FULL_PANEL_DESCRIPTIVE_ONLY',n_gobp=1691,n_controls=51,n_common_references=max(0,len(common)-1),
        n_cases=len(op['cases']),joint_1691_diagnostic_rank=float(B.diagnostic_rank(T[:,:1691].max(axis=1)[:,None])[0]) if T is not None else None))


def calibrate(sid,family):
    d=OUT/'halo_edge'/sid/family
    if (d/'calibration_complete.json').exists():return
    meta=json.loads((d/'operators_complete.json').read_text())
    if meta['n_common_references']<19:
        M.js(d/'calibration_complete.json',dict(status='NOT_TESTABLE_COMMON_REFERENCE_RESOLUTION' if meta['n_bins'] else meta['status'],n_common_references=meta['n_common_references'],formal_p_allowed=False));return
    op=load(sid,family);ids=op['valid_cases'];positions=op['common_positions'];g=np.load(OUT/'geometry'/f'{sid}.npz');xy=g['xy'];labels=g['labels']
    dtls=spatial.cKDTree(xy[labels=='TLS']).query(xy)[0];dist=dtls if family=='halo' else spatial.cKDTree(xy[g['edge']]).query(xy)[0]
    patterns=[];templates=[]
    for kind in (['halo'] if family=='halo' else ['edge_peak','edge_jump']):
        sign=np.where(labels=='TUM',1,np.where(np.isin(labels,['NOR','INFL','LN']),-1,0)) if kind=='edge_jump' else np.ones(len(xy))
        for length in LENGTHS:patterns.append(np.exp(-dist/length)*sign);templates.append(dict(kind=kind,length_um=length))
    H=np.column_stack(patterns);raw=(op['W']@H).reshape(len(op['offsets']),len(op['rows']),-1)
    contrasts=np.array([raw[0,int(op['cases'][i]['a'])]-raw[0,int(op['cases'][i]['b'])] for i in ids])
    target=np.argmax(np.abs(contrasts),axis=0);normalizer=np.max(np.abs(contrasts),axis=0);usable=normalizer>1e-12
    H[:,usable]/=normalizer[usable]
    # Keep exact raw linear responses; visible SD is recomputed after each injection.
    null=[];power=[]
    for model in B.MODELS:
        receipt=d/f'{model}_complete.json'
        if receipt.exists():null.extend(M.tab(d/f'{model}_null.tsv'));power.extend(M.tab(d/f'{model}_power.tsv'));continue
        seed=int.from_bytes(hashlib.sha256((sid+'|'+family+'|'+model+'|D179').encode()).digest()[:8],'little')
        S=B.simulate(xy,model,seed);_,E,SD=evaluate(S,op);T,_=search(E[ids][:,positions]);p=B.diagnostic_rank(T);k=int(sum(p<=.05))
        nr=dict(section_id=sid,family=family,model=model,n_replicates=400,n_false_positive=k,rate=k/400,upper95=B.bound(k,400,True),passed=B.bound(k,400,True)<=.08,n_common_references=len(positions)-1)
        null.append(nr);tsv(d/f'{model}_null.tsv',[nr]);np.savez_compressed(d/f'{model}_null.npz',max_by_position=T,diagnostic_rank=p)
        pr=[]
        for j,t in enumerate(templates):
            if not usable[j]:continue
            case=op['cases'][ids[target[j]]];radius=RADII.index(int(case['mask_um']));background_sd=SD[0,radius,:200]
            for amplitude in [.5,1.]:
                injected=S[:,:200]+H[:,j,None]*(amplitude*background_sd)
                _,e,_=evaluate(injected,op);selected,own=search(e[ids][:,positions],target=int(target[j]));rank=B.diagnostic_rank(selected,own)
                hits=int(sum((rank<=.05)&(own>0)));lower=B.bound(hits,200,False)
                pr.append(dict(section_id=sid,family=family,model=model,template=t['kind'],length_um=t['length_um'],target_case=case['case_id'],target_mask_um=case['mask_um'],amplitude_sd=amplitude,n_replicates=200,n_detected=hits,power=hits/200,lower95=lower,passed=lower>=.8))
        tsv(d/f'{model}_power.tsv',pr);power.extend(pr);M.js(receipt,dict(status='COMPLETE',seed=seed,n_null=400,n_power_settings=len(pr)))
        print('halo/edge calibration',sid,family,model,k,'/400',flush=True)
    def yes(value):return value is True or value=='True'
    passed=all(yes(r['passed']) for r in null);half=[r for r in power if float(r['amplitude_sd'])==.5]
    M.js(d/'calibration_complete.json',dict(status='NECESSARY_SINGLE_FIELD_GATE_PASSED_ONLY' if passed else 'TRANSLATION_GATE_FAILED',n_common_references=len(positions)-1,
        null_settings=len(null),failed_null_settings=sum(not yes(r['passed']) for r in null),power_settings=len(power),half_sd_settings=len(half),failed_half_sd_settings=sum(not yes(r['passed']) for r in half),formal_p_allowed=False,joint_GO_null_validated=False))


def run_source(task):
    sid,stage=task
    for family in FAMILIES:
        prepare(sid,family)
        if stage in ['score','all']:score(sid,family)
        if stage in ['calibrate','all']:calibrate(sid,family)
    return sid


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare','score','calibrate','all']);ap.add_argument('--section',default='');ap.add_argument('--workers',type=int,default=1);a=ap.parse_args()
    tasks=[(r['section_id'],a.stage) for r in json.loads((OUT/'contract.json').read_text())['sources'] if not a.section or a.section==r['section_id']]
    if a.workers==1:
        for task in tasks:run_source(task)
    else:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=a.workers) as pool:
            for sid in pool.map(run_source,tasks):print('source finished',sid,a.stage,flush=True)
