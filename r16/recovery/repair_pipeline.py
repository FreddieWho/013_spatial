"""D-140 repair replay. Historical artifacts are input-only. CPU/local data only."""
from __future__ import annotations
import argparse, csv, hashlib, json, time, shutil
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import rankdata
from r16 import axes as A, census as C
from r16.section_io import load_section_symbols, usz_tls_labels
from r16.recovery.corrected import (NotTestable, NaturalSpline, classify_curve,
    classify_window, hex_coordinates, tls_components, signed_tls_distance,
    feature_block, fit_from_moments, fit_linear, masked_predict, neighbor_values,
    patient_interval)
from r16.recovery.mask_tasks import generate_windows, window_mask

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'infra/repair_20260921'
OLD=ROOT/'infra/r16/recovery_20260918'

def digest(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8<<20),b''): h.update(b)
    return h.hexdigest()

def dump(name,obj):
    p=OUT/name; p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def table(name,rows):
    rows=list(rows)
    p=OUT/name; p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,'w') as f:
        if rows:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def readtable(p):
    with open(p) as f: return list(csv.DictReader(f,delimiter='\t'))

def axes():
    d=json.load(open(ROOT/'infra/r04/marker_proxy_combined.json'))
    v={k:list(r['voted_genes']) for k,r in d['classes'].items() if k not in A.RETIRED_AXES}
    v['Plasma']=list(A.PLASMA_GENES)
    return v

def definitions(go=False):
    return json.load(open(OUT/('go_representatives.json' if go else 'six_programs.json')))

def sections(cohorts=None):
    rows=json.load(open(OUT/'sections.json'))
    return [r for r in rows if cohorts is None or r['cohort'] in cohorts]

def load(r):
    p=OUT/'cache'/r['key']
    return np.load(str(p)+'_counts.npy',mmap_mode='r'),np.load(str(p)+'_coords.npy'),np.load(str(p)+'_labels.npy')

def source_coverage():
    """Source metadata coverage, kept distinct from eligibility in the matched panel."""
    import gzip,h5py
    from r04.io_contract import _decode_strings
    defs=definitions();rec={r['section']:r for r in sections()};out=[];hashes={}
    for mf,cohort in [('internal_validation','ST-CRC'),('external_validation','USZ')]:
        for row in json.load(open(ROOT/f'infra/r04/role_manifests/{mf}_manifest.json'))['rows']:
            p=Path(row['matrix_locator']);r=rec[row['section_id']]
            if p.is_dir():
                f=next(p/name for name in ['features.tsv.gz','features.tsv','genes.tsv.gz','genes.tsv'] if (p/name).exists())
                with (gzip.open(f,'rt') if f.suffix=='.gz' else open(f)) as stream:gs=[line.rstrip('\n').split('\t')[1] for line in stream]
                hashes[str(f)]=digest(f)
                matrix_file=next(p/name for name in ['matrix.mtx','matrix.mtx.gz'] if (p/name).exists())
                hashes[str(matrix_file)]=digest(matrix_file)
            else:
                with h5py.File(p) as handle:gs=list(_decode_strings(handle['matrix']['features']['name'][()]))
                hashes[str(p)]=digest(p)
            if len(set(gs))!=len(gs):gs=sorted(set(gs))
            source=OUT/'cache'/f'{r["key"]}_source.npz'
            with np.load(source) as z:assert int(z['shape'][1])==len(gs)
            dump('cache/'+r['key']+'_source_genes.json',gs)
            for g,d in defs.items():
                missing_in=sorted(set(d['input_genes'])-set(gs));missing_ro=sorted(set(d['readout_genes'])-set(gs))
                out.append(dict(cohort=cohort,section=r['section'],program=g,missing_input=';'.join(missing_in),missing_readout=';'.join(missing_ro),source_full_coverage=not(missing_in or missing_ro),matched_panel_full_coverage=not d['missing_genes']))
    table('six_source_coverage.tsv',out)
    cache=ROOT/'infra/r16/census_cache_hvg10k/sections';records=sections();coverage=[]
    panel=set(json.load(open(OUT/'genes.json')))
    for r in records:
        if r['cohort']!='DISCOVERY':continue
        meta_path=cache/(r['section']+'.json');meta=json.load(open(meta_path));p=Path(meta['source_matrix_locator'])
        with h5py.File(p) as f:gs=list(_decode_strings(f['var']['_index'][()]))
        absent=sorted(panel-set(gs));r['unmeasured_genes']=absent
        r['source_var_sha256']=hashlib.sha256(json.dumps(gs).encode()).hexdigest()
        dump('cache/'+r['key']+'_source_genes.json',gs)
        with np.load(cache/(r['section']+'.npz')) as z:
            assert len(set(meta['genes'])-set(gs))==meta['n_genes_missing_in_source']
        hashes[str(meta_path)]=digest(meta_path);hashes[str(cache/(r['section']+'.npz'))]=digest(cache/(r['section']+'.npz'))
        coverage.append(dict(section=r['section'],patient=r['patient'],n_unmeasured_panel_genes=len(absent),genes=';'.join(absent)))
    dump('sections.json',records);table('discovery_source_coverage.tsv',coverage);dump('source_feature_hashes.json',hashes)

def prepare():
    if shutil.disk_usage(ROOT).free < 1_200_000_000_000:
        raise RuntimeError('Storage reserve below 1.2 TB')
    OUT.mkdir(exist_ok=True);(OUT/'cache').mkdir(exist_ok=True)
    records=[]; common=None; hashes={}; t=time.time()
    # Temporary section-specific genes remain on disk so alignment is checked for every section.
    cache=ROOT/'infra/r16/census_cache_hvg10k/sections'
    for i,stem in enumerate(C.list_stems(cache)):
        mat,bc,genes,xy,meta=C.load_section(cache,stem)
        gs=[str(g) for g in genes];key=f'disc_{i:02d}'
        sparse.save_npz(OUT/'cache'/f'{key}_source.npz',sparse.csr_matrix(mat))
        np.save(OUT/'cache'/f'{key}_coords.npy',xy)
        np.save(OUT/'cache'/f'{key}_labels.npy',np.full((len(xy),2),np.nan))
        records.append(dict(key=key,cohort='DISCOVERY',patient=meta['patient_id'],section=stem,genes=gs,barcodes=[str(b) for b in bc]))
        common=set(gs) if common is None else common&set(gs)
    for mf,cohort in [('internal_validation','ST-CRC'),('external_validation','USZ')]:
        p=ROOT/f'infra/r04/role_manifests/{mf}_manifest.json';hashes[str(p.relative_to(ROOT))]=digest(p)
        for i,row in enumerate(json.load(open(p))['rows']):
            sec=load_section_symbols(row);gs=[str(g) for g in sec.gene_id];key=f'{cohort}_{i:02d}'
            if len(gs)!=len(set(gs)):
                # Sum duplicate symbols rather than silently picking the first column.
                unique=sorted(set(gs));gx={g:j for j,g in enumerate(unique)}
                mapping=sparse.csr_matrix((np.ones(len(gs)),(np.arange(len(gs)),[gx[g] for g in gs])),shape=(len(gs),len(unique)))
                mat=sparse.csr_matrix(sec.counts)@mapping;gs=unique
            else: mat=sparse.csr_matrix(sec.counts)
            sparse.save_npz(OUT/'cache'/f'{key}_source.npz',mat)
            xy=hex_coordinates(np.asarray(sec.coords,float))
            np.save(OUT/'cache'/f'{key}_coords.npy',xy)
            labels=np.full((len(xy),2),np.nan)
            if cohort=='USZ':
                labels[:,0]=usz_tls_labels(sec.section_id,tuple(sec.barcode),ROOT)
                import h5py
                from r04.io_contract import _decode_strings
                path=ROOT/'data/other_sources/zenodo_usz_tls_visium/TLS_VISIUM_USZ/h5ad_preprocessed'/f'{sec.section_id.split("::")[-1]}.h5ad'
                with h5py.File(path) as f:
                    bs=_decode_strings(f['obs']['_index'][()]);cats=_decode_strings(f['obs']['ground_truth']['categories'][()]);codes=f['obs']['ground_truth']['codes'][()]
                lookup={b:(cats[int(c)] if 0<=int(c)<len(cats) else 'UNASSIGNED') for b,c in zip(bs,codes)}
                labels[:,1]=[np.nan if lookup.get(b,'UNASSIGNED')=='UNASSIGNED' else float(lookup[b]=='TUM') for b in sec.barcode]
                hashes[str(path.relative_to(ROOT))]=digest(path)
            np.save(OUT/'cache'/f'{key}_labels.npy',labels)
            records.append(dict(key=key,cohort=cohort,patient=sec.patient_id,section=sec.section_id,genes=gs,barcodes=list(sec.barcode)))
            common&=set(gs)
            print('loaded',key,len(xy),flush=True)
    common=sorted(common);dump('genes.json',common)
    for r in records:
        key=r['key'];gx={g:i for i,g in enumerate(r.pop('genes'))}
        mat=sparse.load_npz(OUT/'cache'/f'{key}_source.npz')[:,[gx[g] for g in common]].toarray().astype(np.float32)
        assert np.isfinite(mat).all() and (mat>=0).all()
        np.save(OUT/'cache'/f'{key}_counts.npy',mat)
        r['n_spots']=len(mat);r['counts_sha256']=digest(OUT/'cache'/f'{key}_counts.npy')
    dump('sections.json',records)
    six=json.load(open(OLD/'programs/program_definitions.json'))
    gd=json.load(open(OLD/'programs_go/definitions.json'));im=json.load(open(OLD/'programs_go/id_map.json'))
    reps={};aliases={};stamp=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
    for target,source in [(six,six),(reps,{g:gd[g] for g in im})]:
        for g,d in list(source.items()):
            d=dict(d);d['frozen_for_external_at']=stamp
            d['normalization_rule']='raw mean log1p; QC and composition exclude input+readout; D-140'
            d['nuisance_rule']='intercept + log1p background counts/detected genes + raw mean-log six clean axes; no per-section z-score'
            d['geometry_definition']='external Visium hex pitch; annotated TLS spot-set distance for spline; GT-independent windows for mask'
            d['coverage_policy']='complete program coverage in actual 69-section common gene universe; empty composition axis -> NOT_SEPARABLE'
            d['exposure']='post-audit exploratory reanalysis of previously observed cohorts; not a new independent confirmatory test'
            d['missing_genes']=sorted(set(d['input_genes']+d['readout_genes'])-set(common))
            assert not set(d['input_genes'])&set(d['readout_genes'])
            d.pop('definition_hash',None);d['definition_hash']=hashlib.sha256(json.dumps(d,sort_keys=True).encode()).hexdigest();target[g]=d
    for rep,ids in im.items():
        for alias in ids: aliases[alias]={'representative':rep,'input_genes':reps[rep]['input_genes'],'readout_genes':reps[rep]['readout_genes'],'definition_hash':reps[rep]['definition_hash']}
    dump('six_programs.json',six);dump('go_representatives.json',reps);dump('go_aliases.json',aliases)
    for p in [OLD/'programs/program_definitions.json',OLD/'programs_go/definitions.json',OLD/'programs_go/id_map.json']:
        hashes[str(p.relative_to(ROOT))]=digest(p)
    dump('input_receipt.json',dict(status='PREPARED',elapsed_seconds=time.time()-t,common_genes=len(common),sections=len(records),input_sha256=hashes,definition_hashes={g:d['definition_hash'] for g,d in {**six,**reps}.items()}))
    source_coverage()

def weight_matrix(gx,sets):
    rr=[];cc=[];vv=[]
    for j,gs in enumerate(sets):
        for g in gs:rr.append(gx[g]);cc.append(j);vv.append(1/len(gs))
    return sparse.csc_matrix((vv,(rr,cc)),shape=(len(gx),len(sets)))

def moran_columns(values,w):
    z=values-values.mean(0);den=np.sum(z*z,axis=0)
    return np.divide(len(z)/float(w.sum())*np.sum(z*(w@z),axis=0),den,out=np.zeros_like(den),where=den>0)

def raw_blocks(counts,xy,genes,defs,unmeasured=()):
    """Yield bounded program batches; QC background explicitly excludes both gene sets."""
    gx={g:i for i,g in enumerate(genes)};gene_set=set(gx)-set(unmeasured);ax=axes()
    counts=np.asfortranarray(counts,dtype=float);log=np.asfortranarray(np.log1p(counts));present=np.asfortranarray(counts>0,dtype=float)
    axis_sets=[set(g for g in gs if g in gene_set) for gs in ax.values()]
    axis_sizes=np.array([len(gs) for gs in axis_sets])
    axis_sums=np.column_stack([log[:,[gx[g] for g in gs]].sum(1) for gs in axis_sets])
    total=counts.sum(1);detected=present.sum(1)
    dist,ix=cKDTree(xy).query(xy,k=min(9,len(xy)));dist,ix=dist[:,1:],ix[:,1:]
    nw=1/np.maximum(dist,1e-6);nw/=nw.sum(1,keepdims=True)
    ids=list(defs)
    for start in range(0,len(ids),256):
        batch=ids[start:start+256];valid=[];bad=[];sets=[];exc=[];overlap_sizes=[]
        for g in batch:
            d=defs[g];ex=set(d['input_genes']+d['readout_genes'])
            overlap=[sorted(gs&ex) for gs in axis_sets]
            if not ex<=gene_set:bad.append((g,'MISSING_PROGRAM_GENES'));continue
            if any(len(gs)==len(ov) for gs,ov in zip(axis_sets,overlap)):bad.append((g,'NOT_SEPARABLE'));continue
            valid.append(g);sets.extend(overlap+[d['input_genes'],d['readout_genes']]);exc.append([gx[s] for s in ex]);overlap_sizes.append([len(v) for v in overlap])
        if not valid:yield [],None,None,bad;continue
        scores=(log@weight_matrix(gx,sets)).reshape(len(counts),len(valid),8)
        removed=np.asarray(overlap_sizes)
        scores[:,:,:6]=(axis_sums[:,None,:]-scores[:,:,:6]*removed[None,:,:])/(axis_sizes[None,None,:]-removed[None,:,:])
        x=np.ones((len(counts),len(valid),11));x[:,:,3:9]=scores[:,:,:6];x[:,:,9]=scores[:,:,6]
        er=[i for cols in exc for i in cols];ec=[j for j,cols in enumerate(exc) for _ in cols]
        ew=sparse.csc_matrix((np.ones(len(er)),(er,ec)),shape=(len(genes),len(valid)))
        x[:,:,1]=np.log1p(np.maximum(total[:,None]-counts@ew,0))
        x[:,:,2]=np.log1p(np.maximum(detected[:,None]-present@ew,0))
        x[:,:,10]=np.sum(scores[ix,:,6]*nw[:,:,None],axis=1)
        yield valid,x,scores[:,:,7],bad

def _section_replay(task):
    r,defs,genes,mode=task
    mom=defaultdict(dict);disc=[];invalid=[];molecular=[]
    counts,xy,_=load(r)
    if mode=='raw':blocks=raw_blocks(counts,xy,genes,defs,r.get('unmeasured_genes',()))
    else:
        def normblocks():
            for g,d in defs.items():
                try:
                    f=feature_block(counts,genes,d['input_genes'],d['readout_genes'],axes(),mode)
                    x=np.c_[f['X'],neighbor_values(xy,f['X'][:,-1],exclude_self=True)]
                    yield [g],x[:,None,:],f['y'][:,None],[]
                except NotTestable as e:yield [],None,None,[(g,str(e))]
        blocks=normblocks()
    w=C.spatial_weights(xy,C.KNN_SPATIAL)
    for ids,x,y,bad in blocks:
        invalid.extend(dict(program=g,cohort=r['cohort'],patient=r['patient'],section=r['section'],status=s) for g,s in bad)
        if r['cohort']!='DISCOVERY':
            gx={g:i for i,g in enumerate(genes)}
            for g,status in bad:
                d=defs[g];ex=set(d['input_genes']+d['readout_genes'])
                if not ex<=set(gx):continue
                scale=(1e4/np.maximum(counts[:,[i for i,s in enumerate(genes) if s not in ex]].sum(1,dtype=float),1))[:,None] if mode=='depth_normalized' else 1
                input_score=np.log1p(counts[:,[gx[s] for s in d['input_genes']]]*scale).mean(1)
                ro=np.log1p(counts[:,[gx[s] for s in d['readout_genes']]]*scale).mean(1)
                corr=float(np.corrcoef(input_score,ro)[0,1]) if min(input_score.var(),ro.var())>1e-12 else None
                molecular.append(dict(program=g,cohort=r['cohort'],patient=r['patient'],section=r['section'],input_readout_corr=corr,input_moran=C.moran_i(input_score,w),readout_moran=C.moran_i(ro,w),status='MOLECULAR_ONLY_'+status))
        if not ids:continue
        xb=np.ascontiguousarray(x.transpose(1,0,2))
        gram=(xb.transpose(0,2,1)@xb)/len(x)
        rhs=np.einsum('nbi,nb->bi',x,y,optimize=True)/len(x);yy=np.mean(y*y,axis=0)
        if r['cohort']=='DISCOVERY':
            v=x[:,:,9];bq=fit_from_moments(gram[:,:3,:3],gram[:,:3,9]);bc=fit_from_moments(gram[:,:9,:9],gram[:,:9,9])
            rq=v-np.einsum('nbi,bi->nb',x[:,:,:3],bq);rc=v-np.einsum('nbi,bi->nb',x[:,:,:9],bc)
            mi,mq,mc=moran_columns(v,w),moran_columns(rq,w),moran_columns(rc,w)
            for j,g in enumerate(ids):disc.append(dict(program=g,patient=r['patient'],section=r['section'],raw=mi[j],Q=mq[j],QC=mc[j]))
        else:
            mi,mr=moran_columns(x[:,:,9],w),moran_columns(y,w)
            for j,g in enumerate(ids):
                mom[(r['cohort'],r['patient'])].setdefault(g,[]).append((gram[j],rhs[j],yy[j]))
                vy=np.var(y[:,j]);vx=np.var(x[:,j,9]);corr=float(np.corrcoef(x[:,j,9],y[:,j])[0,1]) if min(vy,vx)>1e-12 else None
                molecular.append(dict(program=g,cohort=r['cohort'],patient=r['patient'],section=r['section'],input_readout_corr=corr,input_moran=mi[j],readout_moran=mr[j],status='COMPUTED' if corr is not None else 'CONSTANT_SCORE'))
    return mom,disc,invalid,molecular


def replay(go=False,mode='raw'):
    from concurrent.futures import ProcessPoolExecutor
    defs=definitions(go);genes=json.load(open(OUT/'genes.json'))
    prefix='go' if go else 'six_'+mode
    rows=sections() if mode=='raw' else sections(['ST-CRC','USZ'])
    mom=defaultdict(dict);disc=[];invalid=[];molecular=[]
    tasks=[(r,defs,genes,mode) for r in rows]
    pool=ProcessPoolExecutor(max_workers=8) if go else None
    results=pool.map(_section_replay,tasks) if pool else map(_section_replay,tasks)
    try:
        for si,(sm,sd,sv,sr) in enumerate(results):
            for key,programs in sm.items():
                for g,v in programs.items():mom[key].setdefault(g,[]).extend(v)
            disc.extend(sd);invalid.extend(sv);molecular.extend(sr)
            print(prefix,si+1,len(rows),rows[si]['key'],flush=True)
            dump(prefix+'_progress.json',dict(status='RUNNING',completed_sections=si+1,total_sections=len(rows),last_section=rows[si]['key']))
    finally:
        if pool:pool.shutdown(wait=True,cancel_futures=True)
    result=[]
    for cohort in ['ST-CRC','USZ']:
        pats=sorted(p for c,p in mom if c==cohort)
        for g in defs:
            available=[p for p in pats if g in mom[(cohort,p)]]
            if len(available)<3:continue
            gs=[];bs=[];ys=[]
            for p in available:
                vals=mom[(cohort,p)][g];gs.append(np.mean([v[0] for v in vals],0));bs.append(np.mean([v[1] for v in vals],0));ys.append(np.mean([v[2] for v in vals]))
            gs=np.array(gs);bs=np.array(bs);ys=np.array(ys)
            for k,p in enumerate(available):
                tr=np.arange(len(available))!=k;gt=gs[tr].mean(0);bt=bs[tr].mean(0)
                mse=[]
                for n in [9,10,11]:
                    beta=fit_from_moments(gt[:n,:n],bt[:n]);val=ys[k]-2*beta@bs[k,:n]+beta@gs[k,:n,:n]@beta;mse.append(float(max(val,0)))
                variance=ys[k]-bs[k,0]**2;ok=variance>1e-12 and mse[0]>1e-12
                result.append(dict(program=g,cohort=cohort,patient=p,status='COMPUTED' if ok else 'NOT_TESTABLE_CONSTANT_TARGET',mse_M0=mse[0],mse_M1=mse[1],mse_M2=mse[2],delta_M1=100*(mse[0]-mse[1])/mse[0] if ok else None,delta_M2=100*(mse[0]-mse[2])/mse[0] if ok else None,n_sections=len(mom[(cohort,p)][g])))
    table(prefix+'_prediction.tsv',result);table(prefix+'_molecular.tsv',molecular);table(prefix+'_not_testable.tsv',invalid)
    if disc:
        grouped=defaultdict(list)
        for r in disc:grouped[(r['program'],r['patient'])].append(r)
        table(prefix+'_discovery.tsv',[dict(program=g,patient=p,n_sections=len(v),**{key:float(np.median([r[key] for r in v])) for key in ['raw','Q','QC']}) for (g,p),v in grouped.items()])
    summary=[];result_groups=defaultdict(list)
    for r in result:
        if r['status']=='COMPUTED':result_groups[(r['cohort'],r['program'])].append(r)
    for cohort in ['ST-CRC','USZ']:
        for g in defs:
            v=result_groups[(cohort,g)]
            a=[r['delta_M1'] for r in v];b=[r['delta_M2'] for r in v];ci=patient_interval(a,draws=1000)
            summary.append(dict(program=g,cohort=cohort,n_patients=len(v),status='COMPUTED' if v else 'NOT_TESTABLE',median_delta_M1=float(np.median(a)) if a else None,median_delta_M2=float(np.median(b)) if b else None,patient_bootstrap_low=ci[0],patient_bootstrap_high=ci[1]))
    table(prefix+'_summary.tsv',summary)
    dump(prefix+'_receipt.json',dict(status='COMPLETED',programs=len(defs),prediction_rows=len(result),not_testable_sections=len(invalid),discovery_sections=len({r['section'] for r in disc}),normalization=mode,claims='Descriptive patient LOPO; not spatial-null-calibrated discovery'))
    dump(prefix+'_progress.json',dict(status='COMPLETED',completed_sections=len(rows),total_sections=len(rows)))

def marker_replay(mode='raw'):
    defs=definitions();genes=json.load(open(OUT/'genes.json'));gx={g:i for i,g in enumerate(genes)};mom=defaultdict(list)
    for r in sections(['ST-CRC','USZ']):
        counts,xy,_=load(r)
        for g,d in defs.items():
            try:f=feature_block(counts,genes,d['input_genes'],d['readout_genes'],axes(),mode)
            except NotTestable:continue
            scale=(1e4/np.maximum(f['background_library'],1))[:,None] if mode=='depth_normalized' else 1
            markers=np.log1p(counts[:,[gx[s] for s in d['input_genes']]]*scale)
            X=np.c_[f['X'][:,:9],markers];y=f['y']
            mom[(r['cohort'],g,r['patient'])].append((X.T@X/len(X),X.T@y/len(X),np.mean(y*y)))
    out=[]
    for c in ['ST-CRC','USZ']:
        for g,d in defs.items():
            pats=sorted(p for cc,gg,p in mom if cc==c and gg==g)
            if len(pats)<3:continue
            vals=[tuple(np.mean([r[j] for r in mom[(c,g,p)]],axis=0) for j in range(3)) for p in pats]
            for k,p in enumerate(pats):
                train=[v for j,v in enumerate(vals) if j!=k];gt=np.mean([v[0] for v in train],0);bt=np.mean([v[1] for v in train],0);yt=np.mean([v[2] for v in train]);gh,bh,yh=vals[k]
                candidates=[]
                for j in range(len(d['input_genes'])):
                    ix=list(range(9))+[9+j];beta=fit_from_moments(gt[np.ix_(ix,ix)],bt[ix]);risk=yt-2*beta@bt[ix]+beta@gt[np.ix_(ix,ix)]@beta
                    candidates.append((risk,j,ix,beta))
                _,j,ix,beta=min(candidates,key=lambda r:(r[0],r[1]));mse=max(0,float(yh-2*beta@bh[ix]+beta@gh[np.ix_(ix,ix)]@beta))
                out.append(dict(cohort=c,program=g,patient=p,selected_marker=d['input_genes'][j],heldout_mse=mse,selection='training_patients_only'))
    table('six_'+mode+'_single_marker.tsv',out)

def pc_replay():
    from r16 import axis_factory as F
    from r16.recovery import pc_matching as M
    artifact=json.load(open(ROOT/'infra/r16/axis_factory_20260915.json'));cand=artifact['candidates']
    z=np.load(ROOT/'infra/r16/axis_factory_loadings_20260915.npz');L=z['loadings'].astype(float);genes=list(z['genes']);gx={g:i for i,g in enumerate(genes)}
    groups=M.match_groups(L,cutoff=C.MATCH_COSINE)
    historical=json.load(open(OLD/'audit/pc_matching_reanalysis.json'))
    selected=[r['new_group'] for r in historical['groups'] if r['grade_new_rule']=='EXPLORATORY_REPRODUCED']
    bystem=defaultdict(list)
    for g in selected:
        mem=np.flatnonzero(groups==g);mean,_,_=M.oriented_mean(L,list(mem))
        for i in mem:bystem[cand[i]['stem']].append((g,i,1 if L[i]@mean>=0 else -1))
    rows=[];cache=ROOT/'infra/r16/census_cache_hvg10k/sections'
    for stem,items in bystem.items():
        counts,bc,gs,xy,meta=C.load_section(cache,stem)
        assert list(gs)==genes
        X=C.log1p_zscore(counts);log=np.log1p(counts.toarray()).astype(float);W=F.spatial_mean_weights(xy,15);neighbor=F.zscore_columns(W@log)
        comp=np.column_stack([np.ones(len(X))]+[X[:,[gx[s] for s in ss if s in gx]].mean(1) for ss in axes().values()])
        for g in selected:
            members=[(i,sign) for gg,i,sign in items if gg==g]
            if not members:continue
            scores=[]
            for i,sign in members:
                lam=cand[i]['lambda'];score=(np.sqrt(1-lam)*(X@L[i])+np.sqrt(lam)*(neighbor@L[i]))*sign
                scores.append((score-score.mean())/max(score.std(),1e-12))
            score=np.mean(scores,axis=0);prediction=comp@fit_linear(comp,score);r2=1-np.mean((score-prediction)**2)/max(score.var(),1e-12)
            rows.append(dict(group=int(g),patient=meta['patient_id'],section=stem,n_members=len(members),spot_joint_composition_R2=float(r2)))
        print('pc',stem,flush=True)
    table('pc_actual_spot_regression.tsv',rows)
    summarize_pc(rows,historical)

def summarize_pc(rows,historical):
    selected=[r['new_group'] for r in historical['groups'] if r['grade_new_rule']=='EXPLORATORY_REPRODUCED']
    assert len({(r['group'],r['section']) for r in rows})==len(rows)
    summaries=[]
    for g in selected:
        expected=next(r['n_members'] for r in historical['groups'] if r['new_group']==g)
        assert sum(r['n_members'] for r in rows if r['group']==g)==expected
        vals=defaultdict(list)
        for r in rows:
            if r['group']==g:vals[r['patient']].append(r['spot_joint_composition_R2'])
        v=[np.median(vv) for vv in vals.values()];summaries.append(dict(group=int(g),n_patients=len(v),median_patient_spot_R2=float(np.median(v)),n_patients_R2_ge_049=int(sum(x>=.49 for x in v)),interpretation='Composition association diagnostic, not evidence of spatial novelty'))
    dump('pc_receipt.json',dict(status='COMPLETED',reproduced_groups=len(selected),historical_loading_projection_new_groups=historical['corrected']['new_reproduced_joint_R2'],groups=summaries,claim='Two loading groups reproduce; loading projection cannot classify them as biologically novel'))

def reconstructed_scores(counts,xy,genes,sets,hidden):
    """Score only visible molecules, then interpolate into the hidden window."""
    gx={g:i for i,g in enumerate(genes)};vis=~hidden
    if vis.sum()<20 or not hidden.any():raise NotTestable('INSUFFICIENT_VISIBLE_OR_HIDDEN')
    vals=np.column_stack([np.log1p(counts[vis][:,[gx[g] for g in gs if g in gx]]).mean(1) for gs in sets])
    if not np.isfinite(vals).all():raise NotTestable('MISSING_SCORE_GENES')
    return neighbor_values(xy[vis],vals,xy[hidden]).mean(0)

def auc(y,score):
    y=np.asarray(y);score=np.asarray(score);a=y==1;b=y==0
    if not a.any() or not b.any():return None
    return float((rankdata(score)[a].sum()-a.sum()*(a.sum()+1)/2)/(a.sum()*b.sum()))

def mask_replay():
    defs=definitions();genes=json.load(open(OUT/'genes.json'));ax=axes();field=[];struct=[];leaks=[];fine=[];focus_meta=[]
    for r in sections(['ST-CRC','USZ']):
        counts,xy,labels=load(r);wins=generate_windows(xy,12,8,20260921)
        base={k:r[k] for k in ['cohort','patient','section']}
        for g,d in defs.items():
            try:f=feature_block(counts,genes,d['input_genes'],d['readout_genes'],ax)
            except NotTestable as e:
                field.append(dict(**base,program=g,window=-1,model='ALL',mse=None,status=str(e)));continue
            probed=False
            for wi,win in enumerate(wins):
                hidden=window_mask(xy,win)
                try:pred=masked_predict(counts,xy,genes,d['input_genes'],d['readout_genes'],ax,hidden)
                except NotTestable as e:
                    field.append(dict(**base,program=g,window=wi,model='ALL',mse=None,status=str(e)));continue
                for model,yhat in pred.items():field.append(dict(**base,program=g,window=wi,model=model,mse=float(np.mean((f['y'][hidden]-yhat)**2)),status='COMPUTED'))
                if not probed:
                    changed=np.array(counts);changed[hidden]=1e6
                    altered=masked_predict(changed,xy,genes,d['input_genes'],d['readout_genes'],ax,hidden)
                    delta=max(float(np.max(abs(pred[k]-altered[k]))) for k in pred)
                    assert delta==0,(r['section'],g,delta)
                    leaks.append(dict(**base,program=g,window=wi,max_prediction_change=delta));probed=True
        if r['cohort']=='USZ':
            sets=[defs['P-mhc2']['input_genes'],ax['B']]
            for wi,win in enumerate(wins):
                m=window_mask(xy,win);state=classify_window(labels[m,0])
                try:sc=reconstructed_scores(counts,xy,genes,sets,m)
                except NotTestable:sc=[None,None]
                struct.append(dict(**base,window=wi,label=state,n_spots=int(m.sum()),n_known=int(np.isfinite(labels[m,0]).sum()),program_score=sc[0],B_score=sc[1]))
            foci=tls_components(xy,labels[:,0]);focus_meta.append(dict(**base,n_hex_foci=len(foci)))
            tree=cKDTree(xy)
            for half in [3,5,8]:
                neg=[]
                # Enumerate all on-tissue centers, requiring wholly known non-TLS windows.
                for ci,center in enumerate(xy):
                    inds=tree.query_ball_point(center,half*np.sqrt(2)+1e-8)
                    inds=np.asarray(inds,int);inds=inds[(abs(xy[inds]-center)<=half).all(1)]
                    if len(inds)>=5 and classify_window(labels[inds,0])=='NEGATIVE':neg.append(ci)
                for fi,fidx in enumerate(foci):
                    center=xy[fidx].mean(0);m=window_mask(xy,dict(center=center,half=half))
                    status=classify_window(labels[m,0]);eligible=status=='POSITIVE' and len(neg)>0
                    if not eligible:
                        fine.append(dict(**base,focus=fi,half=half,seed=None,status='NO_KNOWN_NEGATIVE' if not neg else status,positive_program=None,negative_program=None,positive_B=None,negative_B=None,delta_auc=None));continue
                    pos=reconstructed_scores(counts,xy,genes,sets,m)
                    # Seeds select controls, never increase the biological sample count.
                    chosen=set()
                    for seed in [11,22,33]:
                        pool=[c for c in neg if c not in chosen]
                        if not pool:break
                        ci=int(np.random.default_rng(seed+fi).choice(pool));chosen.add(ci)
                        nm=window_mask(xy,dict(center=xy[ci],half=half));ns=reconstructed_scores(counts,xy,genes,sets,nm)
                        pa=float(pos[0]>ns[0])+.5*float(pos[0]==ns[0]);ba=float(pos[1]>ns[1])+.5*float(pos[1]==ns[1])
                        fine.append(dict(**base,focus=fi,half=half,seed=seed,status='COMPUTED',positive_program=float(pos[0]),negative_program=float(ns[0]),positive_B=float(pos[1]),negative_B=float(ns[1]),delta_auc=pa-ba))
        print('mask',r['key'],flush=True)
    table('masked_field_windows.tsv',field);table('masked_structure_windows.tsv',struct);table('l009_focus_controls.tsv',fine);table('hex_foci.tsv',focus_meta)
    dump('masked_prediction_intervention.json',dict(status='PASS',interventions=leaks,contract='All hidden counts changed to 1e6; complete M0/M1/M2/NN/KNN8 predictions identical'))
    # Average windows within section, then sections within patient, then infer across patients.
    persection=defaultdict(list)
    for r in field:
        if r['status']=='COMPUTED':persection[(r['cohort'],r['program'],r['patient'],r['section'],r['model'])].append(r['mse'])
    perpatient=defaultdict(list)
    for (c,g,p,s,m),v in persection.items():perpatient[(c,g,p,m)].append(float(np.mean(v)))
    patient=[]
    for (c,g,p,m),v in perpatient.items():patient.append(dict(cohort=c,program=g,patient=p,model=m,mse=float(np.mean(v))))
    table('masked_field_patients.tsv',patient)
    comparison=[]
    for c in ['ST-CRC','USZ']:
        for g in defs:
            for model in ['M1','M2']:
                for baseline in ['M0','NN','KNN8']:
                    a={p:np.mean(v) for (cc,gg,p,m),v in perpatient.items() if cc==c and gg==g and m==model}
                    b={p:np.mean(v) for (cc,gg,p,m),v in perpatient.items() if cc==c and gg==g and m==baseline}
                    dif=[100*(b[p]-a[p])/b[p] for p in a.keys()&b.keys() if b[p]>1e-12];ci=patient_interval(dif)
                    comparison.append(dict(cohort=c,program=g,model=model,baseline=baseline,n_patients=len(dif),median_improvement=float(np.median(dif)) if dif else None,ci_low=ci[0],ci_high=ci[1]))
    table('masked_field_summary.tsv',comparison)
    structure_pat=[]
    for p in sorted({r['patient'] for r in struct}):
        v=[r for r in struct if r['patient']==p and r['label'] in ['POSITIVE','NEGATIVE'] and r['program_score'] is not None]
        y=[int(r['label']=='POSITIVE') for r in v];pa=auc(y,[r['program_score'] for r in v]);ba=auc(y,[r['B_score'] for r in v])
        structure_pat.append(dict(patient=p,n_positive=sum(y),n_negative=len(y)-sum(y),program_auc=pa,B_auc=ba,status='COMPUTED' if pa is not None else 'NOT_TESTABLE'))
    table('masked_structure_patients.tsv',structure_pat)
    focus=defaultdict(list)
    for r in fine:
        if r['status']=='COMPUTED':focus[(r['patient'],r['focus'])].append(r['delta_auc'])
    pats=defaultdict(list)
    for (p,f),v in focus.items():pats[p].append(float(np.mean(v)))
    table('l009_patients.tsv',[dict(patient=p,n_foci=len(v),mean_delta_auc=float(np.mean(v))) for p,v in pats.items()])
    vals=[np.mean(v) for v in pats.values()]
    dump('mask_receipt.json',dict(status='COMPLETED',label_unknown_preserved=True,field_interventions=len(leaks),hex_foci=sum(r['n_hex_foci'] for r in focus_meta),l009_n_patients=len(vals),l009_median_delta=float(np.median(vals)) if vals else None,l009_patient_interval=patient_interval(vals),structure_claim='Targeted focus controls; not blind localization or cross-cohort transfer'))

def spline_replay():
    old=readtable(OLD/'programs_go/go_summary_full.tsv')
    oldhits={r['rep_id'] for r in old if float(r['STCRC_dM1'])>20 and float(r['USZ_dM1'])>20}
    new=readtable(OUT/'go_summary.tsv');hitsets=[]
    for c in ['ST-CRC','USZ']:hitsets.append({r['program'] for r in new if r['cohort']==c and r['median_delta_M1'] and float(r['median_delta_M1'])>20 and int(r['n_patients'])==(7 if c=='ST-CRC' else 8)})
    newhits=hitsets[0]&hitsets[1];selected=oldhits|newhits;gd=definitions(True);defs={**definitions(),**{g:gd[g] for g in sorted(selected)}}
    dump('spline_selection.json',dict(original_hits=sorted(oldhits),corrected_hits=sorted(newhits),union=sorted(selected),selection='Exploratory reanalysis; selected on molecular prediction, no independent pathway-discovery claim'))
    genes=json.load(open(OUT/'genes.json'));gx={g:i for i,g in enumerate(genes)};res=[];injections=[]
    for r in sections(['USZ']):
        counts,xy,labels=load(r);tum=labels[:,1]==1;distance=signed_tls_distance(xy,labels[:,0]);d=distance/np.quantile(distance[distance>0],.99)
        t=d[tum];grid=np.linspace(*np.quantile(t,[.05,.95]),101);basis=NaturalSpline(t);B=basis(t);BG=basis(grid)
        injected=1+2*t;beta=np.linalg.lstsq(B,injected,rcond=None)[0];pred=B@beta
        injections.append(dict(section=r['section'],linear_R2=float(1-np.sum((injected-pred)**2)/np.sum((injected-injected.mean())**2)),shape=classify_curve(BG@beta,injected.std())))
        ranks=rankdata(counts,axis=1,method='average').astype(np.float32)
        tx=xy[tum];order=np.argsort(tx[:,0],kind='stable');fold=np.empty(len(t),int)
        for k,indices in enumerate(np.array_split(order,5)):fold[indices]=k
        splits=[]
        for k in range(5):
            te=fold==k;tr=(fold!=k)&(cKDTree(tx[te]).query(tx)[0]>1.01)
            if tr.sum()<30 or te.sum()<5:continue
            nb=NaturalSpline(t[tr]);splits.append((tr,te,nb(t[tr]),nb(t[te])))
        for ids,x,y,bad in raw_blocks(counts,xy,genes,defs):
            for g,status in bad:
                # Raw/rank descriptive curves still exist even if composition separation fails.
                cols=[gx[s] for s in defs[g]['input_genes'] if s in gx]
                for route in ['raw','rank','resid','rank_resid']:
                    if status=='MISSING_PROGRAM_GENES' or route in ['resid','rank_resid'] or not cols:
                        res.append(dict(program=g,section=r['section'],patient=r['patient'],route=route,status=status,shape='NOT_TESTABLE',amplitude_sd=None,R2=None,blocked_cv_improvement=None,p_status='NOT_CALIBRATED'));continue
                    score=(np.log1p(counts[:,cols]).mean(1) if route=='raw' else (ranks[:,cols].sum(1)-len(cols)*(len(cols)+1)/2)/(len(cols)*(len(genes)-len(cols))))[tum]
                    coef=np.linalg.lstsq(B,score,rcond=None)[0];curve=BG@coef
                    res.append(dict(program=g,section=r['section'],patient=r['patient'],route=route,status='DESCRIPTIVE_ONLY',shape=classify_curve(curve,score.std()),amplitude_sd=float(np.ptp(curve)/max(score.std(),1e-12)),R2=float(1-np.mean((score-B@coef)**2)/max(score.var(),1e-12)),blocked_cv_improvement=None,p_status='NOT_CALIBRATED'))
            for j,g in enumerate(ids):
                cols=[gx[s] for s in defs[g]['input_genes']];rk=(ranks[:,cols].sum(1)-len(cols)*(len(cols)+1)/2)/(len(cols)*(len(genes)-len(cols)))
                Q=x[tum,j,:9]
                for route in ['raw','rank','resid','rank_resid']:
                    source=(rk if 'rank' in route else x[:,j,9])[tum].astype(float)
                    residual='resid' in route;score=source-Q@fit_linear(Q,source) if residual else source
                    coef=np.linalg.lstsq(B,score,rcond=None)[0];curve=BG@coef;sd=score.std();err0=[];err1=[]
                    for tr,te,bt,be in splits:
                        q=Q if residual else np.ones((len(t),1));b0=fit_linear(q[tr],source[tr]);b1=fit_linear(np.c_[q[tr],bt],source[tr])
                        err0.append(np.mean((source[te]-q[te]@b0)**2));err1.append(np.mean((source[te]-np.c_[q[te],be]@b1)**2))
                    cv=100*(np.mean(err0)-np.mean(err1))/np.mean(err0) if err0 and np.mean(err0)>1e-12 else None
                    res.append(dict(program=g,section=r['section'],patient=r['patient'],route=route,status='COMPUTED' if sd>1e-12 else 'NOT_TESTABLE_CONSTANT_SCORE',shape=classify_curve(curve,sd),amplitude_sd=float(np.ptp(curve)/sd) if sd>1e-12 else None,R2=float(1-np.mean((score-B@coef)**2)/score.var()) if sd>1e-12 else None,blocked_cv_improvement=float(cv) if cv is not None else None,p_status='NOT_CALIBRATED'))
        print('spline',r['key'],len(defs),flush=True)
    table('spline_results.tsv',res);table('spline_linear_injection.tsv',injections)
    assert all(v['linear_R2']>1-1e-10 and v['shape'] in ['low-to-high','increasing'] for v in injections),injections
    dump('spline_receipt.json',dict(status='COMPLETED',original_go=len(oldhits),new_go=len(newhits),union_go=len(selected),programs=len(defs),rows=len(res),spatial_p='NOT_CALIBRATED',injected_linear_field='PASS',cv='five contiguous x-stripes; >1 pitch buffer; train-only knots and nuisance; descriptive'))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('phase',choices=['prepare','six','go','mask','spline','pc']);args=ap.parse_args()
    if shutil.disk_usage(ROOT).free < 1_200_000_000_000:
        raise RuntimeError('Storage reserve below 1.2 TB')
    if args.phase=='prepare':prepare()
    elif args.phase=='six':
        replay();replay(mode='depth_normalized');marker_replay();marker_replay(mode='depth_normalized')
    elif args.phase=='go':replay(go=True)
    elif args.phase=='mask':mask_replay()
    elif args.phase=='spline':spline_replay()
    elif args.phase=='pc':pc_replay()

if __name__=='__main__':main()
