#!/usr/bin/env python3
"""R-17: metadata qualification, visible-only scoring and bounded spatial calibration.

No GPU; no downloads. All patient identifiers require source metadata evidence.
"""
from pathlib import Path
import argparse, csv, gzip, hashlib, importlib.util, json, time
import h5py
import numpy as np
from scipy import sparse, io, optimize, spatial, stats, ndimage

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'infra/tls_peripheral_20261002'
BASE = ROOT / 'data/other_sources/zenodo_usz_tls_visium/TLS_VISIUM_USZ'
RAW = ROOT / 'data/GEO/GSE175540/raw'
LENGTHS = np.array([100., 300., 800., 1600.])
LAGS = np.array([0., 150., 250., 400., 600., 900., 1300., 1800., 2500., 4000., 8000., np.inf])

def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT/'scripts'/f'{name}.py')
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m

def table(path):
    with open(path) as f: return list(csv.DictReader(f, delimiter='\t'))

def write_table(path, rows):
    if not rows: return
    with open(path,'w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def save_json(path,obj):
    path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def decode(x): return [v.decode() if isinstance(v,bytes) else str(v) for v in x]

def hash_file(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()

def source_metadata(cohort,sid):
    """Read only barcodes, gene names, coordinates, labels; never expression."""
    if cohort=='USZ':
        alias=sid.removeprefix('usz-'); directory=BASE/'10x_Visium'/alias/'filtered_feature_bc_matrix'
        with gzip.open(directory/'barcodes.tsv.gz','rt') as f:barcodes=[l.strip() for l in f]
        with gzip.open(directory/'features.tsv.gz','rt') as f:features=[l.rstrip().split('\t') for l in f]
        genes=[r[1] for r in features]; gene_keep=np.array([r[2]=='Gene Expression' for r in features])
        annotation=BASE/'h5ad_preprocessed'/f'{alias}.h5ad'
        with h5py.File(annotation) as f:
            bc=decode(f['obs/_index'][:]); cat=decode(f['obs/ground_truth/categories'][:]);code=f['obs/ground_truth/codes'][:]
            amap=dict(zip(bc,[cat[c] if c>=0 else '' for c in code]))
        positions=BASE/'10x_Visium'/alias/'spatial/tissue_positions.csv'
        with open(positions) as f:rows=list(csv.DictReader(f))
        pos={r['barcode']:(int(r['array_row']),int(r['array_col'])) for r in rows if r['in_tissue']=='1'}
        sources=[directory/'matrix.mtx.gz',directory/'features.tsv.gz',directory/'barcodes.tsv.gz',positions,annotation]
        Xpath=directory/'matrix.mtx.gz'
    else:
        gsm=sid.removeprefix('gse175540-');Xpath=RAW/f'{gsm}_filtered_feature_bc_matrix.h5'
        with h5py.File(Xpath) as f:
            barcodes=decode(f['matrix/barcodes'][:]);genes=decode(f['matrix/features/name'][:])
            types=decode(f['matrix/features/feature_type'][:]); gene_keep=np.array([v=='Gene Expression' for v in types])
        annotation=RAW/f'{gsm}_TLS_annotation.csv.gz'
        with gzip.open(annotation,'rt',encoding='utf-8-sig') as f:
            amap={r['Barcode'].strip():r['TLS_2_cat'].strip() for r in csv.DictReader(f)}
        positions=RAW/f'{gsm}_tissue_positions_list.csv.gz'
        with gzip.open(positions,'rt') as f:rows=list(csv.reader(f))
        pos={r[0]:(int(r[2]),int(r[3])) for r in rows if r[1]=='1'}
        sources=[Xpath,positions,annotation]
    order=np.array([i for i,b in enumerate(barcodes) if b in pos],dtype=int)
    coords=np.array([pos[barcodes[i]] for i in order],dtype=int)
    labels=np.array([amap.get(barcodes[i],'') for i in order])
    xy=module('tls_pool_expand').hex_xy(coords[:,0],coords[:,1],100)
    lost_tls=sum(v=='TLS' and k not in {barcodes[i] for i in order} for k,v in amap.items())
    return dict(order=order,coords=coords,xy=xy,labels=labels,genes=np.array(genes),gene_keep=gene_keep,
                Xpath=Xpath,sources=sources,lost_tls=lost_tls,n_matrix_rows=len(barcodes))

def geometry(xy,labels,radius):
    core=labels=='TLS'
    d=spatial.cKDTree(xy[core]).query(xy)[0] if core.any() else np.full(len(xy),np.inf)
    visible=d>radius
    near=visible&(d>=200)&(d<400);far=visible&(d>=700)&(d<1000)
    return d,visible,near,far

def prepare():
    OUT.mkdir(exist_ok=True)
    if (OUT/'contract.json').exists(): raise FileExistsError('frozen run already exists')
    (OUT/'geometry').mkdir(exist_ok=True)
    diag=[r for r in table(ROOT/'infra/tls_pool_expand_20260930/pool_diagnostics.tsv') if r['status']=='ok' and r['cohort'] in ['USZ','GSE175540']]
    assert len(diag)==26
    h=module('tls_field_hallmark').read_gmt(ROOT/'data/geneset/h.all.v2026.1.Hs.symbols.gmt')
    proxy=json.loads((ROOT/'infra/r04/marker_proxy_combined.json').read_text())
    h['B_AXIS']=sorted(set(proxy['classes']['B']['voted_genes']))
    assert len(h)==51
    rows=[];identities=[];files={};coverage=[]
    for r in diag:
        sid,cohort=r['section_id'],r['cohort'];m=source_metadata(cohort,sid)
        for f in m['sources']:
            files[str(f.relative_to(ROOT))]=dict(size=f.stat().st_size,sha256=hash_file(f))
        present=set(m['genes'][m['gene_keep']])
        for name,genes in h.items():
            overlap=sorted(set(genes)&present)
            coverage.append(dict(section_id=sid,set_id=name,n_defined=len(set(genes)),n_present=len(overlap),fraction=len(overlap)/len(set(genes)),missing=','.join(sorted(set(genes)-present))))
        np.savez_compressed(OUT/'geometry'/f'{sid}.npz',xy=m['xy'],labels=m['labels'],coords=m['coords'],order=m['order'])
        identities.append(dict(section_id=sid,cohort=cohort,patient_id='',patient_status='UNKNOWN_NO_EXPLICIT_CROSSWALK',
            evidence='SAMPLE only; source record gives 3 kidney/5 lung tumors' if cohort=='USZ' else 'GEO specimen; atlas Table S2 corroborates 24 patients but no GSM-to-R_P mapping',
            aggregation_unit='section/specimen; not asserted patient'))
        for radius in [130,230]:
            d,v,n,f=geometry(m['xy'],m['labels'],radius)
            reason=[]
            if sum(m['labels']=='TLS')<20:reason.append('TLS_LT20')
            if n.sum()<30:reason.append('NEAR_LT30')
            if f.sum()<30:reason.append('FAR_LT30')
            if m['lost_tls']:reason.append('KNOWN_TLS_MISSING_FROM_MATRIX')
            rows.append(dict(section_id=sid,cohort=cohort,mask_um=radius,n_spots=len(d),n_tls=int(sum(m['labels']=='TLS')),n_visible=int(v.sum()),n_near=int(n.sum()),n_far=int(f.sum()),n_unknown=int(sum(np.isin(m['labels'],['','UNASSIGNED']))),lost_tls=m['lost_tls'],status='ELIGIBLE' if not reason else 'NOT_TESTABLE',reason=';'.join(reason)))
        print('qualified',sid,flush=True)
    write_table(OUT/'eligibility.tsv',rows);write_table(OUT/'identity.tsv',identities);write_table(OUT/'coverage.tsv',coverage)
    save_json(OUT/'source_hashes.json',files);save_json(OUT/'panel.json',h)
    contract=dict(status='FROZEN_BEFORE_NEW_SCORING',decision='D-169',sections=[{k:r[k] for k in ['section_id','cohort']} for r in diag],
       mask_um=[130,230],near_um=[200,400],far_um=[700,1000],min_spots=30,min_TLS=20,
       unknown_policy='Unknown is not TLS-negative. Distance refers only to labelled TLS; missing known TLS excludes geometry.',
       score='visible-only per-spot AUCell top 5% after raw count sum for duplicated gene symbols',
       coverage_min=.8,min_genes=10,panel_sha256=hash_file(OUT/'panel.json'),
       statistic='near mean minus far mean, raw and divided by visible sample SD; no composition/B/trend residualization',
       useful_effect_sd=.2,patient_inference='BLOCKED_EXPLICIT_IDENTITY_CROSSWALK',
       spatial_candidate='NNLS empirical variogram: nugget plus Gaussian correlations ell=100,300,800,1600um; plug-in two-sided normal statistic',
       calibration=dict(seed=20261002,replicates=400,power_replicates=200,alpha=.05,
          families=['white','short150','long800','mixture','non_gaussian','random_trend'],
          generators='Gaussian-filtered iid rectangular lattice padded 4 kernel SD; ell/sqrt(2) convolution; no thinning of tissue spots',
          injected_effects_sd=[.2,.5],injection_widths_um=[200,600],
          gate='Every eligible section/mask/family Wilson upper 95% FPR <=0.10. Power lower95%>=0.80 at 0.5SD for both widths; report 0.2SD separately.',
          fit_pairs=16000,lag_edges_um=[float(x) for x in LAGS[:-1]],max_fit_lag='unbounded final bin',
          real_gate='Simulation gate + per-readout variogram fit relative RMSE <=0.25; even then model-conditional only, not proven empirical exchangeability'),
       multiplicity='Holm over51 readouts per eligible section/mask, only if formal gates pass; no patient inference without IDs.',
       stop='Failed calibration keeps all real scores descriptive; no new null tuning, radius search, inversion or field-absence claim.',
       figure_contract='Python existing workflow; evidence=all-section mask grid, complete readout effects, null calibration. SVG/PDF editable text plus PNG; no unsupported patient CI.')
    save_json(OUT/'contract.json',contract)
    print('FROZEN',hash_file(OUT/'contract.json'),flush=True)


def counts(cohort,sid,m):
    if cohort=='USZ':
        with gzip.open(m['Xpath'],'rb') as f:X=io.mmread(f).T.tocsr()
    else:X=module('tls_pool_expand').read_10x_h5_matrix(m['Xpath'])[0]
    X=X[m['order']][:,m['gene_keep']].astype(float)
    genes=m['genes'][m['gene_keep']];uniq,inverse=np.unique(genes,return_inverse=True)
    mapping=sparse.csr_matrix((np.ones(len(genes)),(np.arange(len(genes)),inverse)),shape=(len(genes),len(uniq)))
    return (X@mapping).tocsr(),uniq


def visible_scores(X,genes,visible,panel):
    idx={g:i for i,g in enumerate(genes)}
    members=[np.array([idx[g] for g in sorted(set(gl)) if g in idx],dtype=int) for gl in panel.values()]
    # Slicing precedes every score/QC calculation. No cross-spot normalization.
    return module('gobp_halo_screen').score_sets(X[visible],members)[0]


def fit_setup(xy,weight,seed):
    rng=np.random.default_rng(seed); n=len(xy)
    i=rng.integers(0,n,20000);j=rng.integers(0,n,20000);keep=i!=j;i=i[keep][:16000];j=j[keep][:16000]
    dist=np.linalg.norm(xy[i]-xy[j],axis=1);band=np.digitize(dist,LAGS)-1
    parts=[np.flatnonzero(band==k) for k in range(len(LAGS)-1)];parts=[p for p in parts if len(p)>=40]
    kernel=1-np.exp(-dist[:,None]**2/(2*LENGTHS[None,:]**2))
    design=np.array([np.r_[1,kernel[p].mean(axis=0)] for p in parts])
    nz=np.flatnonzero(weight);ww=weight[nz]; dd=spatial.distance.cdist(xy[nz],xy[nz])**2
    variance=np.r_[weight@weight,[ww@np.exp(-dd/(2*l*l))@ww for l in LENGTHS]]
    return i,j,parts,design,np.maximum(variance,0)


def fit_null(Y,setup,weight):
    if Y.ndim==1:Y=Y[:,None]
    i,j,parts,A,V=setup
    dif=(Y[i]-Y[j])**2/2
    gamma=np.array([dif[p].mean(axis=0) for p in parts])
    coeff=np.column_stack([optimize.nnls(A,gamma[:,k])[0] for k in range(Y.shape[1])])
    var=V@coeff;err=np.sqrt(np.mean((A@coeff-gamma)**2,axis=0))/(gamma.mean(axis=0)+1e-15)
    observed=weight@Y;se=np.sqrt(var)
    p=np.full(Y.shape[1],np.nan);ok=(se>1e-12)&np.isfinite(se)
    p[ok]=2*stats.norm.sf(np.abs(observed[ok]/se[ok]))
    return p,err


def gaussian_lattice(coords,ell,n_draws,rng):
    """Full support, padded filtered white noise; convolution variance normalized analytically."""
    c=coords-coords.min(axis=0);r,cx=c[:,0],c[:,1]
    if ell==0:return rng.standard_normal((len(coords),n_draws))
    sy=ell/np.sqrt(2)/(100*np.sqrt(3)/2);sx=ell/np.sqrt(2)/50
    py=int(4*sy+.5);px=int(4*sx+.5)
    yy=np.arange(-py,py+1);xx=np.arange(-px,px+1)
    ky=np.exp(-yy**2/(2*sy**2));ky/=ky.sum();kx=np.exp(-xx**2/(2*sx**2));kx/=kx.sum()
    norm=np.sqrt(np.sum(ky**2)*np.sum(kx**2))
    out=np.empty((len(coords),n_draws))
    for a in range(0,n_draws,25):
        count=min(25,n_draws-a)
        z=rng.standard_normal((int(r.max())+1+2*py,int(cx.max())+1+2*px,count))
        z=ndimage.gaussian_filter(z,sigma=(sy,sx,0),mode='constant',truncate=4)
        out[:,a:a+count]=z[r+py,cx+px]/norm
    return out


def wilson(k,n):
    z=stats.norm.ppf(.975);p=k/n;den=1+z*z/n
    centre=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return float(centre-half),float(centre+half)


def calibrate():
    contract=json.loads((OUT/'contract.json').read_text());elig=table(OUT/'eligibility.tsv')
    target=OUT/'calibration.tsv'
    if target.exists():raise FileExistsError(target)
    rng=np.random.default_rng(20261002);records=[]
    for info in contract['sections']:
        sid=info['section_id'];g=np.load(OUT/'geometry'/f'{sid}.npz');xy=g['xy'];co=g['coords'];lab=g['labels']
        masks=[]
        for rad in [130,230]:
            q=next(r for r in elig if r['section_id']==sid and int(r['mask_um'])==rad)
            if q['status']!='ELIGIBLE':continue
            d,v,near,far=geometry(xy,lab,rad);weight=near[v]/near.sum()-far[v]/far.sum()
            setup=fit_setup(xy[v],weight,20261002+rad);masks.append((rad,d[v],v,near[v],far[v],weight,setup))
        if not masks:continue
        short=gaussian_lattice(co,150,400,rng);long=gaussian_lattice(co,800,400,rng);white=gaussian_lattice(co,0,400,rng)
        mix=np.sqrt(.35)*short+np.sqrt(.5)*long+np.sqrt(.15)*white
        # Random slope orientation is independent of TLS; stress-tests nonstationarity.
        coordz=(xy-xy.mean(axis=0))/np.maximum(xy.std(axis=0),1)
        slope=coordz@rng.standard_normal((2,400))/np.sqrt(2)
        families={'white':white,'short150':short,'long800':long,'mixture':mix,
                  'non_gaussian':(mix*mix-1)/np.sqrt(2),'random_trend':mix+0.7*slope}
        for family,Yall in families.items():
            for rad,d,v,near,far,w,setup in masks:
                Y=Yall[v];p,err=fit_null(Y,setup,w)
                valid=np.isfinite(p);k=int(sum(p[valid]<=.05));n=int(sum(valid));lo,hi=wilson(k,n)
                records.append(dict(section_id=sid,cohort=info['cohort'],mask_um=rad,family=family,kind='null',effect_sd=0,width_um=0,n=n,rejections=k,rate=k/n,lo=lo,hi=hi,pass_gate=int(n==400 and hi<=.1),fit_rmse_median=float(np.median(err))))
                for width in [200,600]:
                    bump=np.exp(-d/width);contrast=w@bump;bump=bump/contrast
                    for effect in [.2,.5]:
                        base=Y[:,:200];injected=base+bump[:,None]*effect*base.std(axis=0,ddof=1)[None,:]
                        pp,ee=fit_null(injected,setup,w);ok=np.isfinite(pp);kk=int(sum(pp[ok]<=.05));nn=int(sum(ok));ll,hh=wilson(kk,nn)
                        records.append(dict(section_id=sid,cohort=info['cohort'],mask_um=rad,family=family,kind='power',effect_sd=effect,width_um=width,n=nn,rejections=kk,rate=kk/nn,lo=ll,hi=hh,pass_gate=int(nn==200 and ll>=.8),fit_rmse_median=float(np.median(ee))))
        write_table(target,records);print('calibrated',sid,len(records),flush=True)
    null=[r for r in records if r['kind']=='null'];power=[r for r in records if r['kind']=='power' and r['effect_sd']==.5]
    save_json(OUT/'calibration_gate.json',dict(status='PASS_MODEL_ONLY' if all(r['pass_gate'] for r in null+power) else 'FAILED',n_null_cells=len(null),n_null_fail=sum(not r['pass_gate'] for r in null),n_power_cells=len(power),n_power_fail=sum(not r['pass_gate'] for r in power),note='Real-data/patient confirmation remains separately gated.'))


def score():
    contract=json.loads((OUT/'contract.json').read_text());panel=json.loads((OUT/'panel.json').read_text());elig=table(OUT/'eligibility.tsv')
    if (OUT/'effects.tsv').exists():raise FileExistsError('effects already exist')
    effects=[];curves=[];checks=[];raw_cov=[]
    for info in contract['sections']:
        sid,cohort=info['section_id'],info['cohort'];m=source_metadata(cohort,sid);X,genes=counts(cohort,sid,m)
        for rad in [130,230]:
            q=next(r for r in elig if r['section_id']==sid and int(r['mask_um'])==rad)
            if q['status']!='ELIGIBLE':continue
            d,v,near,far=geometry(m['xy'],m['labels'],rad);S=visible_scores(X,genes,v,panel)
            lib=np.asarray(X[v].sum(axis=1)).ravel();ndetect=np.diff(X[v].indptr)
            near=near[v];far=far[v];dist=d[v];w=near/near.sum()-far/far.sum()
            p,err=fit_null(np.nan_to_num(S),fit_setup(m['xy'][v],w,20261002+rad),w)
            # Actual hidden-value intervention: perturb one feature by 1e9 in every hidden spot.
            if rad==130:
                changed=X.copy().tolil(); hidden=np.flatnonzero(~v)
                # Sparse changed values in every hidden row suffice to catch any cross-spot dependency.
                changed[hidden,0]=1e9; changed=changed.tocsr()
                S2=visible_scores(changed,genes,v,panel)
                delta=float(np.nanmax(np.abs(S2-S)));assert delta==0
                checks.append(dict(section_id=sid,mask_um=rad,hidden_rows=len(hidden),max_abs_change=delta))
            for j,(name,gl) in enumerate(panel.items()):
                k=len(set(gl)&set(genes));coverage=k/len(set(gl));values=S[:,j]
                sd=float(np.std(values,ddof=1));status='DESCRIPTIVE_ONLY'
                if coverage<.8 or k<10:status='NOT_TESTABLE_COVERAGE'
                elif not np.isfinite(sd) or sd<1e-12:status='NOT_TESTABLE_CONSTANT'
                delta=float(values[near].mean()-values[far].mean())
                eff=delta/sd if sd>1e-12 else None
                rho=float(stats.spearmanr(values,lib).statistic) if sd>1e-12 else None
                effects.append(dict(section_id=sid,cohort=cohort,mask_um=rad,set_id=name,status=status,n_near=int(near.sum()),n_far=int(far.sum()),genes_present=k,genes_defined=len(set(gl)),coverage=coverage,raw_delta=delta,effect_sd=eff if eff is not None else '',depth_spearman=rho if rho is not None else '',candidate_p=float(p[j]) if np.isfinite(p[j]) else '',variogram_relative_rmse=float(err[j]),formal_p='',claim='NOT_CALIBRATED_PATIENT_ID_UNKNOWN'))
                for lo in range(0,1000,100):
                    ix=(dist>=lo)&(dist<lo+100)
                    curves.append(dict(section_id=sid,cohort=cohort,mask_um=rad,set_id=name,lo_um=lo,hi_um=lo+100,n_spots=int(ix.sum()),mean_score=float(values[ix].mean()) if ix.any() else '',readout_status=status))
            # Depth is a technical sensitivity; it is not regressed out of the primary signal.
            raw_cov.append(dict(section_id=sid,mask_um=rad,depth_near=float(lib[near].mean()),depth_far=float(lib[far].mean()),detected_near=float(ndetect[near].mean()),detected_far=float(ndetect[far].mean())))
        print('scored',sid,flush=True)
    write_table(OUT/'effects.tsv',effects);write_table(OUT/'curves.tsv',curves);write_table(OUT/'hidden_intervention.tsv',checks);write_table(OUT/'depth_sensitivity.tsv',raw_cov)
    summary=[]
    for cohort in ['USZ','GSE175540']:
        for rad in [130,230]:
            for name in panel:
                rows=[r for r in effects if r['cohort']==cohort and r['mask_um']==rad and r['set_id']==name and r['status']=='DESCRIPTIVE_ONLY']
                e=np.array([r['effect_sd'] for r in rows],float)
                summary.append(dict(cohort=cohort,mask_um=rad,set_id=name,n_sections=len(e),n_confirmed_patients='',median_effect_sd=float(np.median(e)) if len(e) else '',n_positive=int(sum(e>0)),n_useful_positive=int(sum(e>=.2)),formal_p='',status='SECTION_DESCRIPTIVE_NOT_PATIENT_INFERENCE'))
    write_table(OUT/'summary.tsv',summary)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare','calibrate','score']);a=ap.parse_args()
    {'prepare':prepare,'calibrate':calibrate,'score':score}[a.stage]()
