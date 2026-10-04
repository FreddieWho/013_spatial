#!/usr/bin/env python3
"""Full-panel multiscale TLS study. Staged immutable artifacts, CPU/local only."""
import argparse,csv,hashlib,importlib.util,json,time
from pathlib import Path
import numpy as np
from scipy import sparse,spatial
from scipy.sparse import csgraph

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'infra/tls_multiscale_gobp_20261002'
SLIM=ROOT/'infra/gobp_compress_20261001/survivors_specific_floor20_jac70.tsv'
GMT=Path('/home/huyudi/006/data/pathway/c5.go.v2025.1.Hs.symbols.gmt')

def module(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def js(path,obj):path.write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def sha(path):return module('tls_peripheral').hash_file(path)
def tab(path):return module('tls_peripheral').table(path)
def tsv(path,rows):module('tls_peripheral').write_table(path,rows)


def tissue_graph(xy):
    pairs=spatial.cKDTree(xy).query_pairs(130.01,output_type='ndarray')
    i=np.r_[pairs[:,0],pairs[:,1]];j=np.r_[pairs[:,1],pairs[:,0]]
    dist=np.linalg.norm(xy[i]-xy[j],axis=1)
    return sparse.csr_matrix((dist,(i,j)),shape=(len(xy),len(xy)))


def components(graph,mask,min_size):
    ids=np.flatnonzero(mask)
    if not len(ids):return []
    n,labels=csgraph.connected_components(graph[ids][:,ids],directed=False)
    return [ids[labels==k] for k in range(n) if sum(labels==k)>=min_size]


def representative(xy,ids):
    return int(ids[np.argmin(np.linalg.norm(xy[ids]-xy[ids].mean(axis=0),axis=1))])


def build_endpoints(xy,labels):
    graph=tissue_graph(xy);endpoints=[];members=[]
    for family,minsize,mask in [('TLS',2,labels=='TLS'),('TUMOR_COMPONENT',10,labels=='TUM')]:
        for ids in components(graph,mask,minsize):
            idx=representative(xy,ids)
            diam=float(spatial.distance.pdist(xy[ids]).max()) if len(ids)>1 else 0.
            endpoints.append(dict(endpoint_id=len(endpoints),family=family,spot_index=idx,n_spots=len(ids),diameter_um=diam,x_um=float(xy[idx,0]),y_um=float(xy[idx,1])))
            members.append(ids)
    # Explicit TUM-versus-NOR/INFL/LN interface; unknown is never non-tumor.
    explicit_non=np.isin(labels,['NOR','INFL','LN'])
    edge=(labels=='TUM')&(np.asarray((graph>0)@explicit_non.astype(int)).ravel()>0)
    for ids in components(graph,edge,2):
        idx=representative(xy,ids)
        endpoints.append(dict(endpoint_id=len(endpoints),family='TUMOR_EDGE',spot_index=idx,n_spots=len(ids),diameter_um=float(spatial.distance.pdist(xy[ids]).max()),x_um=float(xy[idx,0]),y_um=float(xy[idx,1])))
        members.append(ids)
    # Location-defined tumor subregions at two scales, without expression selection.
    tumor=np.flatnonzero(labels=='TUM')
    for width in [500,1000]:
        tiles=np.floor(xy[tumor]/width).astype(int)
        for tile in np.unique(tiles,axis=0):
            ids=tumor[np.all(tiles==tile,axis=1)]
            if len(ids)<5:continue
            idx=representative(xy,ids)
            endpoints.append(dict(endpoint_id=len(endpoints),family=f'TUMOR_TILE_{width}',spot_index=idx,n_spots=len(ids),diameter_um=float(spatial.distance.pdist(xy[ids]).max()),x_um=float(xy[idx,0]),y_um=float(xy[idx,1])))
            members.append(ids)
    return graph,endpoints,members,edge


def prepare():
    OUT.mkdir(exist_ok=True)
    if (OUT/'contract.json').exists():raise FileExistsError('Existing contract')
    for d in ['geometry','scores','receipts']:(OUT/d).mkdir(exist_ok=True)
    r17=module('tls_peripheral');r17c=json.loads((ROOT/'infra/tls_peripheral_20261002/contract.json').read_text())
    panelrows=tab(SLIM);gmt=module('tls_field_hallmark').read_gmt(GMT)
    panel={r['program_id']:sorted(set(gmt[r['program_id']])) for r in panelrows}
    assert len(panel)==1691
    # Side-by-side fixed comparators; all GO terms remain separate in every table.
    controls=json.loads((ROOT/'infra/tls_peripheral_20261002/panel.json').read_text())
    js(OUT/'panel_gobp.json',panel);js(OUT/'panel_controls.json',controls)
    sections=[];source_receipts=[];coverage=[]
    for info in r17c['sections']:
        sid,cohort=info['section_id'],info['cohort'];m=r17.source_metadata(cohort,sid)
        graph,ends,members,edge=build_endpoints(m['xy'],m['labels'])
        tsv(OUT/'geometry'/f'{sid}_endpoints.tsv',ends)
        np.savez_compressed(OUT/'geometry'/f'{sid}.npz',xy=m['xy'],labels=m['labels'],order=m['order'],edge=edge,
                            endpoint_members=np.concatenate(members) if members else np.array([],int),
                            endpoint_indptr=np.r_[0,np.cumsum([len(v) for v in members])])
        sparse.save_npz(OUT/'geometry'/f'{sid}_graph.npz',graph)
        present=set(m['genes'][m['gene_keep']]);tls=[e for e in ends if e['family']=='TLS'];tum=[e for e in ends if e['family'].startswith('TUMOR') and e['family']!='TUMOR_EDGE'];edges=[e for e in ends if e['family']=='TUMOR_EDGE']
        sections.append(dict(section_id=sid,cohort=cohort,n_spots=len(m['labels']),n_TLS_spots=int(sum(m['labels']=='TLS')),n_TLS_components=len(tls),n_tumor_endpoints=len(tum),n_edge_endpoints=len(edges),n_known_edge_spots=int(edge.sum()),n_TT_pairs=len(tls)*(len(tls)-1)//2,n_TR_pairs=len(tls)*len(tum),n_TE_pairs=len(tls)*len(edges),lost_TLS=m['lost_tls'],geometry_status='KNOWN_TLS_MISSING' if m['lost_tls'] else 'SOURCE_GEOMETRY_AVAILABLE',tumor_status='EXPLICIT_TUM' if sum(m['labels']=='TUM') else 'NOT_TESTABLE_NO_TUMOR_LABEL'))
        for name,genes in panel.items():
            used=set(genes)&present;fraction=len(used)/len(genes)
            coverage.append(dict(section_id=sid,set_id=name,defined=len(genes),present=len(used),fraction=fraction,status='COVERED' if fraction>=.8 and len(used)>=10 else 'NOT_TESTABLE_COVERAGE'))
        source_receipts.extend([{'section_id':sid,'source':str(f.relative_to(ROOT)),'size':f.stat().st_size,'sha256':r17.hash_file(f)} for f in m['sources']])
        print('metadata',sid,'TLS components',len(tls),'TR/TE pairs',len(tls)*len(tum),len(tls)*len(edges),flush=True)
    tsv(OUT/'sections.tsv',sections);tsv(OUT/'coverage.tsv',coverage);tsv(OUT/'source_hashes.tsv',source_receipts)
    contract=dict(status='FROZEN_BEFORE_GOBP_SCORING',decision='D-171',scope='FULL_1691_MULTISCALE_HALO_EDGE_BRIDGE_AND_REVERSE',
      panel=dict(n_gobp=1691,go_sha256=sha(OUT/'panel_gobp.json'),slim_sha256=sha(SLIM),gmt_sha256=sha(GMT),controls=51),
      sources=r17c['sections'],score='raw count duplicate-symbol sum; per-spot AUCell top5%; no cross-spot normalization/smoothing; mask applied by each operator',
      coverage={'min_fraction':.8,'min_genes':10,'all_terms_retained':True},
      radial_edges_um=[130,230,400,800,1600,3200,6400,'section_max'],
      length_strata_um=[0,500,1000,2000,4000,8000,'section_max'],corridor_halfwidth_um=[100,250,500,1000],mask_um=[130,230],
      endpoint_definitions={'TLS':'source label connected components, >=2 spots, graph edges<=130.01um; singleton labels remain in exclusion mask',
        'tumor_regions':'explicit TUM connected components >=10 spots and fixed 500/1000um tiles >=5 TUM spots; medoid-nearest-centroid',
        'tumor_edge':'TUM neighboring explicit NOR/INFL/LN, interface connected components >=2 spots; UNKNOWN and TLS do not define non-tumor'},
      geometry={'paths':['straight','tissue_graph_shortest'],'graph':'<=130.01um spatial neighbor graph; no expression cost',
        'all_pairs':'all unordered TLS-TLS and all TLS-to-region/interface endpoints; no signal-based selection',
        'endpoint_masks':'all known TLS cores plus radius ring for forward analyses; target-end exclusion disk; uniform query mask in reverse',
        'support_min':'>=5 spots per longitudinal segment in corridor and each of two matched-width sidebands; 5 segments; report partial/not_testable',
        'continuity':'all 5 segments measured; require at least4 same-direction contrasts for a descriptive bridge candidate; formal gate separate',
        'scale_policy':'no minimum-mm assumption, full continuous lengths and length/TLS-component-diameter retained'},
      inference={'primary':'exploratory complete screen; no reuse of failed R17 p','formal':'pending selection-aware spatial calibration and sample identity',
        'reverse':'uniform query-centered masks; no TLS labels/distances/GT-shaped holes in features; training-fold-only selection; all1691 programs attempted; section/leave-cohort-out internal validation'},
      resource={'GPU':False,'local_only':True,'min_free_bytes':1200000000000})
    js(OUT/'contract.json',contract)
    js(OUT/'goal_progress.json',{'goal_complete':False,'requirements':{'full_panel_source_qualification':'COMPLETE','full_panel_scoring':'NOT_RUN','multiscale_halo':'NOT_RUN','multiscale_edge':'NOT_RUN','TLS_TLS_bridge':'NOT_RUN','TLS_tumor_region_bridge':'NOT_RUN','TLS_tumor_edge_bridge':'NOT_RUN','spatial_selection_calibration':'NOT_RUN','molecular_to_TLS_prediction':'NOT_RUN','full_scope_acceptance':'NOT_RUN'}})
    print('FROZEN',sha(OUT/'contract.json'),flush=True)


def score():
    import shutil
    assert shutil.disk_usage(ROOT).free>=1200000000000
    contract=json.loads((OUT/'contract.json').read_text());panel=json.loads((OUT/'panel_gobp.json').read_text());control=json.loads((OUT/'panel_controls.json').read_text())
    assert sha(OUT/'panel_gobp.json')==contract['panel']['go_sha256']
    allsets=panel|control;r17=module('tls_peripheral');scorer=module('gobp_halo_screen');all_receipts=[]
    for info in contract['sources']:
        sid,cohort=info['section_id'],info['cohort'];receipt=OUT/'receipts'/f'{sid}_scores.json';out=OUT/'scores'/f'{sid}.npz'
        if receipt.exists():
            rr=json.loads(receipt.read_text());assert sha(out)==rr['sha256'];all_receipts.append(rr);print('verified completed cache',sid,flush=True);continue
        t=time.time();m=r17.source_metadata(cohort,sid);X,genes=r17.counts(cohort,sid,m);index={g:i for i,g in enumerate(genes)}
        members=[np.array([index[g] for g in gl if g in index],dtype=int) for gl in allsets.values()]
        S,k=scorer.score_sets(X,members)
        # No pathway/readout is omitted; coverage and constant flags travel with cache.
        fraction=k/np.array([len(v) for v in allsets.values()]);valid=(fraction>=.8)&(k>=10)&(np.nanstd(S,axis=0)>1e-12)
        lib=np.asarray(X.sum(axis=1)).ravel().astype(np.float32)
        with (out.with_suffix('.partial')).open('wb') as f:
            np.savez_compressed(f,scores=S.astype(np.float32),set_ids=np.array(list(allsets)),coverage=fraction,n_genes=k,valid=valid,libsize=lib,detected=np.diff(X.indptr))
        out.with_suffix('.partial').replace(out)
        rr=dict(section_id=sid,cohort=cohort,status='FULL_PANEL_SCORED',n_spots=X.shape[0],n_gobp=1691,n_controls=51,n_covered_variable_GO=int(valid[:1691].sum()),elapsed_seconds=time.time()-t,sha256=sha(out),contract_sha256=sha(OUT/'contract.json'))
        js(receipt,rr);all_receipts.append(rr);print('scored',sid,rr['n_covered_variable_GO'],'/1691',round(time.time()-t,1),'s',flush=True)
        assert shutil.disk_usage(ROOT).free>=1200000000000
    js(OUT/'scoring_complete.json',dict(status='COMPLETE',sections=len(all_receipts),gobp_per_section=1691,attempted_section_gobp=len(all_receipts)*1691,all_receipts=all_receipts))
    p=json.loads((OUT/'goal_progress.json').read_text());p['requirements']['full_panel_scoring']='COMPLETE';js(OUT/'goal_progress.json',p)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare','score']);a=ap.parse_args();{'prepare':prepare,'score':score}[a.stage]()
