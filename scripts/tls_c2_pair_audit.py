#!/usr/bin/env python3
"""D-186: fixed c_2 audit; no fitting new rules or changing failed predictions."""
from pathlib import Path
import importlib.util,json,csv,hashlib,collections
import numpy as np
from scipy import sparse
from scipy.sparse import csgraph
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'infra/tls_multiscale_gobp_20261002';OUT=ROOT/'infra/tls_c2_pair_audit_20261003'
def mod(name):
    s=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
R=mod('tls_multiscale_reverse_models');V=mod('tls_multiscale_reverse');P=mod('tls_peripheral');S=mod('gobp_halo_screen')
PAIR=['gse175540-GSM5924030_ffpe_c_2','gse175540-GSM5924050_frozen_c_2']
PROGRAMS=['B_AXIS','GOBP_GAMMA_DELTA_T_CELL_ACTIVATION','GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY']
def read(p):return list(csv.DictReader(Path(p).open(),delimiter='\t'))
def write(p,rows):
    with Path(p).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def sha(p):return P.hash_file(p)
def main():
    OUT.mkdir(exist_ok=True)
    identity=ROOT/'infra/bioinf-data-index/raw/user_gse175540_20261003/current_26_source_identity_confirmed.tsv'
    panel=json.loads((BASE/'panel_gobp.json').read_text())|json.loads((BASE/'panel_controls.json').read_text())
    contract={'decision':'D-186','pair':PAIR,'fixed_programs':PROGRAMS,'identity_sha256':sha(identity),'scoring':'recompute original raw-count rank scores for every1742 program','reverse':'reconstruct labels and fixed-readout features at230/500/1000; decompose fixed1000 group-held-out model only','source_sha256':{},'prediction_sha256':{}}
    for sid in PAIR:
        m=P.source_metadata('GSE175540',sid)
        for f in m['sources']:contract['source_sha256'][str(f)]=sha(f)
        for f in [BASE/'geometry'/f'{sid}.npz',BASE/'scores'/f'{sid}.npz']:contract['source_sha256'][str(f)]=sha(f)
    predpath=ROOT/'infra/tls_patient_validation_20261003/models/1000_leave_patient_GSE175540_c_2.npz';contract['prediction_sha256'][str(predpath)]=sha(predpath)
    (OUT/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    data=R.load_sources(1000);train=[s for s in data if s not in PAIR];fit=R.fit_moments(*R.combine_stats(data,train));saved=np.load(predpath);summaries=[];supports=[];qrows=[];contrib=[];classmetrics=[];checks=[];scorestats=[]
    for sid in PAIR:
        g=np.load(BASE/'geometry'/f'{sid}.npz');z={k:v for k,v in np.load(BASE/'scores'/f'{sid}.npz').items()};names=z['set_ids'];ix=np.array([int(np.flatnonzero(names==x)[0]) for x in PROGRAMS]);bidx=ix[0]
        source=P.source_metadata('GSE175540',sid);np.testing.assert_array_equal(g['labels'],source['labels']);np.testing.assert_array_equal(g['order'],source['order']);np.testing.assert_allclose(g['xy'],source['xy'],rtol=0,atol=0);assert source['lost_tls']==0
        X,genes=P.counts('GSE175540',sid,source);index={g:i for i,g in enumerate(genes)};members=[np.array([index[g] for g in gl if g in index],dtype=int) for gl in panel.values()]
        raw,k=S.score_sets(X,members);np.testing.assert_allclose(raw,z['scores'],rtol=1e-6,atol=1e-7,equal_nan=True);np.testing.assert_array_equal(k,z['n_genes']);np.testing.assert_allclose(np.asarray(X.sum(axis=1)).ravel(),z['libsize']);np.testing.assert_array_equal(np.diff(X.indptr),z['detected'])
        graph=sparse.load_npz(BASE/'geometry'/f'{sid}_graph.npz');tls=np.flatnonzero(g['labels']=='TLS');nc,cc=csgraph.connected_components(graph[tls][:,tls],directed=False);component=np.full(len(g['labels']),-1,int);component[tls]=cc
        summaries.append(dict(section_id=sid,n_spots=len(g['labels']),n_TLS=len(tls),n_TLS_components=nc,n_source_genes=len(genes),n_valid_GO=int(z['valid'][:1691].sum()),median_library=float(np.median(z['libsize'])),median_detected=float(np.median(z['detected'])),x_span_um=float(np.ptp(g['xy'][:,0])),y_span_um=float(np.ptp(g['xy'][:,1])),lost_TLS=source['lost_tls']))
        for j,name in zip(ix,PROGRAMS):
            for label in ['TLS','NO_TLS']:
                v=z['scores'][g['labels']==label,j]
                scorestats.append(dict(section_id=sid,program=name,label=label,n_spots=len(v),mean=float(np.mean(v)),median=float(np.median(v)),q25=float(np.quantile(v,.25)),q75=float(np.quantile(v,.75)),gene_coverage=float(z['coverage'][j]),n_member_genes=int(z['n_genes'][j])))
        for radius in [230,500,1000]:
            qs=V.query_labels(g['xy'],g['labels'],graph,radius);old=read(BASE/'reverse'/f'{sid}_{radius}_queries.tsv');assert len(old)==len(qs)
            for a,b in zip(qs,old):
                for key in ['spot_index','label','n_first_ring','n_second_ring']:assert a[key]==int(b[key])
                assert a['status']==b['status']
            elig=[q for q in qs if q['label']>=0];F,B=V.pool_query_features(g['xy'],np.nan_to_num(z['scores'][:,ix]),z['libsize'],elig,radius);cached={k:v for k,v in np.load(BASE/'reverse'/f'{sid}_{radius}_features.npz').items()}
            np.testing.assert_array_equal(cached['spot_index'],[q['spot_index'] for q in elig]);np.testing.assert_array_equal(cached['y'],[q['label'] for q in elig]);np.testing.assert_allclose(F,cached['features'][:,ix,:],rtol=2e-5,atol=2e-5);np.testing.assert_allclose(B,cached['baseline'],rtol=1e-6,atol=1e-6)
            counts=collections.Counter(q['status'] for q in qs)
            for status,n in counts.items():supports.append(dict(section_id=sid,mask_um=radius,status=status,n_queries=n,n_positive=sum(q['label']==1 for q in qs),n_negative=sum(q['label']==0 for q in qs)))
            if radius!=1000:continue
            pr=R.predict(cached['features'],fit);mask=saved['section_ids']==sid;np.testing.assert_allclose(pr[:,ix],saved['prediction'][mask][:,ix],rtol=2e-5,atol=2e-5)
            y=cached['y'];positives=y==1;negatives=y==0;assert positives.any() and negatives.any()
            for pi,(j,name) in enumerate(zip(ix,PROGRAMS)):
                auc,ap=R.metrics(y,pr[:,j]);classmetrics.append(dict(section_id=sid,program=name,n_positive=int(positives.sum()),n_negative=int(negatives.sum()),n_positive_TLS_components=len(set(component[cached['spot_index'][positives]])),auc=float(auc[0]),ap=float(ap[0]),mean_positive_prediction=float(pr[positives,j].mean()),mean_negative_prediction=float(pr[negatives,j].mean())))
                for dim in range(15):
                    posmean=float(F[positives,pi,dim].mean());negmean=float(F[negatives,pi,dim].mean());coef=float(fit[0][j,dim]);fa,_=R.metrics(y,F[:,pi,dim])
                    contrib.append(dict(section_id=sid,program=name,feature=('radial' if dim<5 else 'angular_order1' if dim<10 else 'angular_order2')+f'_ring{dim%5}',coefficient=coef,mean_positive=posmean,mean_negative=negmean,mean_prediction_gap_contribution=coef*(posmean-negmean),single_feature_auc=float(fa[0])))
            for qi,q in enumerate(elig):
                spot=q['spot_index'];dist=np.linalg.norm(g['xy']-g['xy'][spot],axis=1);visible=dist>radius;vtls=visible&(g['labels']=='TLS');component_here=component[spot]
                # Source annotations below are audit-only, never fed to the model.
                qrows.append(dict(section_id=sid,spot_index=spot,barcode_order=int(g['order'][spot]),y=q['label'],target_TLS_component=int(component_here),n_visible=int(visible.sum()),n_visible_TLS=int(vtls.sum()),visible_TLS_fraction=float(vtls.sum()/visible.sum()),visible_B_mean=float(z['scores'][visible,bidx].mean()),visible_B_sd=float(z['scores'][visible,bidx].std(ddof=1)),first_ring_support=q['n_first_ring'],second_ring_support=q['n_second_ring'],n_supported_radial_bins=int(sum(B[qi,:5]>=np.log1p(20)-1e-5)),B_prediction=float(pr[qi,bidx]),geometry_prediction=float(saved['baseline'][mask][qi]),combined_prediction=float(saved['combined'][mask][qi])))
        checks.append(dict(section_id=sid,original_metadata_match=True,all1742_scores_recomputed_match=True,raw_depth_detected_match=True,query_labels_all3_masks_match=True,fixed_features_all3_masks_match=True,grouped_B_and_fixed_GO_predictions_recomputed_match=True))
        print('audited',sid,flush=True)
    write(OUT/'source_comparison.tsv',summaries);write(OUT/'query_support.tsv',supports);write(OUT/'raw_score_by_annotation.tsv',scorestats);write(OUT/'R1000_fixed_readout_metrics.tsv',classmetrics);write(OUT/'R1000_feature_contributions.tsv',contrib);write(OUT/'R1000_query_diagnostics.tsv',qrows)
    byclass=[]
    for sid in PAIR:
        for y in [0,1]:
            rows=[r for r in qrows if r['section_id']==sid and r['y']==y];r={'section_id':sid,'y':y,'n_queries':len(rows)}
            for key in ['n_visible','n_visible_TLS','visible_TLS_fraction','visible_B_mean','visible_B_sd','first_ring_support','second_ring_support','n_supported_radial_bins','B_prediction','combined_prediction']:
                r[key+'_median']=float(np.median([x[key] for x in rows]));r[key+'_mean']=float(np.mean([x[key] for x in rows]))
            byclass.append(r)
    write(OUT/'R1000_query_class_summary.tsv',byclass)
    (OUT/'complete.json').write_text(json.dumps(dict(status='COMPLETE_NO_REPRODUCTION_MISMATCH_FOUND',checks=checks,causal_preparation_effect_identified=False,model_or_labels_changed=False,scope='Source-label agreement is not independent pathology validation; diagnostic association is not cause'),indent=2)+'\n')
if __name__=='__main__':main()
