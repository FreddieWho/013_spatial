#!/usr/bin/env python3
"""D-188: same query population, paired saved R500/R1000 predictions."""
from pathlib import Path
import importlib.util,csv,json,hashlib
import numpy as np
from scipy.sparse import csgraph
from scipy import sparse
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'infra/tls_multiscale_gobp_20261002';OUT=ROOT/'infra/tls_paired_masks_20261003';VALID=ROOT/'infra/tls_patient_validation_20261003'
s=importlib.util.spec_from_file_location('reverse_models',ROOT/'scripts/tls_multiscale_reverse_models.py');R=importlib.util.module_from_spec(s);s.loader.exec_module(R)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return list(csv.DictReader(Path(p).open(),delimiter='\t'))
def write(p,rows):
    with Path(p).open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def align_queries(a,ya,b,yb):
    # Empty feature archives used NumPy's default float dtype; indices remain integer.
    assert np.all(np.asarray(a)==np.asarray(a,dtype=np.int64)) and np.all(np.asarray(b)==np.asarray(b,dtype=np.int64))
    a=np.asarray(a,dtype=np.int64);b=np.asarray(b,dtype=np.int64)
    assert len(np.unique(a))==len(a) and len(np.unique(b))==len(b)
    common,ia,ib=np.intersect1d(a,b,return_indices=True);np.testing.assert_array_equal(ya[ia],yb[ib]);return common,ia,ib,ya[ia]

def main():
    OUT.mkdir(exist_ok=True);c=json.loads((VALID/'contract.json').read_text());mapping=c['group_map'];jobs={(j['mask_um'],j['split'],j['heldout_group']):j for j in c['jobs']};pairs=sorted({(j['split'],j['heldout_group']) for j in c['jobs']});rows=[];support=[];querymap=[];hashes={}
    for split,group in pairs:
        js=[jobs[(rad,split,group)] for rad in [500,1000]];assert js[0]['test_sections']==js[1]['test_sections']
        zs=[];qual=[]
        for job in js:
            path=ROOT/job['directory']/(job['stem']+'.npz');receipt=json.loads(path.with_suffix('.json').read_text());digest=sha(path);assert digest==receipt['sha256'];hashes[str(path.relative_to(ROOT))]=digest
            zs.append({k:v for k,v in np.load(path).items()})
            qual.append({(r['section_id'],r['program']):r['status']=='INTERNAL_SECTION_EVALUATION' for r in read(path.parent/(job['stem']+'_metrics.tsv')) if r['arm'] in ['full','geometry_depth','train_selected_GO']})
        np.testing.assert_array_equal(zs[0]['set_ids'],zs[1]['set_ids']);names=list(zs[0]['set_ids'])+['geometry_depth','train_selected_GO'];assert len(names)==1744
        for sid in js[0]['test_sections']:
            cohort='USZ' if sid.startswith('usz-') else 'GSE175540';features=[];predictions=[]
            for radius,z in zip([500,1000],zs):
                cache=np.load(BASE/'reverse'/f'{sid}_{radius}_features.npz');mask=z['section_ids']==sid;np.testing.assert_array_equal(z['y'][mask],cache['y']);features.append((cache['spot_index'],cache['y']));predictions.append(np.column_stack([z['prediction'][mask],z['baseline'][mask],z['combined'][mask]]))
            common,ia,ib,y=align_queries(*features[0],*features[1]);g=np.load(BASE/'geometry'/f'{sid}.npz');tls=np.flatnonzero(g['labels']=='TLS');graph=sparse.load_npz(BASE/'geometry'/f'{sid}_graph.npz');nc,cc=csgraph.connected_components(graph[tls][:,tls],directed=False);lookup=np.full(len(g['labels']),-1,int);lookup[tls]=cc
            npos=int(sum(y==1));nneg=int(sum(y==0));components=set(lookup[common[y==1]]);assert -1 not in components
            support.append(dict(split=split,heldout_group=group,section_id=sid,cohort=cohort,evaluation_unit=mapping[sid],n_queries_500=len(features[0][0]),n_queries_1000=len(features[1][0]),n_common=len(common),n_common_positive=npos,n_common_negative=nneg,n_common_TLS_instances=len(components),both_classes=npos>0 and nneg>0,n_positive_500=int(sum(features[0][1]==1)),n_positive_1000=int(sum(features[1][1]==1))))
            for spot,label in zip(common,y):querymap.append(dict(split=split,section_id=sid,evaluation_unit=mapping[sid],spot_index=int(spot),label=int(label),component_id=int(lookup[spot]) if label==1 else '',matching_rule='same source spot, same label, eligible at both masks'))
            a,ap=R.metrics(y,predictions[0][ia]);b,bp=R.metrics(y,predictions[1][ib])
            for j,name in enumerate(names):
                ok=npos>0 and nneg>0 and qual[0][(sid,name)] and qual[1][(sid,name)]
                rows.append(dict(split=split,heldout_group=group,section_id=sid,cohort=cohort,evaluation_unit=mapping[sid],arm='full' if j<1742 else name,program=str(name),n_common_positive=npos,n_common_negative=nneg,n_common_TLS_instances=len(components),auc_500=float(a[j]) if ok else '',auc_1000=float(b[j]) if ok else '',ap_500=float(ap[j]) if ok else '',ap_1000=float(bp[j]) if ok else '',delta_auc=float(b[j]-a[j]) if ok else '',status='PAIRED_INTERNAL_ONLY' if ok else 'NOT_TESTABLE_COMMON_LABELS_OR_PROGRAM'))
        print('paired',split,group,flush=True)
    # Load project SciPy before pandas chooses its C++ runtime.
    import pandas as pd
    f=pd.DataFrame(rows);num=['auc_500','auc_1000','ap_500','ap_1000','delta_auc']
    for key in num:f[key]=pd.to_numeric(f[key],errors='coerce')
    keys=['split','cohort','arm','program']
    u=f.groupby(keys+['evaluation_unit']).agg(**{k:(k,'mean') for k in num},n_evaluable_sections=('delta_auc','count'),n_sections=('section_id','size')).reset_index()
    summary=u.groupby(keys).agg(**{k:(k,'mean') for k in num},n_evaluable_units=('delta_auc','count'),n_units=('evaluation_unit','size'),n_positive_delta=('delta_auc',lambda a:int((a>0).sum())),n_negative_delta=('delta_auc',lambda a:int((a<0).sum()))).reset_index()
    assert np.allclose(summary.delta_auc,summary.auc_1000-summary.auc_500,equal_nan=True)
    old=pd.read_csv(VALID/'all_program_summary.tsv',sep='\t')
    for radius in [500,1000]:
        subset=old.loc[old.mask_um.eq(radius),keys+['mean_auc','n_evaluable_units']].rename(columns={'mean_auc':f'original_full_population_auc_{radius}','n_evaluable_units':f'original_full_population_n_{radius}'})
        summary=summary.merge(subset,on=keys,how='left',validate='one_to_one')
    summary['unit_type']=np.where(summary.cohort.eq('GSE175540'),'patient','tumor_sample_unknown_patient');summary['status']=np.where(summary.n_evaluable_units>0,'PAIRED_INTERNAL_ONLY','NOT_TESTABLE');summary['formal_p']=''
    f.to_csv(OUT/'all_section_paired_metrics.tsv.gz',sep='\t',index=False);u.to_csv(OUT/'all_unit_paired_metrics.tsv.gz',sep='\t',index=False);summary.to_csv(OUT/'all_program_paired_summary.tsv',sep='\t',index=False)
    fixed=['B_AXIS','geometry_depth','train_selected_GO','GOBP_GAMMA_DELTA_T_CELL_ACTIVATION','GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY'];summary[summary.program.isin(fixed)].to_csv(OUT/'fixed_readout_paired_summary.tsv',sep='\t',index=False)
    write(OUT/'common_query_support.tsv',support);write(OUT/'common_query_identity.tsv',querymap)
    contract=dict(decision='D-188',patient_validation_contract_sha256=sha(VALID/'contract.json'),predictions_sha256=hashes,arms='full1742 plus geometry/depth and training-selected GO; radial/angular not recomputed',match='source spot_index intersection with exact label agreement',new_training=False,formal_confirmation=False)
    (OUT/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    (OUT/'complete.json').write_text(json.dumps(dict(status='COMPLETE_PAIRED_COMMON_QUERY_EVALUATION',paired_model_jobs=len(pairs),section_comparisons=len(support),section_metric_rows=len(f),unit_metric_rows=len(u),summary_rows=len(summary),programs_per_comparison=1744,matched_query_rows=len(querymap),source_predictions_sha_verified=True,query_label_agreement=True,unit_mean_delta_consistency=True,original_full_population_preserved=True,scientific_confirmation=False),indent=2)+'\n')
if __name__=='__main__':main()
