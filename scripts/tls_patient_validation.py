#!/usr/bin/env python3
"""D-185: exact-fold reuse and grouped c_2 holdout, with patient-macro summaries."""
import argparse,csv,json,hashlib,importlib.util,shutil
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('reverse_models',ROOT/'scripts/tls_multiscale_reverse_models.py')
R=importlib.util.module_from_spec(spec);spec.loader.exec_module(R)
OLD=R.OUT;OUT=ROOT/'infra/tls_patient_validation_20261003'
IDENTITY=ROOT/'infra/bioinf-data-index/raw/user_gse175540_20261003/current_26_source_identity_confirmed.tsv'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def groups():
    rows=list(csv.DictReader(IDENTITY.open(),delimiter='\t'))
    return {r['section_id']:(r['patient_id'] if r['cohort']=='GSE175540' else 'sample:'+r['section_id']) for r in rows}
def split_group(mapping,group):
    test=[s for s,g in mapping.items() if g==group];train=[s for s,g in mapping.items() if g!=group]
    assert test and not ({mapping[s] for s in test}&{mapping[s] for s in train})
    return train,test

def prepare():
    OUT.mkdir(exist_ok=True);mapping=groups();allids=set(mapping)
    assert len(mapping)==26 and len({g for g in mapping.values() if not g.startswith('sample:')})==17
    jobs=[]
    for radius in R.RADII:
        for group in sorted(set(mapping.values())):
            if group.startswith('sample:'):continue
            train,test=split_group(mapping,group)
            if len(test)==1:
                stem=f'{radius}_leave_section_{test[0]}';directory=OLD/'models';action='REUSE_EXACT_FOLD'
            else:
                assert group=='GSE175540:c_2' and len(test)==2
                stem=f'{radius}_leave_patient_GSE175540_c_2';directory=OUT/'models';action='REFIT_GROUP'
            jobs.append(dict(mask_um=radius,split='leave_patient',heldout_group=group,train_sections=train,test_sections=test,directory=str(directory.relative_to(ROOT)),stem=stem,action=action))
        for cohort in ['USZ','GSE175540']:
            test=[s for s in mapping if (s.startswith('usz-'))==(cohort=='USZ')];train=sorted(allids-set(test))
            jobs.append(dict(mask_um=radius,split='leave_cohort',heldout_group=cohort,train_sections=train,test_sections=test,directory=str((OLD/'models').relative_to(ROOT)),stem=f'{radius}_leave_cohort_{cohort}',action='REUSE_EXACT_FOLD'))
    assert len(jobs)==57
    for j in jobs:
        assert not ({mapping[s] for s in j['train_sections']}&{mapping[s] for s in j['test_sections']})
        if j['action']=='REUSE_EXACT_FOLD':
            p=ROOT/j['directory'];receipt=json.loads((p/(j['stem']+'.json')).read_text())
            assert set(receipt['train_sections'])==set(j['train_sections']) and set(receipt['test_sections'])==set(j['test_sections'])
            assert sha(p/(j['stem']+'.npz'))==receipt['sha256']
        j['training_units']=[mapping[s] for s in j['train_sections']]
        j['test_units']=[mapping[s] for s in j['test_sections']]
    contract=dict(decision='D-185',identity_sha256=sha(IDENTITY),base_contract_sha256=sha(OLD.parent/'contract.json'),source_model_script_sha256=sha(ROOT/'scripts/tls_multiscale_reverse_models.py'),group_map=mapping,jobs=jobs,training_weights='original equal section/class weights; not patient-balanced',evaluation='mean evaluable section metrics within patient, then equal patients; USZ equal samples',scientific_status='INTERNAL_ONLY',disk_free_bytes=shutil.disk_usage(OUT).free)
    assert contract['disk_free_bytes']>=1.2e12
    target=OUT/'contract.json'
    if target.exists():
        prior=json.loads(target.read_text());assert all(prior[k]==contract[k] for k in contract if k!='disk_free_bytes')
    else:R.M.js(target,contract)
    print('57 folds: 54 exact reuse, 3 refit',flush=True)

def train():
    prepare();mapping=groups();_,test=split_group(mapping,'GSE175540:c_2')
    R.train(jobs_override=[('leave_patient','GSE175540_c_2',test)],dest_override=OUT/'models',completion_path=OUT/'refit_complete.json')

def summarize():
    import pandas as pd
    c=json.loads((OUT/'contract.json').read_text());mapping=c['group_map'];frames=[];audit=[]
    assert sha(IDENTITY)==c['identity_sha256']
    for j in c['jobs']:
        d=ROOT/j['directory'];receipt=json.loads((d/(j['stem']+'.json')).read_text());prediction=d/(j['stem']+'.npz');metrics=d/(j['stem']+'_metrics.tsv')
        assert sha(prediction)==receipt['sha256']
        assert set(receipt['test_sections'])==set(j['test_sections']) and set(receipt['train_sections'])==set(j['train_sections'])
        assert not ({mapping[s] for s in receipt['train_sections']}&{mapping[s] for s in receipt['test_sections']})
        f=pd.read_csv(metrics,sep='\t');assert set(f.section_id)==set(j['test_sections']);f['split']=j['split'];f['heldout_job']=j['heldout_group'];f['evaluation_unit']=f.section_id.map(mapping);f['fold_action']=j['action'];frames.append(f)
        audit.append(dict(stem=j['stem'],action=j['action'],prediction_sha256=receipt['sha256'],metrics_sha256=sha(metrics),train_test_group_overlap=False))
    f=pd.concat(frames,ignore_index=True);f.loc[~f.status.eq('INTERNAL_SECTION_EVALUATION'),['auc','ap']]=np.nan
    keys=['mask_um','split','cohort','arm','program']
    assert len(f)==3*(18+26)*(1742*3+2)
    u=f.groupby(keys+['evaluation_unit'],dropna=False).agg(auc=('auc','mean'),ap=('ap','mean'),n_evaluable_sections=('auc','count'),n_sections=('section_id','size'),n_positive=('n_positive','sum'),n_negative=('n_negative','sum')).reset_index()
    u['all_sections_evaluable']=u.n_evaluable_sections.eq(u.n_sections)
    a=u.groupby(keys,dropna=False).agg(mean_auc=('auc','mean'),mean_ap=('ap','mean'),n_evaluable_units=('auc','count'),n_units=('evaluation_unit','size'),n_all_sections_evaluable_units=('all_sections_evaluable','sum')).reset_index()
    a['unit_type']=np.where(a.cohort.eq('GSE175540'),'patient','tumor_sample_unknown_patient')
    a['status']=np.where(a.n_evaluable_units>0,'INTERNAL_VALIDATION','NOT_TESTABLE')
    # Same source/program/arm old predictions; compare equal-section vs corrected folds and equal-patient aggregation.
    old=[]
    for radius in R.RADII:
        for split,ids in [('leave_section',[s for s in mapping if s.startswith('gse175540-')]),('leave_cohort',['USZ','GSE175540'])]:
            for sid in ids:
                x=pd.read_csv(OLD/'models'/f'{radius}_{split}_{sid}_metrics.tsv',sep='\t');x['split']='leave_patient' if split=='leave_section' else split;old.append(x)
    old=pd.concat(old,ignore_index=True);old.loc[~old.status.eq('INTERNAL_SECTION_EVALUATION'),['auc','ap']]=np.nan
    old_summary=old.groupby(keys).agg(old_section_mean_auc=('auc','mean'),old_evaluable_sections=('auc','count')).reset_index()
    new_section=f.groupby(keys).agg(new_section_mean_auc=('auc','mean')).reset_index()
    a=a.merge(old_summary,on=keys,validate='one_to_one').merge(new_section,on=keys,validate='one_to_one')
    a['delta_total']=a.mean_auc-a.old_section_mean_auc;a['delta_group_holdout']=a.new_section_mean_auc-a.old_section_mean_auc;a['delta_patient_weighting']=a.mean_auc-a.new_section_mean_auc
    f.to_csv(OUT/'all_section_metrics.tsv.gz',sep='\t',index=False);u.to_csv(OUT/'all_unit_metrics.tsv.gz',sep='\t',index=False);a.to_csv(OUT/'all_program_summary.tsv',sep='\t',index=False)
    selected=['B_AXIS','train_selected_GO','geometry_depth','GOBP_GAMMA_DELTA_T_CELL_ACTIVATION','GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY']
    a[a.program.isin(selected)&a.arm.isin(['full','train_selected_GO','geometry_depth'])].to_csv(OUT/'fixed_readout_comparison.tsv',sep='\t',index=False)
    # c_2 pair's actual old-v-new fold metrics, without pooling query points.
    pair=f[f.evaluation_unit.eq('GSE175540:c_2')&f.split.eq('leave_patient')]
    old_pair=old[old.section_id.isin(pair.section_id)&old.split.eq('leave_patient')]
    pair=pair.merge(old_pair[keys+['section_id','auc']].rename(columns={'auc':'old_auc'}),on=keys+['section_id'],validate='one_to_one');pair['delta_auc']=pair.auc-pair.old_auc;pair.to_csv(OUT/'c2_pair_all_program_comparison.tsv.gz',sep='\t',index=False)
    R.M.js(OUT/'complete.json',dict(status='COMPLETE_GROUPED_INTERNAL_VALIDATION',jobs=len(audit),refit_jobs=3,reused_jobs=54,metric_rows=len(f),summary_rows=len(a),GSE_patients=17,USZ_unit='tumor_sample_unknown_patient',scientific_confirmation=False,fold_audit=audit))
    print(a[a.program.isin(['B_AXIS','train_selected_GO'])&a.arm.isin(['full','train_selected_GO'])][['mask_um','split','cohort','program','mean_auc','n_evaluable_units','delta_total']].to_string(index=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','train','summarize']);args=p.parse_args();globals()[args.stage]()
