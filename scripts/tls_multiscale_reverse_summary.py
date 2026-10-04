#!/usr/bin/env python3
"""All-program held-out summaries with paired section baselines and no p-values."""
from pathlib import Path
import importlib.util,json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'infra/tls_multiscale_gobp_20261002'
spec=importlib.util.spec_from_file_location('m',ROOT/'scripts/tls_multiscale_gobp.py')
M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M)
# Load the project's SciPy stack before pandas/Arrow selects a C++ runtime.
import pandas as pd


def main():
    directory=OUT/'reverse';done=json.loads((directory/'model_complete.json').read_text())
    assert done['folds']==84
    paths=sorted((directory/'models').glob('*_metrics.tsv'));assert len(paths)==84
    frame=pd.concat([pd.read_csv(p,sep='\t') for p in paths],ignore_index=True)
    keys=['mask_um','split','heldout_job','section_id','cohort']
    assert not frame.duplicated(keys+['arm','program']).any()
    assert len(frame)==3*2*26*(1742*3+2)
    assert set(frame.loc[frame.arm=='full','program'])==set(np.load(OUT/'scores/usz-KC1.npz')['set_ids'])
    frame['valid']=frame.status.eq('INTERNAL_SECTION_EVALUATION')
    frame.loc[~frame.valid,['auc','ap']]=np.nan
    # Pair comparisons within the same held-out section and operating point.
    base=frame.loc[frame.arm.eq('geometry_depth'),keys+['auc','ap']].rename(columns={'auc':'geometry_auc','ap':'geometry_ap'})
    b=frame.loc[frame.arm.eq('full')&frame.program.eq('B_AXIS'),keys+['auc','ap']].rename(columns={'auc':'B_axis_auc','ap':'B_axis_ap'})
    frame=frame.merge(base,on=keys,validate='many_to_one').merge(b,on=keys,validate='many_to_one')
    frame['delta_auc_geometry']=frame.auc-frame.geometry_auc
    frame['delta_auc_B_axis']=frame.auc-frame.B_axis_auc
    frame['delta_auc_chance']=frame.auc-.5
    # A descriptive conservative envelope, not a test-selected deployable model.
    frame['fixed_comparator_envelope']=np.maximum(.5,frame[['geometry_auc','B_axis_auc']].max(axis=1).fillna(.5))
    frame['delta_auc_comparator_envelope']=frame.auc-frame.fixed_comparator_envelope
    frame['prevalence']=frame.n_positive/(frame.n_positive+frame.n_negative).replace(0,np.nan)
    frame['ap_above_prevalence']=frame.ap-frame.prevalence
    dest=directory/'summary';dest.mkdir(exist_ok=True)
    groups=['mask_um','split','cohort','arm','program']
    summary=frame.groupby(groups,dropna=False).agg(
        attempted_sections=('section_id','nunique'),valid_sections=('valid','sum'),
        mean_auc=('auc','mean'),median_auc=('auc','median'),min_auc=('auc','min'),max_auc=('auc','max'),
        mean_ap=('ap','mean'),mean_ap_above_prevalence=('ap_above_prevalence','mean'),
        mean_delta_auc_geometry=('delta_auc_geometry','mean'),median_delta_auc_geometry=('delta_auc_geometry','median'),
        n_paired_geometry=('delta_auc_geometry','count'),mean_delta_auc_B_axis=('delta_auc_B_axis','mean'),
        n_paired_B_axis=('delta_auc_B_axis','count'),mean_delta_auc_chance=('delta_auc_chance','mean'),
        mean_delta_auc_comparator_envelope=('delta_auc_comparator_envelope','mean')).reset_index()
    summary['status']=np.where(summary.valid_sections>0,'DESCRIPTIVE_INTERNAL_SECTION_VALIDATION','NOT_TESTABLE')
    summary.to_csv(dest/'all_programs.tsv',sep='\t',index=False)
    frame.to_csv(dest/'paired_section_metrics.tsv.gz',sep='\t',index=False,compression='gzip')
    statuses=frame.groupby(groups+['status']).size().rename('n_sections').reset_index()
    statuses.to_csv(dest/'all_program_statuses.tsv',sep='\t',index=False)
    # The directional increment is paired, not a difference of differently supported aggregates.
    arms=frame[frame.arm.isin(['full','radial','angular'])].pivot(index=keys+['program'],columns='arm',values='auc').reset_index()
    arms['full_minus_radial']=arms['full']-arms['radial'];arms['angular_minus_radial']=arms['angular']-arms['radial']
    arms.to_csv(dest/'paired_directional_increment.tsv.gz',sep='\t',index=False,compression='gzip')
    arm_summary=arms.groupby(['mask_um','split','cohort','program']).agg(
        n_paired=('full_minus_radial','count'),mean_full_minus_radial=('full_minus_radial','mean'),
        median_full_minus_radial=('full_minus_radial','median'),mean_angular_minus_radial=('angular_minus_radial','mean')).reset_index()
    arm_summary.to_csv(dest/'all_program_directional_increment.tsv',sep='\t',index=False)
    support=frame[frame.arm.eq('geometry_depth')][keys+['n_positive','n_negative','status']]
    support.to_csv(dest/'heldout_label_support.tsv',sep='\t',index=False)
    queries=pd.concat([pd.read_csv(p,sep='\t') for p in sorted(directory.glob('*_queries.tsv'))],ignore_index=True)
    queries.groupby(['section_id','cohort','mask_um','centre_annotation','status'],dropna=False).size().rename('n_queries').reset_index().to_csv(dest/'query_qualification.tsv',sep='\t',index=False)
    selection=[]
    for job in done['jobs']:
        for rank,name in enumerate(job['selected_GO'],1):
            selection.append(dict(mask_um=job['mask_um'],split=job['split'],job=job['job'],training_rank=rank,program=name))
    M.tsv(dest/'training_selected_GO.tsv',selection)
    M.js(dest/'receipt.json',dict(status='FULL_INTERNAL_DESCRIPTIVE_SUMMARY',folds=84,metric_rows=len(frame),
         summary_rows=len(summary),n_gobp=1691,n_controls=51,patient_confirmation=False,spatial_p='NOT_CLAIMED',
         location_probability_calibration='NOT_ESTABLISHED',test_set_selection=False,
         interpretation='Report all masks/arms/cohorts. Inspecting this table does not create an independent validation set.'))
    print('full reverse summary',len(frame),'metric rows',len(summary),'aggregate rows',flush=True)


if __name__=='__main__':main()
