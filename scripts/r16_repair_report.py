#!/usr/bin/env python3
"""Collect completed D-140/D-141 receipts; fail if any required replay is incomplete."""
import json
from collections import Counter,defaultdict
import numpy as np
from r16.recovery.repair_pipeline import OUT,readtable,dump,table,digest,patient_interval


def collect():
    receipts={}
    for name in ['six_raw','six_depth_normalized','go','mask','spline','pc']:
        receipt=json.load(open(OUT/(name+'_receipt.json')))
        assert receipt['status']=='COMPLETED',name
        receipts[name]=receipt
        if name in ['six_raw','six_depth_normalized','go']:
            count=22 if name=='six_depth_normalized' else 69
            dump(name+'_progress.json',dict(status='COMPLETED',completed_sections=count,total_sections=count))
    assert receipts['go']['programs']==6870
    defs=json.load(open(OUT/'go_representatives.json'));aliases=json.load(open(OUT/'go_aliases.json'))
    for a,d in aliases.items():
        rep=defs[d['representative']]
        assert (d['input_genes'],d['readout_genes'],d['definition_hash'])==(rep['input_genes'],rep['readout_genes'],rep['definition_hash']),a
    pred=readtable(OUT/'go_prediction.tsv');missing=readtable(OUT/'go_not_testable.tsv')
    discovery=readtable(OUT/'go_discovery.tsv');molecular=readtable(OUT/'go_molecular.tsv')
    assert sum(int(r['n_sections']) for r in discovery)+sum(r['cohort']=='DISCOVERY' for r in missing)==6870*47
    assert len(molecular)==6870*22
    missing_external={r['program'] for r in missing if r['cohort']!='DISCOVERY'}
    assert len(pred)==(6870-len(missing_external))*15
    selection=json.load(open(OUT/'spline_selection.json'))
    result={'receipts':receipts,'alias_definitions_checked':len(aliases),
            'go':{'representatives':6870,'molecular_section_rows':len(molecular),'prediction_patient_rows':len(pred),
                  'prediction_status_counts':dict(Counter(r['status'] for r in pred)),
                  'section_ineligibility':dict(Counter(r['cohort']+':'+r['status'] for r in missing)),
                  'original_both_cohort_hits':len(selection['original_hits']),
                  'corrected_both_cohort_hits':len(selection['corrected_hits'])}}
    result['six']={mode:readtable(OUT/f'six_{mode}_summary.tsv') for mode in ['raw','depth_normalized']}
    markers=[]
    for mode in ['raw','depth_normalized']:
        p={(r['cohort'],r['program'],r['patient']):float(r['mse_M1']) for r in readtable(OUT/f'six_{mode}_prediction.tsv') if r['status']=='COMPUTED'}
        groups=defaultdict(list)
        for r in readtable(OUT/f'six_{mode}_single_marker.tsv'):
            key=(r['cohort'],r['program'],r['patient']);mse=float(r['heldout_mse'])
            if key in p and mse>1e-12:groups[key[:2]].append(100*(mse-p[key])/mse)
        for (c,g),v in groups.items():
            ci=patient_interval(v);markers.append(dict(mode=mode,cohort=c,program=g,n_patients=len(v),median_program_improvement_over_selected_marker=float(np.median(v)),ci_low=ci[0],ci_high=ci[1]))
    result['single_marker']=markers;table('single_marker_comparison.tsv',markers)
    result['masked_field']=readtable(OUT/'masked_field_summary.tsv')
    result['structure']={'patient_results':readtable(OUT/'masked_structure_patients.tsv'),
                         'window_label_counts':dict(Counter(r['label'] for r in readtable(OUT/'masked_structure_windows.tsv'))),
                         'fine_control_status_counts':dict(Counter(r['status'] for r in readtable(OUT/'l009_focus_controls.tsv')))}
    curves=readtable(OUT/'spline_results.tsv');shape_counts=Counter();curves_by=defaultdict(list)
    for r in curves:
        shape_counts[r['route']+':'+r['shape']]+=1;curves_by[(r['program'],r['route'])].append(r)
    for g,d in json.load(open(OUT/'six_programs.json')).items():
        if d['missing_genes']:
            assert all(r['status']=='MISSING_PROGRAM_GENES' for route in ['raw','rank','resid','rank_resid'] for r in curves_by[(g,route)])
    consistency=[]
    for route in ['raw','rank','resid','rank_resid']:
        for group,ids in [('original820',set(selection['original_hits'])),('union',set(selection['union']))]:
            same=0;all_computed=0
            for g in ids:
                v=curves_by[(g,route)];assert len(v)==8
                if all(r['status']=='COMPUTED' for r in v):all_computed+=1
                shapes={r['shape'] for r in v}
                if len(shapes)==1 and shapes<= {'low-to-high','high-to-low'}:same+=1
            consistency.append(dict(route=route,selection=group,n_programs=len(ids),all8_computed=all_computed,all8_same_monotonic_direction=same))
    result['spline']={'shape_counts':dict(shape_counts),'direction_consistency_descriptive':consistency}
    table('spline_direction_summary.tsv',consistency)
    result['r04']=json.load(open(OUT/'r04_final_gate.json'))
    dump('final_results.json',result)
    return result


if __name__=='__main__':
    result=collect()
    print(json.dumps({k:v for k,v in result.items() if k in ['go','spline','structure','single_marker']},indent=2))
