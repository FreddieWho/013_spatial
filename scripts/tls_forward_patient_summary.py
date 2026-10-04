#!/usr/bin/env python3
"""D-186: descriptive section -> patient mean -> cohort median, fixed conditions."""
from pathlib import Path
import importlib.util,json,csv,hashlib,warnings,shutil
import numpy as np
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'infra/tls_multiscale_gobp_20261002/integrated';OUT=ROOT/'infra/tls_forward_patient_20261003'
IDENTITY=ROOT/'infra/bioinf-data-index/raw/user_gse175540_20261003/current_26_source_identity_confirmed.tsv'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return list(csv.DictReader(Path(p).open(),delimiter='\t'))
def write(p,rows):
    with Path(p).open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def mean_finite(x):
    n=np.isfinite(x).sum(axis=0);return np.divide(np.nansum(x,axis=0),n,out=np.full(x.shape[1:],np.nan),where=n>0),n
def aggregate(meta,values,keys,mapping):
    """No case borrowing: every key, including scale and family, stays fixed."""
    groups={};seen=set();units_total={}
    for sid,g in mapping.items():units_total[g]=units_total.get(g,0)+1
    for i,row in enumerate(meta):
        sid=row['section_id'];cond=tuple(row[k] for k in ['cohort']+keys)
        assert (sid,cond) not in seen;seen.add((sid,cond))
        groups.setdefault(cond+(mapping[sid],),[]).append(i)
    unitrows=[];unitvalues={k:[] for k in values};counts={k:[] for k in values}
    for group,ix in groups.items():
        uid=group[-1];row=dict(zip(['cohort']+keys,group[:-1]));row.update(unit_row=len(unitrows),evaluation_unit=uid,unit_type='patient' if row['cohort']=='GSE175540' else 'tumor_sample_unknown_patient',n_sections_with_condition=len(ix),n_sections_in_unit=units_total[uid],source_rows=','.join(map(str,ix)),section_ids=';'.join(meta[i]['section_id'] for i in ix))
        unitrows.append(row)
        for k,v in values.items():
            mean,n=mean_finite(v[ix]);unitvalues[k].append(mean);counts[k].append(n)
    unitvalues={k:np.asarray(v) for k,v in unitvalues.items()};counts={k:np.asarray(v,dtype=np.int16) for k,v in counts.items()}
    cases={}
    for i,row in enumerate(unitrows):cases.setdefault(tuple(row[k] for k in ['cohort']+keys),[]).append(i)
    summary=[];arrays={}
    for k in values:
        for suffix in ['median','old_section_median','n_units','n_positive_units','n_negative_units','n_complete_units']:arrays[k+'_'+suffix]=[]
    for condition,ix in cases.items():
        originals=[i for u in ix for i in map(int,unitrows[u]['source_rows'].split(','))]
        row=dict(zip(['cohort']+keys,condition));row.update(summary_row=len(summary),n_units_with_condition=len(ix),n_sections_with_condition=len(originals),unit_type='patient' if row['cohort']=='GSE175540' else 'tumor_sample_unknown_patient',status='DESCRIPTIVE_ONLY_NOT_SPATIAL_CONFIRMATION');summary.append(row)
        for k,v in unitvalues.items():
            a=v[ix];total=np.array([unitrows[i]['n_sections_in_unit'] for i in ix])[:,None]
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning)
                arrays[k+'_median'].append(np.nanmedian(a,axis=0));arrays[k+'_old_section_median'].append(np.nanmedian(values[k][originals],axis=0))
            arrays[k+'_n_units'].append(np.isfinite(a).sum(axis=0));arrays[k+'_n_positive_units'].append((a>0).sum(axis=0));arrays[k+'_n_negative_units'].append((a<0).sum(axis=0));arrays[k+'_n_complete_units'].append(((counts[k][ix]==total)&np.isfinite(a)).sum(axis=0))
    return unitrows,unitvalues,counts,summary,{k:np.asarray(v) for k,v in arrays.items()}

def bidirectional_update():
    reverse=read(ROOT/'infra/tls_patient_validation_20261003/all_program_summary.tsv')
    reverse={(r['cohort'],r['program'],r['mask_um'],r['arm']):r for r in reverse if r['split']=='leave_cohort'}
    index=read(OUT/'contrasts_summary_index.tsv');z={k:v for k,v in np.load(OUT/'contrasts_summary.npz').items()};programs={str(p):i for i,p in enumerate(z['set_ids'])}
    selected={}
    for i,r in enumerate(index):
        if r['family']!='halo' or r['kind']!='adjacent_distance':continue
        first,second=('130.0','230.0') if r['mask_um']=='130' else ('230.0','400.0')
        if r['a_lo_um']==first and r['b_lo_um']==second:selected[(r['cohort'],r['mask_um'])]=i
    rows=read(BASE/'all_program_bidirectional_screen.tsv')
    for r in rows:
        j=programs[r['program']]
        for cohort in ['USZ','GSE175540']:
            i=selected[(cohort,r['forward_core_margin_um'])]
            for field,key in [('n','n_units'),('positive','n_positive_units'),('negative','n_negative_units'),('median_sd','median')]:
                value=z['effect_sd_'+key][i,j];r[cohort+'_forward_'+field]=float(value) if np.isfinite(value) else ''
            for radius in [230,500,1000]:
                for arm in ['radial','angular','full']:
                    v=reverse[(cohort,r['program'],str(radius),arm)]
                    r[f'{cohort}_inverse_R{radius}_{arm}_auc']=v['mean_auc'];r[f'{cohort}_inverse_R{radius}_{arm}_n']=v['n_evaluable_units']
        r['GSE_evaluation_unit']='patient';r['USZ_evaluation_unit']='tumor_sample_unknown_patient';r['formal_p']='';r['inference']='DESCRIPTIVE_FORWARD_AND_GROUPED_INTERNAL_REVERSE_NOT_CONFIRMATION'
    assert len(rows)==1742*2
    write(OUT/'all_program_bidirectional_patient_screen.tsv',rows)

def main():
    OUT.mkdir(exist_ok=True);assert shutil.disk_usage(OUT).free>=1.2e12
    identity=read(IDENTITY);mapping={r['section_id']:r['patient_id'] if r['cohort']=='GSE175540' else 'sample:'+r['section_id'] for r in identity};assert all(mapping.values())
    specs=[('bins','all_bins.tsv','halo_edge_bins.npz',['family','mask_um','lo_um','hi_um','compartment'],['means']),('contrasts','all_contrasts.tsv','halo_edge_contrasts.npz',['family','kind','mask_um','compartment','a_lo_um','a_hi_um','b_lo_um','b_hi_um'],['effect_sd','reference_mean','reference_sd','native_quantile']),('bridges','bridge_groups.tsv','bridge_all_programs.npz',['family','mode','mask_um','halfwidth_um','length_stratum','target_family'],['raw_median','adjusted_median','positive_continuity','negative_continuity','side_asymmetry','native_z_median','native_quantile_median','tissue_z_median'])]
    contract={'decision':'D-186','identity_sha256':sha(IDENTITY),'input_sha256':{str(BASE/n):sha(BASE/n) for _,a,b,_,_ in specs for n in [a,b]},'rule':'Within-patient finite section mean at identical condition; cohort median over available patient units. USZ tumor samples. Missing is not zero.','formal_p':None}
    target=OUT/'contract.json'
    if target.exists():assert json.loads(target.read_text())==contract
    else:target.write_text(json.dumps(contract,indent=2)+'\n')
    receipts=[];examples=[]
    for family,mfile,afile,keys,metrics in specs:
        meta=read(BASE/mfile);z=np.load(BASE/afile);names=z['set_ids'];assert len(names)==1742
        u,v,n,s,a=aggregate(meta,{k:z[k] for k in metrics},keys,mapping)
        write(OUT/f'{family}_unit_index.tsv',u);write(OUT/f'{family}_summary_index.tsv',s)
        np.savez_compressed(OUT/f'{family}_unit_values.npz',set_ids=names,**v,**{k+'_n_sections':x for k,x in n.items()})
        np.savez_compressed(OUT/f'{family}_summary.npz',set_ids=names,**a)
        # USZ must be numerically unchanged; each section is one explicit sample unit.
        ix=[i for i,r in enumerate(s) if r['cohort']=='USZ']
        for k in metrics:np.testing.assert_allclose(a[k+'_median'][ix],a[k+'_old_section_median'][ix],equal_nan=True)
        for i,r in enumerate(s):
            for name in ['B_AXIS','GOBP_GAMMA_DELTA_T_CELL_ACTIVATION','GOBP_B_CELL_RECEPTOR_SIGNALING_PATHWAY']:
                j=int(np.flatnonzero(names==name)[0])
                for k in metrics:
                    if k not in ['means','effect_sd','raw_median','adjusted_median']:continue
                    examples.append(dict(output_family=family,summary_row=i,cohort=r['cohort'],condition=json.dumps({x:r[x] for x in keys},sort_keys=True),program=name,metric=k,median=float(a[k+'_median'][i,j]) if np.isfinite(a[k+'_median'][i,j]) else '',old_section_median=float(a[k+'_old_section_median'][i,j]) if np.isfinite(a[k+'_old_section_median'][i,j]) else '',n_units=int(a[k+'_n_units'][i,j]),n_positive=int(a[k+'_n_positive_units'][i,j]),n_negative=int(a[k+'_n_negative_units'][i,j]),n_complete_units=int(a[k+'_n_complete_units'][i,j])))
        receipts.append(dict(family=family,input_section_rows=len(meta),unit_rows=len(u),summary_rows=len(s),programs=len(names),metrics=metrics));print(receipts[-1],flush=True)
    write(OUT/'fixed_readout_examples.tsv',examples)
    bidirectional_update()
    (OUT/'complete.json').write_text(json.dumps({'status':'COMPLETE_DESCRIPTIVE_PATIENT_AGGREGATION','receipts':receipts,'GSE_patients':17,'USZ_tumor_samples':8,'formal_confirmation':False,'calibration_status':'UNCHANGED_NOT_PASSED','missing_cases':'original source eligibility retained; index lists observed conditions only'},indent=2)+'\n')
if __name__=='__main__':main()
