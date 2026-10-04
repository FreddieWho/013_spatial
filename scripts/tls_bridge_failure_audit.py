#!/usr/bin/env python3
"""D-189: diagnose frozen bridge controls, without rescuing failed tests."""
from pathlib import Path
import csv,json,hashlib,importlib.util
import numpy as np
from scipy import sparse
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'infra/tls_multiscale_gobp_20261002';CAL=BASE/'bridge_calibration';OUT=ROOT/'infra/tls_bridge_failure_audit_20261003'
s=importlib.util.spec_from_file_location('bridge',ROOT/'scripts/tls_multiscale_bridge_calibration.py');C=importlib.util.module_from_spec(s);s.loader.exec_module(C)
def read(p):return list(csv.DictReader(Path(p).open(),delimiter='\t'))
def write(p,rows):
    with Path(p).open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    OUT.mkdir(exist_ok=True);null=read(CAL/'null_calibration.tsv');power=read(CAL/'power_calibration.tsv');sources=json.loads((BASE/'contract.json').read_text())['sources'];source_rows=[];angles=[];normrows=[];template_rows=[];checked=0;hashes={}
    allpass={sid:all(r['passed']=='True' for r in null if r['section_id']==sid) for sid in {r['section_id'] for r in null}};cellpass={(r['section_id'],r['model']):r['passed']=='True' for r in null}
    for info in sources:
        sid=info['section_id'];d=CAL/sid;receipt=json.loads((d/'operators_complete.json').read_text());done=json.loads((d/'calibration_complete.json').read_text());nr=[r for r in null if r['section_id']==sid]
        source_rows.append(dict(**info,status=done['status'],n_common_geometry=receipt['n_geometry'],n_attempted_geometry=receipt.get('attempted_geometry',0),n_null_settings=len(nr),n_failed_null=sum(r['passed']!='True' for r in nr),min_FP_rate=min([float(r['rate']) for r in nr],default=''),max_FP_rate=max([float(r['rate']) for r in nr],default='')))
        if not receipt['n_geometry']:continue
        for r in nr:
            path=d/f"{r['model']}_null.npz";hashes[str(path.relative_to(ROOT))]=sha(path);z=np.load(path);T=z['max_by_direction'];assert T.shape==(40,400)
            ranks=np.array([(T>=T[i]).sum(axis=0)/40 for i in range(40)])
            np.testing.assert_array_equal(ranks[0],z['diagnostic_rank']);assert int((ranks[0]<=.05).sum())==int(r['n_false_positive'])
            rejection=(ranks<=.05).mean(axis=1);assert rejection.mean()<=.05+1e-12
            for i in range(40):angles.append(dict(section_id=sid,model=r['model'],direction_index=i,angle_degrees=4.5*i,identity_direction=i==0,rate=float(rejection[i]),mean_max_statistic=float(T[i].mean()),median_max_statistic=float(np.median(T[i])),direction_average_rate=float(rejection.mean()),original_setting_passed=r['passed'],n_simulations=400))
        op=C.load_operators(sid);n=op['n'];nref=len(op['ref_slots']);meanop=sparse.kron(sparse.eye(nref,format='csr'),np.ones((1,5))/5,format='csr');U=(meanop@op['RW']+op['L']).tocsr();var=np.empty(n*40);var[np.arange(n)*40]=np.asarray(op['A'].multiply(op['A']).sum(axis=1)).ravel()
        for group in np.unique(op['ref_groups']):
            indices=np.flatnonzero(op['ref_groups']==group);Z=op['Z'][group*13:(group+1)*13];gram=Z@Z.T
            for start in range(0,len(indices),256):
                ix=indices[start:start+256];a=U[ix];coef=op['coefficients'][ix];cross=np.asarray(a@Z.T);v=np.asarray(a.multiply(a).sum(axis=1)).ravel()-2*np.sum(cross*coef,axis=1)+np.einsum('ni,ij,nj->n',coef,gram,coef)
                assert np.all(np.isfinite(v)) and np.min(v)>-1e-6
                var[op['ref_slots'][ix]]=np.maximum(v,0)
            # Independent direct-vector check, without materializing all references.
            i=indices[0];dense=U[i].toarray().ravel()-op['coefficients'][i]@Z;expected=float(dense@dense);actual=var[op['ref_slots'][i]];assert np.isclose(actual,expected,rtol=2e-5,atol=1e-8);checked+=1
        v=var.reshape(n,40);geoms=read(d/'common_geometries.tsv');assert len(geoms)==n
        for i,r in enumerate(geoms):
            med=float(np.median(v[i,1:]));normrows.append(dict(section_id=sid,effect_row=r['effect_row'],family=r['family'],mode=r['mode'],mask_um=r['mask_um'],halfwidth_um=r['width_um'],length_um=r['length_um'],length_stratum=r['length_stratum'],n_visible=r['n_visible'],actual_operator_variance=float(v[i,0]),median_reference_operator_variance=med,actual_to_reference_variance_ratio=float(v[i,0]/med) if med>1e-15 else '',min_direction_variance=float(v[i].min()),max_direction_variance=float(v[i].max()),diagnostic_model='unit_independent_noise_not_real_spatial_covariance'))
        reps=np.load(d/'power_representatives.npy')
        H=np.column_stack([(np.asarray(op['W'][g*5:(g+1)*5].sum(axis=0)).ravel()>0).astype(float) for g in reps])
        for j,g in enumerate(reps):H[:,j]/=float(np.mean(op['W'][g*5:(g+1)*5]@H[:,j]))
        dR,dE=C.evaluate(H,op)
        for j,g in enumerate(reps):
            maximum,own=C.selected_statistics(dE[:,:,j:j+1],(dR[:,:,:,j:j+1]>0).sum(axis=2),(dR[:,:,:,j:j+1]<0).sum(axis=2),target=int(g))
            rank=float(C.diagnostic_rank(maximum,own)[0]);r=geoms[g]
            template_rows.append(dict(section_id=sid,effect_row=r['effect_row'],family=r['family'],mode=r['mode'],halfwidth_um=r['width_um'],length_um=r['length_um'],length_stratum=r['length_stratum'],target_raw_gain=float(dR[g,0,:,j].mean()),target_adjusted_gain=float(dE[g,0,j]),signal_only_diagnostic_rank=rank,signal_only_target_selected=rank<=.05 and float(own[0])>0,all_source_null_settings_passed=allpass[sid],interpretation='deterministic original template, not real-data power or p'))
        print('audited directions and operators',sid,n,flush=True)
    write(OUT/'signal_only_template_diagnostics.tsv',template_rows)
    for r in power:
        r['matched_null_setting_passed']=cellpass[(r['section_id'],r['model'])];r['all_source_null_settings_passed']=allpass[r['section_id']];r['power_and_matching_null_passed']=r['passed']=='True' and r['matched_null_setting_passed'];r['power_and_whole_source_null_passed']=r['passed']=='True' and r['all_source_null_settings_passed']
    write(OUT/'source_coverage.tsv',source_rows);write(OUT/'all_direction_rejection_rates.tsv',angles);write(OUT/'all_geometry_operator_variance.tsv',normrows);write(OUT/'power_with_null_qualification.tsv',power)
    import pandas as pd
    f=pd.DataFrame(null);f['rate']=f.rate.astype(float);f['failed']=f.passed.ne('True')
    f.groupby('model').agg(n_settings=('rate','size'),n_failed=('failed','sum'),mean_setting_rate=('rate','mean'),min_rate=('rate','min'),max_rate=('rate','max')).reset_index().to_csv(OUT/'null_by_field.tsv',sep='\t',index=False)
    pf=pd.DataFrame(power);pf['power']=pf.power.astype(float);pf['power_passed']=pf.passed.eq('True')
    groupcols=['family','mode','mask_um','length_stratum','halfwidth_um','amplitude_sd','model']
    pf.groupby(groupcols).agg(n_settings=('power','size'),mean_setting_power=('power','mean'),min_power=('power','min'),max_power=('power','max'),n_power_passed=('power_passed','sum'),n_power_and_matching_null_passed=('power_and_matching_null_passed','sum'),n_power_and_whole_source_null_passed=('power_and_whole_source_null_passed','sum')).reset_index().to_csv(OUT/'power_by_geometry_and_field.tsv',sep='\t',index=False)
    for cols,name in [(['amplitude_sd','model'],'power_by_field'),(['amplitude_sd','family'],'power_by_family'),(['amplitude_sd','halfwidth_um'],'power_by_width'),(['amplitude_sd','length_stratum'],'power_by_length')]:
        pf.groupby(cols).agg(n_settings=('power','size'),mean_setting_power=('power','mean'),n_power_passed=('power_passed','sum'),n_power_and_matching_null_passed=('power_and_matching_null_passed','sum'),n_power_and_whole_source_null_passed=('power_and_whole_source_null_passed','sum')).reset_index().to_csv(OUT/(name+'.tsv'),sep='\t',index=False)
    # Conditions simultaneously passing all six pre-existing fields, not posthoc selected fields.
    g=pf.groupby(['section_id','effect_row','amplitude_sd']).agg(n_fields=('model','nunique'),all_power_fields_passed=('power_passed','all'),source_null_passed=('all_source_null_settings_passed','all')).reset_index();assert (g.n_fields==6).all();g['both_necessary_gates_passed']=g.all_power_fields_passed & g.source_null_passed;g.to_csv(OUT/'six_field_joint_necessary_gate.tsv',sep='\t',index=False)
    hf=pf[pf.amplitude_sd.eq('0.5')];ratios=np.array([r['actual_to_reference_variance_ratio'] for r in normrows if r['actual_to_reference_variance_ratio']!=''])
    complete=dict(status='COMPLETE_FROZEN_FAILURE_DIAGNOSIS',sources=26,null_settings=len(null),null_failed=int(f.failed.sum()),direction_rows=len(angles),geometry_norm_rows=len(normrows),direct_vector_checks=checked,half_sd_settings=len(hf),half_sd_power_passed=int(hf.power_passed.sum()),half_sd_power_and_matching_null_passed=int(hf.power_and_matching_null_passed.sum()),half_sd_power_and_source_null_passed=int(hf.power_and_whole_source_null_passed.sum()),half_sd_representatives_passing_all6_power_and_source_null=int(g[g.amplitude_sd.eq('0.5')].both_necessary_gates_passed.sum()),operator_variance_ratio_quantiles=dict(zip(['min','q25','median','q75','max'],np.quantile(ratios,[0,.25,.5,.75,1]).tolist())),signal_only_templates=len(template_rows),signal_only_templates_selected=sum(r['signal_only_target_selected'] for r in template_rows),formal_p_allowed=False,method_or_threshold_changed=False)
    (OUT/'complete.json').write_text(json.dumps(complete,indent=2)+'\n');(OUT/'contract.json').write_text(json.dumps(dict(decision='D-189',null_table_sha256=sha(CAL/'null_calibration.tsv'),power_table_sha256=sha(CAL/'power_calibration.tsv'),saved_simulation_sha256=hashes,real_GO_expression_used=False,source_patient_status='GSE known patients; USZ tumor samples only; no pooled patient inference'),indent=2)+'\n');print(complete,flush=True)
if __name__=='__main__':main()
