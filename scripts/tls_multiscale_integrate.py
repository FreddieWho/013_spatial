#!/usr/bin/env python3
"""All-program evidence matrices; preserve dimensions and failed inference gates."""
from pathlib import Path
import argparse,csv,gzip,importlib.util,json,warnings
import numpy as np

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'infra/tls_multiscale_gobp_20261002';DEST=OUT/'integrated'
def mod(name):
    s=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
M=mod('tls_multiscale_gobp');H=mod('tls_multiscale_halo_edge')


def halo_edge():
    DEST.mkdir(exist_ok=True);sources=json.loads((OUT/'contract.json').read_text())['sources']
    names=np.load(OUT/'scores/usz-KC1.npz')['set_ids'];bins=[];means=[];cases=[];effects=[];refmeans=[];refsd=[];quantiles=[];counts=[]
    reports=[];null=[];power=[];diagnostic=[]
    for source in sources:
        sid=source['section_id']
        for family in ['halo','edge']:
            d=OUT/'halo_edge'/sid/family;meta=json.loads((d/'operators_complete.json').read_text());score=json.loads((d/'scoring_complete.json').read_text());cal=json.loads((d/'calibration_complete.json').read_text())
            reports.append(dict(**source,family=family,**meta,scoring_status=score['status'],calibration_status=cal['status']))
            for f in sorted(d.glob('*_null.tsv')):null.extend(M.tab(f))
            for f in sorted(d.glob('*_power.tsv')):power.extend(M.tab(f))
            diagnostic.extend(dict(**source,family=family,**r,calibration_status=cal['status']) for r in M.tab(d/'all_program_diagnostics.tsv'))
            if not meta['n_bins']:continue
            z=np.load(d/'full_panel.npz');geometry=np.load(d/'geometry.npz');np.testing.assert_array_equal(z['set_ids'],names)
            binrows=M.tab(d/'bins.tsv');C=M.tab(d/'cases.tsv');E=z['effect_sd'];R=E[:,1:];valid=np.isfinite(R);n=valid.sum(axis=1)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning)
                mu=np.nanmean(R,axis=1);sd=np.nanstd(R,axis=1,ddof=1)
            observed=E[:,0];q=np.divide(np.sum(R<observed[:,None,:],axis=1)+.5*np.sum(R==observed[:,None,:],axis=1),n,out=np.full_like(observed,np.nan),where=(n>0)&np.isfinite(observed))
            for i,r in enumerate(binrows):
                bins.append(dict(bin_row=len(bins),**source,family=family,local_bin=i,**r,n_geometry_references=int(sum(geometry['status'][1:,i]==0))))
                means.append(z['bin_means'][0,i])
            for i,r in enumerate(C):
                cases.append(dict(contrast_row=len(cases),**source,family=family,**r,calibration_status=cal['status'],n_common_references=meta['n_common_references']))
                effects.append(observed[i]);refmeans.append(mu[i]);refsd.append(sd[i]);quantiles.append(q[i]);counts.append(n[i])
            print('integrated profiles',sid,family,len(binrows),len(C),flush=True)
    H.tsv(DEST/'all_bins.tsv',bins);H.tsv(DEST/'all_contrasts.tsv',cases)
    np.savez_compressed(DEST/'halo_edge_bins.npz',means=np.array(means,dtype=np.float32),set_ids=names)
    np.savez_compressed(DEST/'halo_edge_contrasts.npz',effect_sd=np.array(effects,dtype=np.float32),reference_mean=np.array(refmeans,dtype=np.float32),reference_sd=np.array(refsd,dtype=np.float32),native_quantile=np.array(quantiles,dtype=np.float32),n_references=np.array(counts,dtype=np.int16),set_ids=names)
    H.tsv(DEST/'halo_edge_source_completion.tsv',reports);H.tsv(DEST/'halo_edge_null_calibration.tsv',null);H.tsv(DEST/'halo_edge_power_calibration.tsv',power);H.tsv(DEST/'halo_edge_all_program_diagnostics.tsv',diagnostic)
    with gzip.open(DEST/'all_source_program_contrasts.tsv.gz','wt') as f:
        fields=list(cases[0])+['program','effect_sd','reference_mean','reference_sd','native_quantile','n_references','program_status','formal_p']
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
        for i,c in enumerate(cases):
            for j,name in enumerate(names):
                row=dict(c,program=str(name),n_references=int(counts[i][j]),formal_p='')
                for field,values in [('effect_sd',effects),('reference_mean',refmeans),('reference_sd',refsd),('native_quantile',quantiles)]:
                    row[field]=float(values[i][j]) if np.isfinite(values[i][j]) else ''
                row['program_status']='DESCRIPTIVE_ONLY' if np.isfinite(effects[i][j]) else ('NOT_TESTABLE_BINS' if c['status']!='EVALUABLE' else 'NOT_TESTABLE_PROGRAM_COVERAGE_OR_VARIANCE')
                w.writerow(row)
    M.js(DEST/'halo_edge_complete.json',dict(status='ALL26_ALL1742_FULL_DESCRIPTIVE_INTEGRATION',sources=26,source_families=len(reports),n_bins=len(bins),n_contrasts=len(cases),n_long_contrast_rows=len(cases)*len(names),n_gobp=1691,n_controls=51,formal_p_allowed=False))


def bridge():
    DEST.mkdir(exist_ok=True);contract=json.loads((OUT/'contract.json').read_text());names=np.load(OUT/'scores/usz-KC1.npz')['set_ids']
    metadata=[];results={key:[] for key in ['raw_median','adjusted_median','positive_continuity','negative_continuity','side_asymmetry','native_z_median','native_quantile_median','tissue_z_median','n_valid_geometry','n_native_geometry']}
    for source in contract['sources']:
        sid=source['section_id'];d=OUT/'forward'/sid;ctl=OUT/'controls'/sid;meta=json.loads((d/'operator_complete.json').read_text())
        if not meta['evaluable_geometries']:continue
        rows=M.tab(d/'effect_rows.tsv');pairs={(r['pair_id'],r['mode']):r for r in M.tab(d/'pairs.tsv')};qual=M.tab(ctl/'reference_qualification.tsv')
        arrays={k:np.load(d/f'{file}.npy',mmap_mode='r') for k,file in [('raw','raw_effect'),('adjusted','endpoint_residual'),('positive','positive_segments'),('negative','negative_segments'),('side','side_asymmetry')]}
        for key in ['native_z','native_quantile','tissue_z']:arrays[key]=np.load(ctl/f'{key}.npy',mmap_mode='r')
        q=np.array([r['reference_qualified']=='True' for r in qual]);tq=np.array([r['tissue_reference_qualified']=='True' for r in qual]);aspect=np.array([r['bridge_aspect_qualified']=='True' for r in qual])
        groups={}
        for row in rows:
            pair=pairs[row['pair_id'],row['mode']];key=tuple(row[k] for k in ['family','mode','mask_um','width_um','length_stratum'])+(pair['target_family'],)
            groups.setdefault(key,[]).append((int(row['effect_row']),pair['source_endpoint']))
        for key,items in sorted(groups.items()):
            components=sorted({c for _,c in items});values={k:[] for k in results}
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning)
                for comp in components:
                    ids=np.array([i for i,c in items if c==comp]);r=np.array(arrays['raw'][ids]);valid=np.isfinite(r);den=valid.sum(axis=0)
                    values['raw_median'].append(np.nanmedian(r,axis=0));values['adjusted_median'].append(np.nanmedian(arrays['adjusted'][ids],axis=0))
                    values['side_asymmetry'].append(np.nanmedian(np.abs(arrays['side'][ids]),axis=0));values['n_valid_geometry'].append(den)
                    for direction in ['positive','negative']:
                        values[direction+'_continuity'].append(np.divide(np.sum(valid&(arrays[direction][ids]>=4),axis=0),den,out=np.full(len(names),np.nan),where=den>0))
                    native=np.array(arrays['native_z'][ids]);native[~q[ids]]=np.nan
                    quantile=np.array(arrays['native_quantile'][ids]);quantile[~q[ids]]=np.nan
                    tissue=np.array(arrays['tissue_z'][ids]);tissue[~tq[ids]]=np.nan
                    values['native_z_median'].append(np.nanmedian(native,axis=0));values['native_quantile_median'].append(np.nanmedian(quantile,axis=0));values['tissue_z_median'].append(np.nanmedian(tissue,axis=0));values['n_native_geometry'].append(np.isfinite(native).sum(axis=0))
                for name in results:results[name].append(np.sum(values[name],axis=0) if name.startswith('n_') else np.nanmedian(values[name],axis=0))
            indices=[i for i,_ in items];metadata.append(dict(bridge_group=len(metadata),**source,family=key[0],mode=key[1],mask_um=key[2],halfwidth_um=key[3],length_stratum=key[4],target_family=key[5],n_source_components=len(components),n_geometry_rows=len(items),n_reference_qualified=int(q[indices].sum()),n_bridge_aspect=int(aspect[indices].sum()),status='DESCRIPTIVE_ONLY_SPATIAL_GATE_NOT_PASSED'))
        print('integrated bridge',sid,len(groups),flush=True)
    H.tsv(DEST/'bridge_groups.tsv',metadata)
    np.savez_compressed(DEST/'bridge_all_programs.npz',set_ids=names,**{k:np.array(v,dtype=np.int32 if k.startswith('n_') else np.float32) for k,v in results.items()})
    # The existing full forward summary must agree exactly in group order and values.
    old=np.load(OUT/'summary/bridge_all_programs.npz');np.testing.assert_allclose(results['raw_median'],old['raw_median'],rtol=1e-6,atol=1e-7,equal_nan=True)
    M.js(DEST/'bridge_complete.json',dict(status='ALL_FORWARD_GROUPS_WITH_NATIVE_REFERENCES',groups=len(metadata),n_gobp=1691,n_controls=51,formal_p_allowed=False))


def rank_and_index():
    names=np.load(OUT/'scores/usz-KC1.npz')['set_ids'];rows=M.tab(OUT/'reverse/summary/all_programs.tsv')
    lookup={(r['mask_um'],r['split'],r['cohort'],r['arm'],r['program']):r for r in rows};ranking=[]
    for radius in ['230','500','1000']:
        for arm in ['radial','angular','full']:
            group=[]
            for j,name in enumerate(names):
                source_rows=[lookup[radius,'leave_cohort',cohort,arm,str(name)] for cohort in ['USZ','GSE175540']]
                n=[int(r['valid_sections']) for r in source_rows]
                expected=[int(lookup[radius,'leave_cohort',cohort,'geometry_depth','geometry_depth']['valid_sections']) for cohort in ['USZ','GSE175540']]
                values=[float(r['mean_auc']) if r['mean_auc'] else np.nan for r in source_rows]
                valid=all(nn>0 and nn==ee for nn,ee in zip(n,expected)) and all(np.isfinite(values))
                group.append(dict(mask_um=radius,arm=arm,program=str(name),panel='GO_BP' if j<1691 else 'FIXED_CONTROL',USZ_auc=values[0] if np.isfinite(values[0]) else '',GSE175540_auc=values[1] if np.isfinite(values[1]) else '',USZ_sections=n[0],GSE175540_sections=n[1],worst_source_auc=min(values) if valid else '',rank='',status='EXPLORATORY_POST_VALIDATION' if valid else 'NOT_RANKED_INCOMPLETE_EVALUATION'))
            candidates=sorted([r for r in group if r['panel']=='GO_BP' and r['worst_source_auc']!=''],key=lambda r:(-r['worst_source_auc'],r['program']))
            for i,r in enumerate(candidates,1):r['rank']=i
            ranking.extend(group)
    H.tsv(DEST/'all_program_exploratory_ranking.tsv',ranking)
    curves=np.load(DEST/'halo_edge_bins.npz');contrasts=np.load(DEST/'halo_edge_contrasts.npz');bridges=np.load(DEST/'bridge_all_programs.npz')
    for z in [curves,contrasts,bridges]:np.testing.assert_array_equal(z['set_ids'],names)
    evidence=[];inverse_counts={}
    for r in rows:
        if int(r['valid_sections'])>0:inverse_counts[r['program']]=inverse_counts.get(r['program'],0)+1
    for j,name in enumerate(names):
        evidence.append(dict(program_index=j,program=str(name),panel='GO_BP' if j<1691 else 'FIXED_CONTROL',n_valid_distance_bins=int(np.isfinite(curves['means'][:,j]).sum()),n_valid_distance_contrasts=int(np.isfinite(contrasts['effect_sd'][:,j]).sum()),n_valid_bridge_groups=int(np.isfinite(bridges['raw_median'][:,j]).sum()),n_bridge_groups_with_qualified_references=int(np.isfinite(bridges['native_z_median'][:,j]).sum()),n_valid_inverse_aggregate_rows=inverse_counts.get(name,0),inference_status='INTERNAL_AND_DESCRIPTIVE_NOT_SPATIAL_CONFIRMATION',formal_p=''))
    H.tsv(DEST/'program_evidence_index.tsv',evidence)
    representatives=['B_AXIS']
    for radius in ['500','1000']:
        best=[r['program'] for r in ranking if r['mask_um']==radius and r['arm']=='full' and r['rank']==1]
        if best and best[0] not in representatives:representatives.append(best[0])
    M.js(DEST/'representatives.json',dict(programs=representatives,rule='D-180: B_AXIS plus each mask500/1000 full-arm worst-source-AUC rank1, illustration only',independent_validation=False))
    M.js(DEST/'complete.json',dict(status='ALL1742_INDEXED_WITH_UNCERTAINTY',n_gobp=1691,n_controls=51,programs=len(evidence),ranking_rows=len(ranking),formal_spatial_confirmation=False,figures='SEPARATE'))
    print('full program index',len(evidence),'ranking rows',len(ranking),flush=True)


def bidirectional_screen():
    z=np.load(DEST/'halo_edge_contrasts.npz');names=z['set_ids'];meta=M.tab(DEST/'all_contrasts.tsv');E=z['effect_sd'];inverse=M.tab(OUT/'reverse/summary/all_programs.tsv')
    inv={(r['mask_um'],r['split'],r['cohort'],r['arm'],r['program']):r for r in inverse};output=[]
    for radius,a_hi,b_hi in [(130,230,400),(230,400,800)]:
        selected={cohort:[i for i,r in enumerate(meta) if r['family']=='halo' and r['kind']=='adjacent_distance' and int(r['mask_um'])==radius and float(r['a_lo_um'])==radius and float(r['a_hi_um'])==a_hi and float(r['b_lo_um'])==a_hi and float(r['b_hi_um'])==b_hi and r['cohort']==cohort] for cohort in ['USZ','GSE175540']}
        for j,name in enumerate(names):
            row=dict(program_index=j,program=str(name),panel='GO_BP' if j<1691 else 'FIXED_CONTROL',forward_core_margin_um=radius,inner_bin=f'{radius}-{a_hi}',outer_bin=f'{a_hi}-{b_hi}',formal_p='',inference='DESCRIPTIVE_ASSOCIATION_AND_INTERNAL_PREDICTION_ONLY')
            for cohort,ids in selected.items():
                effects=E[ids,j];effects=effects[np.isfinite(effects)]
                row.update({f'{cohort}_forward_n':len(effects),f'{cohort}_forward_positive':int(sum(effects>0)),f'{cohort}_forward_negative':int(sum(effects<0)),f'{cohort}_forward_median_sd':float(np.median(effects)) if len(effects) else ''})
                for mask in ['230','500','1000']:
                    for arm in ['radial','angular','full']:
                        r=inv[mask,'leave_cohort',cohort,arm,str(name)]
                        row[f'{cohort}_inverse_R{mask}_{arm}_auc']=r['mean_auc'];row[f'{cohort}_inverse_R{mask}_{arm}_n']=r['valid_sections']
            output.append(row)
    H.tsv(DEST/'all_program_bidirectional_screen.tsv',output)
    scope=[]
    for source in json.loads((OUT/'contract.json').read_text())['sources']:
        d=OUT/'forward'/source['section_id'];rows=M.tab(d/'eligibility.tsv') if (d/'eligibility.tsv').exists() else []
        for family in ['TLS_TLS','TLS_REGION','TLS_EDGE']:
            r=[x for x in rows if x['family']==family];scope.append(dict(**source,family=family,n_endpoint_pairs=len(set(x['pair_id'] for x in r)),n_geometry_attempted=len(r),n_geometry_evaluable=sum(x['status']=='ELIGIBLE' for x in r),status='DESCRIPTIVE_ONLY' if any(x['status']=='ELIGIBLE' for x in r) else 'NOT_TESTABLE_FROM_AVAILABLE_GEOMETRY'))
    H.tsv(DEST/'bridge_family_scope.tsv',scope)
    M.js(DEST/'bidirectional_screen_complete.json',dict(status='ALL1742_BOTH_FORWARD_MARGINS_ALL_INVERSE_ARMS',rows=len(output),programs=len(names),same_distance_definition=False,spatial_confirmation=False))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['halo_edge','bridge','index','screen','all']);a=ap.parse_args()
    for name,fn in [('halo_edge',halo_edge),('bridge',bridge),('index',rank_and_index),('screen',bidirectional_screen)]:
        if a.stage in [name,'all']:fn()
