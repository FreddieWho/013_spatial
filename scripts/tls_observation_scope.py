#!/usr/bin/env python3
"""D-187: eligibility of the frozen masked-query design; metadata only."""
from pathlib import Path
import csv,json,hashlib,collections
import numpy as np
from scipy import sparse,spatial
from scipy.sparse import csgraph
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'infra/tls_multiscale_gobp_20261002';OUT=ROOT/'infra/tls_observation_scope_20261003'
IDENTITY=ROOT/'infra/bioinf-data-index/raw/user_gse175540_20261003/current_26_source_identity_confirmed.tsv'
def read(p):return list(csv.DictReader(Path(p).open(),delimiter='\t'))
def write(p,rows):
    with Path(p).open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def component_extent(xy):
    maxima=[]
    for start in range(0,len(xy),256):maxima.extend(spatial.distance.cdist(xy[start:start+256],xy).max(axis=1))
    a=np.asarray(maxima);return float(a.max()),a

def main():
    OUT.mkdir(exist_ok=True);identity=read(IDENTITY);meta={r['section_id']:r for r in read(BASE/'sections.tsv')};records=[];sections=[];querycounts=[];sourcehash={};classsupport=[]
    for info in identity:
        sid=info['section_id'];cohort=info['cohort'];unit=info['patient_id'] if cohort=='GSE175540' else 'sample:'+sid
        geometry=BASE/'geometry'/f'{sid}.npz';sourcehash[str(geometry.relative_to(ROOT))]=sha(geometry);z=np.load(geometry);xy=z['xy'];labels=z['labels'];tls=np.flatnonzero(labels=='TLS');graph=sparse.load_npz(BASE/'geometry'/f'{sid}_graph.npz');nc,cc=csgraph.connected_components(graph[tls][:,tls],directed=False)
        comp={}
        for c in range(nc):
            ids=tls[cc==c];diameter,maximum=component_extent(xy[ids]);comp[c]=(ids,diameter,maximum)
        for radius in [230,500,1000]:
            path=BASE/'reverse'/f'{sid}_{radius}_queries.tsv';sourcehash[str(path.relative_to(ROOT))]=sha(path);qs=read(path);byindex={int(q['spot_index']):q for q in qs};reason=collections.Counter(q['status'] for q in qs)
            for status,n in reason.items():querycounts.append(dict(section_id=sid,cohort=cohort,evaluation_unit=unit,mask_um=radius,query_status=status,n_queries=n))
            source_ok=meta[sid]['geometry_status']=='SOURCE_GEOMETRY_AVAILABLE';local=[]
            for c,(ids,diameter,maximum) in comp.items():
                candidates=[(i,byindex[int(spot)]) for i,spot in enumerate(ids) if int(spot) in byindex];feasible=[q for i,q in candidates if len(ids)>=2 and maximum[i]+130<=radius+1e-6];qualified=[q for _,q in candidates if int(q['label'])==1]
                if not source_ok:status='KNOWN_TLS_MISSING_IN_SOURCE'
                elif len(ids)<2:status='SINGLETON_TLS_COMPONENT'
                elif not candidates:status='NO_FIXED_GRID_CENTRE'
                elif not feasible:status='TARGET_CORE_NOT_FULLY_HIDDEN'
                elif not qualified:status='INSUFFICIENT_VISIBLE_SUPPORT'
                else:status='EVALUABLE_POSITIVE_COMPONENT'
                if source_ok:assert len(qualified)==sum(q['status']=='LABELLED' and int(q['label'])==1 for _,q in candidates)
                row=dict(section_id=sid,cohort=cohort,evaluation_unit=unit,mask_um=radius,component_id=c,n_TLS_spots=len(ids),diameter_um=diameter,necessary_diameter_ceiling_um=2*(radius-130),passes_diameter_necessary_condition=diameter<=2*(radius-130)+1e-6,min_mask_at_any_member_um=float(maximum.min()+130),min_mask_at_grid_member_um=float(min(maximum[i] for i,_ in candidates)+130) if candidates else '',n_grid_TLS_queries=len(candidates),n_core_hidden_grid_queries=len(feasible),n_evaluable_positive_queries=len(qualified),status=status)
                records.append(row);local.append(row)
            pos=sum(int(q['label'])==1 for q in qs);neg=sum(int(q['label'])==0 for q in qs);contributing=[r for r in local if r['n_evaluable_positive_queries']>0]
            assert sum(r['n_evaluable_positive_queries'] for r in local)==pos
            sections.append(dict(section_id=sid,cohort=cohort,evaluation_unit=unit,unit_type='patient' if cohort=='GSE175540' else 'tumor_sample_unknown_patient',mask_um=radius,n_spots=len(xy),n_TLS_spots=len(tls),n_explicit_non_TLS_spots=int(np.isin(labels,['NOR','TUM','INFL','LN','NO_TLS']).sum()),n_unknown_label_spots=int((~np.isin(labels,['TLS','NOR','TUM','INFL','LN','NO_TLS'])).sum()),n_TLS_components=nc,n_nonsingleton_components=sum(r['n_TLS_spots']>=2 for r in local),n_components_with_positive_queries=len(contributing),n_positive_queries=pos,n_negative_queries=neg,n_excluded_queries=len(qs)-pos-neg,geometry_allows_section_AUC=pos>0 and neg>0,n_components_core_hideable_at_grid=sum(r['n_core_hidden_grid_queries']>0 for r in local),largest_qualified_component_diameter_um=max((r['diameter_um'] for r in contributing),default=''),source_geometry_status=meta[sid]['geometry_status']))
            for label in [0,1]:
                selected=[q for q in qs if int(q['label'])==label]
                if not selected:continue
                nvisible=[];supported=[]
                for q in selected:
                    dist=np.linalg.norm(xy-xy[int(q['spot_index'])],axis=1);vis=dist>radius;b=np.searchsorted([radius,radius+250,radius+750,radius+1750,radius+3750],dist,side='right')-1
                    counts=np.bincount(b[vis],minlength=5);nvisible.append(int(vis.sum()));supported.append(int((counts>=20).sum()))
                classsupport.append(dict(section_id=sid,cohort=cohort,evaluation_unit=unit,mask_um=radius,label=label,n_queries=len(selected),median_visible_spots=float(np.median(nvisible)),min_visible_spots=min(nvisible),median_first_ring=float(np.median([int(q['n_first_ring']) for q in selected])),median_second_ring=float(np.median([int(q['n_second_ring']) for q in selected])),median_supported_radial_bins=float(np.median(supported)),min_supported_radial_bins=min(supported)))
        print('scope',sid,nc,'TLS components',flush=True)
    unitrows=[];buckets={}
    for r in sections:buckets.setdefault((r['cohort'],r['mask_um'],r['evaluation_unit']),[]).append(r)
    sumfields=['n_TLS_components','n_nonsingleton_components','n_components_with_positive_queries','n_positive_queries','n_negative_queries','n_excluded_queries','n_components_core_hideable_at_grid']
    for (cohort,radius,unit),rs in buckets.items():
        row=dict(cohort=cohort,mask_um=radius,evaluation_unit=unit,unit_type=rs[0]['unit_type'],n_sections=len(rs),n_sections_with_both_classes=sum(r['geometry_allows_section_AUC'] for r in rs),n_sections_with_positive_queries=sum(r['n_positive_queries']>0 for r in rs),**{k:sum(r[k] for r in rs) for k in sumfields},component_count_scope='section-specific annotated instances; not independent patients or cross-section deduplicated TLS')
        unitrows.append(row)
    summary=[];component_reasons=[]
    for cohort in ['USZ','GSE175540']:
        for radius in [230,500,1000]:
            ur=[r for r in unitrows if r['cohort']==cohort and r['mask_um']==radius];cr=[r for r in records if r['cohort']==cohort and r['mask_um']==radius];valid=[r for r in cr if r['n_evaluable_positive_queries']>0]
            summary.append(dict(cohort=cohort,mask_um=radius,unit_type=ur[0]['unit_type'],n_total_units=len(ur),n_units_with_evaluable_section=sum(r['n_sections_with_both_classes']>0 for r in ur),n_units_with_positive_queries=sum(r['n_positive_queries']>0 for r in ur),n_annotated_component_instances=len(cr),n_nonsingleton_instances=sum(r['n_TLS_spots']>=2 for r in cr),n_contributing_component_instances=len(valid),n_positive_queries=sum(r['n_positive_queries'] for r in ur),n_negative_queries=sum(r['n_negative_queries'] for r in ur),necessary_diameter_ceiling_um=2*(radius-130),min_evaluable_component_diameter_um=min((r['diameter_um'] for r in valid),default=''),max_evaluable_component_diameter_um=max((r['diameter_um'] for r in valid),default=''),n_instances_diameter_too_large=sum(r['diameter_um']>2*(radius-130)+1e-6 for r in cr),n_units_single_contributing_instance=sum(r['n_components_with_positive_queries']==1 for r in ur)))
            for status,n in collections.Counter(r['status'] for r in cr).items():component_reasons.append(dict(cohort=cohort,mask_um=radius,status=status,n_instances=n))
    overlap=[]
    for cohort in ['USZ','GSE175540']:
        for low,high in [(230,500),(500,1000),(230,1000)]:
            sets={rad:{(r['section_id'],r['component_id']) for r in records if r['cohort']==cohort and r['mask_um']==rad and r['n_evaluable_positive_queries']>0} for rad in [low,high]}
            overlap.append(dict(cohort=cohort,smaller_mask_um=low,larger_mask_um=high,n_common_instances=len(sets[low]&sets[high]),n_new_at_larger_mask=len(sets[high]-sets[low]),n_lost_at_larger_mask=len(sets[low]-sets[high])))
    write(OUT/'cross_mask_component_overlap.tsv',overlap)
    for r in summary:
        good={s['section_id'] for s in sections if s['cohort']==r['cohort'] and s['mask_um']==r['mask_um'] and s['geometry_allows_section_AUC']}
        chosen=[c for c in records if c['section_id'] in good and c['mask_um']==r['mask_um'] and c['n_evaluable_positive_queries']>0]
        r['n_component_instances_in_AUC_evaluable_sections']=len(chosen)
        r['n_positive_queries_in_AUC_evaluable_sections']=sum(c['n_evaluable_positive_queries'] for c in chosen)
    write(OUT/'all_component_scope.tsv',records);write(OUT/'all_section_scope.tsv',sections);write(OUT/'all_unit_scope.tsv',unitrows);write(OUT/'cohort_scope.tsv',summary);write(OUT/'component_exclusion_reasons.tsv',component_reasons);write(OUT/'query_exclusion_reasons.tsv',querycounts);write(OUT/'visible_support_by_section_class.tsv',classsupport)
    contract=dict(decision='D-187',identity_sha256=sha(IDENTITY),base_contract_sha256=sha(BASE/'contract.json'),reverse_contract_sha256=sha(BASE/'reverse/contract.json'),source_sha256=sourcehash,expression_used=False,changes_to_model_or_geometry=False)
    (OUT/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    # Geometry is an upper bound; absent training positives can still prevent a fitted model.
    model_support=[]
    old=read(ROOT/'infra/tls_patient_validation_20261003/all_program_summary.tsv')
    for r in summary:
        match=next(x for x in old if x['cohort']==r['cohort'] and int(x['mask_um'])==r['mask_um'] and x['split']=='leave_cohort' and x['arm']=='full' and x['program']=='B_AXIS')
        fitted=int(match['n_evaluable_units']);assert r['n_units_with_evaluable_section']>=fitted
        reason='MATCHED' if r['n_units_with_evaluable_section']==fitted else 'NO_GSE_TRAINING_POSITIVES_AT_R230'
        if reason!='MATCHED':assert r['cohort']=='USZ' and r['mask_um']==230 and fitted==0
        model_support.append(dict(cohort=r['cohort'],mask_um=r['mask_um'],geometry_evaluable_units=r['n_units_with_evaluable_section'],source_holdout_B_evaluable_units=fitted,reason=reason))
    write(OUT/'geometry_vs_model_support.tsv',model_support)
    (OUT/'complete.json').write_text(json.dumps(dict(status='COMPLETE_FROZEN_DESIGN_SCOPE_AUDIT',component_mask_rows=len(records),section_mask_rows=len(sections),unit_mask_rows=len(unitrows),summary=summary,query_counts_reconcile=True,B_axis_qualified_units_within_geometry_upper_bound=True,training_qualification_distinguished=True,independent_structure_claim=False,scientific_confirmation=False),indent=2)+'\n')
if __name__=='__main__':main()
