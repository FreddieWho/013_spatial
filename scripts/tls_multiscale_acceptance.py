#!/usr/bin/env python3
"""Execution/evidence audit; never equates completed computation with confirmation."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import csv,hashlib,json,shutil
import numpy as np

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'infra/tls_multiscale_gobp_20261002'
def tab(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8<<20),b''):h.update(block)
    return h.hexdigest()
def js(path):return json.loads(Path(path).read_text())


def main():
    contract=js(OUT/'contract.json');sources=contract['sources'];names=list(js(OUT/'panel_gobp.json'))+list(js(OUT/'panel_controls.json'))
    assert len(sources)==26 and len(names)==1742 and len(set(names))==1742
    assert len(js(OUT/'panel_gobp.json'))==1691
    inputs=tab(OUT/'source_hashes.tsv')
    def verify_input(row):
        path=ROOT/row['source'];assert path.stat().st_size==int(row['size']);assert sha(path)==row['sha256'];return row['source']
    with ThreadPoolExecutor(max_workers=4) as pool:verified_inputs=list(pool.map(verify_input,inputs))
    records=[]
    for source in sources:
        sid=source['section_id'];receipt=js(OUT/'receipts'/f'{sid}_scores.json');path=OUT/'scores'/f'{sid}.npz';assert sha(path)==receipt['sha256']
        assert receipt['contract_sha256']==sha(OUT/'contract.json')
        z=np.load(path);np.testing.assert_array_equal(z['set_ids'],names);S=z['scores'];valid=z['valid']
        assert S.shape==(receipt['n_spots'],1742) and np.isfinite(S[:,valid]).all()
        assert np.all(z['coverage'][valid]>=.8) and np.all(z['n_genes'][valid]>=10)
        assert np.min(S[:,valid])>=-1e-6 and np.max(S[:,valid])<=1+1e-6
        assert len(np.load(OUT/'geometry'/f'{sid}.npz')['xy'])==len(S)
        op=js(OUT/'forward'/sid/'operator_complete.json');forward=js(OUT/'forward'/sid/'scoring_complete.json')
        assert forward['go_terms']==1691 and forward['controls']==51
        if op['evaluable_geometries']:
            for file in ['raw_effect','endpoint_residual','positive_segments','negative_segments']:
                assert np.load(OUT/'forward'/sid/f'{file}.npy',mmap_mode='r').shape==(op['evaluable_geometries'],1742)
        ctl=js(OUT/'controls'/sid/'geometry_complete.json');assert (OUT/'controls'/sid/'scoring_complete.json').exists();assert (OUT/'controls'/sid/'summary_complete.json').exists()
        if ctl['matched_references']:
            assert np.load(OUT/'controls'/sid/'adjusted.npy',mmap_mode='r').shape==(ctl['matched_references'],1742)
        for family in ['halo','edge']:
            directory=OUT/'halo_edge'/sid/family
            for marker in ['operators_complete','scoring_complete','calibration_complete']:assert (directory/f'{marker}.json').exists()
            assert len(tab(directory/'all_program_diagnostics.tsv'))==1742
        assert (OUT/'bridge_calibration'/sid/'calibration_complete.json').exists();assert (OUT/'bridge_calibration'/sid/'real_complete.json').exists()
        records.append(dict(section_id=sid,score_sha256=receipt['sha256'],n_covered_variable_GO=int(valid[:1691].sum()),geometry_attempts=op['attempted_geometries'],evaluable_geometry=op['evaluable_geometries'],native_references=ctl['matched_references']))
        print('audited source',sid,flush=True)
    assert sum(r['geometry_attempts'] for r in records)==80768
    assert sum(r['evaluable_geometry'] for r in records)==41058
    assert sum(r['native_references'] for r in records)==404049
    reverse=js(OUT/'reverse/model_complete.json');assert reverse['folds']==84
    for job in reverse['jobs']:
        assert not set(job['train_sections'])&set(job['test_sections'])
        assert len(set(job['train_sections'])|set(job['test_sections']))==26
        target=OUT/'reverse/models'/f"{job['mask_um']}_{job['split']}_{job['job']}.npz"
        assert sha(target)==job['sha256']
    assert js(OUT/'reverse/summary/receipt.json')['metric_rows']==815568
    assert len(tab(OUT/'integrated/program_evidence_index.tsv'))==1742
    assert len(tab(OUT/'integrated/all_program_bidirectional_screen.tsv'))==3484
    assert len(tab(OUT/'integrated/all_program_exploratory_ranking.tsv'))==15678
    for filename in ['halo_edge_bins','halo_edge_contrasts','bridge_all_programs']:
        np.testing.assert_array_equal(np.load(OUT/'integrated'/f'{filename}.npz')['set_ids'],names)
    assert np.load(OUT/'integrated/bridge_all_programs.npz')['raw_median'].shape==(1559,1742)
    for filename in ['integrated/halo_edge_all_program_diagnostics.tsv','bridge_calibration/all_source_program_diagnostics.tsv','integrated/all_program_bidirectional_screen.tsv']:
        assert all(not r['formal_p'] for r in tab(OUT/filename))
    assert js(OUT/'figures/forward_atlas/pdf_content_audit.json')['status']=='PASS_ALL1742_PAGES'
    assert js(OUT/'figures/forward_atlas/complete.json')['pages']==1742
    assert len(js(OUT/'figures/bridge_overviews_complete.json'))==6
    maps=js(OUT/'figures/representative_maps_complete.json');assert len(maps)==6 and sum(len(m['sections']) for m in maps)==78
    assert (OUT/'tests_halo_edge.txt').read_text().find('14 passed')>=0
    free=shutil.disk_usage(ROOT).free;assert free>=1_200_000_000_000
    audit=dict(status='FULL_LOCAL_EXECUTION_AUDITED_CONFIRMATION_NOT_PROVEN',verified_input_files=len(verified_inputs),sources=records,
        gobp_attempted=1691,fixed_controls=51,reverse_folds=84,atlas_pages=1742,bridge_groups=1559,
        raw_inputs_and_score_and_prediction_receipts_sha256_verified=True,failed_tests_preserved=True,
        scientific_confirmation_proven=False,goal_confirmation_requirement='NOT_ESTABLISHED',
        why_not_confirmed=['I-031: bridge null gate failed and halo/edge power inadequate; joint real-GO spatial applicability not proven.',
                           'I-033: no current source-to-patient mapping to establish the required independent repeat units.',
                           'I-034/I-035: short/narrow support and complete tumor-boundary truth unavailable for parts of the requested scope.'],
        free_bytes=free,external_data_acquired=False,GPU_used=False)
    (OUT/'acceptance_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2)+'\n')
    print('FULL_LOCAL_EXECUTION_AUDITED; scientific confirmation NOT PROVEN',flush=True)


if __name__=='__main__':main()
