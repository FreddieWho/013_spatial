#!/usr/bin/env python3
"""Independent stored-output, split and same-section transfer audit for D-259."""
import csv
import hashlib
import json
import shutil
import pickle
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'infra/field_source_holdout_20261018'
OUT=ROOT/'infra/field_source_holdout_review_20261018'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def write(name,rr):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,list(rr[0]),delimiter='\t');w.writeheader();w.writerows(rr)


def main():
    OUT.mkdir(exist_ok=True); checks=Counter()
    c=json.loads((WORK/'contract.json').read_text());rows=json.loads((WORK/'input_manifest.json').read_text());byid={r['section_id']:r for r in rows}
    version=json.loads((WORK/'prepare_version.json').read_text())
    assert version['code_sha256']==sha(ROOT/'scripts/field_source_holdout_20261018.py')
    assert version['contract_sha256']==sha(WORK/'contract.json')
    # Check the frozen inputs actually reused by this run.
    for folder in ['field_discovery_20261018','field_followup_20261018']:
        root=ROOT/'infra'/folder; frozen=json.loads((root/'freeze_manifest.json').read_text())
        for name,entry in frozen['files'].items():
            needed=(name.startswith('cache/') if folder=='field_discovery_20261018' else name.startswith('genes/')) or name in ['input_manifest.json','visible_definitions.json','vocabulary.json','module_definitions.json']
            if needed:
                assert sha(root/name)==entry['sha256'],name
                checks['reused_frozen_inputs']+=1
    vocab=json.loads((WORK/'vocabulary.json').read_text());defs=json.loads((WORK/'definitions.json').read_text())
    panels=json.loads((ROOT/'infra/tls_multiscale_gobp_20261002/panel_controls.json').read_text())
    for t in c['targets']:
        for view in c['views']:
            assert set(vocab[t][view])==set().union(*map(set,defs[t][view].values()))
            assert not set(vocab[t][view])&set(panels[t]);checks['disjoint_visible_budgets']+=1
    records=list(csv.DictReader((WORK/'section_metrics.tsv').open(),delimiter='\t'))
    lookup={(r['arm'],r['section_id'],r['target'],r['view'],r['model']):r for r in records};assert len(lookup)==len(records)==5000
    coverage=defaultdict(list); max_error=0.; old_error=0.; old_rmse_error=0.
    for s in json.loads((WORK/'splits.json').read_text()):
        out=WORK/s['name']; saved=json.loads((out/'split.json').read_text())
        assert saved['code_sha256']==version['code_sha256'] and saved['contract_sha256']==version['contract_sha256']
        assert saved['train']==s['train'] and saved['test']==s['test']
        tr={byid[x]['unit_id'] for x in s['train']};te={byid[x]['unit_id'] for x in s['test']};assert not tr&te
        if s['arm']=='source':
            assert all(byid[x]['source_key'].startswith('http') for x in s['train']+s['test'])
            assert {byid[x]['source_key'] for x in s['test']}=={s['held_source']}
            assert s['held_source'] not in {byid[x]['source_key'] for x in s['train']}
        checks['group_partitions']+=1;coverage[s['arm']]+=s['test']
        indices=np.load(out/'train_indices.npz');assert set(indices.files)==set(s['train'])
        if s['arm']=='internal':
            prior_dir=ROOT/f"infra/field_followup_20261018/folds/fold{s['name'].split('_')[1]}"
            prior_indices=np.load(prior_dir/'train_indices.npz')
            for sid in indices.files:np.testing.assert_array_equal(indices[sid],prior_indices[sid])
            for ti in range(4):
                with (out/f'models_t{ti}_full15.pkl').open('rb') as f:new_model=pickle.load(f)['models']['gene_ridge']
                with (prior_dir/f'models_t{ti}.pkl').open('rb') as f:old_model=pickle.load(f)['models']['gene_ridge']
                np.testing.assert_array_equal(new_model[0].mean_,old_model[0].mean_)
                np.testing.assert_array_equal(new_model[0].scale_,old_model[0].scale_)
                assert new_model[1].alpha==old_model[1].alpha
                checks['prior_train_normalizers_identical']+=1
        for sid in s['test']:
            old=np.load(ROOT/f'infra/field_discovery_20261018/cache/{sid}.npz')
            for ti,t in enumerate(c['targets']):
                for view in c['views']:
                    z=np.load(out/f'{sid}_t{ti}_{view}.npz'); y=z['y'].astype(float)
                    np.testing.assert_array_equal(z['barcodes'],old['barcodes']);np.testing.assert_array_equal(z['y'],old['y'][:,ti])
                    assert set(z.files)=={'y','barcodes'}|set(c['models']);checks['prediction_packages']+=1
                    for model in c['models']:
                        p=z[model].astype(float);assert np.isfinite(p).all()
                        rr=lookup[s['arm'],sid,t,view,model]
                        mse=np.mean((y-p)**2);cmse=np.mean(((y-y.mean())-(p-p.mean()))**2)
                        err=max(abs(mse-float(rr['mse'])),abs(cmse-float(rr['centered_mse'])))
                        assert err<2e-6;max_error=max(max_error,err);checks['recomputed_metrics']+=1
                    if s['arm']=='internal' and view=='full15':
                        prior=np.load(ROOT/f"infra/field_followup_20261018/folds/fold{byid[sid]['fold']}/{sid}_t{ti}.npz")
                        e=float(np.max(np.abs(z['gene_ridge']-prior['gene_ridge'])))
                        old_error=max(old_error,e)
                        old_rmse_error=max(old_rmse_error,abs(float(np.sqrt(np.mean((z['y'].astype(float)-z['gene_ridge'])**2)))-float(np.sqrt(np.mean((prior['y'].astype(float)-prior['gene_ridge'])**2)))))
                        checks['prior_ridge_differences_measured']+=1
    for arm,n in [('internal',72),('source',53)]:assert len(coverage[arm])==len(set(coverage[arm]))==n
    # Transfer comparisons pair exactly the same units/sections across the two arms.
    uu=list(csv.DictReader((WORK/'unit_metrics.tsv').open(),delimiter='\t'))
    ul={(r['arm'],r['target'],r['view'],r['model'],r['unit_id']):r for r in uu}
    transfer=[]
    for t in c['targets']:
        for view in c['views']:
            for m in c['models']:
                source=[r for r in uu if r['arm']=='source' and r['target']==t and r['view']==view and r['model']==m]
                for src in sorted({r['source_key'] for r in source}):
                    ss=[r for r in source if r['source_key']==src];ii=[ul['internal',t,view,m,r['unit_id']] for r in ss]
                    si=float(np.mean([float(r['rmse']) for r in ss]));ins=float(np.mean([float(r['rmse']) for r in ii]))
                    transfer.append(dict(target=t,view=view,model=m,source=src,n_units=len(ss),internal_rmse=ins,source_rmse=si,relative_error_change=si/ins-1))
    write('paired_transfer.tsv',transfer)
    # Source-macro summaries weight each source once, independent of cohort size.
    macro=[]
    for t in c['targets']:
        for view in c['views']:
            for m in c['models']:
                rr=[r for r in transfer if r['target']==t and r['view']==view and r['model']==m]
                macro.append(dict(target=t,view=view,model=m,n_sources=len(rr),internal_rmse=float(np.mean([r['internal_rmse'] for r in rr])),source_rmse=float(np.mean([r['source_rmse'] for r in rr]))))
    write('source_macro.tsv',macro)
    for name in ['prepare']+[s['name'] for s in json.loads((WORK/'splits.json').read_text())]:
        receipt=json.loads((WORK/('input_receipt_prepare.json' if name=='prepare' else f'{name}/input_receipt.json')).read_text())
        assert receipt['human_tls_files_opened']==[]
        for path in receipt['python_read_events']:
            assert '/annotations/' not in path and '/incoming/' not in path
        checks['guard_receipts']+=1
    free=shutil.disk_usage(WORK).free;assert free>=c['min_free_bytes']
    result=dict(status='PASS',checks=dict(checks),max_metric_difference=max_error,max_prior_ridge_difference=old_error,max_prior_section_rmse_difference=old_rmse_error,prior_bitwise_reproduction=False,free_bytes=free,limitations=c['not_run'],audit_scope='stored prediction recomputation, source/unit partitions, target exclusion, identical prior training indices/scalers and quantified Ridge differences (not exact reproduction); not a full independent training rerun')
    (OUT/'audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
