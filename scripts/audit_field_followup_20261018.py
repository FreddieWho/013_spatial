#!/usr/bin/env python3
"""Recompute stored predictions and verify train-only representation provenance."""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
WORK=ROOT/'infra/field_followup_20261018'
OLD=ROOT/'infra/field_discovery_20261018'
OUT=ROOT/'infra/field_followup_review_20261018'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def table(path, rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)


def main():
    OUT.mkdir(exist_ok=True);checks=Counter()
    c=json.loads((WORK/'contract.json').read_text());freeze=json.loads((WORK/'freeze_manifest.json').read_text())
    assert freeze['code_sha256']==sha(ROOT/'scripts/field_followup_20261018.py')
    assert freeze['prior_freeze_sha256']==sha(OLD/'freeze_manifest.json')
    version=json.loads((WORK/'prepare_version.json').read_text())
    assert version['code_sha256']==freeze['code_sha256'] and version['contract_sha256']==sha(WORK/'contract.json')
    for name,entry in freeze['files'].items():
        p=WORK/name;assert p.stat().st_size==entry['bytes'] and sha(p)==entry['sha256'],name
        checks['frozen_hashes']+=1
    rows=json.loads((WORK/'input_manifest.json').read_text());byid={r['section_id']:r for r in rows}
    vocab=json.loads((WORK/'vocabulary.json').read_text());defs=json.loads((WORK/'module_definitions.json').read_text())
    panel=json.loads((ROOT/'infra/tls_multiscale_gobp_20261002/panel_controls.json').read_text())
    for target in c['targets']:
        assert set(vocab[target])==set().union(*(set(g) for g in defs[target].values()))
        assert not set(vocab[target]) & set(panel[target]);checks['matched_disjoint_gene_budgets']+=1
    stored=list(csv.DictReader((WORK/'section_metrics.tsv').open(),delimiter='\t'))
    lookup={(r['section_id'],r['target'],r['model']):r for r in stored};assert len(lookup)==len(stored)
    max_metric=0.;max_latent=0.
    with threadpool_limits(limits=2):
        for fold in c['folds']:
            folder=WORK/f'folds/fold{fold}';split=json.loads((folder/'split_receipt.json').read_text())
            assert set(split['test'])=={r['section_id'] for r in rows if r['fold']==fold}
            assert set(split['train'])==set(byid)-set(split['test'])
            train_units={byid[s]['unit_id'] for s in split['train']};test_units={byid[s]['unit_id'] for s in split['test']}
            assert not train_units & test_units;checks['outer_group_partitions']+=1
            di=dict(np.load(folder/'dictionary_indices.npz',allow_pickle=False));ii=dict(np.load(folder/'train_indices.npz',allow_pickle=False))
            assert set(di)==set(ii)==set(split['train'])
            count=Counter()
            for sid,index in di.items():
                assert sid in split['train'] and len(set(index.tolist()))==len(index)
                assert np.all(index>=0) and np.all(index<byid[sid]['n_spots']);count[byid[sid]['unit_id']]+=len(index)
            assert set(count.values())=={c['dictionary_spots_per_unit']};checks['unit_balanced_dictionary_samples']+=1
            for ti,target in enumerate(c['targets']):
                with np.load(folder/f'residual_oof_t{ti}.npz',allow_pickle=False) as residual:
                    expected=np.concatenate([np.full(len(ii[sid]),split['inner_fold_by_unit'][byid[sid]['unit_id']]) for sid in split['train']])
                    assert np.array_equal(residual['inner_fold'],expected)
                    assert np.isfinite(residual['oof']).all();checks['unit_crossfit_residual_partitions']+=1
                d=dict(np.load(folder/f'dictionary_t{ti}.npz',allow_pickle=False))
                assert np.array_equal(d['genes'],vocab[target])
                train=np.concatenate([np.load(WORK/f'genes/{sid}_t{ti}.npy',mmap_mode='r')[di[sid]] for sid in split['train']])
                np.testing.assert_allclose(d['mean'],train.mean(axis=0,dtype=np.float64),atol=1e-7)
                scale=np.sqrt(train.var(axis=0,dtype=np.float64));scale[scale==0]=1
                np.testing.assert_allclose(d['scale'],scale,atol=1e-7)
                np.testing.assert_allclose(d['components']@d['components'].T,np.eye(8),atol=3e-5)
                checks['training_only_dictionary_statistics']+=1
                for sid in split['test']:
                    with np.load(folder/f'{sid}_t{ti}.npz',allow_pickle=False) as p, np.load(OLD/f'cache/{sid}.npz',allow_pickle=False) as old:
                        assert np.array_equal(p['barcodes'],old['barcodes']) and np.array_equal(p['xy_source_pixels'],old['xy_source_pixels'])
                        assert np.array_equal(p['y'],old['y'][:,ti]);checks['aligned_prediction_packages']+=1
                        gene=np.load(WORK/f'genes/{sid}_t{ti}.npy',mmap_mode='r')
                        latent=((gene-d['mean'])/d['scale']-d['pca_mean'])@d['components'].T
                        err=float(np.max(np.abs(latent-p['latent_fields'])));max_latent=max(max_latent,err)
                        assert err<3e-4,(sid,target,err);checks['latent_field_recomputations']+=1
                        y=p['y'].astype(float)
                        for model in p.files:
                            if model in ['y','barcodes','xy_source_pixels','latent_fields']:continue
                            pred=p[model].astype(float);assert pred.shape==y.shape and np.isfinite(pred).all()
                            mse=np.mean((y-pred)**2);cmse=np.mean(((y-y.mean())-(pred-pred.mean()))**2)
                            row=lookup[sid,target,model]
                            err=max(abs(mse-float(row['mse'])),abs(cmse-float(row['centered_mse'])))
                            assert err<2e-5;max_metric=max(max_metric,float(err));checks['independent_metric_recomputations']+=1
    assert checks['independent_metric_recomputations']==len(stored)
    errors=list(csv.DictReader((WORK/'train_test_errors.tsv').open(),delimiter='\t'));group=defaultdict(list)
    for r in errors:group[r['target'],r['model'],r['split'],r['unit_id'],r['fold']].append(float(r['mse']))
    within=defaultdict(list)
    for (target,model,split,unit,fold),values in group.items():within[target,model,split,unit].append(float(np.sqrt(np.mean(values))))
    totals=defaultdict(list)
    for (target,model,split,unit),values in within.items():totals[target,model,split].append(np.mean(values))
    diagnostic=[{'target':t,'model':m,'split':s,'macro_unit_rmse':float(np.mean(v)),'n_units':len(v)} for (t,m,s),v in totals.items()]
    table(OUT/'train_test_diagnosis.tsv',diagnostic)
    result={'status':'PASS','checks':dict(checks),'max_metric_absolute_difference':max_metric,'max_latent_absolute_difference':max_latent,
            'n_sections':len(rows),'n_units':len({r['unit_id'] for r in rows}),
            'known_patient_units':len({r['unit_id'] for r in rows if r['patient_known']}),
            'human_tls_labels_read':False,'run_freeze_sha256':sha(WORK/'freeze_manifest.json'),'audit_code_sha256':sha(Path(__file__)),
            'scope':'artifact recomputation and representation provenance; no independent retraining or biological validation'}
    (OUT/'audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
