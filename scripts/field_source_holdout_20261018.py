#!/usr/bin/env python3
"""D-259: fixed reduced-input and metadata-study-heldout comparisons."""
import argparse
import csv
import json
import multiprocessing
import pickle
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import h5py
import numpy as np
from threadpoolctl import threadpool_limits
import field_discovery_20261018 as base
import field_followup_20261018 as follow

WORK = base.ROOT/'infra/field_source_holdout_20261018'
PREV = follow.WORK


def splits(rows):
    result = []
    for fold in range(5):
        result.append(dict(name=f'internal_{fold}', arm='internal',
            train=[r['section_id'] for r in rows if r['fold'] != fold],
            test=[r['section_id'] for r in rows if r['fold'] == fold]))
    known = [r for r in rows if r['source_key'].startswith('http')]
    for i, source in enumerate(sorted({r['source_key'] for r in known})):
        train = [r for r in known if r['source_key'] != source]
        test = [r for r in known if r['source_key'] == source]
        result.append(dict(name=f'source_{i}', arm='source', held_source=source,
            train=[r['section_id'] for r in train], test=[r['section_id'] for r in test],
            same_organ_train_sections=sum(r['organ']==test[0]['organ'] for r in train)))
    byid = {r['section_id']:r for r in rows}
    for s in result:
        assert not {byid[x]['unit_id'] for x in s['train']} & {byid[x]['unit_id'] for x in s['test']}
    return result


def setup():
    c = json.loads((WORK/'contract.json').read_text())
    rows = json.loads((PREV/'input_manifest.json').read_text())
    allowed = {Path(__file__),Path(base.__file__),Path(follow.__file__)}
    for root in [base.WORK, PREV]:
        freeze = json.loads((root/'freeze_manifest.json').read_text())
        allowed.add(root/'freeze_manifest.json')
        allowed.update(root/p for p in freeze['files'])
    allowed.update(base.HEST/f"st/{r['section_id']}.h5ad" for r in rows)
    guard = base.InputGuard(base.ROOT,allowed,WORK); guard.install()
    base.reserve(c)
    return c, rows, guard


def prepare():
    c, rows, guard = setup()
    if (WORK/'PREPARED').exists(): raise RuntimeError('Preserve prepared version')
    full = json.loads((PREV/'vocabulary.json').read_text())
    olddefs = json.loads((base.WORK/'visible_definitions.json').read_text())
    definitions = {}; vocab = {}
    for target in c['targets']:
        definitions[target] = {}; vocab[target] = {}
        for view in c['views']:
            definitions[target][view] = {m:sorted(set(genes)&set(full[target])) for m,genes in olddefs[target][view].items()}
            vocab[target][view] = sorted(set().union(*map(set,definitions[target][view].values())))
            assert min(map(len,definitions[target][view].values())) >= 10
        assert vocab[target]['full15'] == full[target]
    base.json_write(WORK/'vocabulary.json',vocab)
    base.json_write(WORK/'definitions.json',definitions)
    base.json_write(WORK/'input_manifest.json',rows)
    base.json_write(WORK/'splits.json',splits(rows))
    base.tsv_write(WORK/'source_eligibility.tsv', [dict(section_id=r['section_id'],unit_id=r['unit_id'],organ=r['organ'],source_key=r['source_key'],status='ELIGIBLE' if r['source_key'].startswith('http') else 'NOT_TESTABLE_SOURCE_IDENTITY') for r in rows])
    (WORK/'genes').mkdir(exist_ok=True)
    frozen = json.loads((PREV/'freeze_manifest.json').read_text())
    for r in rows:
        base.reserve(c); sid = r['section_id']; path = base.HEST/f'st/{sid}.h5ad'
        guard.record(path,['X','var/_index','obs/_index','obs/in_tissue'],'reduced visible gene input')
        assert base.sha256(path) == r['expression_sha256']
        with h5py.File(guard.check(path)) as f:
            genes = base.strings(f['var/_index'][:]); barcodes=base.strings(f['obs/_index'][:])
            keep=np.flatnonzero(f['obs/in_tissue'][:]==1)
            x=base.read_matrix(f['X'],len(barcodes),len(genes))[keep]
        cache=base.load_cache(r)
        assert np.array_equal(barcodes[keep],cache['barcodes'])
        for ti,target in enumerate(c['targets']):
            z,qc=follow.gene_input(x,genes,vocab[target]['small5'])
            np.save(WORK/f'genes/{sid}_t{ti}.npy',z)
            np.save(WORK/f'genes/{sid}_t{ti}_qc.npy',qc)
            for suffix in ['', '_qc']:
                rel=f'genes/{sid}_t{ti}{suffix}.npy'
                assert base.sha256(PREV/rel)==frozen['files'][rel]['sha256']
        print('prepared '+sid,flush=True)
    base.json_write(WORK/'input_receipt_prepare.json',guard.receipt())
    base.json_write(WORK/'prepare_version.json',dict(code_sha256=base.sha256(__file__),contract_sha256=base.sha256(WORK/'contract.json'),prior_freeze_sha256=base.sha256(PREV/'freeze_manifest.json')))
    (WORK/'PREPARED').write_text('READY\n')


def run_split(s):
    c, rows, guard=setup(); start=time.time()
    out=WORK/s['name']; out.mkdir(exist_ok=True)
    if (out/'COMPLETE').exists(): return s['name']+' already complete'
    byid={r['section_id']:r for r in rows}; train=[byid[x] for x in s['train']]; test=[byid[x] for x in s['test']]
    active=train+test; cache={r['section_id']:base.load_cache(r) for r in active}
    idx={r['section_id']:base.sample_indices(r['n_spots'],r['section_id'],c['seed'],512) for r in train}
    w=follow.weights(train,idx)
    base.json_write(out/'split.json',s|dict(code_sha256=base.sha256(__file__),contract_sha256=base.sha256(WORK/'contract.json')))
    np.savez_compressed(out/'train_indices.npz',**idx)
    vocab=json.loads((WORK/'vocabulary.json').read_text()); defs=json.loads((WORK/'definitions.json').read_text()); metrics=[]
    with threadpool_limits(limits=c['threads']):
        for ti,target in enumerate(c['targets']):
            ys={sid:a['y'][:,ti] for sid,a in cache.items()}; y=follow.stack(train,ys,idx)
            for view in c['views']:
                base.reserve(c); root=WORK if view=='small5' else PREV
                raw={}; point={}; local={}; gidx={g:i for i,g in enumerate(vocab[target][view])}
                for r in active:
                    sid=r['section_id']; z=np.load(root/f'genes/{sid}_t{ti}.npy',mmap_mode='r'); qc=np.load(root/f'genes/{sid}_t{ti}_qc.npy',mmap_mode='r')
                    scores=np.column_stack([z[:,[gidx[g] for g in members]].mean(axis=1) for members in defs[target][view].values()])
                    point[sid]=np.c_[scores,qc,follow.context(r,len(z))].astype(np.float32)
                    local[sid]=follow.local_features(point[sid],cache[sid])
                    raw[sid]=np.c_[z,qc,follow.context(r,len(z))].astype(np.float32)
                designs={'gene_ridge':raw,'fixed_point_hgb':point,'fixed_local_hgb':local}
                models={name:(follow.ridge if name=='gene_ridge' else follow.hgb)(c,follow.stack(train,xx,idx),y,w) for name,xx in designs.items()}
                mean=float(np.average(y,weights=w))
                with (out/f'models_t{ti}_{view}.pkl').open('wb') as f: pickle.dump(dict(models=models,train_mean=mean,train_sections=s['train']),f,protocol=5)
                for r in test:
                    sid=r['section_id']; pred={name:model.predict(designs[name][sid]).astype(np.float32) for name,model in models.items()}
                    pred['gene_ridge_smooth_half']=.5*(pred['gene_ridge']+base.neighbor_mean(cache[sid]['near'],pred['gene_ridge']))
                    pred['train_mean']=np.full(len(ys[sid]),mean,dtype=np.float32)
                    for name,p in pred.items():
                        assert np.isfinite(p).all()
                        metrics.append({k:r[k] for k in ['section_id','unit_id','organ','source_key','patient_known']}|dict(arm=s['arm'],split=s['name'],target=target,view=view,model=name,**base.metrics(ys[sid],p)))
                    np.savez_compressed(out/f'{sid}_t{ti}_{view}.npz',y=ys[sid],barcodes=cache[sid]['barcodes'],**pred)
                print(s['name']+' '+target+' '+view,flush=True)
    base.tsv_write(out/'metrics.tsv',metrics);base.json_write(out/'input_receipt.json',guard.receipt())
    (out/'COMPLETE').write_text('COMPLETE\n')
    return dict(split=s['name'],seconds=time.time()-start)


def summarize():
    c,rows,guard=setup(); records=[]
    for s in json.loads((WORK/'splits.json').read_text()):
        assert (WORK/s['name']/'COMPLETE').exists()
        records.extend(csv.DictReader((WORK/s['name']/'metrics.tsv').open(),delimiter='\t'))
    groups=defaultdict(list)
    for r in records: groups[tuple(r[k] for k in ['arm','target','view','model','unit_id'])].append(r)
    units=[]
    for key,rr in groups.items():
        units.append(dict(zip(['arm','target','view','model','unit_id'],key))|dict(source_key=rr[0]['source_key'],organ=rr[0]['organ'],split=rr[0]['split'],rmse=float(np.sqrt(np.mean([float(r['mse']) for r in rr]))),centered_rmse=float(np.sqrt(np.mean([float(r['centered_mse']) for r in rr])))))
    summary=[]; bysource=[]; decisions=[]
    for arm in ['internal','source']:
        for target in c['targets']:
            for view in c['views']:
                group=[r for r in units if r['arm']==arm and r['target']==target and r['view']==view]
                for model in c['models']:
                    mm=[r for r in group if r['model']==model]
                    summary.append(dict(arm=arm,target=target,view=view,model=model,n_units=len(mm),rmse=float(np.mean([r['rmse'] for r in mm])),centered_rmse=float(np.mean([r['centered_rmse'] for r in mm]))))
                    if arm=='source':
                        for src in sorted({r['source_key'] for r in mm}):
                            ss=[r for r in mm if r['source_key']==src]
                            bysource.append(dict(target=target,view=view,model=model,source=src,organ=ss[0]['organ'],n_units=len(ss),rmse=float(np.mean([r['rmse'] for r in ss])),centered_rmse=float(np.mean([r['centered_rmse'] for r in ss]))))
                loc={r['unit_id']:r for r in group if r['model']=='fixed_local_hgb'}
                for model in c['models']:
                    if model=='fixed_local_hgb':continue
                    bb={r['unit_id']:r for r in group if r['model']==model}
                    gain=1-np.mean([r['rmse'] for r in loc.values()])/np.mean([r['rmse'] for r in bb.values()])
                    cg=1-np.mean([r['centered_rmse'] for r in loc.values()])/np.mean([r['centered_rmse'] for r in bb.values()])
                    frac=np.mean([loc[u]['rmse']<bb[u]['rmse'] for u in loc])
                    decisions.append(dict(arm=arm,target=target,view=view,baseline=model,gain=float(gain),centered_gain=float(cg),positive_fraction=float(frac),passes=bool(gain>=.03 and cg>0 and frac>=.6)))
    base.tsv_write(WORK/'section_metrics.tsv',records);base.tsv_write(WORK/'unit_metrics.tsv',units)
    base.tsv_write(WORK/'model_comparison.tsv',summary);base.tsv_write(WORK/'source_comparison.tsv',bysource);base.tsv_write(WORK/'candidate_comparisons.tsv',decisions)
    base.json_write(WORK/'summary.json',dict(n_sections=len(rows),n_units=len({r['unit_id'] for r in rows}),metrics=len(records),not_run=c['not_run']))
    (WORK/'STATUS').write_text('COMPLETE_PENDING_AUDIT\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','run','summarize']);a=p.parse_args()
    if a.command=='prepare':prepare()
    elif a.command=='summarize':summarize()
    else:
        c=json.loads((WORK/'contract.json').read_text()); assert (WORK/'PREPARED').exists()
        with ProcessPoolExecutor(max_workers=c['workers'],mp_context=multiprocessing.get_context('spawn')) as pool:
            futures=[pool.submit(run_split,s) for s in json.loads((WORK/'splits.json').read_text())]
            for f in as_completed(futures):print(f.result(),flush=True)
        summarize()
