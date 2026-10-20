#!/usr/bin/env python3
"""Bounded background/residual diagnosis and an equal-gene-budget PCA probe."""
import argparse
import csv
import hashlib
import json
import multiprocessing
import pickle
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

import field_discovery_20261018 as base

ROOT = base.ROOT
OLD = base.WORK
WORK = ROOT / 'infra/field_followup_20261018'


def common_vocab(visible, gene_sets, hidden):
    return sorted((set(visible) - set(hidden)).intersection(*gene_sets))


def gene_input(x, genes, vocabulary):
    counts_by_name = Counter(genes)
    if any(counts_by_name[g] != 1 for g in vocabulary):
        raise ValueError('Visible gene identity must occur exactly once')
    lookup = {g: i for i, g in enumerate(genes)}
    columns = [lookup[g] for g in vocabulary]
    counts = x[:, columns].tocsr()
    library = np.asarray(counts.sum(axis=1)).ravel()
    z = base.normalized_subset(x, columns, library).toarray().astype(np.float32)
    qc = np.c_[np.log1p(library), np.asarray((counts > 0).sum(axis=1)).ravel()/len(columns)].astype(np.float32)
    return z, qc


def inner_fold(unit):
    return int.from_bytes(hashlib.sha256(('residual|' + unit).encode()).digest()[:8], 'little') % 3


def dictionary_indices(rows, limit, seed):
    counts = Counter(r['unit_id'] for r in rows)
    return {r['section_id']: base.sample_indices(r['n_spots'], r['section_id'], seed,
             limit // counts[r['unit_id']]) for r in rows}


def setup():
    c = json.loads((WORK/'contract.json').read_text())
    frozen = json.loads((OLD/'freeze_manifest.json').read_text())
    assert base.sha256(Path(base.__file__)) == frozen['code_sha256']
    rows = [r for r in json.loads((OLD/'input_manifest.json').read_text()) if r['status']=='QUALIFIED']
    identity = {r['section_id']: r for r in csv.DictReader((ROOT/c['identity_source']).open(), delimiter='\t')}
    for r in rows:
        r['patient_known'] = identity[r['section_id']]['patient_confirmed_within_study'] == 'True'
    base.assert_grouped(rows)
    oldc = json.loads((OLD/'contract.json').read_text())
    allowed = {OLD/p for p in frozen['files']} | {OLD/'freeze_manifest.json', Path(base.__file__), Path(__file__),
              ROOT/c['identity_source'], ROOT/oldc['panel_file'], WORK/'contract.json'}
    allowed.update(base.HEST/f"st/{r['section_id']}.h5ad" for r in rows)
    guard = base.InputGuard(ROOT, allowed, WORK); guard.install()
    for p in [OLD/'freeze_manifest.json', ROOT/c['identity_source'], ROOT/oldc['panel_file'], Path(__file__), WORK/'contract.json']:
        guard.record(p, ['frozen definition'], 'input provenance')
    base.reserve(c)
    return c, rows, oldc, guard


def prepare():
    c, rows, oldc, guard = setup()
    if (WORK/'PREPARED').exists():
        raise RuntimeError('Prepared output exists; preserve this version')
    (WORK/'genes').mkdir(exist_ok=True)
    definitions = json.loads((OLD/'visible_definitions.json').read_text())
    panels = json.loads((ROOT/oldc['panel_file']).read_text())
    source_genes = {}; gene_sets = []; duplicates = {}
    for r in rows:
        sid = r['section_id']; p = guard.check(base.HEST/f'st/{sid}.h5ad')
        with h5py.File(p) as f:
            gg = base.strings(f['var/_index'][:])
        duplicates[sid] = {g:n for g,n in Counter(gg).items() if n > 1}
        source_genes[sid] = gg; gene_sets.append(set(gg))
    vocab = {}; modules = {}
    for ti, target in enumerate(c['targets']):
        visible = set().union(*(set(v) for v in definitions[target]['full15'].values()))
        vocab[target] = common_vocab(visible, gene_sets, panels[target])
        modules[target] = {name: [g for g in members if g in set(vocab[target])]
                           for name, members in definitions[target]['full15'].items()}
        if min(map(len, modules[target].values())) < 10:
            raise RuntimeError(f'Common module coverage too small: {target}')
        assert not set(vocab[target]) & set(panels[target])
        for sid, repeated in duplicates.items():
            if set(repeated) & (set(vocab[target]) | set(panels[target])):
                raise ValueError(f'Ambiguous measured/visible gene identity: {sid}/{target}')
    base.json_write(WORK/'source_duplicate_symbols.json', duplicates)
    base.json_write(WORK/'vocabulary.json', vocab)
    base.json_write(WORK/'module_definitions.json', modules)
    base.json_write(WORK/'input_manifest.json', rows)
    old_inputs = json.loads((OLD/'input_receipt_prepare.json').read_text())['inputs']
    for r in rows:
        base.reserve(c); sid = r['section_id']; p = base.HEST/f'st/{sid}.h5ad'
        guard.record(p, ['X', 'var/_index', 'obs/_index', 'obs/in_tissue'], 'visible gene matrix')
        assert guard.records[str(p.relative_to(ROOT))]['sha256'] == old_inputs[str(p.relative_to(ROOT))]['sha256']
        with h5py.File(p) as f:
            barcodes = base.strings(f['obs/_index'][:]); keep = np.flatnonzero(f['obs/in_tissue'][:] == 1)
            x = base.read_matrix(f['X'], len(barcodes), len(source_genes[sid]))[keep]
        with np.load(OLD/f'cache/{sid}.npz', allow_pickle=False) as z:
            assert np.array_equal(barcodes[keep], z['barcodes'])
        for ti, target in enumerate(c['targets']):
            z, qc = gene_input(x, source_genes[sid], vocab[target])
            np.save(WORK/f'genes/{sid}_t{ti}.npy', z)
            np.save(WORK/f'genes/{sid}_t{ti}_qc.npy', qc)
        print(json.dumps({'prepared':sid}), flush=True)
    base.json_write(WORK/'input_receipt.json', guard.receipt())
    base.json_write(WORK/'prepare_version.json', {'code_sha256':base.sha256(__file__), 'contract_sha256':base.sha256(WORK/'contract.json')})
    (WORK/'PREPARED').write_text('ALL_72_INPUTS_READY\n')


def context(row, n):
    vector = [float(row['organ']==o) for o in base.ORGANS] + [float(row['metadata']['st_technology']==p) for p in base.PLATFORMS]
    return np.broadcast_to(np.array(vector, dtype=np.float32), (n, len(vector)))


def local_features(point, cache):
    data = point[:, :-6]
    return np.c_[point, base.neighbor_mean(cache['near'], data), base.neighbor_mean(cache['outer'], data),
                 np.log1p(np.diff(cache['near'].indptr)), np.log1p(np.diff(cache['outer'].indptr))].astype(np.float32)


def stack(rows, values, indices):
    return np.concatenate([values[r['section_id']][indices[r['section_id']]] for r in rows])


def weights(rows, indices):
    counts = Counter(r['unit_id'] for r in rows)
    w = np.concatenate([np.full(len(indices[r['section_id']]), 1/(counts[r['unit_id']]*len(indices[r['section_id']]))) for r in rows])
    return w/w.mean()


def hgb(c, x, y, w):
    model = HistGradientBoostingRegressor(**c['hgb'])
    model.fit(x, y, sample_weight=w)
    return model


def ridge(c, x, y, w):
    model = make_pipeline(StandardScaler(), Ridge(alpha=c['ridge_alpha']))
    model.fit(x, y, standardscaler__sample_weight=w, ridge__sample_weight=w)
    return model


def record_metrics(row, target, model, y, pred):
    return {k:row[k] for k in ['section_id','unit_id','fold','organ','patient_known','source_key']} | {
        'target':target, 'model':model, **base.metrics(y, pred)}


def run_fold(fold):
    c, rows, oldc, guard = setup(); start = time.time()
    assert (WORK/'PREPARED').exists()
    out = WORK/f'folds/fold{fold}'; out.mkdir(parents=True, exist_ok=True)
    if (out/'COMPLETE').exists():
        return {'fold':fold,'status':'already_complete'}
    train = [r for r in rows if r['fold'] != fold]; test = [r for r in rows if r['fold']==fold]
    assert not {r['unit_id'] for r in train} & {r['unit_id'] for r in test}
    cache = {r['section_id']:base.load_cache(r) for r in rows}
    idx = {r['section_id']:base.sample_indices(r['n_spots'],r['section_id'],oldc['seed'],512) for r in train}
    di = dictionary_indices(train, c['dictionary_spots_per_unit'], c['seed']); w = weights(train, idx)
    inner = np.concatenate([np.full(len(idx[r['section_id']]),inner_fold(r['unit_id'])) for r in train])
    vocab = json.loads((WORK/'vocabulary.json').read_text()); definitions = json.loads((WORK/'module_definitions.json').read_text())
    base.json_write(out/'split_receipt.json', {'train':[r['section_id'] for r in train], 'test':[r['section_id'] for r in test],
        'inner_fold_by_unit':{r['unit_id']:inner_fold(r['unit_id']) for r in train}, 'code_sha256':base.sha256(__file__),
        'contract_sha256':base.sha256(WORK/'contract.json')})
    np.savez_compressed(out/'train_indices.npz', **idx); np.savez_compressed(out/'dictionary_indices.npz', **di)
    metrics=[]; fit_errors=[]
    with threadpool_limits(limits=c['threads']):
        for ti,target in enumerate(c['targets']):
            ys = {sid:a['y'][:,ti] for sid,a in cache.items()}; y = stack(train,ys,idx)
            original = {r['section_id']:base.feature_views(r,cache[r['section_id']],ti,'full15') for r in rows}
            xp = {sid:xx['point'] for sid,xx in original.items()}
            xl = {sid:np.c_[xx['point'],xx['space'][:,xx['bag'].shape[1]:]] for sid,xx in original.items()}
            txp, txl = stack(train,xp,idx), stack(train,xl,idx)
            with (OLD/f'folds/fold{fold}/models_t{ti}_full15.pkl').open('rb') as f:
                old_models=pickle.load(f)['models']
            for r in train:
                sid=r['section_id']; ii=idx[sid]
                for level in ['point','bag','space']:
                    pred=old_models[level+'_hgb'].predict(original[sid][level][ii])
                    fit_errors.append({**record_metrics(r,target,level+'_hgb',ys[sid][ii],pred),'split':'training_fit'})
            local = hgb(c,txl,y,w)
            oof = np.empty_like(y)
            for k in range(c['inner_folds']):
                fit = inner != k; held = ~fit
                if not fit.any() or not held.any(): raise RuntimeError('Empty inner fold')
                wi=w[fit]/w[fit].mean()
                model=hgb(c,txp[fit],y[fit],wi);oof[held]=model.predict(txp[held])
            residual_point=hgb(c,txp,y-oof,w);residual_local=hgb(c,txl,y-oof,w)
            np.savez_compressed(out/f'residual_oof_t{ti}.npz', y=y, oof=oof, inner_fold=inner)
            gene={r['section_id']:np.load(WORK/f"genes/{r['section_id']}_t{ti}.npy",mmap_mode='r') for r in rows}
            qc={r['section_id']:np.load(WORK/f"genes/{r['section_id']}_t{ti}_qc.npy",mmap_mode='r') for r in rows}
            dictionary_train=stack(train,gene,di)
            scaler=StandardScaler(); zz=scaler.fit_transform(dictionary_train)
            pca=PCA(n_components=c['dictionary_components'],svd_solver='randomized',iterated_power=4,random_state=c['seed'])
            pca.fit(zz)
            gidx={g:i for i,g in enumerate(vocab[target])}
            fixed={}; learned={}; raw={}; latent={}
            for r in rows:
                sid=r['section_id']; z=gene[sid]
                scores=np.column_stack([z[:,[gidx[g] for g in members]].mean(axis=1) for members in definitions[target].values()])
                latent[sid]=pca.transform(scaler.transform(z)).astype(np.float32)
                fixed[sid]=np.c_[scores,qc[sid],context(r,len(z))].astype(np.float32)
                learned[sid]=np.c_[latent[sid],qc[sid],context(r,len(z))].astype(np.float32)
                raw[sid]=np.c_[z,qc[sid],context(r,len(z))].astype(np.float32)
            fixed_local={sid:local_features(x,cache[sid]) for sid,x in fixed.items()}
            pca_local={sid:local_features(x,cache[sid]) for sid,x in learned.items()}
            designs={'fixed_point_hgb':fixed,'fixed_point_ridge':fixed,'fixed_local_hgb':fixed_local,
                     'pca_point_hgb':learned,'pca_point_ridge':learned,'pca_local_hgb':pca_local,'gene_ridge':raw}
            models={}
            for name,values in designs.items():
                xx=stack(train,values,idx)
                models[name]=(ridge if name.endswith('ridge') else hgb)(c,xx,y,w)
            np.savez_compressed(out/f'dictionary_t{ti}.npz',genes=np.array(vocab[target]),components=pca.components_,
                mean=scaler.mean_,scale=scaler.scale_,pca_mean=pca.mean_,explained_variance_ratio=pca.explained_variance_ratio_)
            top=[]
            for ci,loading in enumerate(pca.components_/scaler.scale_):
                order=np.argsort(-np.abs(loading))[:20]
                top.append({'component':ci,'genes':[{'gene':vocab[target][j],'weight':float(loading[j])} for j in order]})
            base.json_write(out/f'top_loadings_t{ti}.json',top)
            with (out/f'models_t{ti}.pkl').open('wb') as f:
                pickle.dump({'local':local,'residual_point':residual_point,'residual_local':residual_local,
                             'scaler':scaler,'pca':pca,'models':models,'train_sections':[r['section_id'] for r in train]},f,protocol=5)
            for r in test:
                sid=r['section_id']; yy=ys[sid]
                with np.load(OLD/f'folds/fold{fold}/{sid}_t{ti}_full15.npz',allow_pickle=False) as z:
                    pred={'old_'+name:z[name] for name in oldc['baselines']+['space_hgb']}
                point=old_models['point_hgb'].predict(xp[sid])
                np.testing.assert_allclose(point,pred['old_point_hgb'],atol=1e-6,rtol=1e-6)
                pred['diag_local_hgb']=local.predict(xl[sid])
                pred['diag_residual_point']=point+c['residual_multiplier']*residual_point.predict(xp[sid])
                pred['diag_residual_local']=point+c['residual_multiplier']*residual_local.predict(xl[sid])
                for name,model in models.items():pred[name]=model.predict(designs[name][sid])
                for name in ['fixed_point_hgb','pca_point_hgb','gene_ridge']:
                    pred[name+'_smooth_half']=(pred[name]+base.neighbor_mean(cache[sid]['near'],pred[name]))*.5
                pred={name:np.asarray(value,dtype=np.float32) for name,value in pred.items()}
                for name,value in pred.items():
                    if not np.all(np.isfinite(value)):raise ValueError(f'Nonfinite prediction: {sid}/{name}')
                    metrics.append(record_metrics(r,target,name,yy,value))
                for level in ['point','bag','space']:
                    name=level+'_hgb';fit_errors.append({**record_metrics(r,target,name,yy,pred['old_'+name]),'split':'heldout'})
                np.savez_compressed(out/f'{sid}_t{ti}.npz',y=yy,barcodes=cache[sid]['barcodes'],
                    xy_source_pixels=cache[sid]['xy_source_pixels'],latent_fields=latent[sid],**pred)
            print(json.dumps({'fold':fold,'target':target,'seconds':round(time.time()-start,1)}),flush=True)
            del original, xp, xl, txp, txl, gene, raw, fixed, learned, fixed_local, pca_local, designs, dictionary_train, zz
    base.tsv_write(out/'metrics.tsv',metrics);base.tsv_write(out/'train_test_errors.tsv',fit_errors)
    base.json_write(out/'input_receipt.json',guard.receipt());(out/'COMPLETE').write_text('ALL_FOUR_TARGETS_COMPLETE\n')
    return {'fold':fold,'seconds':round(time.time()-start,1)}


def aggregate(rows):
    grouped=defaultdict(list)
    for r in rows: grouped[(r['target'],r['unit_id'],r['model'])].append(r)
    return {key:{'rmse':float(np.sqrt(np.mean([float(r['mse']) for r in rr]))),
                 'centered_rmse':float(np.sqrt(np.mean([float(r['centered_mse']) for r in rr]))),
                 'fold':int(rr[0]['fold'])} for key,rr in grouped.items()}


def summarize():
    c=json.loads((WORK/'contract.json').read_text());oldc=json.loads((OLD/'contract.json').read_text())
    rows=[];fit_errors=[]
    for fold in c['folds']:
        out=WORK/f'folds/fold{fold}'; assert (out/'COMPLETE').exists()
        receipt=json.loads((out/'split_receipt.json').read_text())
        assert receipt['code_sha256']==base.sha256(__file__)
        assert receipt['contract_sha256']==base.sha256(WORK/'contract.json')
        rows+=list(csv.DictReader((out/'metrics.tsv').open(),delimiter='\t'))
        fit_errors+=list(csv.DictReader((out/'train_test_errors.tsv').open(),delimiter='\t'))
    base.tsv_write(WORK/'section_metrics.tsv',rows);base.tsv_write(WORK/'train_test_errors.tsv',fit_errors)
    unit=aggregate(rows); matching=[]; stability={}
    for ti,target in enumerate(c['targets']):
        load=[]
        for fold in c['folds']:
            with np.load(WORK/f'folds/fold{fold}/dictionary_t{ti}.npz',allow_pickle=False) as z:
                a=z['components']/z['scale']; a/=np.linalg.norm(a,axis=1)[:,None];load.append(a)
        similarities=[]
        for i in range(5):
            for j in range(i+1,5):
                cosine=np.abs(load[i]@load[j].T);ii,jj=linear_sum_assignment(-cosine)
                qa=np.linalg.qr(load[i].T)[0];qb=np.linalg.qr(load[j].T)[0]
                subspace=np.linalg.svd(qa.T@qb,compute_uv=False)
                similarities.extend(cosine[ii,jj].tolist())
                for a,b in zip(ii,jj):matching.append({'target':target,'fold_a':i,'fold_b':j,'component_a':int(a),
                    'component_b':int(b),'absolute_cosine':float(cosine[a,b]),'min_subspace_cosine':float(subspace.min())})
        stability[target]=float(np.median(similarities))
    base.tsv_write(WORK/'dictionary_stability.tsv',matching)
    oldbaseline=['old_'+n for n in oldc['baselines']]+['diag_residual_point']
    fixedbase=['fixed_point_hgb','fixed_point_ridge','gene_ridge','fixed_point_hgb_smooth_half','gene_ridge_smooth_half']
    candidates={'diag_local_hgb':oldbaseline,'diag_residual_local':oldbaseline,
        'fixed_local_hgb':fixedbase,'pca_point_hgb':fixedbase,
        'pca_local_hgb':fixedbase+['pca_point_hgb','pca_point_ridge','pca_point_hgb_smooth_half']}
    comparisons=[];summary=[];models=[]
    for target in c['targets']:
        units=sorted({u for t,u,m in unit if t==target})
        for model in sorted({m for t,u,m in unit if t==target}):
            models.append({'target':target,'model':model,'n_units':len(units),
                           'macro_unit_rmse':float(np.mean([unit[target,u,model]['rmse'] for u in units]))})
        for candidate,baselines in candidates.items():
            selected=[]
            for u in units:
                p=unit[target,u,candidate];best=min(baselines,key=lambda m:unit[target,u,m]['rmse'])
                b=unit[target,u,best]['rmse'];cb=min(unit[target,u,m]['centered_rmse'] for m in baselines)
                r={'target':target,'model':candidate,'unit_id':u,'fold':p['fold'],'best_ordinary':best,
                   'relative_gain':1-p['rmse']/max(b,1e-12),'centered_gain':1-p['centered_rmse']/max(cb,1e-12)}
                selected.append(r);comparisons.append(r)
            gain=float(np.mean([r['relative_gain'] for r in selected]));frac=float(np.mean([r['relative_gain']>0 for r in selected]))
            center=float(np.mean([r['centered_gain'] for r in selected]));positive=sum(np.mean([r['relative_gain'] for r in selected if r['fold']==f])>0 for f in c['folds'])
            stable=not candidate.startswith('pca') or stability[target]>=c['screen']['dictionary_median_matched_cosine']
            passed=gain>=.03 and frac>=.6 and positive>=4 and center>0 and stable
            summary.append({'target':target,'model':candidate,'macro_relative_gain':gain,'positive_fraction':frac,
                 'positive_folds':int(positive),'centered_gain':center,'median_loading_cosine':stability[target] if candidate.startswith('pca') else None,
                 'status':'INTERNAL_CANDIDATE' if passed else 'NOT_RETAINED_BY_FIXED_RULE'})
    base.tsv_write(WORK/'unit_comparisons.tsv',comparisons);base.tsv_write(WORK/'model_comparison.tsv',models)
    base.tsv_write(WORK/'relations.tsv',summary);base.json_write(WORK/'summary.json',{'results':summary,'not_run':c['not_run']})
    (WORK/'STATUS').write_text('FOLLOWUP_COMPLETE_FROZEN_NO_HUMAN_TLS\n')
    files={str(p.relative_to(WORK)):{'sha256':base.sha256(p),'bytes':p.stat().st_size} for p in sorted(WORK.rglob('*'))
           if p.is_file() and p.name!='freeze_manifest.json' and p.suffix!='.log'}
    base.json_write(WORK/'freeze_manifest.json',{'files':files,'code_sha256':base.sha256(__file__),
                    'prior_freeze_sha256':base.sha256(OLD/'freeze_manifest.json'),'human_tls_labels_read':False})
    print(json.dumps(summary,indent=2),flush=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('command',choices=['prepare','run','summarize']);a=p.parse_args()
    if a.command=='prepare':prepare()
    elif a.command=='summarize':summarize()
    else:
        c=json.loads((WORK/'contract.json').read_text());base.reserve(c)
        (WORK/'STATUS').write_text('RUNNING\n')
        with ProcessPoolExecutor(max_workers=c['workers'],mp_context=multiprocessing.get_context('spawn')) as pool:
            jobs=[pool.submit(run_fold,f) for f in c['folds']]
            for job in as_completed(jobs):print(json.dumps(job.result()),flush=True)
        summarize()


if __name__=='__main__':main()
