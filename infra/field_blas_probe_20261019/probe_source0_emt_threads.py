#!/usr/bin/env python3
"""Controlled BLAS-thread probe for the frozen-baseline reproduction gap (review N-01).

Scope fixed by the numerical review's required_fix.minimal_verification:
  * only source_0 / HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION (target index 3)
  * one training design shared by every setting: same sampled rows, row order, target,
    sample weights, training-only constant-column mask, dtype, alpha and solver
  * compare effective BLAS thread counts 1 and 8 (recorded), each fitted twice
  * report point and derived fixed-smoothing predictions for NCBI642 and NCBI643 against
    (a) each other, (b) the frozen gene_ridge / gene_ridge_smooth_half package,
    (c) the currently stored run package
  * record whether the design facts match the committed fit_diagnostics row, and record the
    design / package paths actually used, so a design or save-path difference is not blamed
    on the environment

Writes only into this probe directory. Never modifies the accepted implementation, its machine
contract or any stored run artifact. No evaluate / summarize, no tolerance change, no solver
change, no broad scan.
"""
import argparse
import csv
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_info, threadpool_limits

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import field_raw_gene_spatial_20261019 as R          # noqa: E402  (frozen implementation)
import field_followup_20261018 as follow             # noqa: E402
import field_discovery_20261018 as base              # noqa: E402

RUN = ROOT / 'infra/field_raw_gene_spatial_20261019'
PROBE = ROOT / 'infra/field_blas_probe_20261019'
PRED = ROOT / 'infra/field_source_holdout_20261018'
SPLIT = 'source_0'
TARGET_INDEX = 3
SECTIONS = ['NCBI642', 'NCBI643']
THREADS = [1, 8]
REPEATS = 2


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def max_abs(a, b):
    return float(np.max(np.abs(np.asarray(a, np.float64) - np.asarray(b, np.float64))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe-dir', default=str(PROBE))
    args = parser.parse_args()
    probe = Path(args.probe_dir)
    probe.mkdir(parents=True, exist_ok=True)
    # the frozen module reads its machine contract and guard output from R.WORK
    (probe / 'contract.json').write_bytes((RUN / 'contract.json').read_bytes())
    R.WORK = probe
    ctx = R.setup(verify='none', sids=SECTIONS)
    # the probe must read the accepted run's own stored artifacts; allowlist them explicitly
    accepted_outputs = [p for p in (RUN / 'splits').rglob('*') if p.is_file()]
    accepted_outputs += [RUN / 'splits.json', RUN / 'contract.json', RUN / 'input_checks.json']
    ctx['guard'].allowed.update(p.resolve() for p in accepted_outputs if p.exists())
    c = ctx['c']
    split = next(s for s in ctx['splits'] if s['name'] == SPLIT)
    target = c['targets'][TARGET_INDEX]

    started = time.time()
    train = [ctx['byid'][s] for s in split['train']]
    idx = {r['section_id']: base.sample_indices(r['n_spots'], r['section_id'], c['seed'],
                                                c['train_spots_per_section']) for r in train}
    weights = follow.weights(train, idx)
    perms = {s: R.permutation_for(s, c['seed'], ctx['byid'][s]['n_spots']) for s in ctx['byid']}
    y_train = R.target_vector(train, idx, TARGET_INDEX)
    stacks, widths = R.stack_training(ctx, train, idx, TARGET_INDEX, perms, ['point_ridge'])
    design = stacks['point_ridge']
    scaled_width = widths['point_ridge']
    std = design.std(axis=0)
    keep = std > 0
    dropped = int(np.count_nonzero(~keep))

    # design identity and the committed fit_diagnostics row must agree
    diag_path = RUN / 'splits' / SPLIT / 'fit_diagnostics.tsv'
    diag = [r for r in csv.DictReader(diag_path.open(), delimiter='\t')
            if r['target'] == target and r['model'] == 'point_ridge'][0]
    design_identity = dict(
        split=SPLIT, target=target, target_index=TARGET_INDEX,
        sections=list(SECTIONS),
        n_train_rows=int(design.shape[0]), n_features=int(design.shape[1]),
        n_columns_dropped_train_constant=dropped,
        n_context_columns=R.CONTEXT_DIM, scaled_width=scaled_width,
        dtype=str(design.dtype), y_dtype=str(y_train.dtype),
        weights_sum=float(weights.sum()), weights_len=int(len(weights)),
        design_sha256=array_digest(design), y_sha256=array_digest(y_train),
        weights_sha256=array_digest(weights),
        keep_mask_sha256=array_digest(keep.astype(np.uint8)),
        diagnostic_row={k: diag[k] for k in
                        ['n_train_rows', 'n_features', 'n_scaled_columns', 'n_context_columns',
                         'n_columns_dropped_train_constant', 'solver', 'alpha', 'block_scale',
                         'train_rmse', 'converged']},
        matches_diagnostics=dict(
            rows=int(diag['n_train_rows']) == int(design.shape[0]),
            features=int(diag['n_features']) == int(design.shape[1]),
            dropped=int(diag['n_columns_dropped_train_constant']) == dropped,
            scaled_columns=int(diag['n_scaled_columns']) == scaled_width,
            context_columns=int(diag['n_context_columns']) == R.CONTEXT_DIM),
        sources=dict(
            genes={f'{sid}_t{ti}': sha256(ROOT / f'infra/field_followup_20261018/genes/{sid}_t{ti}.npy')
                   for sid in [SECTIONS[0]] for ti in [TARGET_INDEX]},
            splits_run=sha256(RUN / 'splits.json'), splits_frozen=sha256(PRED / 'splits.json'),
            splits_identical=(sha256(RUN / 'splits.json') == sha256(PRED / 'splits.json')),
            contract_sha256=sha256(RUN / 'contract.json'),
            implementation_sha256=sha256(ROOT / 'scripts/field_raw_gene_spatial_20261019.py'),
            frozen_package_templates=str(PRED / SPLIT / '{sid}_t3_full15.npz'),
            stored_package_templates=str(RUN / 'splits' / SPLIT / '{sid}_t3.npz')),
        env={k: os.environ.get(k) for k in
             ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS',
              'CUDA_VISIBLE_DEVICES']},
    )
    np.savez_compressed(probe / 'sampled_rows.npz', **idx)

    # frozen and stored packages, and the frozen graphs for the derived smoothing
    reference = {}
    for sid in SECTIONS:
        with np.load(PRED / SPLIT / f'{sid}_t{TARGET_INDEX}_full15.npz', allow_pickle=False) as z:
            frozen = {k: z[k].astype(np.float64)
                      for k in ['gene_ridge', 'gene_ridge_smooth_half', 'y']}
            frozen['barcodes'] = z['barcodes']
        with np.load(RUN / 'splits' / SPLIT / f'{sid}_t{TARGET_INDEX}.npz', allow_pickle=False) as z:
            stored = {k: z[k].astype(np.float64)
                      for k in ['point_ridge', 'point_ridge_smooth_half', 'y']}
            stored['barcodes'] = z['barcodes']
        cache = base.load_cache(ctx['byid'][sid])
        reference[sid] = dict(frozen=frozen, stored=stored, near=cache['near'])
        assert np.array_equal(frozen['barcodes'], stored['barcodes']), sid

    fits = {}
    for threads in THREADS:
        for rep in range(REPEATS):
            with threadpool_limits(limits=threads):
                info = threadpool_info()
                effective = sorted({int(lib.get('num_threads', -1)) for lib in info
                                    if lib.get('num_threads')})
                t0 = time.time()
                model = R.fit_linear(c, design, y_train, weights, scaled_width)
                seconds = time.time() - t0
            predictions = {}
            for sid in SECTIONS:
                blocks = R.section_designs(ctx['byid'][sid], TARGET_INDEX,
                                           base.load_cache(ctx['byid'][sid]),
                                           perms[sid], ['point_ridge'])['point_ridge']['design']
                point = R.model_predict(model, blocks)
                smooth = 0.5 * (point + base.neighbor_mean(reference[sid]['near'], point))
                predictions[sid] = dict(point=point.astype(np.float64),
                                        smooth=smooth.astype(np.float64))
            fits[(threads, rep)] = dict(
                seconds=seconds, effective_threads=effective,
                train_rmse=model['train_rmse'],
                coef_l2=float(np.linalg.norm(model['ridge'].coef_)),
                intercept=float(model['ridge'].intercept_),
                scaler=model['scaler'], keep=model['keep'], scaled_width=model['scaled_width'],
                block_scale=model['block_scale'],
                threadpool_info=info, predictions=predictions)

    # condition numbers of the design actually used (script for the earlier probe numbers)
    model0 = fits[(THREADS[0], 0)]
    scaled = R.prepare_scaled(dict(scaler=model0['scaler'], keep=model0['keep'],
                                   scaled_width=model0['scaled_width'],
                                   block_scale=model0['block_scale']), design)
    gram = scaled.T @ scaled
    eig = np.linalg.eigvalsh(gram)
    conditioning = dict(gram_min_eig=float(eig.min()), gram_max_eig=float(eig.max()),
                        raw_rcond=float(eig.min() / eig.max()),
                        ridge_rcond=float((eig.min() + c['alpha']) / (eig.max() + c['alpha'])))

    rows = []
    for threads in THREADS:
        other = 8 if threads == 1 else 1
        for rep in range(REPEATS):
            fit = fits[(threads, rep)]
            base_fit = fits[(threads, 0)]
            for sid in SECTIONS:
                for kind in ['point', 'smooth']:
                    pred = fit['predictions'][sid][kind]
                    frozen_key = 'gene_ridge_smooth_half' if kind == 'smooth' else 'gene_ridge'
                    stored_key = 'point_ridge_smooth_half' if kind == 'smooth' else 'point_ridge'
                    rows.append(dict(
                        threads=threads, repeat=rep, section=sid, kind=kind,
                        effective_threads=';'.join(str(t) for t in fit['effective_threads']),
                        vs_repeat0=max_abs(pred, base_fit['predictions'][sid][kind]),
                        vs_other_thread=max_abs(pred, fits[(other, rep)]['predictions'][sid][kind]),
                        vs_frozen=max_abs(pred, reference[sid]['frozen'][frozen_key]),
                        vs_stored=max_abs(pred, reference[sid]['stored'][stored_key]),
                        train_rmse=fit['train_rmse'], coef_l2=fit['coef_l2'],
                        fit_seconds=round(fit['seconds'], 2)))

    summary = dict(
        scope=dict(split=SPLIT, target=target, target_index=TARGET_INDEX, sections=SECTIONS,
                   threads=THREADS, repeats=REPEATS,
                   solver=c['solver'], alpha=c['alpha'], seed=c['seed'],
                   train_spots_per_section=c['train_spots_per_section']),
        design_identity=design_identity, conditioning=conditioning,
        within_setting_floor=max(row['vs_repeat0'] for row in rows),
        across_threads={f'{sid}/{kind}': max(row['vs_other_thread'] for row in rows
                                            if row['section'] == sid and row['kind'] == kind)
                        for sid in SECTIONS for kind in ['point', 'smooth']},
        vs_frozen={f'{sid}/{kind}': max(row['vs_frozen'] for row in rows
                                        if row['section'] == sid and row['kind'] == kind)
                   for sid in SECTIONS for kind in ['point', 'smooth']},
        vs_stored={f'{sid}/{kind}': max(row['vs_stored'] for row in rows
                                        if row['section'] == sid and row['kind'] == kind)
                   for sid in SECTIONS for kind in ['point', 'smooth']},
        threadpool_effective={f'{t}': sorted({','.join(str(x) for x in f['effective_threads'])
                                              for (tt, _), f in fits.items() if tt == t})
                              for t in THREADS},
        elapsed_seconds=round(time.time() - started, 1),
        note='within_setting_floor is the same-design repeat floor; across_threads is the controlled '
             'thread contrast; vs_frozen and vs_stored compare with the frozen and stored packages')
    with (probe / 'probe_results.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    (probe / 'probe_results.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'design_identity'}, indent=2))
    print(json.dumps(design_identity, indent=2), flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
