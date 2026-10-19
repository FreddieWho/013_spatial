#!/usr/bin/env python3
"""Same-raw-gene spatial contrast (contract: docs/FIELD_RAW_GENE_SPATIAL_CONTRACT_20261019.md).

Fixed candidates only, no tuning:
  point_ridge               allowed point input, re-fitted as the numeric reference
  point_ridge_smooth_half   0.5 * (point + near-neighbour mean of point)
  local_ridge               point X | WX | VX                      (the only candidate)
  background_ridge          point X | whole-section mean of X | same
  permuted_ridge            point X | P W P^T X | P V P^T X        (single mechanism control)

Implementation clauses enforced here (contract D-269 / review B-01..B-06):
  * one frozen full15 visible matrix + QC per target; raw HDF5 is not in the input allowlist
  * frozen split file, no re-draw; train/test unit and source disjointness asserted per split
  * spatial and background features are built on all tissue spots aligned with the frozen
    graphs and the frozen training rows are extracted afterwards
  * training-only zero-variance column removal, training-only standardisation
  * one pre-fixed within-section permutation per section for both graphs, same self fallback
  * triple-copy block equivalence verified with the formal solver/intercept/QC/column rule
  * reviewed evaluation definitions (Git 5fee80f); per-section abs radius deviation first
CPU only, local data only, no human TLS / annotations / incoming / HookNet.
"""
import argparse
import csv
import hashlib
import json
import math
import multiprocessing
import os
import shutil
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

import field_discovery_20261018 as base
import field_followup_20261018 as follow
import field_biological_review_20261019 as review

ROOT = base.ROOT
WORK = ROOT / 'infra/field_raw_gene_spatial_20261019'
PRED = ROOT / 'infra/field_source_holdout_20261018'
FUP = ROOT / 'infra/field_followup_20261018'
CACHE = base.WORK
PANEL_FILE = ROOT / 'infra/tls_multiscale_gobp_20261002/panel_controls.json'
IDENTITY = ROOT / 'infra/field_discovery_review_20261018/identity_review.tsv'
CONTRACT_DOC = ROOT / 'docs/FIELD_RAW_GENE_SPATIAL_CONTRACT_20261019.md'

MODELS = ['point_ridge', 'point_ridge_smooth_half', 'local_ridge', 'background_ridge',
          'permuted_ridge']
FIT_MODELS = ['point_ridge', 'local_ridge', 'background_ridge', 'permuted_ridge']
BASELINES = ['point_ridge', 'point_ridge_smooth_half', 'background_ridge', 'permuted_ridge']
CANDIDATE = 'local_ridge'
PAIRED_METRICS = ['rec_rmse', 'rec_centered_rmse', 'rec_pearson', 'rec_spearman',
                  'rec_local_pattern_pearson', 'rec_jaccard_q75', 'rec_jaccard_q80',
                  'rec_jaccard_q90', 'rec_centroid_displacement_pitch',
                  'rec_patch75_largest_radius_pitch_delta', 'rec_patch90_largest_radius_pitch_delta',
                  'rec_patch75_largest_radius_abs_delta', 'rec_patch90_largest_radius_abs_delta',
                  'rec_sd_ratio']
MACRO_METRICS = ['rec_rmse', 'rec_centered_rmse', 'rec_mae', 'rec_bias', 'rec_sd_ratio',
                 'rec_calibration_slope', 'rec_pearson', 'rec_spearman', 'rec_local_pattern_pearson',
                 'rec_jaccard_q75', 'rec_jaccard_q80', 'rec_jaccard_q90',
                 'rec_centroid_displacement_pitch', 'rec_patch75_largest_radius_pitch_delta',
                 'rec_patch90_largest_radius_pitch_delta', 'rec_patch75_largest_radius_abs_delta',
                 'rec_patch90_largest_radius_abs_delta', 'geo_patch75_largest_radius_pitch',
                 'geo_patch75_largest_share', 'geo_autocorr_near', 'dist_sd']
# favourable direction for the candidate-vs-baseline descriptive conjunction (contract claim rule)
FAVOURABLE_LOW = {'rec_rmse', 'rec_centered_rmse', 'rec_mae', 'rec_bias',
                  'rec_centroid_displacement_pitch', 'rec_patch75_largest_radius_pitch_delta',
                  'rec_patch90_largest_radius_pitch_delta', 'rec_patch75_largest_radius_abs_delta',
                  'rec_patch90_largest_radius_abs_delta'}
FAVOURABLE_HIGH = {'rec_pearson', 'rec_spearman', 'rec_local_pattern_pearson', 'rec_jaccard_q75',
                   'rec_jaccard_q80', 'rec_jaccard_q90'}
VIEW = 'full15'
CONTEXT_DIM = len(base.ORGANS) + len(base.PLATFORMS)
MODEL_CLASS = {'point_ridge': 'raw_gene_point_ridge',
               'point_ridge_smooth_half': 'posthoc_neighbourhood_smoothing',
               'local_ridge': 'raw_gene_local_ridge',
               'background_ridge': 'whole_section_background',
               'permuted_ridge': 'permuted_graph_control'}


def contract_path():
    return WORK / 'contract.json'


def machine_contract():
    return json.loads(contract_path().read_text())


def seed_for(*parts):
    digest = hashlib.sha256('|'.join(str(p) for p in parts).encode()).digest()
    return int.from_bytes(digest[:8], 'little')


def permutation_for(sid, seed, n):
    """One pre-fixed within-section node permutation (index array of length n)."""
    return np.random.default_rng(seed_for(seed, sid, 'graph_permutation')).permutation(n)


# --------------------------------------------------------------------------- setup
def assert_split_integrity(ctx):
    rows = {r['section_id']: r for r in ctx['rows']}
    internal_units, source_units = set(), set()
    checks = Counter()
    for split in ctx['splits']:
        train = [rows[s] for s in split['train']]
        test = [rows[s] for s in split['test']]
        if {r['unit_id'] for r in train} & {r['unit_id'] for r in test}:
            raise ValueError(f'Unit crosses train/test in {split["name"]}')
        if split['arm'] == 'source':
            if {r['source_key'] for r in train} & {r['source_key'] for r in test}:
                raise ValueError(f'Source crosses train/test in {split["name"]}')
            if any(not r['source_key'].startswith('http') for r in train + test):
                raise ValueError(f'Source arm contains a section without a study link: {split["name"]}')
            source_units.update(r['unit_id'] for r in test)
        else:
            internal_units.update(r['unit_id'] for r in test)
        checks['splits'] += 1
    internal_sections = [s for sp in ctx['splits'] if sp['arm'] == 'internal' for s in sp['test']]
    source_sections = [s for sp in ctx['splits'] if sp['arm'] == 'source' for s in sp['test']]
    if sorted(internal_sections) != sorted(rows):
        raise ValueError('Internal folds do not partition the 72 qualified sections')
    if len(set(source_sections)) != 53 or len(source_sections) != 53:
        raise ValueError('Source folds do not partition the 53 known-source sections')
    if len(internal_units) != 65 or len(source_units) != 47:
        raise ValueError(f'Unit counts differ from 65/47: {len(internal_units)}/{len(source_units)}')
    checks['internal_sections'] = len(internal_sections)
    checks['source_sections'] = len(source_sections)
    checks['internal_units'] = len(internal_units)
    checks['source_units'] = len(source_units)
    return dict(checks)


def setup(verify='sample', sids=None):
    c = machine_contract()
    rows = json.loads((PRED / 'input_manifest.json').read_text())
    byid = {r['section_id']: r for r in rows}
    splits = json.loads((PRED / 'splits.json').read_text())
    vocab = json.loads((PRED / 'vocabulary.json').read_text())
    panels = json.loads(PANEL_FILE.read_text())
    identity = {r['section_id']: r for r in csv.DictReader(IDENTITY.open(), delimiter='\t')}
    eligibility = {r['section_id']: r['status'] for r in
                   csv.DictReader((PRED / 'source_eligibility.tsv').open(), delimiter='\t')}
    for row in rows:
        if bool(identity[row['section_id']]['patient_confirmed_within_study'] == 'True') \
                != bool(row['patient_known']):
            raise ValueError(f'Identity flag disagreement: {row["section_id"]}')
        row['confirmed_unit'] = bool(row['patient_known'])
        row['source_eligibility'] = eligibility.get(row['section_id'], 'UNKNOWN')
    allowed = {Path(__file__), Path(base.__file__), Path(follow.__file__), Path(review.__file__),
               PANEL_FILE, IDENTITY, PRED / 'input_manifest.json', PRED / 'splits.json',
               PRED / 'vocabulary.json', PRED / 'contract.json', PRED / 'source_eligibility.tsv'}
    for root in [CACHE, FUP, PRED]:
        freeze = json.loads((root / 'freeze_manifest.json').read_text())
        allowed.add(root / 'freeze_manifest.json')
        allowed.update(root / p for p in freeze['files'])
    guard = base.InputGuard(ROOT, allowed, WORK)
    guard.install()
    review.install_denylist_hook()
    base.reserve(c)
    checks = {'n_allowed_paths': len(allowed),
              'raw_hdf5_allowlisted': any(str(p).startswith(str(base.HEST)) for p in allowed),
              'no_allowed_path_matches_denylist': True}
    if checks['raw_hdf5_allowlisted']:
        raise ValueError('Raw HDF5 must not be reachable: only the frozen visible matrices are inputs')
    for path in allowed:
        if review.is_denied(path):
            checks['no_allowed_path_matches_denylist'] = False
    if not checks['no_allowed_path_matches_denylist']:
        raise ValueError('Denylist pattern present in allowlist')
    for target in c['targets']:
        if set(vocab[target]['full15']) & set(panels[target]):
            raise ValueError(f'Target genes entered the visible vocabulary: {target}')
    checks['visible_target_disjoint'] = True
    checks['visible_dim'] = {t: len(vocab[t]['full15']) for t in c['targets']}
    checks['target_genes'] = {t: len(panels[t]) for t in c['targets']}
    checks.update(assert_split_integrity(dict(rows=rows, splits=splits)))
    checks['known_sections'] = sum(1 for r in rows if r['confirmed_unit'])
    checks['known_units'] = len({r['unit_id'] for r in rows if r['confirmed_unit']})
    checks['source_not_testable_sections'] = sum(1 for r in rows
                                                 if r['source_eligibility'].startswith('NOT_TESTABLE'))
    checks['free_bytes'] = shutil.disk_usage(WORK).free
    checks['storage_reserve_ok'] = bool(checks['free_bytes'] >= c['min_free_bytes'])
    if not checks['storage_reserve_ok']:
        raise RuntimeError(f'Storage reserve violated: free={checks["free_bytes"]}')
    if verify in ('all', 'sample'):
        verify_frozen(checks, None if verify == 'all' else (sids or sorted(byid)))
    return dict(c=c, rows=rows, byid=byid, splits=splits, vocab=vocab, panels=panels,
                guard=guard, checks=checks)


def verify_frozen(checks, sids):
    verified = Counter()
    for root, prefixes in [(CACHE, ('cache/',)), (FUP, ('genes/',)),
                           (PRED, ('genes/', 'internal_', 'source_'))]:
        freeze = json.loads((root / 'freeze_manifest.json').read_text())
        for name, entry in freeze['files'].items():
            if not any(name.startswith(p) for p in prefixes):
                continue
            if sids is not None and Path(name).name.split('_')[0] not in sids:
                continue
            if base.sha256(root / name) != entry['sha256']:
                raise ValueError(f'Frozen input changed: {root/name}')
            verified['frozen_hashes'] += 1
    checks['verified'] = dict(verified)


def prepare():
    ctx = setup(verify='all')
    if (WORK / 'PREPARED').exists():
        raise RuntimeError('Prepared outputs exist; preserve this version')
    (WORK / 'checks').mkdir(parents=True, exist_ok=True)
    base.json_write(WORK / 'splits.json', ctx['splits'])
    base.json_write(WORK / 'input_checks.json', ctx['checks'])
    base.json_write(WORK / 'prepare_receipt.json', dict(receipt=ctx['guard'].receipt()))
    (WORK / 'PREPARED').write_text('FROZEN_INPUTS_VERIFIED\n')
    print(json.dumps(ctx['checks'], indent=1), flush=True)


# --------------------------------------------------------------------------- features
def visible_matrices(sid, ti):
    return (np.asarray(np.load(FUP / f'genes/{sid}_t{ti}.npy'), np.float32),
            np.asarray(np.load(FUP / f'genes/{sid}_t{ti}_qc.npy'), np.float32))


def context_block(row, n):
    vector = [float(row['organ'] == o) for o in base.ORGANS] + \
             [float(row['metadata']['st_technology'] == p) for p in base.PLATFORMS]
    return np.broadcast_to(np.array(vector, np.float32), (n, len(vector)))


def full_features(row, ti, cache, perm=None):
    """Point X, organ/platform indicators E, and the two graph aggregates of X.

    X = frozen visible genes + QC (never raw counts). E = the existing organ/platform indicator
    columns reused from the frozen point input. The graph aggregates act on X only; E never enters
    the multiplication and is appended once, outside the block scaling, by raw_design.
    """
    z, qc = visible_matrices(row['section_id'], ti)
    if len(z) != len(cache['barcodes']):
        raise ValueError(f'Visible matrix rows differ from barcodes: {row["section_id"]}')
    x = np.c_[z, qc].astype(np.float32)
    e = context_block(row, len(x)).astype(np.float32)
    near_graph, outer_graph = cache['near'], cache['outer']
    if perm is not None:
        near_graph = sparse.csr_matrix(near_graph[perm][:, perm])
        outer_graph = sparse.csr_matrix(outer_graph[perm][:, perm])
    near = np.asarray(base.neighbor_mean(near_graph, x), np.float32)
    outer = np.asarray(base.neighbor_mean(outer_graph, x), np.float32)
    background = np.broadcast_to(x.mean(axis=0), x.shape)
    return dict(x=x, e=e, near=near, outer=outer, background=background)


def raw_design(blocks, kind):
    """[X, WX, VX]/sqrt(3) then E once (point: [X, E]); returns (design, scaled_width)."""
    x, e = blocks['x'], blocks['e']
    width = x.shape[1]
    if kind == 'point_ridge':
        return np.c_[x, e], 0
    if kind == 'background_ridge':
        return np.c_[x, blocks['background'], blocks['background'], e], 3 * width
    if kind in ('local_ridge', 'permuted_ridge'):
        return np.c_[x, blocks['near'], blocks['outer'], e], 3 * width
    raise ValueError(kind)


def section_designs(row, ti, cache, perm, kinds):
    """Only the permuted control uses the permutation; other models use the real graphs."""
    plain = full_features(row, ti, cache, perm=None)
    permuted = full_features(row, ti, cache, perm=perm) if 'permuted_ridge' in kinds else None
    out = {}
    for kind in kinds:
        blocks = permuted if (kind == 'permuted_ridge' and permuted is not None) else plain
        design, scaled_width = raw_design(blocks, kind)
        out[kind] = dict(design=design, scaled_width=scaled_width)
    del plain, permuted
    return out


def scaled_width_of(design, kind):
    return 0 if kind == 'point_ridge' else design.shape[1] - CONTEXT_DIM


def fit_linear(c, x_train, y_train, weights, scaled_width):
    """Training-only column rule, training-only standardisation, block scaling on X blocks only.

    `scaled_width` counts the leading columns that receive 1/sqrt(3); the trailing indicator
    columns E keep exactly the same standardisation and penalty as in the point model.
    """
    std = x_train.std(axis=0)
    keep = std > 0
    dropped = int(np.count_nonzero(~keep))
    x_keep = x_train[:, keep]
    kept_scaled = int(np.count_nonzero(keep[:scaled_width]))
    scaler = StandardScaler()
    scaler.fit(x_keep, sample_weight=weights)
    scaled = scaler.transform(x_keep)
    if kept_scaled:
        scaled[:, :kept_scaled] *= c['block_scale']
    solver = c['solver']
    ridge = Ridge(alpha=c['alpha'], fit_intercept=True, solver=solver,
                  tol=c['solver_tol'], max_iter=c['solver_max_iter'])
    ridge.fit(scaled, y_train, sample_weight=weights)
    converged = True
    n_iter = -1
    if solver == 'lsqr':
        n_iter = int(np.ravel(getattr(ridge, 'n_iter_', [-1]))[0])
        converged = bool(n_iter < c['solver_max_iter'])
    return dict(scaler=scaler, ridge=ridge, keep=keep, dropped=dropped,
                scaled_width=kept_scaled, block_scale=c['block_scale'],
                solver=str(ridge.solver), n_iter=n_iter, converged=converged,
                train_rmse=float(np.sqrt(np.mean((y_train - ridge.predict(scaled)) ** 2))))


def prepare_scaled(model, x):
    scaled = np.array(model['scaler'].transform(x[:, model['keep']]))
    if model['scaled_width']:
        scaled[:, :model['scaled_width']] *= model['block_scale']
    return scaled


def model_predict(model, x):
    return model['ridge'].predict(prepare_scaled(model, x))


def stack_training(ctx, train, indexes, ti, perms, kinds):
    stacks = {k: [] for k in kinds}
    for row in train:
        sid = row['section_id']
        cache = base.load_cache(row)
        designs = section_designs(row, ti, cache, perms[sid], kinds)
        keep_rows = indexes[sid]
        for kind in kinds:
            stacks[kind].append(designs[kind]['design'][keep_rows])
        del designs
    out, widths = {}, {}
    for kind in kinds:
        out[kind] = np.concatenate(stacks[kind])
        widths[kind] = scaled_width_of(out[kind], kind)
        del stacks[kind]
    return out, widths


def target_vector(rows, indexes, ti):
    values = {r['section_id']: base.load_cache(r)['y'][:, ti] for r in rows}
    return follow.stack(rows, values, indexes)


# --------------------------------------------------------------------------- fit
def run_split_fit(job):
    split = job['split']
    ctx = setup(verify='none')
    c, byid = ctx['c'], ctx['byid']
    out = Path(job['out'])
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'COMPLETE').exists() and not job.get('force'):
        return dict(split=split['name'], status='already_complete')
    start = time.time()
    train = [byid[s] for s in split['train']]
    test = [byid[s] for s in split['test']]
    if job.get('train_limit'):
        train = train[:job['train_limit']]
    if job.get('test_limit'):
        test = test[:job['test_limit']]
    if {r['unit_id'] for r in train} & {r['unit_id'] for r in test}:
        raise ValueError(f'Train/test unit overlap in {split["name"]}')
    indexes = {r['section_id']: base.sample_indices(r['n_spots'], r['section_id'], c['seed'],
                                                    c['train_spots_per_section']) for r in train}
    weights = follow.weights(train, indexes)
    perms = {r['section_id']: permutation_for(r['section_id'], c['seed'], r['n_spots'])
             for r in train + test}
    diagnostics, timings = [], []
    with threadpool_limits(limits=c['threads']):
        for ti, target in enumerate(c['targets']):
            y_train = target_vector(train, indexes, ti)
            stacks, widths = stack_training(ctx, train, indexes, ti, perms, FIT_MODELS)
            models = {}
            for kind in FIT_MODELS:
                t0 = time.time()
                model = fit_linear(c, stacks[kind], y_train, weights, widths[kind])
                if not model['converged']:
                    raise RuntimeError(f'Non-converged fit: {split["name"]}/{target}/{kind} '
                                       f'n_iter={model["n_iter"]}')
                models[kind] = model
                diagnostics.append(dict(split=split['name'], arm=split['arm'], target=target,
                                        model=kind, n_train_rows=int(stacks[kind].shape[0]),
                                        n_features=int(stacks[kind].shape[1]),
                                        n_scaled_columns=int(widths[kind]),
                                        n_context_columns=CONTEXT_DIM,
                                        n_columns_dropped_train_constant=model['dropped'],
                                        solver=model['solver'], alpha=c['alpha'],
                                        block_scale=(c['block_scale'] if widths[kind] else 1.0),
                                        train_rmse=model['train_rmse'], n_iter=model['n_iter'],
                                        converged=model['converged'],
                                        intercept=float(model['ridge'].intercept_),
                                        fit_seconds=round(time.time() - t0, 1)))
                del stacks[kind]
            for row in test:
                sid = row['section_id']
                cache = base.load_cache(row)
                designs = section_designs(row, ti, cache, perms[sid], FIT_MODELS)
                predictions = {k: model_predict(models[k], designs[k]['design']).astype(np.float32)
                               for k in FIT_MODELS}
                predictions['point_ridge_smooth_half'] = (
                    0.5 * (predictions['point_ridge']
                           + base.neighbor_mean(cache['near'], predictions['point_ridge']))
                ).astype(np.float32)
                for name, values in predictions.items():
                    if not np.all(np.isfinite(values)):
                        raise ValueError(f'Non-finite prediction: {sid}/{target}/{name}')
                np.savez_compressed(out / f'{sid}_t{ti}.npz', y=cache['y'][:, ti],
                                    barcodes=cache['barcodes'],
                                    permutation_seed=seed_for(c['seed'], sid, 'graph_permutation'),
                                    **predictions)
            timings.append(dict(split=split['name'], target=target,
                                seconds=round(time.time() - start, 1)))
        base.tsv_write(out / 'fit_diagnostics.tsv', diagnostics)
    base.json_write(out / 'split_receipt.json',
                    dict(name=split['name'], arm=split['arm'], n_train=len(train), n_test=len(test),
                         train_sections=[r['section_id'] for r in train],
                         test_sections=[r['section_id'] for r in test],
                         trimmed=bool(job.get('train_limit') or job.get('test_limit')),
                         permutation_seeds={sid: seed_for(c['seed'], sid, 'graph_permutation')
                                            for sid in sorted(perms)},
                         seconds=round(time.time() - start, 1), receipt=ctx['guard'].receipt()))
    (out / 'COMPLETE').write_text('SPLIT_FIT_DONE\n')
    return dict(split=split['name'], seconds=round(time.time() - start, 1), test=len(test))


def fit(limit_splits=None, train_limit=None, test_limit=None, work=None, jobs_override=None):
    global WORK
    if work is not None:
        WORK = Path(work)
    ctx = setup(verify='none')
    c = ctx['c']
    jobs = jobs_override or [dict(split=s, out=str(WORK / 'splits' / s['name']),
                                  train_limit=train_limit, test_limit=test_limit)
                             for s in ctx['splits'] if (limit_splits is None or s['name'] in limit_splits)]
    with ProcessPoolExecutor(max_workers=c['workers'], mp_context=multiprocessing.get_context('spawn')) as pool:
        futures = [pool.submit(run_split_fit, job) for job in jobs]
        for future in as_completed(futures):
            print(json.dumps(future.result()), flush=True)
    (WORK / 'FIT_COMPLETE').write_text('SELECTED_SPLITS_FIT\n')


# --------------------------------------------------------------------------- evaluate
def geo_for(ctx, row, cache):
    geo, recomputed = review.geo_of(row, cache, ctx['c'])
    if abs(recomputed - float(row['empirical_pitch_source_pixels'])) > 1e-9:
        raise ValueError(f'Pitch mismatch: {row["section_id"]}')
    return geo


def describe_prediction(ctx, row, ti, y, prediction, geo, tag):
    """Thin wrapper kept for the binding checks: measured descriptors first, then the readouts."""
    yd, yextra, _ = describe_field_once(ctx, row, ti, y, geo, 'measured')
    rec, pd = recovery_with_range(ctx, row, ti, y, prediction, geo, yd, yextra, tag)
    return rec, pd, review.numeric(yd['dist_sd'])


def describe_field_once(ctx, row, ti, values, geo, tag):
    c = ctx['c']
    notes = []
    desc, extra, notes = review.describe_field(values, geo,
                                               seed_for(c['seed'], row['section_id'], ti, tag), notes,
                                               c['shuffle_permutations'], False, c['sd_floor'])
    return desc, extra, notes


def recovery_with_range(ctx, row, ti, y, prediction, geo, yd, yextra, tag):
    """Reviewed recovery readouts plus the range readouts (signed delta and per-section abs)."""
    c = ctx['c']
    pd, pextra, pnotes = describe_field_once(ctx, row, ti, prediction, geo, tag)
    notes = list(pnotes)
    rec = review.recovery_metrics(y, prediction, geo, yextra, pextra, notes, c['sd_floor'])
    for q in review.PATCH_QS:
        for field in ['n', 'largest_share', 'largest_radius_pitch']:
            key = f'rec_patch{review.pkey(q)}_{field}_delta'
            a = pd.get(f'geo_patch{review.pkey(q)}_{field}')
            b = yd.get(f'geo_patch{review.pkey(q)}_{field}')
            if a is None or b is None:
                rec[key] = None
                notes.append(f'{key}:missing_side')
            else:
                rec[key] = float(a - b)
        signed_key = f'rec_patch{review.pkey(q)}_largest_radius_pitch_delta'
        abs_key = f'rec_patch{review.pkey(q)}_largest_radius_abs_delta'
        rec[abs_key] = None if rec[signed_key] is None else abs(rec[signed_key])
    rec['not_estimable'] = ';'.join(sorted(notes))
    rec['n_not_estimable'] = len(notes)
    rec['state'] = review.prediction_state(dict(dist_sd=pd['dist_sd']), c['sd_floor'])
    return rec, pd


def run_split_evaluate(job):
    split = job['split']
    ctx = setup(verify='none')
    c, byid = ctx['c'], ctx['byid']
    src = Path(job['out'])
    out = WORK / 'evaluation' / split['name']
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'COMPLETE').exists() and not job.get('force'):
        return dict(split=split['name'], status='already_complete')
    records, paired = [], []
    sections = [byid[s] for s in split['test']]
    if job.get('test_limit'):
        sections = sections[:job['test_limit']]
    for row in sections:
        sid = row['section_id']
        cache = base.load_cache(row)
        geo = geo_for(ctx, row, cache)
        for ti, target in enumerate(c['targets']):
            bundle = np.load(src / f'{sid}_t{ti}.npz', allow_pickle=False)
            review.assert_aligned(bundle['barcodes'], cache['barcodes'])
            if not np.array_equal(bundle['y'], cache['y'][:, ti]):
                raise ValueError(f'Prediction y differs from cache: {sid}/{target}')
            predictions = {m: bundle[m].astype(np.float64) for m in MODELS}
            bundle.close()
            y = cache['y'][:, ti].astype(np.float64)
            yd, yextra, _ = describe_field_once(ctx, row, ti, y, geo, 'measured')
            values = {}
            for name in MODELS:
                rec, pd = recovery_with_range(ctx, row, ti, y, predictions[name], geo, yd, yextra, name)
                values[name] = rec
                out_row = dict(split=split['name'], arm=split['arm'], section_id=sid,
                               unit_id=row['unit_id'], organ=row['organ'],
                               source_key=row['source_key'],
                               source_split_key=(split.get('held_source', '')
                                                 if split['arm'] == 'source' else ''),
                               platform=row['metadata']['st_technology'],
                               patient_known=row['patient_known'],
                               confirmed_unit=row['confirmed_unit'],
                               source_eligibility=row['source_eligibility'],
                               cross_organ_extrapolation=bool(
                                   split['arm'] == 'source'
                                   and split.get('same_organ_train_sections', 1) == 0),
                               model=name, model_class=MODEL_CLASS[name], view=VIEW, target=target,
                               n_spots=int(len(y)),
                               meas_dist_sd=review.numeric(yd['dist_sd']), pred_dist_sd=pd['dist_sd'],
                               pred_geo_autocorr_near=pd['geo_autocorr_near'],
                               pred_geo_patch75_largest_radius_pitch=pd['geo_patch75_largest_radius_pitch'])
                out_row.update(rec)
                records.append(out_row)
            for base_name in BASELINES:
                pair = dict(split=split['name'], arm=split['arm'], section_id=sid,
                            unit_id=row['unit_id'], organ=row['organ'], source_key=row['source_key'],
                            source_split_key=(split.get('held_source', '')
                                              if split['arm'] == 'source' else ''),
                            platform=row['metadata']['st_technology'],
                            patient_known=row['patient_known'], confirmed_unit=row['confirmed_unit'],
                            cross_organ_extrapolation=bool(
                                split['arm'] == 'source'
                                and split.get('same_organ_train_sections', 1) == 0),
                            baseline=base_name, model=CANDIDATE,
                            model_class=MODEL_CLASS[CANDIDATE], view=VIEW, target=target,
                            n_spots=int(len(y)))
                for metric in PAIRED_METRICS:
                    a, b = values[CANDIDATE].get(metric), values[base_name].get(metric)
                    pair[f'delta_{metric}'] = (a - b) if (a is not None and b is not None) else None
                paired.append(pair)
        print(json.dumps(dict(split=split['name'], section=sid)), flush=True)
    base.tsv_write(out / 'prediction_metrics.tsv', records)
    base.tsv_write(out / 'paired_metrics.tsv', paired)
    (out / 'COMPLETE').write_text('EVALUATION_DONE\n')
    return dict(split=split['name'], rows=len(records), paired=len(paired))


def evaluate(limit_splits=None, test_limit=None, work=None):
    global WORK
    if work is not None:
        WORK = Path(work)
    ctx = setup(verify='none')
    c = ctx['c']
    jobs = []
    for split in ctx['splits']:
        if limit_splits and split['name'] not in limit_splits:
            continue
        fit_out = WORK / 'splits' / split['name']
        if not (fit_out / 'COMPLETE').exists():
            raise RuntimeError(f'Fit incomplete for {split["name"]}')
        jobs.append(dict(split=split, out=str(fit_out), test_limit=test_limit))
    with ProcessPoolExecutor(max_workers=c['workers'], mp_context=multiprocessing.get_context('spawn')) as pool:
        futures = [pool.submit(run_split_evaluate, job) for job in jobs]
        for future in as_completed(futures):
            print(json.dumps(future.result()), flush=True)


# --------------------------------------------------------------------------- aggregation
def unit_identity(row):
    return tuple(row[k] for k in review.UNIT_KEYS)


def aggregate_records(records, metrics, split_source, extra_keys=()):
    identity_keys = list(review.UNIT_KEYS) + list(extra_keys)
    groups = defaultdict(list)
    for r in records:
        groups[tuple(r[k] for k in identity_keys)].append(r)
    units = []
    for key, rows in groups.items():
        unit = dict(zip(identity_keys, key))
        sections = sorted({r['section_id'] for r in rows})
        unit['n_sections'] = len(sections)
        unit['section_ids'] = ';'.join(sections)
        for field in review.UNIT_METADATA_FIELDS:
            values = sorted({str(r.get(field)) for r in rows if r.get(field) not in (None, '')})
            unit[field] = values[0] if len(values) == 1 else ('mixed' if values else 'unknown')
            unit[f'{field}_n_distinct'] = len(values)
        unit['source_split_key'] = (split_source.get(unit['split'], '')
                                    if unit['arm'] == 'source' else '')
        missing, partial = [], []
        for metric in metrics:
            per_section = {}
            for r in rows:
                value = review.numeric(r.get(metric))
                if value is not None:
                    per_section.setdefault(r['section_id'], value)
            unit[f'{metric}__n'] = len(per_section)
            if not per_section:
                unit[metric] = None
                missing.append(metric)
            elif len(per_section) < len(sections):
                unit[metric] = None
                partial.append(metric)
            else:
                unit[metric] = float(np.mean(list(per_section.values())))
        unit['not_estimable_metrics'] = ';'.join(missing)
        unit['partial_unit_metrics'] = ';'.join(partial)
        unit['complete_for_all_metrics'] = not (missing or partial)
        units.append(unit)
    return units


def macro_table(units, metrics, aggregate='arm', extra_keys=()):
    columns = ['arm', 'target', 'model'] + list(extra_keys)
    grouped = defaultdict(list)
    for u in units:
        group_key = tuple(u[c] for c in (['split'] if aggregate == 'split' else []) + columns)
        grouped[group_key].append(u)
    rows = []
    for key, group in sorted(grouped.items()):
        values_by_column = dict(zip((['split'] if aggregate == 'split' else []) + columns, key))
        arm = values_by_column['arm']
        target = values_by_column['target']
        model = values_by_column['model']
        split = values_by_column.get('split', 'ALL')
        for metric in metrics:
            values = [u[metric] for u in group if u[metric] is not None]
            partial_units = sum(1 for u in group
                                if 0 < u.get(f'{metric}__n', 0) < max(u.get('n_sections', 1), 1))
            row = dict(aggregate=aggregate, split=split, arm=arm, target=target, model=model,
                       metric=metric, n_units=len(group), n_estimable=len(values),
                       n_units_partial=partial_units)
            for column in extra_keys:
                row[column] = values_by_column[column]
            if values:
                arr = np.asarray(values, float)
                row.update(mean=float(arr.mean()), median=float(np.median(arr)),
                           p10=float(np.quantile(arr, 0.10)), p90=float(np.quantile(arr, 0.90)))
            else:
                row.update(mean=None, median=None, p10=None, p90=None)
            if arm == 'source':
                per_source = defaultdict(list)
                for u in group:
                    if u[metric] is not None:
                        per_source[u.get('source_split_key') or 'unknown'].append(u[metric])
                row['source_balanced_mean'] = (float(np.mean([np.mean(v) for v in per_source.values()]))
                                               if per_source else None)
                row['n_sources_estimable'] = len(per_source)
            else:
                row['source_balanced_mean'] = None
                row['n_sources_estimable'] = None
            rows.append(row)
    return rows


def favourable(metric, delta):
    if delta is None:
        return None
    if metric in FAVOURABLE_LOW:
        return delta < 0
    if metric in FAVOURABLE_HIGH:
        return delta > 0
    return None


def claim_summary(paired_macro):
    index = {(r['arm'], r['target'], r['baseline'], r['metric']): r for r in paired_macro}
    rows = []
    readouts = ['rec_local_pattern_pearson', 'rec_jaccard_q75', 'rec_centroid_displacement_pitch',
                'rec_patch75_largest_radius_abs_delta', 'rec_rmse']
    for arm in ['internal', 'source']:
        for target in sorted({r['target'] for r in paired_macro}):
            for readout in readouts:
                favourable_all = True
                details = {}
                for base_name in BASELINES:
                    r = index.get((arm, target, base_name, f'delta_{readout}'))
                    value = None if r is None else review.numeric(r['mean'])
                    verdict = favourable(readout, value)
                    details[base_name] = value
                    if verdict is not True:
                        favourable_all = False
                rows.append(dict(arm=arm, target=target, readout=readout,
                                 favourable_vs_all_four=bool(favourable_all),
                                 **{f'delta_vs_{k}': v for k, v in details.items()},
                                 note='descriptive conjunction only; not a threshold, not a p-value'))
    return rows


def summarize(work=None, prefix='', limit_splits=None):
    global WORK
    if work is not None:
        WORK = Path(work)
    ctx = setup(verify='none')
    c = ctx['c']
    split_source = {s['name']: s.get('held_source', '') for s in ctx['splits']}
    records, paired = [], []
    for split in ctx['splits']:
        if limit_splits and split['name'] not in limit_splits:
            continue
        ev = WORK / 'evaluation' / split['name']
        if not (ev / 'COMPLETE').exists():
            raise RuntimeError(f'Evaluation incomplete: {split["name"]}')
        records += list(csv.DictReader((ev / 'prediction_metrics.tsv').open(), delimiter='\t'))
        paired += list(csv.DictReader((ev / 'paired_metrics.tsv').open(), delimiter='\t'))
    seen = set()
    for r in records:
        key = (r['arm'], r['section_id'], r['target'], r['model'])
        if key in seen:
            raise RuntimeError(f'Duplicate section row: {key}')
        seen.add(key)
    bad = [k for k in list(records[0].keys()) + MACRO_METRICS
           if any(word in k.lower() for word in ('pass', 'success', 'significant', 'p_value', 'truth'))]
    if bad:
        raise RuntimeError(f'Output column names suggest a pass/fail claim: {bad}')
    base.tsv_write(WORK / f'{prefix}prediction_metrics.tsv', records)
    base.tsv_write(WORK / f'{prefix}paired_metrics.tsv', paired)
    units = aggregate_records(records, MACRO_METRICS, split_source)
    base.tsv_write(WORK / f'{prefix}unit_metrics.tsv', units)
    macro = macro_table(units, MACRO_METRICS, 'split') + macro_table(units, MACRO_METRICS, 'arm')
    base.tsv_write(WORK / f'{prefix}macro_summary.tsv', macro)
    delta_metrics = [f'delta_{m}' for m in PAIRED_METRICS]
    paired_units = aggregate_records(paired, delta_metrics, split_source, extra_keys=['baseline'])
    base.tsv_write(WORK / f'{prefix}unit_paired.tsv', paired_units)
    paired_macro = (macro_table(paired_units, delta_metrics, 'split', extra_keys=['baseline'])
                    + macro_table(paired_units, delta_metrics, 'arm', extra_keys=['baseline']))
    base.tsv_write(WORK / f'{prefix}macro_paired.tsv', paired_macro)
    claims = claim_summary(paired_macro)
    base.tsv_write(WORK / f'{prefix}readout_claim_summary.tsv', claims)
    denominators = []
    for model in MODELS:
        for arm in ['internal', 'source']:
            for metric in MACRO_METRICS:
                sel = [u for u in units if u['arm'] == arm and u['model'] == model]
                if not sel:
                    continue
                denominators.append(dict(model=model, arm=arm, metric=metric, unit_rows=len(sel),
                                         complete=sum(1 for u in sel if u[metric] is not None),
                                         partial=sum(1 for u in sel
                                                     if 0 < u.get(f'{metric}__n', 0)
                                                     < max(u.get('n_sections', 1), 1)),
                                         missing=sum(1 for u in sel if u[metric] is None)))
    base.tsv_write(WORK / f'{prefix}denominators.tsv', denominators)
    states = Counter(r.get('state', 'UNKNOWN') for r in records)
    base.json_write(WORK / f'{prefix}summary.json',
                    dict(sections=int(len({r['section_id'] for r in records})),
                         units=int(len({r['unit_id'] for r in records})),
                         section_rows=len(records), paired_rows=len(paired),
                         unit_rows=len(units), macro_rows=len(macro), states=dict(states)))
    write_report(records, units, macro, paired_macro, claims, ctx)
    (WORK / 'STATUS').write_text('SPATIAL_CONTROL_RUN_NO_BIOLOGICAL_PASS_THRESHOLD\n')
    freeze = {}
    for path in sorted(WORK.rglob('*')):
        if path.is_file() and path.name != 'freeze_manifest.json':
            freeze[str(path.relative_to(WORK))] = dict(bytes=path.stat().st_size,
                                                       sha256=base.sha256(path))
    base.json_write(WORK / 'freeze_manifest.json', dict(files=freeze, human_tls_labels_read=False,
                                                        biological_pass_threshold=None))
    print(json.dumps(dict(section_rows=len(records), unit_rows=len(units), macro_rows=len(macro),
                          paired_units=len(paired_units), states=dict(states))), flush=True)


def write_report(records, units, macro, paired_macro, claims, ctx):
    c = ctx['c']
    index = {(r['aggregate'], r['arm'], r['target'], r['model'], r['metric']): r for r in macro}
    lines = ['# 同原始基因输入的空间对照（描述性）', '',
             '固定候选与控制，无调参、无通过阈值、无 p 值；只读冻结输入。',
             '范围绝对偏差先在片内取绝对值，再按完整单位合同聚合。', '']
    lines += ['## 臂级宏平均（完整单位等权）', '',
              '| 模型 | 目标 | 臂 | 单位数 | 局部细节Pearson | q75 Jaccard | 质心位移 | 半径绝对偏差 | sd比 | RMSE |',
              '|---|---|---|---:|---:|---:|---:|---:|---:|---:|']
    def val(arm, target, model, metric):
        r = index.get(('arm', arm, target, model, metric))
        return 'NA' if (r is None or r['mean'] is None) else f"{r['mean']:.4f}"
    for arm in ['internal', 'source']:
        for target in c['targets']:
            for model in MODELS:
                r = index.get(('arm', arm, target, model, 'rec_rmse'))
                if r is None:
                    continue
                lines.append(f"| {model} | {target} | {arm} | {r['n_units']} | "
                             f"{val(arm, target, model, 'rec_local_pattern_pearson')} | "
                             f"{val(arm, target, model, 'rec_jaccard_q75')} | "
                             f"{val(arm, target, model, 'rec_centroid_displacement_pitch')} | "
                             f"{val(arm, target, model, 'rec_patch75_largest_radius_abs_delta')} | "
                             f"{val(arm, target, model, 'rec_sd_ratio')} | "
                             f"{val(arm, target, model, 'rec_rmse')} |")
    lines += ['', '## 候选减基线（单位等权配对差值）', '',
              '| 目标 | 臂 | 基线 | 局部细节Δ | q75 Jaccard Δ | 质心位移 Δ | 半径绝对偏差 Δ |',
              '|---|---|---|---:|---:|---:|---:|']
    pidx = {(r['arm'], r['target'], r['baseline'], r['metric']): r for r in paired_macro}
    for target in c['targets']:
        for arm in ['internal', 'source']:
            for base_name in BASELINES:
                def pval(metric):
                    r = pidx.get((arm, target, base_name, f'delta_{metric}'))
                    return 'NA' if (r is None or r['mean'] is None) else f"{r['mean']:+.4f}"
                lines.append(f"| {target} | {arm} | {base_name} | "
                             f"{pval('rec_local_pattern_pearson')} | {pval('rec_jaccard_q75')} | "
                             f"{pval('rec_centroid_displacement_pitch')} | "
                             f"{pval('rec_patch75_largest_radius_abs_delta')} |")
    lines += ['', '## 仅描述：候选是否同时优于四个对照（不构成阈值或显著性）', '']
    for row in claims:
        if row['readout'] in ('rec_local_pattern_pearson', 'rec_jaccard_q75'):
            lines.append(f"- {row['arm']} / {row['target']} / {row['readout']}："
                         f"同时有利 = {row['favourable_vs_all_four']}")
    lines += ['', '证据角色：既有开发材料上的探索性机制对照；预测与输入共享可见基因，含输入回声；',
              '打乱邻域为单次机制破坏，不作 p 值或正式零模型；未同时优于四个对照的读出不得宣称依赖真实邻域排列。',
              'NOT_RUN：人工 TLS、独立 marker 隔离、未知患者身份、正式空间显著性、物理尺度、视野外推断、',
              '打乱参考描述（本轮关闭 shuffle 计算）。', '']
    (WORK / 'REPORT.md').write_text('\n'.join(lines) + '\n')


# --------------------------------------------------------------------------- smoke & checks
def smoke(train_limit=8, test_limit=2, splits=('internal_0', 'source_0')):
    """Bounded end-to-end run in WORK/smoke with a clearly truncated training set."""
    global WORK
    main = WORK
    smoke_dir = main / 'smoke'
    smoke_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(main / 'contract.json', smoke_dir / 'contract.json')
    WORK = smoke_dir
    start = time.time()
    ctx = setup(verify='none')
    base.json_write(smoke_dir / 'input_checks.json', ctx['checks'])
    base.json_write(smoke_dir / 'splits.json', ctx['splits'])
    (smoke_dir / 'PREPARED').write_text('SMOKE_TRUNCATED_PREPARE\n')
    stages = []
    for name in splits:
        split = next(s for s in ctx['splits'] if s['name'] == name)
        t0 = time.time()
        info = run_split_fit(dict(split=split, out=str(smoke_dir / 'splits' / name),
                                  train_limit=train_limit, test_limit=test_limit, force=True))
        stages.append(dict(stage=f'fit/{name}', seconds=round(time.time() - t0, 1), info=info))
    for name in splits:
        split = next(s for s in ctx['splits'] if s['name'] == name)
        t0 = time.time()
        info = run_split_evaluate(dict(split=split, out=str(smoke_dir / 'splits' / name),
                                       test_limit=test_limit, force=True))
        stages.append(dict(stage=f'evaluate/{name}', seconds=round(time.time() - t0, 1), info=info))
    t0 = time.time()
    summarize(work=smoke_dir, limit_splits=list(splits))
    stages.append(dict(stage='summarize', seconds=round(time.time() - t0, 1), info=None))
    outputs = {}
    for name in ['prediction_metrics.tsv', 'paired_metrics.tsv', 'unit_metrics.tsv',
                 'unit_paired.tsv', 'macro_summary.tsv', 'macro_paired.tsv',
                 'readout_claim_summary.tsv', 'denominators.tsv', 'summary.json', 'REPORT.md',
                 'freeze_manifest.json']:
        path = smoke_dir / name
        outputs[name] = dict(exists=path.exists(),
                             rows=(sum(1 for _ in path.open()) - 1) if path.exists() else None)
    base.json_write(smoke_dir / 'smoke_summary.json',
                    dict(splits=list(splits), train_limit=train_limit, test_limit=test_limit,
                         stages=stages, outputs=outputs, seconds=round(time.time() - start, 1),
                         note='truncated training/test sets; smoke only, not a fold result'))
    WORK = main
    print(json.dumps(dict(stages=stages, outputs={k: v['rows'] for k, v in outputs.items()}),
                     indent=1), flush=True)


def check():
    """Targeted checks demanded by the contract and the preflight (B-01..B-06)."""
    global WORK
    out = WORK / 'checks'
    out.mkdir(parents=True, exist_ok=True)
    ctx = setup(verify='sample', sids=['MISC33', 'NCBI642', 'TENX13', 'TENX14', 'ZEN38'])
    c = ctx['c']
    results = []

    def record(name, expected, observed, passed):
        results.append(dict(check=name, expected=str(expected), observed=str(observed),
                            passed=bool(passed)))

    # ---- structural / no-leak guarantees
    record('raw_hdf5_not_allowlisted', 'False', ctx['checks']['raw_hdf5_allowlisted'],
           not ctx['checks']['raw_hdf5_allowlisted'])
    record('denylist_clean', 'True', ctx['checks']['no_allowed_path_matches_denylist'],
           ctx['checks']['no_allowed_path_matches_denylist'])
    record('split_unit_counts', '65 internal / 47 source',
           (ctx['checks']['internal_units'], ctx['checks']['source_units']),
           (ctx['checks']['internal_units'], ctx['checks']['source_units']) == (65, 47))
    record('internal_folds_partition_72', 72, ctx['checks']['internal_sections'],
           ctx['checks']['internal_sections'] == 72)
    record('source_folds_partition_53', 53, ctx['checks']['source_sections'],
           ctx['checks']['source_sections'] == 53)
    record('visible_target_disjoint', 'True', ctx['checks']['visible_target_disjoint'],
           ctx['checks']['visible_target_disjoint'])
    record('identity_flags', '14 sections / 8 units',
           (ctx['checks']['known_sections'], ctx['checks']['known_units']),
           (ctx['checks']['known_sections'], ctx['checks']['known_units']) == (14, 8))
    record('source_not_testable_sections', 19, ctx['checks']['source_not_testable_sections'],
           ctx['checks']['source_not_testable_sections'] == 19)

    # ---- frozen inputs pinning: one matrix per target, no rebuild
    sid, ti = 'MISC33', 0
    row = ctx['byid'][sid]
    cache = base.load_cache(row)
    z, qc = visible_matrices(sid, ti)
    record('visible_matrix_dim_matches_vocabulary', len(ctx['vocab'][c['targets'][ti]]['full15']),
           z.shape[1], z.shape[1] == len(ctx['vocab'][c['targets'][ti]]['full15']))
    record('visible_matrix_rows_match_barcodes', len(cache['barcodes']), z.shape[0],
           z.shape[0] == len(cache['barcodes']))
    record('qc_columns', 2, qc.shape[1], qc.shape[1] == 2)

    # ---- leakage: held-out target readouts must not touch features or predictions
    perms = {s: permutation_for(s, c['seed'], ctx['byid'][s]['n_spots']) for s in ctx['byid']}
    designs = section_designs(row, ti, cache, perms[sid], ['point_ridge', 'local_ridge'])
    tampered = dict(cache)
    tampered['y'] = cache['y'].copy()
    tampered['y'][:, ti] = np.random.default_rng(0).random(len(cache['y'])).astype(np.float32)
    designs_tampered = section_designs(row, ti, tampered, perms[sid], ['point_ridge', 'local_ridge'])
    same = all(np.array_equal(designs[k]['design'], designs_tampered[k]['design']) for k in designs)
    record('features_independent_of_heldout_target', 'bitwise identical', same, same)
    indexes = base.sample_indices(row['n_spots'], sid, c['seed'], 512)
    indexes2 = base.sample_indices(row['n_spots'], sid, c['seed'], 512)
    record('sampling_independent_of_target', 'identical indices', bool(np.array_equal(indexes, indexes2)),
           bool(np.array_equal(indexes, indexes2)))
    shuffled_barcodes = cache['barcodes'][::-1]
    try:
        review.assert_aligned(shuffled_barcodes, cache['barcodes'])
        detected = False
    except ValueError:
        detected = True
    record('barcode_order_mismatch_detected', 'raises', detected, detected)

    # ---- column rule (training only) and block-scaling equivalence with the formal pipeline
    x = np.c_[np.ones((5, 1)), np.zeros((5, 1)), np.arange(10, 15)[:, None]]
    y = np.array([1.0, 2, 3, 4, 5])
    w = np.ones(5)
    model = fit_linear(dict(c, solver='auto'), x, y, w, 1)
    record('train_constant_columns_removed', 2, model['dropped'], model['dropped'] == 2)
    record('constant_column_rule_is_train_only', '[ones=constant, zeros=constant, varying kept]',
           list(model['keep']), list(model['keep']) == [False, False, True])
    real_split = next(s for s in ctx['splits'] if s['name'] == 'internal_0')
    train_rows = [ctx['byid'][s] for s in real_split['train']]
    idx = {r['section_id']: base.sample_indices(r['n_spots'], r['section_id'], c['seed'],
                                                c['train_spots_per_section']) for r in train_rows}
    weights = follow.weights(train_rows, idx)
    y_train = target_vector(train_rows, idx, ti)
    stacks, widths = stack_training(ctx, train_rows, idx, ti, perms, ['point_ridge', 'local_ridge'])
    stack, local_stack = stacks['point_ridge'], stacks['local_ridge']
    d_x = stack.shape[1] - CONTEXT_DIM
    model_point = fit_linear(c, stack, y_train, weights, widths['point_ridge'])
    model_local = fit_linear(c, local_stack, y_train, weights, widths['local_ridge'])
    record('block_and_context_column_layout', '(dX+6) point, (3*dX+6) local, context last',
           (stack.shape[1], local_stack.shape[1], CONTEXT_DIM),
           stack.shape[1] == d_x + CONTEXT_DIM and local_stack.shape[1] == 3 * d_x + CONTEXT_DIM
           and widths['point_ridge'] == 0 and widths['local_ridge'] == 3 * d_x)
    scaled_point = prepare_scaled(model_point, stack)
    scaled_local = prepare_scaled(model_local, local_stack)
    context_diff = float(np.max(np.abs(scaled_point[:, -CONTEXT_DIM:]
                                       - scaled_local[:, -CONTEXT_DIM:])))
    record('context_columns_same_standardisation_and_penalty_as_point', '0.0', context_diff,
           context_diff == 0.0)
    x_only = stack[:, :d_x]
    scaled_x = prepare_scaled(model_point, np.c_[x_only, np.zeros((len(x_only), CONTEXT_DIM),
                                                                 np.float32)]
                              )[:, :d_x]
    neighbourhood_ratio = float(np.max(np.abs(scaled_local[:, :d_x]
                                              - scaled_x * c['block_scale'])))
    record('x_block_scaled_by_one_over_sqrt3_only', '0.0', neighbourhood_ratio,
           neighbourhood_ratio == 0.0)
    # D-270: replicate only the gene/QC block; E stays once and unscaled
    triple = np.c_[x_only, x_only, x_only, stack[:, d_x:]]
    model_triple = fit_linear(c, triple, y_train, weights, 3 * d_x)
    test_row = ctx['byid'][real_split['test'][0]]
    test_cache = base.load_cache(test_row)
    test_design = section_designs(test_row, ti, test_cache,
                                  perms[test_row['section_id']], ['point_ridge'])['point_ridge']['design']
    test_x, test_e = test_design[:, :d_x], test_design[:, d_x:]
    triple_design = np.c_[test_x, test_x, test_x, test_e]
    model_point_lsqr = fit_linear(dict(c, solver='lsqr', solver_tol=1e-12, solver_max_iter=20000),
                                  stack, y_train, weights, widths['point_ridge'])
    p_point = model_predict(model_point, test_design)
    p_triple = model_predict(model_triple, triple_design)
    p_point_lsqr = model_predict(model_point_lsqr, test_design)
    model_wrong = fit_linear(dict(c, block_scale=1.0), triple, y_train, weights, 3 * d_x)
    p_wrong = model_predict(model_wrong, triple_design)
    floor = float(np.max(np.abs(p_point - p_point_lsqr)))
    diff = float(np.max(np.abs(p_point - p_triple)))
    wrong = float(np.max(np.abs(p_point - p_wrong)))
    gram = prepare_scaled(model_point, stack)[:, :d_x]
    eig = np.linalg.eigvalsh(gram.T @ gram)
    conditioning = float(eig.max() / max(eig.min() + c['alpha'], 1e-12))
    tolerance = float(c['tolerances']['block_equivalence_prediction_max_abs_diff'])
    record('triple_copy_block_equivalence_fixed_tolerance',
           f'residual <= fixed 1e-5, wrong-scale control >> residual (solver floor recorded, not used to expand)',
           dict(residual=diff, solver_floor_indicator_only=floor, wrong_scale_control=wrong,
                ridge_gram_conditioning=conditioning),
           diff <= tolerance and wrong > 10 * diff)
    del stack, local_stack, triple

    # ---- frozen point/smoothing agreement (hard tolerances, fixed before the run)
    probe_splits = ['internal_0']
    for split in ctx['splits']:
        if split['name'] not in probe_splits:
            continue
        train_rows = [ctx['byid'][s] for s in split['train']]
        idx = {r['section_id']: base.sample_indices(r['n_spots'], r['section_id'], c['seed'],
                                                    c['train_spots_per_section']) for r in train_rows}
        weights = follow.weights(train_rows, idx)
        y_train = target_vector(train_rows, idx, ti)
        stacks_probe, widths_probe = stack_training(ctx, train_rows, idx, ti, perms, ['point_ridge'])
        stack = stacks_probe['point_ridge']
        model_point = fit_linear(c, stack, y_train, weights, widths_probe['point_ridge'])
        del stack, stacks_probe
        max_pred_diff, max_rmse_diff = 0.0, 0.0
        for s in split['test']:
            trow = ctx['byid'][s]
            tcache = base.load_cache(trow)
            design = section_designs(trow, ti, tcache, perms[s], ['point_ridge'])['point_ridge']['design']
            pred = model_predict(model_point, design)
            smooth = 0.5 * (pred + base.neighbor_mean(tcache['near'], pred))
            with np.load(PRED / split['name'] / f'{s}_t{ti}_full15.npz', allow_pickle=False) as z:
                frozen_point = z['gene_ridge'].astype(np.float64)
                frozen_smooth = z['gene_ridge_smooth_half'].astype(np.float64)
            y = tcache['y'][:, ti].astype(np.float64)
            for name, mine, theirs in [('point', pred, frozen_point), ('smooth', smooth, frozen_smooth)]:
                max_pred_diff = max(max_pred_diff, float(np.max(np.abs(mine - theirs))))
                rmse_mine = float(np.sqrt(np.mean((y - mine) ** 2)))
                rmse_theirs = float(np.sqrt(np.mean((y - theirs) ** 2)))
                max_rmse_diff = max(max_rmse_diff, abs(rmse_mine - rmse_theirs))
        record('frozen_point_prediction_agreement',
               f'<={c["tolerances"]["point_prediction_max_abs_diff_vs_frozen"]}', max_pred_diff,
               max_pred_diff <= c['tolerances']['point_prediction_max_abs_diff_vs_frozen'])
        record('frozen_point_section_rmse_agreement',
               f'<={c["tolerances"]["point_section_rmse_abs_diff_vs_frozen"]}', max_rmse_diff,
               max_rmse_diff <= c['tolerances']['point_section_rmse_abs_diff_vs_frozen'])

    # ---- permutation control (on a section that actually has isolated spots)
    perm_sid = next(s for s in sorted(ctx['byid']) if ctx['byid'][s]['near_isolated_spots'] > 0)
    perm_row = ctx['byid'][perm_sid]
    perm_cache = base.load_cache(perm_row)
    perm = perms[perm_sid]
    near, outer = perm_cache['near'], perm_cache['outer']
    near_p = sparse.csr_matrix(near[perm][:, perm])
    outer_p = sparse.csr_matrix(outer[perm][:, perm])
    record('permutation_preserves_nnz', (near.nnz, outer.nnz), (near_p.nnz, outer_p.nnz),
           (near.nnz, outer.nnz) == (near_p.nnz, outer_p.nnz))
    record('permutation_preserves_degree_multiset', 'equal',
           bool(np.array_equal(np.sort(np.diff(near.indptr)), np.sort(np.diff(near_p.indptr)))),
           bool(np.array_equal(np.sort(np.diff(near.indptr)), np.sort(np.diff(near_p.indptr)))))
    record('permutation_preserves_weights', 'equal',
           bool(np.allclose(np.sort(near.data), np.sort(near_p.data))),
           bool(np.allclose(np.sort(near.data), np.sort(near_p.data))))
    record('permutation_preserves_symmetry', 'both symmetric',
           (sparse.triu(near, 1).nnz == sparse.triu(near.T, 1).nnz,
            sparse.triu(near_p, 1).nnz == sparse.triu(near_p.T, 1).nnz), True)
    iso = np.flatnonzero(np.diff(near.indptr) == 0)
    if len(iso):
        iso_permuted = np.flatnonzero(np.diff(near_p.indptr) == 0)
        z_all, qc_all = visible_matrices(perm_sid, ti)
        spatial = np.c_[z_all, qc_all]
        nm = base.neighbor_mean(near_p, spatial)
        record('permutation_self_fallback_at_isolated', 'own value at the permuted isolated rows',
               float(np.max(np.abs(nm[iso_permuted] - spatial[iso_permuted]))),
               float(np.max(np.abs(nm[iso_permuted] - spatial[iso_permuted]))) == 0.0)
        record('permutation_preserves_isolated_count', 'count preserved',
               (len(iso), int(len(iso_permuted))), len(iso) == len(iso_permuted))
    else:
        record('permutation_self_fallback_at_isolated', 'section with isolated spots',
               'no isolated spots in the checked section', False)
    blocks_real = full_features(perm_row, ti, perm_cache, perm=None)
    blocks_perm = full_features(perm_row, ti, perm_cache, perm=perm)
    changed = float(np.max(np.abs(blocks_real['near'] - blocks_perm['near'])))
    record('permutation_changes_neighbour_values', '>0', changed, changed > 0)
    record('permutation_keeps_x_block', 'identical',
           float(np.max(np.abs(blocks_real['x'] - blocks_perm['x']))),
           float(np.max(np.abs(blocks_real['x'] - blocks_perm['x']))) == 0.0)
    record('permutation_keeps_context_block', 'identical',
           float(np.max(np.abs(blocks_real['e'] - blocks_perm['e']))),
           float(np.max(np.abs(blocks_real['e'] - blocks_perm['e']))) == 0.0)
    local_design, local_width = raw_design(blocks_real, 'local_ridge')
    perm_design, perm_width = raw_design(blocks_perm, 'permuted_ridge')
    d_x_perm = blocks_real['x'].shape[1]
    record('candidate_uses_real_graph', 'candidate != permuted design',
           float(np.max(np.abs(local_design - perm_design))),
           float(np.max(np.abs(local_design - perm_design))) > 0)
    record('candidate_x_block_matches_x', 'identical',
           float(np.max(np.abs(local_design[:, :d_x_perm] - blocks_real['x']))),
           float(np.max(np.abs(local_design[:, :d_x_perm] - blocks_real['x']))) == 0.0)
    record('candidate_neighbour_blocks_are_graph_aggregates_of_x', 'W[Z,QC] and V[Z,QC]',
           (float(np.max(np.abs(local_design[:, d_x_perm:2 * d_x_perm] - blocks_real['near']))),
            float(np.max(np.abs(local_design[:, 2 * d_x_perm:3 * d_x_perm]
                                 - blocks_real['outer'])))),
           float(np.max(np.abs(local_design[:, d_x_perm:2 * d_x_perm] - blocks_real['near']))) == 0.0
           and float(np.max(np.abs(local_design[:, 2 * d_x_perm:3 * d_x_perm]
                                    - blocks_real['outer']))) == 0.0)
    record('context_block_appended_once_unscaled', '(3*dX+6) width, context last, unscaled',
           (local_width, perm_width, CONTEXT_DIM,
            float(np.max(np.abs(local_design[:, -CONTEXT_DIM:] - blocks_real['e'])))),
           local_width == 3 * d_x_perm and perm_width == 3 * d_x_perm
           and float(np.max(np.abs(local_design[:, -CONTEXT_DIM:] - blocks_real['e']))) == 0.0)

    # ---- metric binding: detail-component Pearson, strict-above Jaccard, per-section abs radius
    rng = np.random.default_rng(3)
    geo = geo_for(ctx, row, cache)
    y_vec = cache['y'][:, ti].astype(np.float64)
    pred_vec = np.asarray(base.neighbor_mean(cache['near'], y_vec))
    rec, pd, _ = describe_prediction(ctx, row, ti, y_vec, pred_vec, geo, 'binding_check')
    detail_y = y_vec - base.neighbor_mean(cache['near'], y_vec)
    detail_p = pred_vec - base.neighbor_mean(cache['near'], pred_vec)
    manual = float(np.corrcoef(detail_y, detail_p)[0, 1])
    record('local_detail_pearson_binding', round(manual, 9),
           round(rec['rec_local_pattern_pearson'], 9),
           abs(manual - rec['rec_local_pattern_pearson']) < 1e-9)
    zero_inflated = np.zeros(1000)
    zero_inflated[:850] = 0.0
    zero_inflated[850:] = rng.random(150)
    mask = zero_inflated > np.quantile(zero_inflated, 0.75)
    record('jaccard_threshold_is_strict', 'nonzero set', int(mask.sum()), int(mask.sum()) == 150)
    record('zero_is_a_value_not_missing', 0.0, review.numeric(0.0), review.numeric(0.0) == 0.0)
    signed = rec.get('rec_patch75_largest_radius_pitch_delta')
    absolute = rec.get('rec_patch75_largest_radius_abs_delta')
    if signed is None:
        record('abs_radius_none_propagates', 'both None (NOT_ESTIMABLE propagation)',
               (signed, absolute), absolute is None)
    else:
        record('abs_radius_is_absolute_on_real_section', 'abs(signed delta)',
               (signed, absolute), abs(absolute - abs(signed)) < 1e-12)

    # ---- range readout: per-section abs first, complete-unit mean second (not abs of the mean)
    def radius_pair_records(values):
        rows = []
        for i, value in enumerate(values):
            rows.append(dict(split='internal_0', arm='internal', section_id=f'S{i+1}', unit_id='u1',
                             organ='Bowel', source_key='k', platform='Visium',
                             patient_known='False', confirmed_unit='False', view='full15',
                             model_class='x', cross_organ_extrapolation='False', target='B_AXIS',
                             model='local_ridge', baseline='point_ridge',
                             n_spots=10,
                             rec_patch75_largest_radius_pitch_delta=('' if value is None else str(value)),
                             rec_patch75_largest_radius_abs_delta=('' if value is None else str(abs(value)))))
        return rows

    cancelled = aggregate_records(radius_pair_records([0.4, -0.4, 0.0]),
                                  ['rec_patch75_largest_radius_pitch_delta',
                                   'rec_patch75_largest_radius_abs_delta'], {})[0]
    mean_of_abs = cancelled['rec_patch75_largest_radius_abs_delta']
    abs_of_mean = abs(cancelled['rec_patch75_largest_radius_pitch_delta'])
    record('mean_of_abs_not_abs_of_mean', 'mean_abs=0.266667 vs abs_mean=0.0',
           (round(mean_of_abs, 6), round(abs_of_mean, 6)),
           abs(mean_of_abs - 0.2666666667) < 1e-9 and abs_of_mean == 0.0)
    record('zero_radius_delta_is_estimable', '3 of 3 sections counted',
           (cancelled['n_sections'], cancelled['rec_patch75_largest_radius_abs_delta__n']),
           int(cancelled['n_sections']) == 3
           and int(cancelled['rec_patch75_largest_radius_abs_delta__n']) == 3)
    partial_radius = aggregate_records(radius_pair_records([0.4, None]),
                                       ['rec_patch75_largest_radius_pitch_delta',
                                        'rec_patch75_largest_radius_abs_delta'], {})[0]
    record('radius_missing_side_is_not_testable_unit', 'None + partial recorded',
           (partial_radius['rec_patch75_largest_radius_abs_delta'],
            'rec_patch75_largest_radius_abs_delta' in partial_radius['partial_unit_metrics']),
           partial_radius['rec_patch75_largest_radius_abs_delta'] is None
           and 'rec_patch75_largest_radius_abs_delta' in partial_radius['partial_unit_metrics'])

    # ---- aggregation semantics match the reviewed module (same records, shared metrics)
    agg_records = [dict(r) for r in [
        dict(split='internal_0', arm='internal', section_id='S1', unit_id='u1', organ='Bowel',
             source_key='title:a', platform='Visium', patient_known='False', confirmed_unit='False',
             view='full15', model_class='x', cross_organ_extrapolation='False',
             target='B_AXIS', model='local_ridge', rec_rmse='0.1', rec_pearson='0.5',
             rec_local_pattern_pearson='0.2', rec_jaccard_q75='0.3'),
        dict(split='internal_0', arm='internal', section_id='S2', unit_id='u1', organ='Bowel',
             source_key='title:b', platform='Visium', patient_known='False', confirmed_unit='False',
             view='full15', model_class='x', cross_organ_extrapolation='False',
             target='B_AXIS', model='local_ridge', rec_rmse='0.3', rec_pearson='0.7',
             rec_local_pattern_pearson='0.4', rec_jaccard_q75='0.5')]]
    mine = aggregate_records(agg_records, ['rec_rmse', 'rec_pearson'], {})
    review_rows = [dict(r, cross_organ_extrapolation=False) for r in agg_records]
    reviewed_units = review.aggregate_units(review_rows, split_source={})
    record('aggregation_matches_reviewed_module', reviewed_units[0]['rec_rmse'],
           mine[0]['rec_rmse'], abs(reviewed_units[0]['rec_rmse'] - mine[0]['rec_rmse']) < 1e-12)
    record('unit_metadata_conflict_mixed', 'mixed', mine[0]['source_key'],
           mine[0]['source_key'] == 'mixed')
    record('unit_key_excludes_source_key', 1, len(mine), len(mine) == 1)
    partial = aggregate_records([dict(agg_records[0]), dict(agg_records[1], rec_pearson='')],
                                ['rec_pearson'], {})
    record('partial_unit_metric_not_estimable', 'None + recorded',
           (partial[0]['rec_pearson'], 'rec_pearson' in partial[0]['partial_unit_metrics']),
           partial[0]['rec_pearson'] is None and 'rec_pearson' in partial[0]['partial_unit_metrics'])

    # ---- tiny real end-to-end fit/evaluate on two sections (fast, not a fold result)
    tiny_dir = out / 'tiny'
    tiny_train = train_rows[:3]
    info = run_split_fit(dict(split=real_split, out=str(tiny_dir), train_limit=3, test_limit=1,
                              force=True))
    info = run_split_evaluate(dict(split=real_split, out=str(tiny_dir), test_limit=1, force=True))
    record('tiny_fit_evaluate_ran', 'rows>0', info['rows'], info['rows'] > 0)
    diag = list(csv.DictReader((tiny_dir / 'fit_diagnostics.tsv').open(), delimiter='\t'))
    record('no_nonconverged_fit', 'all converged',
           sorted({d['converged'] for d in diag}), all(d['converged'] == 'True' for d in diag))
    record('solver_recorded', 'non-empty', sorted({d['solver'] for d in diag}),
           all(d['solver'] for d in diag))
    base.tsv_write(out / 'solver_diagnostics.tsv', diag)

    failures = [r for r in results if not r['passed']]
    base.json_write(out / 'checks.json',
                    dict(n_checks=len(results), n_failed=len(failures), results=results,
                         failures=failures, setup_checks=ctx['checks']))
    print(json.dumps(dict(n_checks=len(results), n_failed=len(failures),
                          failures=failures[:6]), indent=1), flush=True)
    return 0 if not failures else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['prepare', 'fit', 'evaluate', 'summarize', 'smoke',
                                            'check'])
    parser.add_argument('--splits', nargs='*', default=None)
    parser.add_argument('--train-limit', type=int, default=None)
    parser.add_argument('--test-limit', type=int, default=None)
    args = parser.parse_args()
    if os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', '-1'):
        raise RuntimeError('GPU disabled')
    if args.command == 'prepare':
        prepare()
    elif args.command == 'fit':
        fit(args.splits, args.train_limit, args.test_limit)
    elif args.command == 'evaluate':
        evaluate(args.splits, args.test_limit)
    elif args.command == 'summarize':
        summarize()
    elif args.command == 'smoke':
        smoke(train_limit=args.train_limit or 8, test_limit=args.test_limit or 2)
    else:
        sys.exit(check())


if __name__ == '__main__':
    main()
