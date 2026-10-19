#!/usr/bin/env python3
"""Two-component candidate-field probe (D-279; contract docs/FIELD_CANDIDATE_SUBFIELD_PROBE_20261019.md).

Frozen rules implemented here (release: review/preflight.json FIT_RELEASED):
  * A/B split by sorted gene names and literal numpy seed 0; A = ceil(2n/3), remainder to B
  * every A and B gene map uses log1p(1e4 * count / L_external), L_external = point total counts
    minus ALL panel genes of that target (A union B); points with L_external <= 0 are dropped on
    both sides; this denominator is NOT the old proxy whole-library scale
  * a gene is measurable in a training unit only when EVERY section of that unit has >= 20
    non-zero tissue spots and the unit-level training map variance > 0; usable A = measurable in
    at least half of the training units; usable A < 6 makes the fold NOT_TESTABLE
  * training sample: common fixed count n_sample = min(512, min available points over training
    units), equal per unit, two-section units split as evenly as possible, seed 0; test units are
    never subsampled
  * scaling: per-gene positive-value mean over the fold's training sample, applied to test
  * NMF K=1 and K=2 fixed (init=nndsvda, solver=cd, random_state=0, max_iter=500, tol=1e-4)
  * in-fold pairing of resample loadings only (5 fixed 90% row resamples, seeds 0-4); cross-fold
    loadings are descriptive only
  * references share the same fold/A/training sample: A mean, K=1, K=2, single gene (chosen by the
    unscaled training map mean) and all-A Ridge(alpha=10) after training StandardScaler
  * B maps are target values only, evaluated after projection; test B needs >= 20 non-zero spots
    and variance > 0; q75 uses strict greater-than
  * technical covariates exclude the whole target panel (A union B)
Commands: cache | smoke | run | summarize | verify
"""
import argparse
import json
import math
import multiprocessing
import sys
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np
from sklearn.decomposition import NMF
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import field_discovery_20261018 as base                       # noqa: E402
import field_proxy_audit_20261019 as proxy                    # noqa: E402

OUT = ROOT / 'infra/field_candidate_subfield_20261019'
CACHE = OUT / 'cache/sections'
RESULTS = OUT / 'results'
PREFLIGHT = ROOT / 'infra/paseo_candidate_subfield_20261019/review/preflight.json'
SPLITS = ROOT / 'infra/field_source_holdout_20261018/splits.json'
MANIFEST = ROOT / 'infra/field_source_holdout_20261018/input_manifest.json'
PANEL = ROOT / 'infra/tls_multiscale_gobp_20261002/panel_controls.json'
TARGETS = ['B_AXIS', 'HALLMARK_INTERFERON_GAMMA_RESPONSE',
           'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
           'HALLMARK_OXIDATIVE_PHOSPHORYLATION']
K_LIST = (1, 2)
NMF_PARAMS = dict(init='nndsvda', solver='cd', beta_loss='frobenius',
                  random_state=0, max_iter=500, tol=1e-4)
N_SAMPLE = 512
MIN_USABLE_A = 6
RESAMPLE_COUNT = 5
RESAMPLE_FRACTION = 0.9
MIN_SPOTS = 20


def json_read(path):
    return json.loads(Path(path).read_text())


def read_tsv(path):
    import csv
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def json_write(path, payload):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(payload, indent=2, ensure_ascii=False) + '\n')


def split_genes(target, panels=None):
    """D-279 A/B split: names only, literal seed 0, A takes ceil(2n/3), remainder to B."""
    panels = panels or json_read(PANEL)
    genes = sorted(set(panels[target]))
    order = [genes[i] for i in np.random.default_rng(0).permutation(len(genes))]
    n_a = math.ceil(2 * len(genes) / 3)
    return order[:n_a], order[n_a:]


def require_release():
    state = json_read(PREFLIGHT)
    if state.get('status') != 'FIT_RELEASED' or state.get('computation_authorized') is not True:
        raise SystemExit(f'refusing to fit: preflight is {state.get("status")!r} '
                         f'(computation_authorized={state.get("computation_authorized")!r})')


def manifest():
    return {r['section_id']: r for r in json_read(MANIFEST)}


def mode_splits():
    out = {}
    for path in sorted(SPLITS.parent.glob('*/split.json')):
        spec = json_read(path)
        out[spec['name']] = spec
    return out


# ---------------------------------------------------------------- cache
def cache_section(args):
    sid, meta, panels = args
    row = meta[sid]
    path = base.HEST / f'st/{sid}.h5ad'
    with h5py.File(path) as handle:
        genes = base.strings(handle['var/_index'][:])
        n_genes = len(genes)
        keep = np.flatnonzero(handle['obs/in_tissue'][:] == 1)
        matrix = base.read_matrix(handle['X'], len(base.strings(handle['obs/_index'][:])), n_genes)
    matrix = matrix[keep].tocsr()
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    out = dict(sid=sid, n_spots=int(matrix.shape[0]))
    lookup = {g: i for i, g in enumerate(genes)}
    nnz_all = np.diff(matrix.indptr)
    for target in TARGETS:
        A, B = split_genes(target, panels)
        a_cols = np.array([lookup.get(g, -1) for g in A])
        b_cols = np.array([lookup.get(g, -1) for g in B])
        a_present = a_cols >= 0
        b_present = b_cols >= 0
        panel_cols = np.concatenate([a_cols[a_present], b_cols[b_present]]).astype(int)
        panel_sum = np.asarray(matrix[:, panel_cols].sum(axis=1)).ravel() if len(panel_cols) else np.zeros_like(totals)
        outside = totals - panel_sum
        ok = outside > 0
        ok_idx = np.flatnonzero(ok)
        nnz_panel = np.asarray((matrix[:, panel_cols] != 0).sum(axis=1)).ravel() if len(panel_cols) else np.zeros_like(totals)
        detection = np.where(ok, (nnz_all - nnz_panel) / max(n_genes - len(panel_cols), 1), np.nan)
        a_map = np.zeros((len(ok_idx), len(A)), np.float32)
        b_map = np.zeros((len(ok_idx), len(B)), np.float32)
        a_nnz = np.zeros(len(A), np.int32)
        b_nnz = np.zeros(len(B), np.int32)
        if len(ok_idx):
            sub = matrix[ok_idx]
            denom = outside[ok_idx][:, None]
            if a_present.any():
                counts = np.asarray(sub[:, a_cols[a_present]].todense())
                a_map[:, a_present] = np.log1p(1e4 * counts / denom)
                a_nnz[a_present] = (counts > 0).sum(axis=0)
            if b_present.any():
                counts = np.asarray(sub[:, b_cols[b_present]].todense())
                b_map[:, b_present] = np.log1p(1e4 * counts / denom)
                b_nnz[b_present] = (counts > 0).sum(axis=0)
        out[f'A_present_{target}'] = a_present
        out[f'B_present_{target}'] = b_present
        out[f'A_map_{target}'] = a_map
        out[f'B_map_{target}'] = b_map
        out[f'A_nnz_{target}'] = a_nnz
        out[f'B_nnz_{target}'] = b_nnz
        out[f'depth_{target}'] = np.log1p(outside[ok_idx]).astype(np.float32)
        out[f'detection_{target}'] = detection[ok_idx].astype(np.float32)
        out[f'ok_rows_{target}'] = ok_idx.astype(np.int32)
    CACHE.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(CACHE / f'{sid}.npz', **out)
    return dict(sid=sid, n_spots=out['n_spots'])


def load_cache(sid):
    with np.load(CACHE / f'{sid}.npz', allow_pickle=True) as handle:
        return {k: handle[k] for k in handle.files}


def cache(workers=4, sids=None):
    meta = manifest()
    sids = list(sids or sorted(meta))
    CACHE.mkdir(parents=True, exist_ok=True)
    todo = [s for s in sids if not (CACHE / f'{s}.npz').exists()]
    panels = json_read(PANEL)
    if todo:
        with multiprocessing.get_context('spawn').Pool(workers) as pool:
            for i, info in enumerate(pool.imap_unordered(cache_section, [(s, meta, panels) for s in todo])):
                if (i + 1) % 10 == 0 or i + 1 == len(todo):
                    print(f'cache {i + 1}/{len(todo)} {info["sid"]}', flush=True)
    rows = [dict(section_id=s, organ=meta[s]['organ'],
                 platform=str((meta[s].get('metadata') or {}).get('st_technology', 'unknown')),
                 unit_id=meta[s]['unit_id'], source_key=str(meta[s]['source_key']),
                 disease=str((meta[s].get('metadata') or {}).get('disease_state', '')),
                 patient_known=str(meta[s].get('patient_known')),
                 n_spots=int(load_cache(s)['n_spots'])) for s in sids]
    base.tsv_write(OUT / 'cache/sections.tsv', rows)
    return rows


# ---------------------------------------------------------------- assemble
def unit_sections(sections, table):
    groups = defaultdict(list)
    for sid in sections:
        groups[table[sid]['unit_id']].append(sid)
    for unit in groups:
        groups[unit] = sorted(groups[unit])
    return groups


def assemble(target, split_name, splits, table, cache_store=None):
    """Freeze usable A on training units, sample training points, return train/test arrays."""
    spec = splits[split_name]
    groups = unit_sections(spec['train'], table)
    test_groups = unit_sections(spec['test'], table)
    train_units = sorted(groups)
    A_genes, B_genes = split_genes(target)
    store = cache_store if cache_store is not None else {}
    def get(sid):
        if sid not in store:
            store[sid] = load_cache(sid)
        return store[sid]
    # per training unit availability and per-section A counts
    unit_rows = {}
    unit_avail = {}
    for unit in train_units:
        unit_rows[unit] = []
        total = 0
        for sid in groups[unit]:
            data = get(sid)
            n_ok = len(data[f'ok_rows_{target}'])
            unit_rows[unit].append((sid, np.arange(n_ok), data[f'A_nnz_{target}']))
            total += n_ok
        unit_avail[unit] = total
    n_sample = min(N_SAMPLE, min(unit_avail.values())) if unit_avail else 0
    n_a = len(A_genes)
    measurable_units = np.zeros(n_a, int)
    for unit in train_units:
        ok_all = np.ones(n_a, bool)
        unit_maps = []
        for sid, rows, nnz in unit_rows[unit]:
            data = get(sid)
            ok_all &= data[f'A_present_{target}'] & (nnz >= MIN_SPOTS)
            unit_maps.append(data[f'A_map_{target}'][:, :])
        if ok_all.any():
            stacked = np.concatenate(unit_maps, axis=0)
            var = stacked.var(axis=0)
            measurable_units[ok_all] += (var[ok_all] > 0)
    needed = 0.5 * len(train_units)
    usable = np.flatnonzero((measurable_units >= needed) & (measurable_units > 0))
    usable_names = [A_genes[i] for i in usable]
    status = dict(usable_A=int(len(usable)), n_sample=int(n_sample),
                  train_units=len(train_units), test_units=len(test_groups),
                  status='OK')
    if len(usable) < MIN_USABLE_A:
        status['status'] = 'NOT_TESTABLE_INSUFFICIENT_A'
        return dict(status=status, target=target, split=split_name)
    # sample training rows: common fixed count per unit, split as evenly as possible
    train_rows = []
    for unit in train_units:
        shares = [n_sample // len(unit_rows[unit])] * len(unit_rows[unit])
        for i in range(n_sample % len(unit_rows[unit])):
            shares[i] += 1
        for (sid, rows, _), share in zip(unit_rows[unit], shares):
            take = base.sample_indices(len(rows), sid, 0, min(share, len(rows)))
            if len(take) < share:
                raise SystemExit(f'{sid}: only {len(take)} of {share} requested training points')
            train_rows.append((unit, sid, rows[take]))
    train_A_raw, train_B, train_B_pres, train_unit, train_section = [], [], [], [], []
    for unit, sid, rows in train_rows:
        data = get(sid)
        train_A_raw.append(data[f'A_map_{target}'][rows][:, usable])
        train_B.append(data[f'B_map_{target}'][rows])
        train_B_pres.append(np.tile(data[f'B_present_{target}'], (len(rows), 1)))
        train_unit += [unit] * len(rows)
        train_section += [sid] * len(rows)
    train_A_raw = np.concatenate(train_A_raw, axis=0).astype(np.float64)
    train_B = np.concatenate(train_B, axis=0).astype(np.float64)
    train_B_pres = np.concatenate(train_B_pres, axis=0)
    pos = train_A_raw > 0
    scale = np.where(pos.any(axis=0),
                     np.where(pos, train_A_raw, 0).sum(axis=0) / np.maximum(pos.sum(axis=0), 1), np.nan)
    keep = np.isfinite(scale) & (scale > 0)
    status['genes_dropped_no_positive'] = [usable_names[i] for i in np.flatnonzero(~keep)]
    usable_names = [usable_names[i] for i in np.flatnonzero(keep)]
    scale = scale[keep]
    train_A = train_A_raw[:, keep] / scale
    test = {}
    for unit, sids in test_groups.items():
        for sid in sids:
            data = get(sid)
            test[sid] = dict(unit=unit,
                             A=data[f'A_map_{target}'][:, usable[keep]] / scale,
                             B=data[f'B_map_{target}'],
                             B_present=data[f'B_present_{target}'],
                             depth=data[f'depth_{target}'], detection=data[f'detection_{target}'])
    status['usable_A_names'] = usable_names
    return dict(status=status, target=target, split=split_name, B_genes=B_genes,
                train_A=train_A, train_B=train_B, train_B_pres=train_B_pres,
                train_unit=train_unit, train_section=train_section, test=test,
                scale=scale)


# ---------------------------------------------------------------- fit
def fit_nmf(matrix, k):
    model = NMF(n_components=k, **NMF_PARAMS)
    w = model.fit_transform(matrix)
    h = model.components_
    info = dict(k=k, n_iter=int(model.n_iter_), status='OK')
    if not np.all(np.isfinite(h)) or not np.all(np.isfinite(w)):
        info['status'] = 'NOT_TESTABLE_NON_FINITE'
    elif model.n_iter_ >= NMF_PARAMS['max_iter']:
        info['status'] = 'NOT_TESTABLE_NOT_CONVERGED'
    elif np.any(np.all(np.abs(h) < 1e-12, axis=1)):
        info['status'] = 'NOT_TESTABLE_ZERO_LOADING'
    return model, w, h, info


def cosine(a, b):
    den = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.dot(a, b) / den) if den else float('nan')


def pair_sum(h_ref, h_new):
    """Two 1-1 pairings for K=2 (identity, swap); larger cosine sum wins, tie keeps identity."""
    k = h_ref.shape[0]
    if k == 1:
        return [0], cosine(h_ref[0], h_new[0])
    ident = cosine(h_ref[0], h_new[0]) + cosine(h_ref[1], h_new[1])
    swap = cosine(h_ref[0], h_new[1]) + cosine(h_ref[1], h_new[0])
    return ([0, 1], ident) if ident >= swap else ([1, 0], swap)


def q75_high(values):
    if len(values) == 0 or np.nanstd(values) == 0:
        return None
    q = np.nanquantile(values, 0.75)
    mask = values > q
    return mask if mask.any() else None


def correlation(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return float('nan')
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def run_split(args):
    target, split_name, mode = args
    splits = mode_splits()
    table = {r['section_id']: r for r in read_tsv(OUT / 'cache/sections.tsv')}
    store = {}
    built = assemble(target, split_name, splits, table, store)
    status = dict(built['status'])
    status.update(target=target, split=split_name, mode=mode)
    root = RESULTS / ('smoke' if mode == 'smoke' else 'main') / target / split_name
    root.mkdir(parents=True, exist_ok=True)
    if status['status'] != 'OK':
        json_write(root / 'status.json', status)
        return status
    train_A = built['train_A']
    train_B = built['train_B']
    train_B_pres = built['train_B_pres']
    test = built['test']
    B_genes = built['B_genes']
    usable_names = status['usable_A_names']
    # candidate models
    fits = {}
    for k in K_LIST:
        model, w, h, info = fit_nmf(train_A, k)
        fits[k] = dict(model=model, w=w, h=h, info=info)
        status[f'nmf_k{k}'] = info
    # in-fold resample pairing (stability only, never used for prediction)
    stability = []
    rng = np.random.default_rng(0)
    n_train = train_A.shape[0]
    take = max(2, int(round(RESAMPLE_FRACTION * n_train)))
    for r in range(RESAMPLE_COUNT):
        rows = np.sort(np.random.default_rng(r).choice(n_train, take, replace=False))
        for k in (2,):
            _, _, h_res, info_res = fit_nmf(train_A[rows], k)
            if info_res['status'] != 'OK':
                stability.append(dict(resample=r, k=k, cos_sum=float('nan'), status=info_res['status']))
                continue
            order, score = pair_sum(fits[k]['h'], h_res)
            stability.append(dict(resample=r, k=k, cos_sum=float(score), order=str(order),
                                  status='OK'))
    # component similarity and technical correlation
    sim = {}
    for k in K_LIST:
        h = fits[k]['h']
        sim[f'k{k}_component_cosine'] = cosine(h[0], h[1]) if k == 2 else float('nan')
    # training-side representation features
    def features(kind, rows='train'):
        if kind == 'A_mean':
            return train_A.mean(axis=1) if rows == 'train' else None
    # per-gene single gene reference: unscaled training map mean
    raw_means = []
    for i in range(train_A.shape[1]):
        raw = train_A[:, i] * built['scale'][i]
        raw_means.append(float(raw.mean()))
    single_idx = int(np.argmax(raw_means))
    # B ridge per arm per gene
    metrics = []
    arms = ['A_mean', 'K1', 'single_gene', 'K2', 'full_A_ridge']
    def train_features(arm):
        if arm == 'A_mean':
            return train_A.mean(axis=1, keepdims=True)
        if arm == 'K1':
            return fits[1]['w']
        if arm == 'single_gene':
            return train_A[:, [single_idx]]
        if arm == 'K2':
            return fits[2]['w']
        return train_A
    test_feat = {}
    for arm in arms:
        test_feat[arm] = {}
        for sid, block in test.items():
            if arm == 'A_mean':
                test_feat[arm][sid] = block['A'].mean(axis=1, keepdims=True)
            elif arm == 'K1':
                test_feat[arm][sid] = fits[1]['model'].transform(block['A'])
            elif arm == 'single_gene':
                test_feat[arm][sid] = block['A'][:, [single_idx]]
            elif arm == 'K2':
                test_feat[arm][sid] = fits[2]['model'].transform(block['A'])
            else:
                test_feat[arm][sid] = block['A']
    for arm in arms:
        X_train = train_features(arm)
        scaler = StandardScaler().fit(X_train)
        Xs = scaler.transform(X_train)
        for j, gene in enumerate(B_genes):
            rows_ok = train_B_pres[:, j]
            if rows_ok.sum() < MIN_SPOTS:
                for sid in test:
                    metrics.append(dict(gene=gene, arm=arm, section=sid, unit=test[sid]['unit'],
                                        pearson='', jaccard='', status='NOT_TESTABLE_TRAIN_ABSENT'))
                continue
            model = Ridge(alpha=10).fit(Xs[rows_ok], train_B[rows_ok, j])
            for sid, block in test.items():
                if not block['B_present'][j]:
                    metrics.append(dict(gene=gene, arm=arm, section=sid, unit=block['unit'],
                                        pearson='', jaccard='', status='NOT_TESTABLE_GENE_ABSENT'))
                    continue
                truth = block['B'][:, j]
                nnz = int((truth > 0).sum())
                if nnz < MIN_SPOTS or np.std(truth) == 0:
                    metrics.append(dict(gene=gene, arm=arm, section=sid, unit=block['unit'],
                                        pearson='', jaccard='', status='NOT_TESTABLE_SPARSE_OR_CONSTANT'))
                    continue
                pred = model.predict(scaler.transform(test_feat[arm][sid]))
                pear = correlation(pred, truth)
                hi_pred = q75_high(pred)
                hi_true = q75_high(truth)
                if hi_pred is None or hi_true is None:
                    jac = float('nan'); jstatus = 'NOT_ESTIMABLE_TIES_OR_CONSTANT'
                else:
                    jac = float((hi_pred & hi_true).sum() / max((hi_pred | hi_true).sum(), 1))
                    jstatus = 'OK'
                metrics.append(dict(gene=gene, arm=arm, section=sid, unit=block['unit'],
                                    pearson=('' if not np.isfinite(pear) else f'{pear:.6f}'),
                                    jaccard=('' if not np.isfinite(jac) else f'{jac:.6f}'),
                                    status='OK' if jstatus == 'OK' else jstatus))
    # technical correlations on test sections
    technical = []
    for sid, block in test.items():
        scores = fits[2]['model'].transform(block['A'])
        for c in range(2):
            technical.append(dict(section=sid, unit=block['unit'], component=c,
                                  depth_corr=correlation(scores[:, c], block['depth']),
                                  detection_corr=correlation(scores[:, c], block['detection'])))
    # training-only B-component association (never used for prediction)
    assoc = []
    for j, gene in enumerate(B_genes):
        rows_ok = train_B_pres[:, j]
        for c in range(2):
            assoc.append(dict(gene=gene, component=c,
                              train_corr=correlation(fits[2]['w'][rows_ok, c], train_B[rows_ok, j])))
    scores = {}
    for sid, block in test.items():
        scores[f'{sid}_k1'] = fits[1]['model'].transform(block['A'])
        scores[f'{sid}_k2'] = fits[2]['model'].transform(block['A'])
    np.savez_compressed(root / 'scores.npz', **scores)
    np.savez_compressed(root / 'loadings.npz',
                        h1=fits[1]['h'], h2=fits[2]['h'], scale=built['scale'],
                        usable=np.array(usable_names, dtype=object), single_gene=usable_names[single_idx],
                        B_genes=np.array(B_genes, dtype=object))
    base.tsv_write(root / 'metrics.tsv', metrics)
    base.tsv_write(root / 'technical.tsv', technical)
    base.tsv_write(root / 'b_association_train.tsv', assoc)
    json_write(root / 'status.json', dict(status, stability=stability, similarity=sim,
                                          single_gene=usable_names[single_idx]))
    return status


def run(mode='main', targets=None, splits=None, workers=3):
    require_release()
    targets = targets or TARGETS
    all_splits = mode_splits()
    names = splits or [k for k in all_splits]
    jobs = [(t, s, mode) for t in targets for s in names]
    with multiprocessing.get_context('spawn').Pool(workers) as pool:
        for i, status in enumerate(pool.imap_unordered(run_split, jobs)):
            print(json.dumps(dict(done=i + 1, of=len(jobs), **{k: status[k] for k in
                                                              ('target', 'split', 'status', 'usable_A')
                                                              if k in status})), flush=True)
    return 0


# ---------------------------------------------------------------- summarize
def summarize(mode='main', out_name='summary'):
    root = RESULTS / mode
    rows = []
    for path in sorted(root.glob('*/*/metrics.tsv')):
        target, split = path.parent.parent.name, path.parent.name
        for row in read_tsv(path):
            row.update(target=target, split=split, mode=mode)
            rows.append(row)
    statuses = {}
    for path in sorted(root.glob('*/*/status.json')):
        statuses[f'{path.parent.parent.name}|{path.parent.name}'] = json_read(path)
    # unit aggregation: per gene/arm/unit, require every section of that unit evaluable
    by_unit = defaultdict(dict)
    for row in rows:
        key = (row['target'], row['split'], row['arm'], row['gene'])
        by_unit[key].setdefault(row['unit'], []).append(row)
    unit_rows = []
    for (target, split, arm, gene), units in by_unit.items():
        for unit, items in units.items():
            if any(i['status'] != 'OK' for i in items):
                unit_rows.append(dict(target=target, split=split, arm=arm, gene=gene, unit=unit,
                                      pearson='', jaccard='', n_sections=len(items),
                                      status='NOT_ESTIMABLE_UNIT'))
                continue
            pears = [float(i['pearson']) for i in items if i['pearson']]
            jacs = [float(i['jaccard']) for i in items if i['jaccard']]
            unit_rows.append(dict(target=target, split=split, arm=arm, gene=gene, unit=unit,
                                  pearson=(f'{np.mean(pears):.6f}' if pears else ''),
                                  jaccard=(f'{np.mean(jacs):.6f}' if jacs else ''),
                                  n_sections=len(items), status='OK'))
    # complete-unit macro average per target/split/arm/gene
    macro = []
    groups = defaultdict(list)
    for row in unit_rows:
        if row['status'] == 'OK':
            groups[(row['target'], row['split'], row['arm'], row['gene'])].append(row)
    for key, items in groups.items():
        pears = [float(i['pearson']) for i in items if i['pearson']]
        jacs = [float(i['jaccard']) for i in items if i['jaccard']]
        macro.append(dict(target=key[0], split=key[1], arm=key[2], gene=key[3],
                          n_units=len(items),
                          pearson=(f'{np.mean(pears):.6f}' if pears else ''),
                          jaccard=(f'{np.mean(jacs):.6f}' if jacs else '')))
    base.tsv_write(OUT / f'{out_name}_per_section.tsv', rows)
    base.tsv_write(OUT / f'{out_name}_per_unit.tsv', unit_rows)
    base.tsv_write(OUT / f'{out_name}_macro.tsv', macro)
    json_write(OUT / f'{out_name}_status.json', statuses)
    counts = Counter((r['target'], r['arm']) for r in macro)
    print(json.dumps({'per_section_rows': len(rows), 'per_unit_rows': len(unit_rows),
                      'macro_rows': len(macro), 'statuses': len(statuses),
                      'macro_by_target_arm': {'|'.join(k): v for k, v in sorted(counts.items())}},
                     ensure_ascii=False)[:1500])
    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['cache', 'smoke', 'run', 'summarize', 'verify'])
    parser.add_argument('--workers', type=int, default=3)
    parser.add_argument('--targets', nargs='*')
    parser.add_argument('--splits', nargs='*')
    args = parser.parse_args()
    if args.command == 'cache':
        cache(args.workers, args.splits)
    elif args.command == 'smoke':
        require_release()
        return run('smoke', args.targets or ['B_AXIS'], args.splits or ['internal_0'], 1)
    elif args.command == 'run':
        return run('main', args.targets, args.splits, args.workers)
    elif args.command == 'summarize':
        return summarize('main')
    elif args.command == 'verify':
        import field_candidate_verify as verify
        return verify.main()
    return 0


if __name__ == '__main__':
    sys.exit(main() or 0)
