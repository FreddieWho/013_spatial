"""CPU full-GOBP stage. Write numerical artifacts during execution; report only at closeout."""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import json
from pathlib import Path
import shutil
import time

import numpy as np
from scipy.spatial import cKDTree

from r16.recovery.corrected import (NaturalSpline, NotTestable, fit_from_moments,
                                    signed_tls_distance, tls_components, classify_window)
from r16.recovery.gobp_stage import FeatureEngine, geometry, mask_predictions, moments, neighbor_operator
from r16.recovery.mask_tasks import generate_windows, window_mask
from r16.recovery.repair_pipeline import axes
from scripts.r16_directional_eligibility import sectors

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'infra/repair_20260921'
OUT = ROOT / 'infra/gobp_stage_20260922'


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


def atomic_json(path, obj):
    tmp = path.with_suffix('.tmp')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def direction_weights(xy, labels):
    records, weights = [], []
    for i, component in enumerate(tls_components(xy, labels[:, 0])):
        distance, bins, status = sectors(xy, labels, component)
        records.append(dict(focus=i, n_tls_spots=len(component), distance=None if not np.isfinite(distance) else distance,
                            n_bands=len(bins), status=status))
        if status != 'ELIGIBLE':
            continue
        w = np.zeros(len(xy))
        for a, b, c in bins:
            w[a] += 1 / len(a) / len(bins)
            w[b] -= .5 / len(b) / len(bins)
            w[c] -= .5 / len(c) / len(bins)
        weights.append(w)
    return records, np.column_stack(weights) if weights else np.empty((len(xy), 0))


def distance_splits(xy, labels):
    tum = labels[:, 1] == 1
    if tum.sum() < 50 or not (labels[:, 0] == 1).any():
        return tum, []
    distance = signed_tls_distance(xy, labels[:, 0])
    t = (distance / np.quantile(distance[distance > 0], .99))[tum]
    tx = xy[tum]
    folds = np.empty(len(tx), int)
    for i, idx in enumerate(np.array_split(np.argsort(tx[:, 0], kind='stable'), 5)):
        folds[idx] = i
    splits = []
    for i in range(5):
        te = folds == i
        tr = (folds != i) & (cKDTree(tx[te]).query(tx)[0] > 1.01)
        if tr.sum() < 30 or te.sum() < 5:
            continue
        try:
            spline = NaturalSpline(t[tr])
        except NotTestable:
            continue
        splits.append((tr, te, spline(t[tr]), spline(t[te])))
    return tum, splits


def distance_errors(x, tum, splits):
    errors = np.full((x.shape[1], 5, 2), np.nan)
    if not splits:
        return errors
    q, target = x[tum, :, :9], x[tum, :, 9]
    for i, (tr, te, bt, be) in enumerate(splits):
        xt = np.concatenate([q[tr], np.broadcast_to(bt[:, None, :], (tr.sum(), x.shape[1], bt.shape[1]))], axis=2)
        xe = np.concatenate([q[te], np.broadcast_to(be[:, None, :], (te.sum(), x.shape[1], be.shape[1]))], axis=2)
        gram, rhs, _ = moments(xt, target[tr])
        for j, n in enumerate([9, xt.shape[2]]):
            beta = fit_from_moments(gram[:, :n, :n], rhs[:, :n])
            pred = np.einsum('nbk,bk->nb', xe[:, :, :n], beta)
            errors[:, i, j] = np.mean((target[te] - pred) ** 2, axis=0)
    return errors


def worker(task):
    section, mode, limit, directory = task
    out = Path(directory)
    key = section['key']
    tag = f'{key}_{mode}'
    receipt_path = out / f'{tag}.json'
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_text())
        if receipt['status'] == 'COMPLETED' and digest(out / f'{tag}.npz') == receipt['output_sha256']:
            return receipt
    started = time.monotonic()
    if shutil.disk_usage(ROOT).free < 1_200_000_000_000:
        raise RuntimeError('Storage reserve below 1.2 TB')
    definitions = json.loads((SOURCE / 'go_representatives.json').read_text())
    definitions = {g: definitions[g] for g in sorted(definitions)[:limit or None]}
    genes = json.loads((SOURCE / 'genes.json').read_text())
    counts_path = SOURCE / 'cache' / f'{key}_counts.npy'
    assert digest(counts_path) == section['counts_sha256']
    counts = np.load(counts_path, mmap_mode='r')
    xy = np.load(SOURCE / 'cache' / f'{key}_coords.npy')
    labels = np.load(SOURCE / 'cache' / f'{key}_labels.npy')
    source_genes = set(json.loads((SOURCE / 'cache' / f'{key}_source_genes.json').read_text()))
    engine = FeatureEngine(counts, genes, axes(), definitions, unavailable=set(genes) - source_genes)
    wins = generate_windows(xy, 12, 8, 20260921)
    geometries, window_states = [], []
    for win in wins:
        hidden = window_mask(xy, win)
        window_states.append(dict(n_hidden=int(hidden.sum()), label=classify_window(labels[hidden, 0])))
        try:
            geometries.append(geometry(xy, hidden))
        except NotTestable:
            geometries.append(None)
    foci, dw = direction_weights(xy, labels) if section['cohort'] == 'USZ' else ([], np.empty((len(xy), 0)))
    tum, splits = distance_splits(xy, labels) if section['cohort'] == 'USZ' else (np.zeros(len(xy), bool), [])
    n = len(definitions)
    ids = list(definitions)
    index = {g: i for i, g in enumerate(ids)}
    grams = np.full((n, 11, 11), np.nan)
    rhs = np.full((n, 11), np.nan)
    yy = np.full(n, np.nan)
    variance = np.full(n, np.nan)
    masked = np.full((n, len(wins), 5), np.nan)
    mask_var = np.full((n, len(wins)), np.nan)
    direction = np.full((n, dw.shape[1]), np.nan)
    cv = np.full((n, 5, 2), np.nan)
    status = np.full(n, 'NOT_RUN', dtype='U40')
    neighbors = neighbor_operator(xy, xy, True)
    for step, (batch, x, y, bad) in enumerate(engine.blocks(mode)):
        for g, reason in bad:
            status[index[g]] = reason
        if not batch:
            continue
        ix = [index[g] for g in batch]
        full = np.concatenate([x, (neighbors @ x[:, :, -1])[:, :, None]], axis=2)
        gg, rr, zz = moments(full, y)
        grams[ix], rhs[ix], yy[ix] = gg, rr, zz
        variance[ix] = y.var(0)
        status[ix] = 'COMPUTED'
        for wi, geo in enumerate(geometries):
            if geo is None:
                continue
            pred = mask_predictions(x, y, geo)
            masked[ix, wi] = np.mean((y[geo['hidden'], :, None] - pred) ** 2, axis=0)
            mask_var[ix, wi] = y[geo['hidden']].var(0)
        if dw.shape[1]:
            beta = fit_from_moments(gg[:, :9, :9], rr[:, :9])
            residual = y - np.einsum('nbk,bk->nb', x[:, :, :9], beta)
            sd = residual.std(0)
            delta = (dw.T @ residual).T
            direction[ix] = np.divide(delta, sd[:, None], out=np.full_like(delta, np.nan), where=sd[:, None] > 1e-12)
        if splits:
            cv[ix] = distance_errors(x, tum, splits)
        if step % 12 == 0:
            atomic_json(out / f'{tag}_progress.json', dict(status='RUNNING', completed=min((step + 1) * 96, n), total=n,
                                                         elapsed_seconds=round(time.monotonic() - started, 2)))
    target = out / f'{tag}.npz'
    tmp = out / f'{tag}.tmp.npz'
    np.savez_compressed(tmp, ids=np.array(ids), status=status, gram=grams, rhs=rhs, yy=yy,
                        target_variance=variance, masked_mse=masked, masked_target_variance=mask_var,
                        directional_delta=direction, distance_cv_mse=cv)
    tmp.replace(target)
    receipt = dict(status='COMPLETED', section=section['section'], patient=section['patient'], cohort=section['cohort'],
                   key=key, mode=mode, n_programs=n, computed=int((status == 'COMPUTED').sum()),
                   not_testable=int((status != 'COMPUTED').sum()), valid_windows=sum(g is not None for g in geometries),
                   windows=window_states, foci=foci, n_direction_foci=dw.shape[1], n_distance_folds=len(splits),
                   output_sha256=digest(target), elapsed_seconds=round(time.monotonic() - started, 2))
    atomic_json(receipt_path, receipt)
    atomic_json(out / f'{tag}_progress.json', dict(status='COMPLETED', completed=n, total=n))
    return receipt


def prepare(out, limit=0):
    out.mkdir(exist_ok=True, parents=True)
    (out / 'sections').mkdir(exist_ok=True)
    input_paths = [SOURCE / name for name in ['sections.json', 'genes.json', 'go_representatives.json', 'go_aliases.json']]
    input_paths += [ROOT / 'infra/r04/marker_proxy_combined.json', ROOT / 'docs/plan.md']
    code_paths = [Path(__file__), ROOT / 'r16/recovery/gobp_stage.py', ROOT / 'r16/recovery/corrected.py',
                  ROOT / 'r16/recovery/mask_tasks.py', ROOT / 'scripts/r16_directional_eligibility.py']
    sections = [r for r in json.loads((SOURCE / 'sections.json').read_text()) if r['cohort'] in {'ST-CRC', 'USZ'}]
    for s in sections:
        input_paths += [SOURCE / 'cache' / f"{s['key']}_{suffix}" for suffix in ['coords.npy', 'labels.npy', 'source_genes.json']]
    contract = dict(version=1, programs=limit or 6870, modes=['raw', 'depth_normalized'], sections=22,
                    scope='Full frozen GOBP; molecular LOPO, masked readout prediction, known-anchor direction, blocked distance prediction',
                    molecular='Separate-cohort LOPO; equal section then patient weight; M0 QC+C6, M1 + input, M2 + input neighbors',
                    mask='12 fixed GT-independent windows/section; M0/M1/M2/NN/KNN8, visible-only; no phenotype labels in prediction',
                    direction='D-143 unchanged; readout QC+C6 residual; equal foci then section then patient; two-sided exact patient sign flips',
                    distance='All 6870, input score; QC+C6 vs +natural spline, five x-stripes, >1 pitch buffer, train-only knots',
                    calibration='Real-geometry smooth Gaussian null/injections; exact patient sign-flip resolution and known-model calibration only; not T2P2 real-data spatial p',
                    inference='Patient sign-flip symmetry assumption; Holm within endpoint families; empirical spatial association p remains NOT_CALIBRATED',
                    exclusions='No new source data/GPU/retraining/general geometry dictionary; structure AUC remains NOT_TESTABLE if no positive windows',
                    predeclared_checks='No threshold tuning; raw/norm whole-panel reporting; constant endpoints explicit; partial results cannot close full stage',
                    input_hashes={str(p.relative_to(ROOT)): digest(p) for p in input_paths},
                    code_hashes={str(p.relative_to(ROOT)): digest(p) for p in code_paths})
    path = out / 'run_contract.json'
    if path.exists():
        assert json.loads(path.read_text()) == contract, 'Run contract changed; use a new output directory'
    else:
        atomic_json(path, contract)
    return sections


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, default=OUT)
    parser.add_argument('--limit', type=int, default=0)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--section', default=None)
    args = parser.parse_args()
    sections = prepare(args.out, args.limit)
    if args.section:
        sections = [s for s in sections if s['key'] == args.section]
    tasks = [(s, mode, args.limit, str(args.out / 'sections')) for s in sections for mode in ['raw', 'depth_normalized']]
    started = time.monotonic()
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(worker, task) for task in tasks]
        for i, future in enumerate(as_completed(futures)):
            r = future.result()
            print(f"{i + 1}/{len(tasks)} {r['key']} {r['mode']} {r['computed']}/{r['n_programs']} {r['elapsed_seconds']:.1f}s", flush=True)
            atomic_json(args.out / 'progress.json', dict(status='RUNNING', completed_tasks=i + 1, total_tasks=len(tasks)))
    atomic_json(args.out / 'compute_receipt.json', dict(status='COMPLETED', tasks=len(tasks), expected_full_tasks=44,
                                                       programs=args.limit or 6870, elapsed_seconds=time.monotonic() - started))
    atomic_json(args.out / 'progress.json', dict(status='COMPLETED', completed_tasks=len(tasks), total_tasks=len(tasks)))


if __name__ == '__main__':
    main()
