#!/usr/bin/env python3
"""FA-13 runner: conditional cross-field increment outside the target neighbourhood.

Commands
  smoke   synthetic positive/negative cases and hidden-entry poisoning checks
  run     the real pilot: 14 sections / 8 units, two rules, two gene-half directions
  assets  read-only inventory of local Xenium / Visium HD archives

Everything is CPU-only, local-data-only and bounded; no human TLS label is read and no file is
downloaded.
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import os
import re
import sys
import tarfile
import time
import zipfile
import zlib
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from field_analyzer import adapters, assets, conditional_relations as cr  # noqa: E402

OUT = ROOT / 'infra/field_analyzer_20261023'
PILOT_ROOT = ROOT / 'data/other_sources/hest1k/st'
IDENTITY = ROOT / 'infra/field_discovery_review_20261018/identity_review.tsv'
ZIP_DATASETS = ROOT / 'data/other_sources/10x_genomics'


def write_tsv(rows: Sequence[dict], path: Path) -> int:
    if not rows:
        path.write_text('')
        return 0
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as handle:
        handle.write('\t'.join(keys) + '\n')
        for row in rows:
            cells = []
            for key in keys:
                value = row.get(key, '')
                if isinstance(value, float):
                    cells.append('' if not np.isfinite(value) else f'{value:.6g}')
                elif isinstance(value, (dict, list)):
                    cells.append(json.dumps(value, ensure_ascii=False))
                else:
                    cells.append(' '.join(str(value).split()))
            handle.write('\t'.join(cells) + '\n')
    return len(rows)


# --------------------------------------------------------------------------- data
def pilot_sections() -> list[tuple[str, str]]:
    identity = assets.read_identity(IDENTITY)
    out: list[tuple[str, str]] = []
    for unit, sections in sorted(identity.units.items()):
        for section in sorted(sections):
            out.append((unit, section))
    return out


def load_designs(panel: dict[str, list[str]]) -> tuple[list[cr.SectionDesign], list[dict]]:
    designs, meta = [], []
    for unit, section in pilot_sections():
        table = adapters.read_h5ad_counts(PILOT_ROOT / f'{section}.h5ad')
        in_tissue = np.asarray(table.aux.get('in_tissue', np.ones(table.n_obs, dtype=bool)),
                               dtype=bool)
        design = cr.build_section(section, unit, table.counts, table.genes, table.xy, in_tissue,
                                  panel, barcodes=table.barcodes)
        assert design.barcodes is not None and len(design.barcodes) == design.n_points
        assert design.barcodes[0] == table.barcodes[0] and design.barcodes[-1] == table.barcodes[-1]
        designs.append(design)
        meta.append(dict(section_id=section, unit_id=unit, n_points=design.n_points,
                         n_in_tissue=int(in_tissue.sum()), n_q=int(design.q_mask.sum()),
                         n_observed=int(design.observed.sum()),
                         n_genes=int(table.counts.shape[1]),
                         genes_a={k[0]: len(design.gene_sets[(k[0], 'A', 'feature')])
                                  for k in design.scores if k[1] == 'A' and k[2] == 'feature'},
                         genes_b={k[0]: len(design.gene_sets[(k[0], 'B', 'feature')])
                                  for k in design.scores if k[1] == 'B' and k[2] == 'feature'},
                         n_background=design.meta['n_background']))
    return designs, meta


# --------------------------------------------------------------------------- smoke
def synthetic_design(*, seed: int, joint: bool) -> cr.SectionDesign:
    """A section where the target depends on its own neighbourhood (always) plus the other field."""
    rng = np.random.default_rng(seed)
    side = 26
    rows, cols = np.meshgrid(np.arange(side), np.arange(side), indexing='ij')
    xy = np.column_stack([cols.ravel(), rows.ravel()]) * 100.0
    n = xy.shape[0]
    in_tissue = np.ones(n, dtype=bool)
    predictor = np.sin(cols.ravel() / 3.0) + 0.4 * rng.normal(size=n)
    own = np.cos(rows.ravel() / 4.0) + 0.4 * rng.normal(size=n)
    from scipy.spatial import cKDTree
    tree = cKDTree(xy)
    _, idx = tree.query(xy, k=9)
    idx = idx[:, 1:]
    own_nb = own[idx].mean(axis=1)
    cross_nb = predictor[idx].mean(axis=1)
    target = 0.9 * own_nb + (0.9 * cross_nb if joint else 0.0) + 0.15 * rng.normal(size=n)
    nil = np.zeros(n)
    scores = {('T', 'A', 'feature'): own, ('T', 'A', 'truth'): target,
              ('T', 'B', 'feature'): own, ('T', 'B', 'truth'): target,
              ('P', 'A', 'feature'): predictor, ('P', 'A', 'truth'): nil,
              ('P', 'B', 'feature'): predictor, ('P', 'B', 'truth'): nil}
    gene_sets = {key: (f'{key[0]}_{i}',) for i, key in enumerate(scores)}
    q_mask = np.zeros(n, dtype=bool)
    q_mask[::10] = True
    return cr.SectionDesign(section_id=f'syn{seed}', unit_id=f'u{seed}', xy=xy,
                            in_tissue=in_tissue, q_mask=q_mask, log_lib=np.zeros(n),
                            scores=scores, gene_sets=gene_sets, meta=dict(synthetic=True))


def run_smoke() -> dict:
    checks: dict[str, Any] = {}
    # positive: cross alone weak, same alone moderate, joint better than same
    sections = [synthetic_design(seed=s, joint=True) for s in range(6)]
    unit_rows, section_rows = cr.run_rule(sections, 'T', 'P', conditioned=False,
                                          protocol='synthetic_positive')
    summary = cr.summarise(unit_rows)
    positive = summary['T<-P|A/B|synthetic_positive']
    checks['synthetic_positive'] = dict(
        n_units=positive['n_units'], n_positive=positive['n_positive'],
        median_delta_r2=positive['median_delta_r2'],
        median_r2_m0=positive['median_r2_m0'], median_r2_m1=positive['median_r2_m1'],
        passes=bool(positive['median_delta_r2'] > 0 and positive['n_positive'] >= 5))
    # negative: target depends only on its own neighbourhood
    sections_neg = [synthetic_design(seed=100 + s, joint=False) for s in range(6)]
    unit_neg, _ = cr.run_rule(sections_neg, 'T', 'P', conditioned=False,
                              protocol='synthetic_negative')
    neg = cr.summarise(unit_neg)['T<-P|A/B|synthetic_negative']
    checks['synthetic_negative'] = dict(
        n_units=neg['n_units'], n_positive=neg['n_positive'],
        median_delta_r2=neg['median_delta_r2'],
        passes=bool(abs(neg['median_delta_r2']) < 0.01 or neg['n_positive'] < 4))
    # separate single-feature references on the positive data
    same_only = cr.run_single_feature(sections, 'T', 'P', feature='same', conditioned=False)
    cross_only = cr.run_single_feature(sections, 'T', 'P', feature='cross', conditioned=False)
    same_med = float(np.median([r['r2'] for r in same_only if r['status'] == 'OK']))
    cross_med = float(np.median([r['r2'] for r in cross_only if r['status'] == 'OK']))
    joint_med = float(np.median([r['r2_m1'] for r in unit_rows if r['status'] == 'OK']))
    checks['separate_models'] = dict(
        same_only_median_r2=same_med, cross_only_median_r2=cross_med,
        joint_median_r2=joint_med,
        cross_weaker_than_same=bool(cross_med < same_med),
        joint_beats_same=bool(joint_med > same_med))
    # poisoning: hidden raw counts and hidden depth must not change features, fits or predictions
    poison = poison_check()
    checks['poison'] = poison
    checks['real_section_poison'] = real_poison_check()
    checks['verdict'] = dict(
        positive_has_increment=checks['synthetic_positive']['passes'],
        negative_has_no_increment=checks['synthetic_negative']['passes'],
        hidden_poison_leaves_inputs_unchanged=poison['features_identical'],
        hidden_poison_changes_truth=poison['truth_changed'],
    )
    return checks


def poison_check() -> dict:
    """Poison Q counts/depth; the design must not move, the offline truth must move."""
    rng = np.random.default_rng(7)
    n, n_genes = 400, 60
    counts = rng.poisson(3.0, size=(n, n_genes)).astype(np.int64)
    genes = np.array([f'G{i}' for i in range(n_genes)], dtype=object)
    panel = {'T': [f'G{i}' for i in range(10)],
             'P': [f'G{i}' for i in range(10, 20)],
             'X': [f'G{i}' for i in range(20, 30)]}
    xy = rng.uniform(0, 400, size=(n, 2))
    in_tissue = np.ones(n, dtype=bool)
    base = cr.build_section('poison', 'unit', sparse.csr_matrix(counts), genes, xy, in_tissue,
                            panel, seed=0)
    poisoned_counts = counts.copy()
    hidden = ~base.observed
    poisoned_counts[hidden, :] = 50_000
    poisoned = cr.build_section('poison', 'unit', sparse.csr_matrix(poisoned_counts), genes, xy,
                                in_tissue, panel, seed=0)
    d0 = cr.fold_design(base, 'T', 'P', conditioned=True)
    d1 = cr.fold_design(poisoned, 'T', 'P', conditioned=True)
    feature_identical = all(
        np.array_equal(np.nan_to_num(d0[key], nan=-12345.0),
                       np.nan_to_num(d1[key], nan=-12345.0))
        for key in ('same', 'cross'))
    q = base.q_mask
    truth_changed = bool(np.any(np.abs(np.nan_to_num(d0['y'][q]) - np.nan_to_num(d1['y'][q])) > 1e-9))
    params_identical = d0['params'] == d1['params']
    # fitted coefficients and predictions must also be identical (fit rows are O rows only)
    design0, design1 = [], []
    for d in (d0, d1):
        sel = d['obs'] & np.isfinite(d['y']) & np.isfinite(d['same']) & np.isfinite(d['cross'])
        design0.append(d) if d is d0 else design1.append(d)
        d['_sel'] = sel
    X = np.column_stack([d0['same'][d0['_sel']], d0['cross'][d0['_sel']]])
    y = d0['y'][d0['_sel']]
    beta0 = cr._ridge(X, y, 1.0)
    X1 = np.column_stack([d1['same'][d1['_sel']], d1['cross'][d1['_sel']]])
    y1 = d1['y'][d1['_sel']]
    beta1 = cr._ridge(X1, y1, 1.0)
    qsel = q & np.isfinite(d0['same']) & np.isfinite(d0['cross'])
    pred0 = beta0[0] + beta0[1] * d0['same'][qsel] + beta0[2] * d0['cross'][qsel]
    pred1 = beta1[0] + beta1[1] * d1['same'][qsel] + beta1[2] * d1['cross'][qsel]
    return dict(features_identical=bool(feature_identical), truth_changed=truth_changed,
                calibration_params_identical=bool(params_identical),
                coefficients_identical=bool(np.array_equal(beta0, beta1)),
                predictions_identical=bool(np.array_equal(pred0, pred1)),
                n_hidden=int(hidden.sum()), n_q=int(q.sum()), n_fit=int(y.size))


def real_poison_check(section_id: str = 'NCBI643') -> dict:
    """Same poisoning contract on one real pilot section (raw counts, hidden Q rows)."""
    panel = cr.load_panel()
    table = adapters.read_h5ad_counts(PILOT_ROOT / f'{section_id}.h5ad')
    in_tissue = np.asarray(table.aux.get('in_tissue', np.ones(table.n_obs, dtype=bool)), dtype=bool)
    counts = table.counts.toarray() if hasattr(table.counts, 'toarray') else np.asarray(table.counts)
    base = cr.build_section(section_id, 'unit', sparse.csr_matrix(counts), table.genes, table.xy,
                            in_tissue, panel)
    poisoned_counts = counts.copy()
    hidden = ~base.observed
    poisoned_counts[hidden, :] = 10_000
    poisoned = cr.build_section(section_id, 'unit', sparse.csr_matrix(poisoned_counts),
                               table.genes, table.xy, in_tissue, panel)
    d0 = cr.fold_design(base, 'EMT', 'IFNG', conditioned=True)
    d1 = cr.fold_design(poisoned, 'EMT', 'IFNG', conditioned=True)
    identical = all(np.array_equal(np.nan_to_num(d0[k], nan=-1.0), np.nan_to_num(d1[k], nan=-1.0))
                    for k in ('same', 'cross'))
    sel = d0['obs'] & np.isfinite(d0['y']) & np.isfinite(d0['same']) & np.isfinite(d0['cross'])
    X = np.column_stack([d0['same'][sel], d0['cross'][sel]])
    y = d0['y'][sel]
    b0 = cr._ridge(X, y, cr.RIDGE_ALPHA)
    X1 = np.column_stack([d1['same'][sel], d1['cross'][sel]])
    b1 = cr._ridge(X1, d1['y'][sel], cr.RIDGE_ALPHA)
    q = base.q_mask
    qsel = q & np.isfinite(d0['same']) & np.isfinite(d0['cross'])
    p0 = b0[0] + b0[1] * d0['same'][qsel] + b0[2] * d0['cross'][qsel]
    p1 = b1[0] + b1[1] * d1['same'][qsel] + b1[2] * d1['cross'][qsel]
    truth0 = d0['y'][q]
    truth1 = d1['y'][q]
    return dict(section_id=section_id, n_points=int(base.n_points), n_q=int(q.sum()),
                n_hidden=int(hidden.sum()),
                features_identical=bool(identical),
                calibration_identical=bool(d0['params'] == d1['params']),
                coefficients_identical=bool(np.array_equal(b0, b1)),
                predictions_identical=bool(np.array_equal(p0, p1)),
                truth_changed=bool(np.any(np.abs(np.nan_to_num(truth0) - np.nan_to_num(truth1)) > 1e-9)))


# --------------------------------------------------------------------------- run
def run_real(out: Path) -> dict:
    start = time.time()
    panel = cr.load_panel()
    designs, meta = load_designs(panel)
    units = sorted({d.unit_id for d in designs})
    results = dict(units=units, sections=meta, settings=dict(
        k=cr.K_NEIGHBOURS, alpha=cr.RIDGE_ALPHA, q_fraction=cr.Q_FRACTION, seed=cr.Q_SEED,
        rules=[list(r) for r in cr.RULES], stable_units=cr.STABLE_UNITS,
        fields=list(cr.FIELDS), panel_path=cr.PANEL_PATH))
    unit_rows: list[dict] = []
    section_rows: list[dict] = []
    for target, predictor in cr.RULES:
        for feature_half, truth_half in (('A', 'B'), ('B', 'A')):
            for protocol, conditioned in (('conditioned', True), ('unconditioned', False)):
                urows, srows = cr.run_rule(designs, target, predictor, feature_half=feature_half,
                                           truth_half=truth_half, conditioned=conditioned,
                                           protocol=protocol)
                for row in urows:
                    row['rule_role'] = 'primary' if (target, predictor) == cr.PRIMARY_RULE \
                        else 'depth_control'
                unit_rows.extend(urows)
                section_rows.extend(srows)
    summary = cr.summarise(unit_rows)
    references = {}
    for target, predictor in cr.RULES:
        for half, truth_half in (('A', 'B'), ('B', 'A')):
            for protocol, conditioned in (('conditioned', True), ('unconditioned', False)):
                same = cr.run_single_feature(designs, target, predictor, feature='same',
                                             feature_half=half, truth_half=truth_half,
                                             conditioned=conditioned)
                cross = cr.run_single_feature(designs, target, predictor, feature='cross',
                                              feature_half=half, truth_half=truth_half,
                                              conditioned=conditioned)
                key = f'{target}<-{predictor}|{half}/{truth_half}|{protocol}'
                references[key] = dict(
                    same_only_median_r2=float(np.median([r['r2'] for r in same
                                                         if r['status'] == 'OK'])),
                    cross_only_median_r2=float(np.median([r['r2'] for r in cross
                                                          if r['status'] == 'OK'])),
                    n_units=int(sum(1 for r in same if r['status'] == 'OK')))
    results['single_feature_references'] = references
    results['summary'] = summary
    results['n_unit_rows'] = len(unit_rows)
    results['n_section_rows'] = len(section_rows)
    results['wall_seconds'] = round(time.time() - start, 1)
    write_tsv(unit_rows, out / 'unit_results.tsv')
    write_tsv(section_rows, out / 'section_results.tsv')
    write_tsv(meta, out / 'section_inputs.tsv')
    (out / 'summary.json').write_text(json.dumps(results, indent=2, ensure_ascii=False,
                                                 default=str) + '\n')
    return results


# --------------------------------------------------------------------------- assets
PANEL_FOR_ASSETS = {}


def _panel_genes() -> dict[str, set[str]]:
    if PANEL_FOR_ASSETS:
        return PANEL_FOR_ASSETS
    panel = json.loads((ROOT / cr.PANEL_PATH).read_text())
    for short, key in cr.FIELDS.items():
        PANEL_FOR_ASSETS[short] = set(panel[key])
    return PANEL_FOR_ASSETS


def _inflate(raw: bytes) -> bytes:
    """Zip method 8 stores raw deflate; a member may itself be a gzip file (``cells.csv.gz``)."""
    for wbits in (-zlib.MAX_WBITS, zlib.MAX_WBITS | 32):
        try:
            out = zlib.decompress(raw, wbits)
            break
        except zlib.error:
            out = None
    if out is None:
        try:
            out = zlib.decompressobj(-zlib.MAX_WBITS).decompress(raw)
        except zlib.error:
            return b''
    if out[:2] == b'\x1f\x8b':
        try:
            out = gzip.decompress(out)
        except Exception:
            pass
    return out


def _zip_tail_has_eocd(path: Path, window: int = 1_000_000) -> bool:
    size = path.stat().st_size
    with path.open('rb') as handle:
        handle.seek(max(0, size - window))
        tail = handle.read()
    return tail.rfind(b'PK\x05\x06') != -1


def _read_zip_feature_axis(path: Path, budget_members: int = 120,
                           budget_bytes: int = 200_000_000) -> dict:
    """Sequential local-header scan (works even when the central directory is missing).

    Members may be written with a data descriptor (sizes unknown up front); those are read by
    feeding chunks to a decompressor until the stream ends or the byte budget is exhausted.
    """
    info = dict(members_scanned=0, gene_axis_n=None, gene_symbols=[], features_member=None,
                experiment={}, scan='sequential_local_headers', scan_error='')
    consumed = 0
    try:
        with path.open('rb') as handle:
            while info['members_scanned'] < budget_members and consumed < budget_bytes:
                header = handle.read(30)
                if len(header) < 30 or header[:4] != b'PK\x03\x04':
                    break
                flag = int.from_bytes(header[6:8], 'little')
                method = int.from_bytes(header[8:10], 'little')
                comp_size = int.from_bytes(header[18:22], 'little')
                name_len = int.from_bytes(header[26:28], 'little')
                extra_len = int.from_bytes(header[28:30], 'little')
                name = handle.read(name_len).decode('utf-8', 'replace')
                handle.read(extra_len)
                info['members_scanned'] += 1
                is_dir = name.endswith('/')
                streamed = bool(flag & 0x08) and not is_dir
                wanted = name.endswith(('features.tsv.gz', 'gene_panel.json', 'experiment.xenium',
                                       'metrics_summary.csv', 'cells.csv.gz'))
                if is_dir:
                    continue
                if not wanted:
                    if streamed:
                        info['scan'] = 'stopped_at_unknown_size_member'
                        break
                    if handle.tell() + comp_size > 8_000_000_000:
                        info['scan'] = 'stopped_by_offset_budget'
                        break
                    handle.seek(comp_size, 1)
                    continue
                raw = b''
                remaining = budget_bytes - consumed
                if streamed:
                    decomp = zlib.decompressobj(-zlib.MAX_WBITS) if method == 8 else None
                    while len(raw) < 40_000_000 and remaining > 0:
                        chunk = handle.read(1_000_000)
                        if not chunk:
                            break
                        consumed += len(chunk)
                        remaining -= len(chunk)
                        if decomp is None:
                            raw += chunk
                        else:
                            try:
                                raw += decomp.decompress(chunk)
                            except zlib.error:
                                break
                            if decomp.eof:
                                break
                else:
                    raw = handle.read(min(comp_size, 120_000_000))
                    consumed += len(raw)
                    if method == 8:
                        raw = _inflate(raw)
                if name.endswith('gene_panel.json') and raw:
                    try:
                        panel = json.loads(raw.decode('utf-8', 'replace'))
                        genes = _collect_gene_symbols(panel)
                        hits = {k: len(set(genes) & v) for k, v in _panel_genes().items()}
                        info.update(gene_axis_n=len(genes), gene_symbols=genes[:5],
                                    features_member=name, panel_hits=hits,
                                    panel_name=panel.get('panel_name', ''),
                                    panel_n_targets=panel.get('n_targets', ''))
                    except Exception as exc:
                        info['scan_error'] = f'gene_panel:{type(exc).__name__}'
                elif name.endswith('features.tsv.gz') and raw:
                    text = raw.decode('utf-8', 'replace')
                    genes = [line.split('\t')[1] for line in text.splitlines() if '\t' in line]
                    hits = {k: len(set(genes) & v) for k, v in _panel_genes().items()}
                    info.update(gene_axis_n=len(genes), gene_symbols=genes[:5],
                                features_member=name, panel_hits=hits)
                elif name.endswith('experiment.xenium') and raw:
                    try:
                        info['experiment'] = json.loads(raw.decode('utf-8', 'replace'))
                    except Exception:
                        info['experiment'] = {}
                elif name.endswith('cells.csv.gz') and raw:
                    lines = raw.decode('utf-8', 'replace').splitlines()
                    if lines:
                        header = [c.strip().strip('"') for c in lines[0].split(',')]
                        info['cells_header'] = header[:12]
                        try:
                            xi = header.index('x_centroid'); yi = header.index('y_centroid')
                            rows_txt = [r.split(',') for r in lines[1:200_001]]
                            xs = np.array([float(r[xi]) for r in rows_txt])
                            ys = np.array([float(r[yi]) for r in rows_txt])
                            info['coord_range'] = dict(n_rows_scanned=int(xs.size),
                                                       x=[float(xs.min()), float(xs.max())],
                                                       y=[float(ys.min()), float(ys.max())])
                        except (ValueError, IndexError):
                            pass
                elif name.endswith('metrics_summary.csv') and raw:
                    info['metrics_summary'] = raw.decode('utf-8', 'replace')[:600]
                if info.get('gene_axis_n') and info.get('experiment') and info.get('coord_range'):
                    break
    except Exception as exc:  # truncated or unsupported archive
        info['scan_error'] = f'{type(exc).__name__}: {exc}'[:160]
    return info


def _collect_gene_symbols(node, depth: int = 0, out=None):
    """Collect every ``gene_name`` / ``name`` leaf from a Xenium gene_panel.json."""
    out = [] if out is None else out
    if depth > 8:
        return out
    if isinstance(node, dict):
        for key, value in node.items():
            if key in ('gene_name', 'gene', 'name') and isinstance(value, str) and value:
                out.append(value)
            elif isinstance(value, (dict, list)):
                _collect_gene_symbols(value, depth + 1, out)
    elif isinstance(node, list):
        for item in node:
            _collect_gene_symbols(item, depth + 1, out)
    return out


def _find_targets(node, depth: int = 0):
    """Locate the target list in a Xenium gene_panel.json (payload nested several levels)."""
    if depth > 4:
        return []
    if isinstance(node, list):
        if node and all(isinstance(x, (dict, str)) for x in node) and any(
                isinstance(x, dict) and (x.get('gene_name') or x.get('name')) for x in node):
            return node
        for item in node:
            found = _find_targets(item, depth + 1)
            if found:
                return found
    elif isinstance(node, dict):
        for key in ('targets', 'genes', 'panel', 'payload'):
            if key in node:
                found = _find_targets(node[key], depth + 1)
                if found:
                    return found
    return []


def _scan_spatial_tar(path: Path, budget_bytes: int = 400_000_000, budget_seconds: float = 25.0
                      ) -> dict:
    info = dict(members=[], coord_range=None, scan='tar_stream')
    start = time.time()
    consumed = 0
    try:
        with tarfile.open(path, 'r|gz') as handle:
            for member in handle:
                info['members'].append(member.name)
                if member.isfile:
                    consumed += member.size
                    if member.name.endswith('tissue_positions.csv'):
                        raw = handle.extractfile(member).read(4_000_000)
                        lines = raw.decode('utf-8', 'replace').splitlines()
                        header = lines[0].split(',')
                        rows = [line.split(',') for line in lines[1:]]
                        if rows and len(rows[0]) >= 6:
                            cols = {name: i for i, name in enumerate(header)}
                            xs = [float(r[cols.get('pxl_col_in_fullres', 4)]) for r in rows if r]
                            ys = [float(r[cols.get('pxl_row_in_fullres', 5)]) for r in rows if r]
                            info['coord_range'] = dict(
                                n_rows=len(rows), x=(min(xs), max(xs)), y=(min(ys), max(ys)))
                if consumed > budget_bytes or time.time() - start > budget_seconds:
                    info['scan'] = 'tar_stream_truncated_by_budget'
                    break
    except Exception as exc:  # truncated archive
        info['scan_error'] = f'{type(exc).__name__}: {exc}'[:160]
    return info


def scan_assets(out: Path) -> list[dict]:
    rows: list[dict] = []
    for directory in sorted(ZIP_DATASETS.iterdir()):
        if not directory.is_dir():
            continue
        platform = ('xenium' if 'Xenium' in directory.name else
                    'visium_hd' if 'Visium_HD' in directory.name else None)
        if platform is None:
            continue
        archives = [p for p in sorted(directory.iterdir()) if p.name.endswith(('.zip', '.tar.gz'))]
        for archive in archives:
            row = dict(dataset=directory.name, platform=platform, archive=archive.name,
                       bytes=archive.stat().st_size)
            if archive.name.endswith('.zip'):
                try:
                    row['zip_has_central_directory'] = _zip_tail_has_eocd(archive)
                except Exception as exc:
                    row['zip_has_central_directory'] = f'ERROR:{type(exc).__name__}'
                info = _read_zip_feature_axis(archive)
                row.update({k: v for k, v in info.items() if k != 'experiment'})
                experiment = info.get('experiment') or {}
                for key in ('run_name', 'sample_id', 'panel_name', 'region_name', 'run_start_time',
                            'experiment_uuid'):
                    row[f'experiment_{key}'] = experiment.get(key, '')
                row['gene_axis'] = ('READ_SEQUENTIAL' if info.get('gene_axis_n')
                                    else 'NOT_READ')
            else:
                spatial = 'spatial' in archive.name
                info = _scan_spatial_tar(archive) if spatial else dict(members=[], scan='skip_big')
                row.update({k: v for k, v in info.items() if k != 'members'})
                row['n_members_listed'] = len(info.get('members', []))
                row['gene_axis'] = 'NOT_READ_BIG_ARCHIVE' if not spatial else 'NOT_IN_SPATIAL_TAR'
                row['coord_range'] = json.dumps(info.get('coord_range')) if info.get('coord_range') \
                    else ''
            rows.append(row)
    write_tsv(rows, out / 'local_asset_fit.tsv')
    return rows


# --------------------------------------------------------------------------- main
def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument('command', choices=['smoke', 'run', 'assets', 'all'])
    parser.add_argument('--out', default=str(OUT))
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {}
    if args.command in ('smoke', 'all'):
        payload['smoke'] = run_smoke()
        (out / 'smoke_checks.json').write_text(
            json.dumps(payload['smoke'], indent=2, ensure_ascii=False) + '\n')
        print(json.dumps(payload['smoke']['verdict'], indent=1))
    if args.command in ('run', 'all'):
        payload['run'] = run_real(out)
        summary = payload['run']['summary']
        print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != 'unit_delta_r2'}
                          for k, v in summary.items()}, indent=1)[:3000])
    if args.command in ('assets', 'all'):
        payload['assets'] = scan_assets(out)
        print(f'assets rows {len(payload["assets"])}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
