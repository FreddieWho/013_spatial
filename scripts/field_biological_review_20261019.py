#!/usr/bin/env python3
"""Descriptive biological review of measured fields and frozen four-target predictions.

Protocol : docs/FIELD_BIOLOGICAL_REVIEW_PROTOCOL_20261019.md
Contract : infra/field_biological_review_20261019/contract.json
Outputs  : infra/field_biological_review_20261019/

Reads only frozen predictions, frozen caches, frozen visible-gene matrices, expression,
gene identity, barcodes and tissue flags. Never reads human TLS labels, annotations or
incoming material (an allowlist plus an explicit denylist is installed before any data
read). CPU only. No training, no model selection, no biological pass threshold, no
p-value. Shuffle references are descriptive perturbations, not null models.

Commands
  prepare   verify frozen inputs, recompute measured support/readout, write measured tables
  run       per frozen split: compare stored predictions with the measured field
  summarize aggregate to unit/arm-level descriptive tables and write the report
  check     synthetic known-answer tests plus a small real-sample end-to-end check
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

import h5py
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from scipy.spatial import ConvexHull, cKDTree
from scipy.stats import rankdata

import field_discovery_20261018 as base
import field_followup_20261018 as follow

ROOT = base.ROOT
WORK = ROOT / 'infra/field_biological_review_20261019'
PRED = ROOT / 'infra/field_source_holdout_20261018'
FUP = ROOT / 'infra/field_followup_20261018'
IDENTITY = ROOT / 'infra/field_discovery_review_20261018/identity_review.tsv'
PROTOCOL = ROOT / 'docs/FIELD_BIOLOGICAL_REVIEW_PROTOCOL_20261019.md'
CACHE = base.WORK
PANEL_FILE = ROOT / 'infra/tls_multiscale_gobp_20261002/panel_controls.json'

QS_OBS = (0.75, 0.80, 0.90)
Q_SHUFFLE = 0.75
PATCH_QS = (0.75, 0.90)
GRID_QS = (0.95,)

# Paths that must never be opened by this workflow, on top of the allowlist.
DENIED_SUBSTRINGS = (
    '/annotations/', '/incoming/', 'tls_annotation_pack', 'tls_label_field_20261001',
    'tls_four_cancer_labels_20261008', 'zenodo_st_crc_cms/Pathology_SpotAnnotations',
    'hooknet', 'HookNet',
)
MODEL_CLASSES = {
    'gene_ridge': ('raw_gene_ridge', False, False),
    'gene_ridge_smooth_half': ('raw_gene_ridge_posthoc_neighbourhood_smoothing', True, True),
    'fixed_point_hgb': ('module_point_hgb', False, False),
    'fixed_local_hgb': ('module_local_hgb', True, False),
    'train_mean': ('constant_prior', False, False),
}


def qkey(q):
    return f'q{int(round(q * 100))}'


def pkey(q):
    return f'{int(round(q * 100))}'


def seed_for(seed, *parts):
    digest = hashlib.sha256('|'.join([str(seed)] + [str(p) for p in parts]).encode()).digest()
    return int.from_bytes(digest[:8], 'little')


def contract():
    return json.loads((WORK / 'contract.json').read_text())


def is_denied(path):
    text = str(path).replace('\\', '/')
    return any(pattern in text for pattern in DENIED_SUBSTRINGS)


def install_denylist_hook():
    def audit(event, args):
        if event != 'open' or not isinstance(args[0], (str, bytes, os.PathLike)):
            return
        path = os.fsdecode(args[0])
        if is_denied(path):
            raise PermissionError(f'Denied path pattern: {path}')
    sys.addaudithook(audit)


# --------------------------------------------------------------------------- setup
def verify_frozen(checks, rows, sids):
    """Verify frozen hashes of every reused input, restricted to the selected sections."""
    verified = Counter()
    for root, prefixes in [(CACHE, ('cache/',)), (PRED, ('genes/',)), (FUP, ('genes/',))]:
        freeze = json.loads((root / 'freeze_manifest.json').read_text())
        for name, entry in freeze['files'].items():
            if not any(name.startswith(p) for p in prefixes):
                continue
            if sids is not None and Path(name).name.split('_')[0] not in sids:
                continue
            if base.sha256(root / name) != entry['sha256']:
                raise ValueError(f'Frozen input changed: {root/name}')
            verified['frozen_hashes'] += 1
    for row in rows:
        if sids is not None and row['section_id'] not in sids:
            continue
        if base.sha256(base.HEST / f"st/{row['section_id']}.h5ad") != row['expression_sha256']:
            raise ValueError(f'Expression hash changed: {row["section_id"]}')
        verified['expression_hashes'] += 1
    if sids is None:
        for root in [PRED, FUP]:
            freeze = json.loads((root / 'freeze_manifest.json').read_text())
            for name, entry in freeze['files'].items():
                if name.startswith('genes/'):
                    continue
                if base.sha256(root / name) != entry['sha256']:
                    raise ValueError(f'Frozen input changed: {root/name}')
                verified['frozen_hashes'] += 1
    checks['verified'] = dict(verified)


def setup(verify='sample', sids=None):
    """Load the frozen manifest, install the guard and denylist, return the run context."""
    c = contract()
    rows = json.loads((PRED / 'input_manifest.json').read_text())
    byid = {r['section_id']: r for r in rows}
    defs = json.loads((PRED / 'definitions.json').read_text())
    vocab = json.loads((PRED / 'vocabulary.json').read_text())
    splits = json.loads((PRED / 'splits.json').read_text())
    identity = {r['section_id']: r for r in csv.DictReader(IDENTITY.open(), delimiter='\t')}
    confirmed_sections = [s for s, r in identity.items() if r['patient_confirmed_within_study'] == 'True']
    for row in rows:
        record = identity[row['section_id']]
        row['identity_source'] = str(IDENTITY.relative_to(ROOT))
        row['confirmed_unit'] = bool(record['patient_confirmed_within_study'] == 'True')
        if bool(record['patient_confirmed_within_study'] == 'True') != bool(row['patient_known']):
            raise ValueError(f'Identity flag disagreement: {row["section_id"]}')
    eligibility = {r['section_id']: r['status'] for r in
                   csv.DictReader((PRED / 'source_eligibility.tsv').open(), delimiter='\t')}
    for row in rows:
        row['source_eligibility'] = eligibility.get(row['section_id'], 'UNKNOWN')
    if set(identity) != set(byid):
        raise ValueError('Identity table and prediction manifest cover different sections')
    allowed = {Path(__file__), Path(base.__file__), Path(follow.__file__),
               WORK / 'contract.json', PANEL_FILE, IDENTITY, PROTOCOL}
    for root in [CACHE, FUP, PRED]:
        freeze = json.loads((root / 'freeze_manifest.json').read_text())
        allowed.add(root / 'freeze_manifest.json')
        allowed.update(root / p for p in freeze['files'])
    allowed.update(base.HEST / f"st/{sid}.h5ad" for sid in byid)
    guard = base.InputGuard(ROOT, allowed, WORK)
    guard.install()
    install_denylist_hook()
    base.reserve(c)
    checks = {'n_allowed_paths': len(allowed),
              'no_allowed_path_matches_denylist': True,
              'allowed_path_matching_denylist': []}
    for path in allowed:
        if is_denied(path):
            checks['no_allowed_path_matches_denylist'] = False
            checks['allowed_path_matching_denylist'].append(str(path))
    if not checks['no_allowed_path_matches_denylist']:
        raise ValueError(f'Denylist pattern present in allowlist: {checks["allowed_path_matching_denylist"]}')
    if verify in ('all', 'sample'):
        verify_frozen(checks, rows, None if verify == 'all' else (sids or sorted(byid)))
    checks['panel_file_sha256'] = base.sha256(PANEL_FILE)
    checks['identity_source_sha256'] = base.sha256(IDENTITY)
    checks['protocol_sha256'] = base.sha256(PROTOCOL)
    receipt = json.loads((CACHE / 'input_receipt_prepare.json').read_text())
    frozen_panel = receipt['inputs'].get(str(PANEL_FILE.relative_to(ROOT)), {}).get('sha256')
    checks['panel_file_matches_frozen_receipt'] = bool(frozen_panel == checks['panel_file_sha256'])
    if not checks['panel_file_matches_frozen_receipt']:
        raise ValueError('Panel file hash differs from the frozen discovery receipt')
    panels = json.loads(PANEL_FILE.read_text())
    for target in c['targets']:
        for view in c['views']:
            assert_visible_target_disjoint(vocab[target][view], panels[target])
            for module, members in defs[target][view].items():
                if set(members) & set(panels[target]):
                    raise ValueError(f'Target genes inside a module: {target}/{view}/{module}')
    checks['visible_target_disjoint'] = True
    checks['module_target_disjoint'] = True
    fup_modules = json.loads((FUP / 'module_definitions.json').read_text())
    for target in c['targets']:
        if fup_modules[target] != defs[target]['full15']:
            raise ValueError(f'Frozen module definitions differ from the followup definitions: {target}')
    checks['module_definitions_match_followup'] = True
    checks['confirmed_sections'] = len(confirmed_sections)
    checks['confirmed_units'] = len({identity[s]['unit_id'] for s in confirmed_sections})
    if (checks['confirmed_sections'], checks['confirmed_units']) != (14, 8):
        raise ValueError('Frozen identity table no longer reports 14 confirmed sections / 8 units')
    checks['source_not_testable_sections'] = sum(
        1 for r in rows if r['source_eligibility'].startswith('NOT_TESTABLE'))
    checks['train_indices_excluded'] = True
    free = shutil.disk_usage(WORK).free
    checks['free_bytes'] = free
    checks['storage_reserve_ok'] = bool(free >= c['min_free_bytes'] + c['estimated_max_new_bytes'])
    if not checks['storage_reserve_ok']:
        raise RuntimeError(f'Storage reserve violated: free={free}')
    return dict(c=c, rows=rows, byid=byid, defs=defs, vocab=vocab, splits=splits,
                guard=guard, panels=panels, identity=identity, checks=checks)


def assert_visible_target_disjoint(visible, target_genes):
    shared = set(visible) & set(target_genes)
    if shared:
        raise ValueError(f'Target genes entered the visible vocabulary: {sorted(shared)[:5]}')


def assert_aligned(barcodes, expected):
    if not np.array_equal(np.asarray(barcodes), np.asarray(expected)):
        raise ValueError('Barcode order mismatch')


# --------------------------------------------------------------------------- metric core
def na(notes, out, name, reason):
    notes.append(f'{name}:{reason}')
    out[name] = None
    return None


def corr(a, b):
    a = np.asarray(a, np.float64)
    b = np.asarray(b, np.float64)
    a = a - a.mean()
    b = b - b.mean()
    den = math.sqrt(float(a @ a) * float(b @ b))
    if not np.isfinite(den) or den <= 0:
        return None
    return float(a @ b / den)


def jaccard(mask_a, mask_b):
    union = int(np.count_nonzero(mask_a | mask_b))
    if union == 0:
        return None
    return int(np.count_nonzero(mask_a & mask_b)) / union


def quantile_summary(values):
    values = [v for v in values if v is not None and np.isfinite(v)]
    if not values:
        return (None, None, None)
    arr = np.asarray(values, np.float64)
    return float(np.median(arr)), float(np.quantile(arr, 0.1)), float(np.quantile(arr, 0.9))


def sd_band(sd, floor):
    if sd is None:
        return 'not_estimable'
    if sd == 0:
        return 'zero'
    if sd < floor:
        return 'below_floor'
    return 'measurable'


def variance_reason(sd, floor):
    if sd == 0:
        return 'constant_field'
    if sd < floor:
        return 'degenerate_variance'
    return None


class Geo:
    """Per-section geometry: pitch, frozen graphs, correlogram bands, convex hull."""

    def __init__(self, xy, pitch, near, outer, band_edges, near_graph_file):
        self.xy = np.asarray(xy, np.float64)
        self.n = len(self.xy)
        self.pitch = float(pitch)
        self.near = near
        self.outer = outer
        self.near_graph_file = near_graph_file
        self.center = self.xy.mean(axis=0)
        degrees = np.diff(near.indptr)
        self.degree_min = int(degrees.min())
        self.degree_median = float(np.median(degrees))
        self.degree_max = int(degrees.max())
        self.isolated = int(np.count_nonzero(degrees == 0))
        outer_degrees = np.diff(outer.indptr)
        self.outer_isolated = int(np.count_nonzero(outer_degrees == 0))
        tri = sparse.triu(near, k=1).tocoo()
        self.edges = np.column_stack([tri.row, tri.col]).astype(np.int64).reshape(-1, 2)
        tree = cKDTree(self.xy)
        pairs = tree.query_pairs(band_edges[-1][1] * self.pitch, output_type='ndarray')
        if len(pairs):
            dd = np.linalg.norm(self.xy[pairs[:, 0]] - self.xy[pairs[:, 1]], axis=1) / self.pitch
        else:
            dd = np.zeros(0)
        self.bands = []
        for lo, hi in band_edges:
            sel = pairs[(dd >= lo) & (dd < hi)] if len(pairs) else pairs
            self.bands.append((sel, float(np.mean(dd[(dd >= lo) & (dd < hi)])) if len(sel) else None))
        self.hull = self._hull()

    def _hull(self):
        if self.n < 3:
            return None
        try:
            hull = ConvexHull(self.xy)
        except Exception:
            return None
        area = float(hull.volume)
        if not np.isfinite(area) or area <= 0:
            return None
        return 0.5 * math.sqrt(area), self.xy[hull.vertices]

    def hull_radial(self, point):
        if self.hull is None:
            return None
        scale, pts = self.hull
        p0 = pts
        p1 = np.roll(pts, -1, axis=0)
        d = p1 - p0
        denom = np.maximum((d * d).sum(axis=1), 1e-12)
        t = np.clip(((point - p0) * d).sum(axis=1) / denom, 0.0, 1.0)
        proj = p0 + t[:, None] * d
        return float(np.linalg.norm(point - proj, axis=1).min() / scale)


def neighbor_mean(geo, values):
    """Exactly the frozen neighbourhood rule: isolated spots keep their own value."""
    return base.neighbor_mean(geo.near, values)


def graph_neighbor_mean(graph, values):
    """Frozen neighbourhood rule for any of the frozen graphs (near or outer)."""
    return base.neighbor_mean(graph, values)


def weighted_centroid(xy, weights):
    total = float(weights.sum())
    if total <= 0 or not np.isfinite(total):
        return None
    return (weights[:, None] * xy).sum(axis=0) / total


def patch_stats(mask, geo, radii=True):
    """Connected components of the thresholded spots on the frozen near graph."""
    sel = np.flatnonzero(mask)
    m = len(sel)
    if m == 0:
        return None
    idx = np.full(geo.n, -1, np.int64)
    idx[sel] = np.arange(m)
    if len(geo.edges):
        keep = (idx[geo.edges[:, 0]] >= 0) & (idx[geo.edges[:, 1]] >= 0)
        a = idx[geo.edges[keep, 0]]
        b = idx[geo.edges[keep, 1]]
    else:
        a = b = np.zeros(0, np.int64)
    if len(a) == 0:
        ncomp = m
        labels = np.arange(m)
    else:
        graph = sparse.csr_matrix((np.ones(len(a), np.int8), (a, b)), shape=(m, m))
        ncomp, labels = connected_components(graph, directed=False)
    sizes = np.bincount(labels, minlength=ncomp)
    out = dict(n_patches=int(ncomp), n_high=int(m), area_frac=m / geo.n,
               largest_share=float(sizes.max() / m))
    if not radii:
        out['largest_radius_pitch'] = None
        out['mean_radius_pitch'] = None
        return out
    xs = geo.xy[sel]
    cx = np.bincount(labels, weights=xs[:, 0], minlength=ncomp) / sizes
    cy = np.bincount(labels, weights=xs[:, 1], minlength=ncomp) / sizes
    dist = np.hypot(xs[:, 0] - cx[labels], xs[:, 1] - cy[labels]) / geo.pitch
    order = np.argsort(labels, kind='stable')
    starts = np.searchsorted(labels[order], np.arange(ncomp))
    radius = np.maximum.reduceat(dist[order], starts)
    out['largest_radius_pitch'] = float(radius.max())
    out['mean_radius_pitch'] = float(radius.mean())
    return out


def band_autocorr(pairs, dv, var):
    if pairs is None or len(pairs) == 0 or var <= 0:
        return None
    i = pairs[:, 0]
    j = pairs[:, 1]
    return float(dv[i] @ dv[j] / (len(i) * var))


def decay_half_pitch(values, distances):
    if not values or values[0] is None or values[0] <= 0:
        return None
    half = values[0] / 2.0
    for k in range(1, len(values)):
        if values[k] is None or distances[k] is None or distances[k - 1] is None:
            return None
        if values[k] <= half:
            v0, v1 = values[k - 1], values[k]
            if v0 - v1 <= 0:
                return None
            frac = (v0 - half) / (v0 - v1)
            return float(distances[k - 1] + frac * (distances[k] - distances[k - 1]))
    return None


def shuffle_pass(v, geo, seed, n_perm):
    """Permute field values across spots: descriptive perturbation, marginals kept."""
    rng = np.random.default_rng(seed)
    threshold = float(np.quantile(v, Q_SHUFFLE))
    masks = np.zeros((n_perm, geo.n), bool)
    stats = []
    for k in range(n_perm):
        permuted = v[rng.permutation(geo.n)]
        masks[k] = permuted > threshold
        stats.append(patch_stats(masks[k], geo, radii=False))
    return masks, stats


def describe_field(v, geo, seed, notes, n_perm, do_shuffle, sd_floor, patch_qs=PATCH_QS):
    """Distribution, geometry, patch structure and (optionally) shuffle perturbation."""
    v = np.asarray(v, np.float64)
    out = {}
    mean = float(v.mean())
    sd = float(v.std(ddof=0))
    q = np.quantile(v, [0.05, 0.25, 0.5, 0.75, 0.95])
    out.update(dist_mean=mean, dist_sd=sd, dist_min=float(v.min()), dist_max=float(v.max()),
               dist_q05=float(q[0]), dist_q25=float(q[1]), dist_median=float(q[2]),
               dist_q75=float(q[3]), dist_q95=float(q[4]),
               dist_frac_zero=float(np.mean(v == 0)), dist_frac_negative=float(np.mean(v < 0)),
               sd_band=sd_band(sd, sd_floor))
    extra = dict(geometry_ok=sd > sd_floor, masks={}, perm_masks=None, shuffle_patch=None)
    reason = variance_reason(sd, sd_floor)
    if reason is None:
        out['dist_skew'] = float(np.mean(((v - mean) / sd) ** 3))
    else:
        na(notes, out, 'dist_skew', reason)
    if not extra['geometry_ok']:
        for name in ['geo_centroid_shift_pitch', 'geo_centroid_box_x', 'geo_centroid_box_y',
                     'geo_centroid_hull_radial', 'geo_autocorr_near', 'geo_autocorr_outer',
                     'geo_decay_half_pitch'] + [f'geo_patch{pkey(p)}_{s}' for p in patch_qs
                                                for s in ['n', 'largest_share', 'largest_radius_pitch',
                                                          'mean_radius_pitch', 'area_frac',
                                                          'at_threshold_frac']]:
            na(notes, out, name, reason)
        return out, extra, notes
    weights = v - v.min()
    centroid = weighted_centroid(geo.xy, weights)
    extra['centroid'] = centroid
    if geo.n < 10:
        na(notes, out, 'geo_centroid_shift_pitch', 'few_spots')
        na(notes, out, 'geo_centroid_hull_radial', 'few_spots')
        out['geo_centroid_box_x'] = None
        out['geo_centroid_box_y'] = None
    else:
        out['geo_centroid_shift_pitch'] = float(np.linalg.norm(centroid - geo.center) / geo.pitch)
        span = geo.xy.max(axis=0) - geo.xy.min(axis=0)
        out['geo_centroid_box_x'] = float((centroid[0] - geo.xy[:, 0].min()) / max(span[0], 1e-12))
        out['geo_centroid_box_y'] = float((centroid[1] - geo.xy[:, 1].min()) / max(span[1], 1e-12))
        radial = geo.hull_radial(centroid)
        if radial is None:
            na(notes, out, 'geo_centroid_hull_radial', 'degenerate_hull')
        else:
            out['geo_centroid_hull_radial'] = radial
    dv = v - neighbor_mean(geo, v)
    extra['dv'] = dv
    dv_energy = float(dv @ dv)
    if dv_energy <= 0:
        na(notes, out, 'geo_autocorr_near', 'constant_detail_component')
        na(notes, out, 'geo_autocorr_outer', 'constant_detail_component')
    else:
        out['geo_autocorr_near'] = float(dv @ (neighbor_mean(geo, dv)) / dv_energy)
        out['geo_autocorr_outer'] = float(dv @ (graph_neighbor_mean(geo.outer, dv)) / dv_energy)
    var = float(np.mean(dv ** 2))
    band_values = []
    band_distances = []
    for k, (pairs, dmean) in enumerate(geo.bands, start=1):
        value = band_autocorr(pairs, dv, var)
        out[f'geo_corr_band{k}'] = value
        out[f'geo_band{k}_mean_dist_pitch'] = dmean
        band_values.append(value)
        band_distances.append(dmean)
        if value is None:
            notes.append(f'geo_corr_band{k}:no_pairs_in_band')
    decay = decay_half_pitch(band_values, band_distances)
    if decay is None:
        reason = ('nonpositive_first_band' if (band_values and band_values[0] is not None
                                               and band_values[0] <= 0) else 'no_decay_within_band_limit')
        na(notes, out, 'geo_decay_half_pitch', reason)
    else:
        out['geo_decay_half_pitch'] = decay
    if geo.n < 10:
        for p in patch_qs:
            for s in ['n', 'largest_share', 'largest_radius_pitch', 'mean_radius_pitch',
                      'area_frac', 'at_threshold_frac']:
                na(notes, out, f'geo_patch{pkey(p)}_{s}', 'few_spots')
        for p in QS_OBS:
            extra['masks'][p] = v > float(np.quantile(v, p))
    else:
        for p in sorted(set(QS_OBS) | set(patch_qs)):
            threshold = float(np.quantile(v, p))
            mask = v > threshold
            extra['masks'][p] = mask
            if p not in patch_qs:
                continue
            out[f'geo_patch{pkey(p)}_area_frac'] = float(mask.mean())
            out[f'geo_patch{pkey(p)}_at_threshold_frac'] = float(np.mean(v == threshold))
            if not mask.any():
                for s in ['n', 'largest_share', 'largest_radius_pitch', 'mean_radius_pitch']:
                    na(notes, out, f'geo_patch{pkey(p)}_{s}', 'empty_threshold_set')
                continue
            stats = patch_stats(mask, geo)
            for name, value in [('n', stats['n_patches']), ('largest_share', stats['largest_share']),
                                ('largest_radius_pitch', stats['largest_radius_pitch']),
                                ('mean_radius_pitch', stats['mean_radius_pitch'])]:
                out[f'geo_patch{pkey(p)}_{name}'] = value
    if do_shuffle:
        perm_masks, patch_quantiles = shuffle_pass(v, geo, seed, n_perm)
        extra['perm_masks'] = perm_masks
        extra['shuffle_patch'] = [s for s in patch_quantiles if s is not None]
        for name in ['n_patches', 'largest_share']:
            med, lo, hi = quantile_summary([p[name] for p in extra['shuffle_patch']])
            suffix = 'n' if name == 'n_patches' else 'largest_share'
            out[f'geo_patch{pkey(Q_SHUFFLE)}_{suffix}_shuffled_median'] = med
            out[f'geo_patch{pkey(Q_SHUFFLE)}_{suffix}_shuffled_p10'] = lo
            out[f'geo_patch{pkey(Q_SHUFFLE)}_{suffix}_shuffled_p90'] = hi
    return out, extra, notes


def shuffled_jaccard(perm_masks, reference_mask):
    if perm_masks is None or reference_mask is None:
        return (None, None, None)
    return quantile_summary([jaccard(m, reference_mask) for m in perm_masks])


# --------------------------------------------------------------------------- relations
def module_matrix(sid, ti, view, defs, vocab, target):
    root = PRED if view == 'small5' else FUP
    z = np.load(root / f'genes/{sid}_t{ti}.npy', mmap_mode='r')
    names = list(defs[target][view].keys())
    index = {g: i for i, g in enumerate(vocab[target][view])}
    columns = [np.array([index[g] for g in defs[target][view][name]], dtype=np.int64) for name in names]
    scores = np.column_stack([np.asarray(z[:, c].mean(axis=1), np.float64) for c in columns])
    counts = {name: len(c) for name, c in zip(names, columns)}
    return names, counts, scores


def relation_rows(field, field_extra, names, counts, scores, geo, sd_floor):
    modules = []
    for j, name in enumerate(names):
        m = scores[:, j]
        sd = float(m.std(ddof=0))
        centroid = weighted_centroid(geo.xy, m - m.min()) if sd > sd_floor else None
        autocorr = None
        if sd > sd_floor:
            dv = m - neighbor_mean(geo, m)
            if float(dv @ dv) > 0:
                autocorr = float(dv @ neighbor_mean(geo, dv) / (dv @ dv))
        modules.append(dict(name=name, sd=sd, centroid=centroid, autocorr=autocorr,
                            mean=float(m.mean()), n_genes=counts[name],
                            mask={q: m > float(np.quantile(m, q)) for q in QS_OBS}))
    rows = []
    for j, name in enumerate(names):
        info = modules[j]
        m = scores[:, j]
        usable = field_extra['geometry_ok'] and info['sd'] > sd_floor
        row = dict(rel_module=name, rel_module_n_genes=info['n_genes'],
                   rel_module_mean=info['mean'], rel_module_sd=info['sd'],
                   rel_module_autocorr_near=info['autocorr'],
                   module_is_view_input=True, module_genes_intersect_target=False,
                   rel_pearson=corr(field, m) if usable else None,
                   rel_spearman=(corr(rankdata(field), rankdata(m)) if usable else None))
        for q in QS_OBS:
            row[f'rel_jaccard_{qkey(q)}'] = (jaccard(field_extra['masks'][q], info['mask'][q])
                                             if usable else None)
        if field_extra.get('perm_masks') is not None and info['sd'] > sd_floor:
            med, lo, hi = shuffled_jaccard(field_extra['perm_masks'], info['mask'][Q_SHUFFLE])
        else:
            med = lo = hi = None
        row['rel_jaccard_q75_shuffled_median'] = med
        row['rel_jaccard_q75_shuffled_p10'] = lo
        row['rel_jaccard_q75_shuffled_p90'] = hi
        row['rel_centroid_distance_pitch'] = (
            float(np.linalg.norm(field_extra['centroid'] - info['centroid']) / geo.pitch)
            if (usable and info['centroid'] is not None) else None)
        rows.append(row)
    return rows


def profile_agreement(measured_rows, predicted_rows):
    a = np.array([r['rel_pearson'] if r['rel_pearson'] is not None else np.nan
                  for r in measured_rows], float)
    b = np.array([r['rel_pearson'] if r['rel_pearson'] is not None else np.nan
                  for r in predicted_rows], float)
    ok = np.isfinite(a) & np.isfinite(b)
    out = dict(prof_relation_pearson=None, prof_relation_l2=None, prof_abs_delta_mean=None,
               prof_top_module_agreement=None, prof_n_modules=int(ok.sum()))
    if ok.sum() < 3:
        return out
    av, bv = a[ok], b[ok]
    out['prof_relation_pearson'] = corr(av, bv)
    out['prof_relation_l2'] = float(np.linalg.norm(av - bv))
    out['prof_abs_delta_mean'] = float(np.mean(np.abs(av - bv)))
    top_a = int(np.argmax(np.abs(av)))
    top_b = int(np.argmax(np.abs(bv)))
    tied = int(np.sum(np.abs(np.abs(av) - np.abs(av[top_a])) < 1e-12)) > 1
    out['prof_top_module_agreement'] = None if tied else int(top_a == top_b)
    return out


# --------------------------------------------------------------------------- recovery
def recovery_metrics(y, p, geo, yd, pd, notes, sd_floor):
    y = np.asarray(y, np.float64)
    p = np.asarray(p, np.float64)
    out = dict(rec_rmse=float(np.sqrt(np.mean((y - p) ** 2))),
               rec_centered_rmse=float(np.sqrt(np.mean(((y - y.mean()) - (p - p.mean())) ** 2))),
               rec_mae=float(np.mean(np.abs(y - p))), rec_bias=float(p.mean() - y.mean()))
    ysd = float(y.std(ddof=0))
    psd = float(p.std(ddof=0))
    if ysd > sd_floor:
        out['rec_sd_ratio'] = float(psd / ysd)
    else:
        na(notes, out, 'rec_sd_ratio', variance_reason(ysd, sd_floor) or 'degenerate_variance')
    if psd > sd_floor:
        out['rec_calibration_slope'] = float(np.mean((p - p.mean()) * (y - y.mean())) / (psd ** 2))
    else:
        na(notes, out, 'rec_calibration_slope', variance_reason(psd, sd_floor) or 'degenerate_variance')
    if yd['geometry_ok'] and pd['geometry_ok']:
        out['rec_pearson'] = corr(y, p)
        out['rec_spearman'] = corr(rankdata(y), rankdata(p))
    else:
        reason = variance_reason(psd, sd_floor) or variance_reason(ysd, sd_floor)
        na(notes, out, 'rec_pearson', reason)
        na(notes, out, 'rec_spearman', reason)
    if not yd['geometry_ok'] or not pd['geometry_ok']:
        na(notes, out, 'rec_local_pattern_pearson',
           variance_reason(psd if not pd['geometry_ok'] else ysd, sd_floor))
    elif float(np.std(pd['dv'])) <= 0:
        na(notes, out, 'rec_local_pattern_pearson', 'constant_detail_component')
    else:
        out['rec_local_pattern_pearson'] = corr(yd['dv'], pd['dv'])
    for q in QS_OBS:
        out[f'rec_jaccard_{qkey(q)}'] = (jaccard(yd['masks'][q], pd['masks'][q])
                                         if (yd['geometry_ok'] and pd['geometry_ok']) else None)
    if pd.get('perm_masks') is not None and yd['geometry_ok']:
        med, lo, hi = shuffled_jaccard(pd['perm_masks'], yd['masks'][Q_SHUFFLE])
    else:
        med = lo = hi = None
    out['rec_jaccard_q75_shuffled_median'] = med
    out['rec_jaccard_q75_shuffled_p10'] = lo
    out['rec_jaccard_q75_shuffled_p90'] = hi
    observed = out.get('rec_jaccard_q75')
    out['rec_jaccard_q75_excess'] = (observed - med) if (observed is not None and med is not None) else None
    if yd['geometry_ok'] and pd['geometry_ok']:
        out['rec_centroid_displacement_pitch'] = float(
            np.linalg.norm(yd['centroid'] - pd['centroid']) / geo.pitch)
    else:
        na(notes, out, 'rec_centroid_displacement_pitch',
           variance_reason(psd, sd_floor) or variance_reason(ysd, sd_floor))
    return out


def delta(out, notes, name, a, b, reason='missing_side'):
    if a is None or b is None:
        na(notes, out, name, reason)
    else:
        out[name] = float(a - b)


# --------------------------------------------------------------------------- data rows
def geo_of(row, cache, c):
    xy = np.asarray(cache['xy_source_pixels'], np.float64)
    pitch = float(row['empirical_pitch_source_pixels'])
    tree = cKDTree(xy)
    distances = tree.query(xy, k=2)[0][:, 1]
    if not np.all(distances > 0):
        raise ValueError(f'Duplicate spot coordinates in {row["section_id"]}')
    geo = Geo(xy, pitch, cache['near'], cache['outer'],
              [tuple(b) for b in c['correlogram_band_edges_pitch']],
              str((CACHE / f'cache/{row["section_id"]}_near.npz').relative_to(ROOT)))
    geo.pitch_recomputed = float(np.median(distances))
    return geo, geo.pitch_recomputed


def support_metrics(ctx, info, cache):
    """Molecular detection coverage and multi-gene support for one section x target."""
    target = info['target']
    with h5py.File(base.HEST / f"st/{info['section_id']}.h5ad") as f:
        genes = base.strings(f['var/_index'][:])
        barcodes = base.strings(f['obs/_index'][:])
        keep = np.flatnonzero(f['obs/in_tissue'][:] == 1)
        x = base.read_matrix(f['X'], len(barcodes), len(genes))[keep]
    assert_aligned(barcodes[keep], cache['barcodes'])
    members = ctx['panels'][target]
    counts = x[:, np.isin(genes, list(members))].tocsr()
    library = np.asarray(x.sum(axis=1)).ravel()
    per_gene = np.asarray((counts > 0).sum(axis=0)).ravel()
    per_spot = np.asarray((counts > 0).sum(axis=1)).ravel()
    det_frac = per_gene / max(len(keep), 1)
    share = np.asarray(counts.sum(axis=1)).ravel() / np.maximum(library, 1.0)
    readout = base.target_readout(x, genes, members)
    n_genes = counts.shape[1]
    n_detected = int(np.count_nonzero(per_gene))
    out = dict(sup_genes_defined=len(members), sup_genes_in_panel=int(n_genes),
               sup_genes_absent_from_panel=int(len(members) - n_genes),
               sup_genes_detected=n_detected,
               sup_genes_in_panel_zero_counts=int(n_genes - n_detected),
               sup_spots_with_any_gene_frac=float(np.mean(per_spot > 0)),
               sup_genes_per_spot_median=float(np.median(per_spot)),
               sup_genes_per_spot_q90=float(np.quantile(per_spot, 0.9)),
               sup_target_count_share_median=float(np.median(share)),
               sup_readout_source='X raw counts; mean log1p(1e4 * target counts / all-gene library)',
               sup_target_gene_source=str(PANEL_FILE.relative_to(ROOT)))
    if n_genes:
        out.update(sup_det_frac_min=float(det_frac.min()),
                   sup_det_frac_q25=float(np.quantile(det_frac, 0.25)),
                   sup_det_frac_median=float(np.median(det_frac)),
                   sup_det_frac_q75=float(np.quantile(det_frac, 0.75)),
                   sup_det_frac_max=float(det_frac.max()))
        weights = np.asarray(counts.sum(axis=0)).ravel().astype(np.float64)
        total = float(weights.sum())
        if total > 0:
            ordered = np.sort(weights)[::-1]
            out['sup_effective_gene_count'] = float(total ** 2 / float((weights ** 2).sum()))
            out['sup_top_gene_share'] = float(ordered[0] / total)
            out['sup_top3_gene_share'] = float(ordered[:3].sum() / total)
        else:
            out.update(sup_effective_gene_count=None, sup_top_gene_share=None, sup_top3_gene_share=None)
    else:
        out.update(sup_det_frac_min=None, sup_det_frac_q25=None, sup_det_frac_median=None,
                   sup_det_frac_q75=None, sup_det_frac_max=None, sup_effective_gene_count=None,
                   sup_top_gene_share=None, sup_top3_gene_share=None)
    detected_spots = per_spot > 0
    out['sup_mean_if_detected'] = (float(readout[detected_spots].mean()) if detected_spots.any() else None)
    out['sup_sd_if_detected'] = (float(readout[detected_spots].std(ddof=0)) if detected_spots.any() else None)
    out['sup_readout_max_abs_diff'] = float(np.max(np.abs(
        readout.astype(np.float64) - cache['y'][:, info['ti']].astype(np.float64))))
    return out


def measured_descriptors(ctx, info, cache, geo, do_shuffle=False):
    c = ctx['c']
    y = cache['y'][:, info['ti']].astype(np.float64)
    notes = []
    desc, extra, notes = describe_field(
        y, geo, seed_for(c['seed'], info['section_id'], info['target'], 'measured'), notes,
        c['shuffle_permutations'], do_shuffle=do_shuffle, sd_floor=c['sd_floor'])
    return y, desc, extra, notes


def identity_row(info, geo, cache):
    return dict(section_id=info['section_id'], unit_id=info['unit_id'], organ=info['organ'],
                source_key=info['source_key'], platform=info['metadata']['st_technology'],
                patient_known=info['patient_known'], confirmed_unit=info['confirmed_unit'],
                identity_source=info['identity_source'],
                source_eligibility=info['source_eligibility'], n_spots=int(len(cache['y'])),
                pitch_source_pixels=geo.pitch, pitch_recomputed=geo.pitch_recomputed,
                empirical_pitch_manifest=info['empirical_pitch_source_pixels'],
                xy_source_field='cache/{sid}.npz:xy_source_pixels',
                near_graph_file=geo.near_graph_file,
                geo_near_degree_min=geo.degree_min, geo_near_degree_median=geo.degree_median,
                geo_near_degree_max=geo.degree_max, geo_near_isolated_spots=geo.isolated,
                geo_outer_isolated_spots=geo.outer_isolated,
                target=info['target'])


def prediction_rows(ctx, split, info, cache, geo, measured_desc, measured_extra, y):
    """All views x models for one held-out section: recovery, geometry and relations."""
    c = ctx['c']
    ti = info['ti']
    organ_extrapolation = bool(split['arm'] == 'source'
                               and split.get('same_organ_train_sections', 1) == 0)
    rows = []
    relation_rows_out = []
    for view in c['views']:
        npz_path = PRED / split['name'] / f'{info["section_id"]}_t{ti}_{view}.npz'
        with np.load(npz_path, allow_pickle=False) as z:
            assert_aligned(z['barcodes'], cache['barcodes'])
            if not np.array_equal(z['y'], cache['y'][:, ti]):
                raise ValueError(f'Prediction y differs from cache: {npz_path}')
            predictions = {name: z[name].astype(np.float64) for name in c['models']}
        names, counts, scores = module_matrix(info['section_id'], ti, view, ctx['defs'],
                                              ctx['vocab'], info['target'])
        measured_relations = relation_rows(y, measured_extra, names, counts, scores, geo, c['sd_floor'])
        for model in c['models']:
            notes = []
            p = predictions[model]
            desc, extra, notes = describe_field(
                p, geo, seed_for(c['seed'], info['section_id'], info['target'], view, model, 'prediction'),
                notes, c['shuffle_permutations'], do_shuffle=True, sd_floor=c['sd_floor'])
            rec = recovery_metrics(y, p, geo, measured_extra, extra, notes, c['sd_floor'])
            delta(rec, notes, 'rec_autocorr_near_delta', desc.get('geo_autocorr_near'),
                  measured_desc.get('geo_autocorr_near'))
            delta(rec, notes, 'rec_autocorr_outer_delta', desc.get('geo_autocorr_outer'),
                  measured_desc.get('geo_autocorr_outer'))
            delta(rec, notes, 'rec_decay_half_pitch_delta', desc.get('geo_decay_half_pitch'),
                  measured_desc.get('geo_decay_half_pitch'))
            for p_q in PATCH_QS:
                for name in ['n', 'largest_share', 'largest_radius_pitch']:
                    delta(rec, notes, f'rec_patch{pkey(p_q)}_{name}_delta',
                          desc.get(f'geo_patch{pkey(p_q)}_{name}'),
                          measured_desc.get(f'geo_patch{pkey(p_q)}_{name}'))
            predicted_relations = relation_rows(p, extra, names, counts, scores, geo, c['sd_floor'])
            model_class, uses_neighbourhood, is_posthoc = MODEL_CLASSES[model]
            row = identity_row(info, geo, cache)
            row.update(split=split['name'], arm=split['arm'], view=view, model=model,
                       model_class=model_class, uses_neighbourhood_input=uses_neighbourhood,
                       is_posthoc_smoothing=is_posthoc, cross_budget_verdict='none',
                       cross_organ_extrapolation=organ_extrapolation,
                       evidence_role='exploratory_reanalysis_of_developed_predictions')
            row.update(desc)
            row.update(rec)
            row.update(profile_agreement(measured_relations, predicted_relations))
            for key in ['dist_sd', 'dist_frac_zero', 'dist_frac_negative', 'geo_autocorr_near',
                        'geo_autocorr_outer', 'geo_centroid_shift_pitch', 'geo_patch75_n',
                        'geo_patch75_largest_share', 'geo_decay_half_pitch']:
                row[f'meas_{key}'] = measured_desc.get(key)
            row['not_estimable'] = ';'.join(sorted(notes))
            row['n_not_estimable'] = len(notes)
            row['state'] = prediction_state(row, c['sd_floor'])
            rows.append(row)
            for a, b in zip(measured_relations, predicted_relations):
                relation_rows_out.append(dict(
                    split=split['name'], arm=split['arm'], section_id=info['section_id'],
                    unit_id=info['unit_id'], organ=info['organ'], source_key=info['source_key'],
                    platform=info['metadata']['st_technology'], patient_known=info['patient_known'],
                    confirmed_unit=info['confirmed_unit'], target=info['target'], view=view, model=model,
                    evidence_role='input_echo_relation_not_spatial_increment',
                    rel_module=a['rel_module'], module_is_view_input=True,
                    module_genes_intersect_target=False,
                    measured_rel_pearson=a['rel_pearson'], predicted_rel_pearson=b['rel_pearson'],
                    rel_pearson_delta=(None if a['rel_pearson'] is None or b['rel_pearson'] is None
                                       else b['rel_pearson'] - a['rel_pearson']),
                    measured_rel_spearman=a['rel_spearman'], predicted_rel_spearman=b['rel_spearman'],
                    measured_rel_jaccard_q75=a['rel_jaccard_q75'],
                    predicted_rel_jaccard_q75=b['rel_jaccard_q75'],
                    measured_rel_jaccard_q75_shuffled_median=a['rel_jaccard_q75_shuffled_median'],
                    predicted_rel_jaccard_q75_shuffled_median=b['rel_jaccard_q75_shuffled_median'],
                    measured_rel_centroid_distance_pitch=a['rel_centroid_distance_pitch'],
                    predicted_rel_centroid_distance_pitch=b['rel_centroid_distance_pitch'],
                    rel_module_n_genes=a['rel_module_n_genes'], rel_module_sd=a['rel_module_sd'],
                    rel_module_autocorr_near=a['rel_module_autocorr_near']))
    return rows, relation_rows_out


def diagnostic_flags(row, sup):
    row['diag_sparse_detection'] = bool((sup['sup_spots_with_any_gene_frac'] or 0) < 0.02)
    row['diag_few_effective_genes'] = bool((sup['sup_effective_gene_count'] or 0) < 2)
    row['diag_dominant_gene'] = bool((sup['sup_top_gene_share'] or 0) > 0.5)
    return row


def measured_state(row, sd_floor):
    """Computability-only label: never a biological pass/fail judgement."""
    if numeric(row.get('sup_genes_in_panel')) == 0:
        return 'NOT_TESTABLE_NO_TARGET_GENES'
    sd = numeric(row.get('dist_sd'))
    if sd is None or sd == 0:
        return 'NOT_TESTABLE_CONSTANT_READOUT'
    if sd < sd_floor:
        return 'NOT_TESTABLE_DEGENERATE_VARIANCE'
    return 'DESCRIBED' if row.get('not_estimable') in (None, '') else 'DESCRIBED_WITH_NOT_ESTIMABLE_METRICS'


def prediction_state(row, sd_floor):
    sd = numeric(row.get('dist_sd'))
    if sd is None or sd == 0:
        return 'NOT_ESTIMABLE_CONSTANT_PREDICTION'
    if sd < sd_floor:
        return 'NOT_ESTIMABLE_DEGENERATE_VARIANCE'
    return 'DESCRIBED' if row.get('not_estimable') in (None, '') else 'DESCRIBED_WITH_NOT_ESTIMABLE_METRICS'


# --------------------------------------------------------------------------- prepare
def prepare_section(job):
    ctx = setup(verify='none')
    c = ctx['c']
    row = ctx['byid'][job['sid']]
    cache = base.load_cache(row)
    geo, recomputed = geo_of(row, cache, c)
    start = time.time()
    measured = []
    relations = []
    integrity = {}
    for ti, target in enumerate(c['targets']):
        info = dict(row, ti=ti, target=target)
        sup = support_metrics(ctx, info, cache)
        y, desc, extra, notes = measured_descriptors(ctx, info, cache, geo, do_shuffle=True)
        mrow = identity_row(info, geo, cache)
        mrow.update(desc)
        mrow.update(sup)
        mrow['evidence_role'] = 'measurement_eligibility'
        diagnostic_flags(mrow, sup)
        for name in ['n_patches', 'largest_share']:
            med, lo, hi = quantile_summary([s[name] for s in extra['shuffle_patch']])
            suffix = 'n' if name == 'n_patches' else 'largest_share'
            mrow[f'geo_patch{pkey(Q_SHUFFLE)}_{suffix}_shuffled_median'] = med
            mrow[f'geo_patch{pkey(Q_SHUFFLE)}_{suffix}_shuffled_p10'] = lo
            mrow[f'geo_patch{pkey(Q_SHUFFLE)}_{suffix}_shuffled_p90'] = hi
        mrow['not_estimable'] = ';'.join(sorted(notes))
        mrow['n_not_estimable'] = len(notes)
        mrow['state'] = measured_state(mrow, c['sd_floor'])
        measured.append(mrow)
        integrity[f'{job["sid"]}|{target}'] = dict(
            readout_max_abs_diff=sup['sup_readout_max_abs_diff'],
            pitch_manifest=row['empirical_pitch_source_pixels'], pitch_recomputed=recomputed,
            isolated_spots_manifest=row['near_isolated_spots'],
            isolated_spots_recomputed=geo.isolated)
        for view in c['views']:
            names, counts, scores = module_matrix(job['sid'], ti, view, ctx['defs'],
                                                  ctx['vocab'], target)
            for r in relation_rows(y, extra, names, counts, scores, geo, c['sd_floor']):
                r.update(identity_row(info, geo, cache), view=view,
                         evidence_role='input_echo_relation_not_spatial_increment')
                relations.append(r)
    return dict(sid=job['sid'], measured=measured, relations=relations, integrity=integrity,
                seconds=time.time() - start, receipt=ctx['guard'].receipt())


def prepare():
    ctx = setup(verify='all')
    c = ctx['c']
    if (WORK / 'PREPARED').exists():
        raise RuntimeError('Prepared outputs exist; preserve this version')
    (WORK / 'checks').mkdir(parents=True, exist_ok=True)
    measured, relations, integrity, receipts, timings = [], [], {}, [], {}
    jobs = [dict(sid=s) for s in sorted(ctx['byid'])]
    with ProcessPoolExecutor(max_workers=c['workers'], mp_context=multiprocessing.get_context('spawn')) as pool:
        futures = [pool.submit(prepare_section, job) for job in jobs]
        for future in as_completed(futures):
            res = future.result()
            measured += res['measured']
            relations += res['relations']
            integrity.update(res['integrity'])
            receipts.append(res['receipt'])
            timings[res['sid']] = res['seconds']
            print(json.dumps(dict(prepared=res['sid'], seconds=round(res['seconds'], 2))), flush=True)
    measured.sort(key=lambda m: (m['section_id'], m['target']))
    relations.sort(key=lambda r: (r['section_id'], r['target'], r['view'], r['rel_module']))
    base.tsv_write(WORK / 'measured_fields.tsv', measured)
    base.tsv_write(WORK / 'measured_module_relations.tsv', relations)
    base.json_write(WORK / 'input_checks.json',
                    dict(ctx['checks'], integrity=integrity, timings_seconds=timings,
                         measured_rows=len(measured), relation_rows=len(relations),
                         human_tls_files_opened=[],
                         guard_scope='repository Python-open allowlist plus explicit HDF5 path/field '
                                     'allowlist; a denylist also rejects annotation/TLS/HookNet paths; '
                                     'not an OS sandbox'))
    base.json_write(WORK / 'prepare_receipt.json',
                    dict(receipts=receipts, code_sha256=base.sha256(__file__),
                         contract_sha256=base.sha256(WORK / 'contract.json')))
    (WORK / 'PREPARED').write_text('MEASURED_FIELDS_FROZEN_HUMAN_TLS_NOT_READ\n')
    print(json.dumps(dict(prepared=len(measured), cpu_seconds=round(sum(timings.values()), 1))), flush=True)


# --------------------------------------------------------------------------- run
def run_split(split):
    ctx = setup(verify='none')
    c = ctx['c']
    out = WORK / 'splits' / split['name']
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'COMPLETE').exists():
        return dict(split=split['name'], status='already_complete')
    start = time.time()
    rows_out = []
    relation_out = []
    for sid in split['test']:
        row = ctx['byid'][sid]
        cache = base.load_cache(row)
        geo, recomputed = geo_of(row, cache, c)
        if abs(recomputed - float(row['empirical_pitch_source_pixels'])) > 1e-9:
            raise ValueError(f'Pitch mismatch for {sid}')
        for ti, target in enumerate(c['targets']):
            info = dict(row, ti=ti, target=target)
            y, mdesc, mex, _ = measured_descriptors(ctx, info, cache, geo)
            prows, prels = prediction_rows(ctx, split, info, cache, geo, mdesc, mex, y)
            rows_out += prows
            relation_out += prels
        print(json.dumps(dict(split=split['name'], section=sid,
                              seconds=round(time.time() - start, 1))), flush=True)
    base.tsv_write(out / 'prediction_metrics.tsv', rows_out)
    base.tsv_write(out / 'prediction_module_relations.tsv', relation_out)
    base.json_write(out / 'split_receipt.json',
                    dict(name=split['name'], arm=split['arm'], n_sections=len(split['test']),
                         rows=len(rows_out), relation_rows=len(relation_out),
                         code_sha256=base.sha256(__file__),
                         contract_sha256=base.sha256(WORK / 'contract.json'),
                         receipt=ctx['guard'].receipt()))
    (out / 'COMPLETE').write_text('SPLIT_DONE\n')
    return dict(split=split['name'], seconds=round(time.time() - start, 1), rows=len(rows_out))


def run():
    ctx = setup(verify='none')
    c = ctx['c']
    if not (WORK / 'PREPARED').exists():
        raise RuntimeError('Run prepare first')
    with ProcessPoolExecutor(max_workers=c['workers'], mp_context=multiprocessing.get_context('spawn')) as pool:
        futures = [pool.submit(run_split, s) for s in ctx['splits']]
        for future in as_completed(futures):
            print(json.dumps(future.result()), flush=True)
    (WORK / 'RUN_COMPLETE').write_text('ALL_SPLITS_DONE_NO_HUMAN_TLS_NO_PASS_THRESHOLD\n')


# --------------------------------------------------------------------------- summarize
MACRO_METRICS = ['rec_rmse', 'rec_centered_rmse', 'rec_mae', 'rec_bias', 'rec_sd_ratio',
                 'rec_calibration_slope', 'rec_pearson', 'rec_spearman', 'rec_local_pattern_pearson',
                 'rec_jaccard_q75', 'rec_jaccard_q80', 'rec_jaccard_q90',
                 'rec_jaccard_q75_shuffled_median', 'rec_jaccard_q75_excess',
                 'rec_centroid_displacement_pitch', 'rec_autocorr_near_delta', 'rec_autocorr_outer_delta',
                 'rec_decay_half_pitch_delta', 'rec_patch75_n_delta', 'rec_patch75_largest_share_delta',
                 'rec_patch75_largest_radius_pitch_delta', 'rec_patch90_n_delta',
                 'rec_patch90_largest_share_delta', 'rec_patch90_largest_radius_pitch_delta',
                 'prof_relation_pearson', 'prof_relation_l2', 'prof_abs_delta_mean',
                 'prof_top_module_agreement', 'dist_sd', 'dist_frac_zero', 'dist_frac_negative',
                 'geo_autocorr_near', 'geo_autocorr_outer', 'geo_decay_half_pitch',
                 'geo_centroid_shift_pitch', 'geo_patch75_n', 'geo_patch75_largest_share']
UNIT_KEYS = ['split', 'arm', 'target', 'view', 'model', 'model_class', 'cross_organ_extrapolation',
             'unit_id', 'source_key', 'organ', 'platform', 'patient_known', 'confirmed_unit']
FORBIDDEN_WORDS = ('pass', 'success', 'significant', 'p_value', 'pvalue', 'truth', 'positive_')
MEASURED_SUMMARY_METRICS = ['sup_genes_detected', 'sup_spots_with_any_gene_frac',
                            'sup_effective_gene_count', 'sup_top_gene_share', 'dist_sd',
                            'dist_frac_zero', 'geo_autocorr_near', 'geo_patch75_n',
                            'geo_patch75_largest_share', 'geo_decay_half_pitch',
                            'geo_centroid_shift_pitch']


def numeric(value):
    if value in (None, '', 'None', 'nan', 'NaN'):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def check_column_names(columns):
    bad = [c for c in columns if any(word in c.lower() for word in FORBIDDEN_WORDS)]
    return bad


def aggregate_units(records):
    """Unit-level aggregation.

    A unit metric is estimable only when every section of that unit produces it; a unit with
    a partial metric keeps the metric as NOT_ESTIMABLE_UNIT (recorded in
    `partial_unit_metrics`) instead of silently up-weighting the remaining section.
    """
    groups = defaultdict(list)
    for r in records:
        groups[tuple(r[k] for k in UNIT_KEYS)].append(r)
    units = []
    for key, rows in groups.items():
        unit = dict(zip(UNIT_KEYS, key))
        sections = sorted({r['section_id'] for r in rows})
        unit['n_sections'] = len(sections)
        missing, partial = [], []
        for metric in MACRO_METRICS:
            per_section = {}
            for r in rows:
                value = numeric(r.get(metric))
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


def section_weighted_index(records, aggregate):
    grouped = defaultdict(list)
    for r in records:
        key = ((r['arm'], r['target'], r['view'], r['model']) if aggregate == 'arm'
               else (r['split'], r['arm'], r['target'], r['view'], r['model']))
        grouped[key].append(r)
    index = {}
    for key, rows in grouped.items():
        for metric in MACRO_METRICS:
            values = [v for v in (numeric(r.get(metric)) for r in rows) if v is not None]
            index[(key, metric)] = (float(np.mean(values)) if values else None, len(values))
    return index


def macro_rows(units, records, aggregate='split'):
    grouped = defaultdict(list)
    for u in units:
        key = ((u['arm'], u['target'], u['view'], u['model']) if aggregate == 'arm'
               else (u['split'], u['arm'], u['target'], u['view'], u['model']))
        grouped[key].append(u)
    index = section_weighted_index(records, aggregate)
    rows = []
    for key, group in sorted(grouped.items()):
        if aggregate == 'arm':
            split = 'ALL'
            arm, target, view, model = key
        else:
            split, arm, target, view, model = key
        arm_key = (arm, target, view, model)
        flags = {str(u['cross_organ_extrapolation']) for u in group}
        if aggregate == 'arm':
            label = ('true' if flags == {'True'} else 'false' if flags == {'False'} else 'mixed')
        else:
            label = group[0]['cross_organ_extrapolation']
        for metric in MACRO_METRICS:
            values = [u[metric] for u in group if u[metric] is not None]
            partial_units = sum(1 for u in group if 0 < u.get(f'{metric}__n', 0)
                                < max(u.get('n_sections', 1), 1))
            row = dict(aggregate=aggregate, split=split, arm=arm, target=target, view=view,
                       model=model, model_class=group[0]['model_class'],
                       cross_organ_extrapolation=label,
                       metric=metric, n_units=len(group), n_estimable=len(values),
                       n_units_partial=partial_units)
            if values:
                arr = np.asarray(values, float)
                row.update(mean=float(arr.mean()), median=float(np.median(arr)),
                           p10=float(np.quantile(arr, 0.1)), p25=float(np.quantile(arr, 0.25)),
                           p75=float(np.quantile(arr, 0.75)), p90=float(np.quantile(arr, 0.90)))
            else:
                row.update(mean=None, median=None, p10=None, p25=None, p75=None, p90=None)
            section_key = (arm_key if aggregate == 'arm' else key)
            value, n_sections = index.get((section_key, metric), (None, 0))
            row['section_weighted_mean'] = value
            row['n_sections_estimable'] = n_sections
            if arm == 'source':
                per_source = defaultdict(list)
                for u in group:
                    if u[metric] is not None:
                        per_source[u['source_key']].append(u[metric])
                row['source_balanced_mean'] = (float(np.mean([np.mean(v) for v in per_source.values()]))
                                               if per_source else None)
                row['n_sources_estimable'] = len(per_source)
            else:
                row['source_balanced_mean'] = None
                row['n_sources_estimable'] = None
            rows.append(row)
    return rows


def strata_rows(units):
    rows = []
    for stratum, key in [('source', 'source_key'), ('organ', 'organ'), ('platform', 'platform'),
                         ('patient_known', 'patient_known'),
                         ('cross_organ_extrapolation', 'cross_organ_extrapolation')]:
        groups = defaultdict(list)
        for u in units:
            groups[(u['arm'], u['target'], u['view'], u['model'], str(u[key]))].append(u)
        for (arm, target, view, model, value), group in sorted(groups.items(), key=lambda kv: str(kv[0])):
            for metric in MACRO_METRICS:
                values = [u[metric] for u in group if u[metric] is not None]
                rows.append(dict(stratum=stratum, value=str(value), arm=arm, target=target,
                                 view=view, model=model, metric=metric, n_units=len(group),
                                 n_estimable=len(values),
                                 mean=float(np.mean(values)) if values else None,
                                 median=float(np.median(values)) if values else None,
                                 reason='' if values else 'no_estimable_units_in_stratum'))
    return rows


def require_columns(row, columns, table):
    """Fail closed when a row does not carry the columns the aggregation needs."""
    missing = [c for c in columns if c not in row]
    if missing:
        raise KeyError(f'{table} row lacks required column(s) {missing}; '
                       f'available relation/profile columns: '
                       f'{[c for c in row if "rel_" in c or "prof_" in c]}')


def within_unit_rows(measured_relations, prediction_relations):
    """Unit-internal consistency of relation profiles, read from the real output columns.

    Measured tables carry `rel_pearson`; prediction tables carry `measured_rel_pearson`
    and `predicted_rel_pearson` (see relation_rows). Missing columns fail closed.
    """
    rows = []
    groups = defaultdict(list)
    for r in measured_relations:
        require_columns(r, ['rel_pearson', 'rel_module', 'unit_id', 'target', 'view'],
                        'measured_module_relations')
        groups[('measured', 'all', r['unit_id'], r['target'], r['view'], 'rel_pearson')].append(r)
    for r in prediction_relations:
        require_columns(r, ['measured_rel_pearson', 'predicted_rel_pearson', 'rel_module',
                            'unit_id', 'target', 'view'], 'prediction_module_relations')
        for column in ['measured_rel_pearson', 'predicted_rel_pearson']:
            groups[(r['model'], r['split'], r['unit_id'], r['target'], r['view'], column)].append(r)
    for (model, split, unit, target, view, column), group in sorted(groups.items()):
        sections = sorted({g['section_id'] for g in group})
        if len(sections) < 2:
            continue
        vectors = [np.array([numeric(g[column]) if numeric(g[column]) is not None else np.nan
                             for g in group if g['section_id'] == sid], float)
                   for sid in sections]
        pairs = [(i, j) for i in range(len(sections)) for j in range(i + 1, len(sections))]
        correlations, l2s = [], []
        for i, j in pairs:
            a, b = vectors[i], vectors[j]
            ok = np.isfinite(a) & np.isfinite(b)
            if ok.sum() < 3:
                continue
            correlations.append(corr(a[ok], b[ok]))
            l2s.append(float(np.linalg.norm(a[ok] - b[ok])))
        correlations = [c for c in correlations if c is not None]
        rows.append(dict(model=model, split=split, unit_id=unit, target=target, view=view,
                         profile_column=column, n_sections=len(sections), n_pairs=len(pairs),
                         n_pairs_estimable=len(correlations),
                         within_unit_relation_pearson_mean=(float(np.mean(correlations))
                                                            if correlations else None),
                         within_unit_relation_l2_mean=(float(np.mean(l2s)) if l2s else None),
                         interpretation='descriptive unit-internal consistency only; not an '
                                        'independent patient or replicate confirmation'))
    return rows


def measured_summary(measured):
    rows = []
    for stratum, key in [('all', None), ('organ', 'organ'), ('platform', 'platform'),
                         ('source', 'source_key'), ('source_eligibility', 'source_eligibility')]:
        groups = defaultdict(list)
        for m in measured:
            groups[(m['target'], 'all' if key is None else m[key])].append(m)
        for (target, value), group in sorted(groups.items(), key=lambda kv: str(kv[0])):
            row = dict(stratum=stratum, value=str(value), target=target, n_sections=len(group),
                       n_not_testable=sum(1 for m in group
                                          if str(m.get('state', '')).startswith('NOT_TESTABLE')))
            for metric in MEASURED_SUMMARY_METRICS:
                values = [v for v in (numeric(m.get(metric)) for m in group) if v is not None]
                row[f'{metric}_median'] = float(np.median(values)) if values else None
                row[f'{metric}_min'] = float(min(values)) if values else None
                row[f'{metric}_max'] = float(max(values)) if values else None
                row[f'{metric}_n'] = len(values)
            rows.append(row)
    return rows


def diagnostic_sensitivity(measured):
    rows = []
    specs = [('diag_sparse_detection', 'sup_spots_with_any_gene_frac', [0.01, 0.02, 0.05], 'lt'),
             ('diag_few_effective_genes', 'sup_effective_gene_count', [1.5, 2.0, 3.0], 'lt'),
             ('diag_dominant_gene', 'sup_top_gene_share', [0.4, 0.5, 0.75], 'gt')]
    for name, field, cutoffs, direction in specs:
        for cutoff in cutoffs:
            hits = 0
            total = 0
            for m in measured:
                value = numeric(m.get(field))
                if value is None:
                    continue
                total += 1
                hits += int(value < cutoff) if direction == 'lt' else int(value > cutoff)
            rows.append(dict(marker=name, field=field, cutoff=cutoff, direction=direction,
                             n_sections_targets=total, n_flagged=hits,
                             nature='diagnostic marker only; gates no metric, state or filter'))
    return rows


def state_counts(records):
    return dict(sorted(Counter(r.get('state', 'UNKNOWN') for r in records).items()))


def write_report(measured, records, macro, ctx):
    c = ctx['c']
    lines = ['# 实测场资格与冻结预测恢复：描述性清单', '',
             '只描述，不判定成功/失败，不设生物学通过线，不给出 p 值。',
             '长度单位一律为该切片经验点间距（pitch）；分位阈值是操作阈值，不是生物真值。',
             '三本账分开：测量资格（coverage/支持）、输入回声关系（与可见模块场）、',
             '已开发预测的探索性再分析（实测 vs 预测）。三者的列名与 evidence_role 不同，不共用“成功”一词。', '',
             '## 0. 语料与身份', '',
             f"- 主清单：来源留出 `{PRED.relative_to(ROOT)}` 的 1000 份冻结预测包（内部 72 片/65 单位、"
             f"来源 53 片/47 单位/6 来源）；`train_indices.npz` 不进入评分。",
             f"- 实测清单覆盖全部 {len(measured) // len(c['targets'])} 片 × {len(c['targets'])} 目标；"
             f"来源身份不足的 {ctx['checks'].get('source_not_testable_sections')} 片保留在内部描述，"
             "不进入来源分母。",
             f"- 研究内确认身份 {ctx['checks'].get('confirmed_sections')} 片 / "
             f"{ctx['checks'].get('confirmed_units')} 个单位；其余为未知单位，不升级为独立患者。",
             f"- 状态计数（实测）：{json.dumps(state_counts(measured), ensure_ascii=False)}",
             f"- 状态计数（预测）：{json.dumps(state_counts(records), ensure_ascii=False)}", '',
             '## 1. 分子检测覆盖（连续报告，不分档删除）', '']
    for target in c['targets']:
        rows = [m for m in measured if m['target'] == target]
        if not rows:
            lines.append(f"- {target}：本表无记录（列表为空，不补值）")
            continue
        any_frac = [v for v in (numeric(m['sup_spots_with_any_gene_frac']) for m in rows) if v is not None]
        eff = [v for v in (numeric(m['sup_effective_gene_count']) for m in rows) if v is not None]
        top = [v for v in (numeric(m['sup_top_gene_share']) for m in rows) if v is not None]
        genes = [v for v in (numeric(m['sup_genes_detected']) for m in rows) if v is not None]
        if not (any_frac and eff and top and genes):
            lines.append(f"- {target}：{len(rows)} 片，但覆盖指标不全（不补值；见 measured_fields.tsv）")
            continue
        lines.append(f"- {target}：{len(rows)} 片；检出基因数中位 {np.median(genes):.0f}；任意目标基因检出点占比 "
                     f"中位 {np.median(any_frac):.3f}（{min(any_frac):.3f}–{max(any_frac):.3f}）；"
                     f"有效基因数中位 {np.median(eff):.2f}（最小 {min(eff):.2f}）；首位基因占比中位 "
                     f"{np.median(top):.3f}（最大 {max(top):.3f}）")
    index = {(r['aggregate'], r['arm'], r['target'], r['view'], r['model'], r['metric']): r
             for r in macro}
    lines += ['', '完整 72 片 × 4 目标逐片数值见 `measured_fields.tsv`；分层中位/极值见 `measured_summary.tsv`。', '',
              '## 2. 恢复描述（单位宏平均，整臂；不跨信息预算判胜负）', '',
              '| 臂 | 目标 | 视图 | 模型 | 模型类别 | 单位 | RMSE | 去均值RMSE | 逐点相关 | 局部图样相关 | '
              'q75 Jaccard | q75 打乱中位 | 质心位移(pitch) | 自相关差 |',
              '|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    combos = sorted({(r['arm'], r['target'], r['view'], r['model'], r['model_class']) for r in macro
                     if r['aggregate'] == 'arm'})
    for arm, target, view, model, model_class in combos:
        def value(metric):
            r = index.get(('arm', arm, target, view, model, metric))
            if r is None or r['mean'] is None:
                return 'NA'
            return f"{r['mean']:.4f}"
        r = index.get(('arm', arm, target, view, model, 'rec_rmse'))
        lines.append(f"| {arm} | {target} | {view} | {model} | {model_class} | "
                     f"{r['n_units'] if r else 'NA'} | {value('rec_rmse')} | {value('rec_centered_rmse')} | "
                     f"{value('rec_pearson')} | {value('rec_local_pattern_pearson')} | "
                     f"{value('rec_jaccard_q75')} | {value('rec_jaccard_q75_shuffled_median')} | "
                     f"{value('rec_centroid_displacement_pitch')} | {value('rec_autocorr_near_delta')} |")
    lines += ['', '打乱中位是保持取值边际、打乱空间排列的描述性扰动参考；不保持空间支撑与密度，未校准，',
              '不构成零模型、不产生 p 值、不用于筛选关系。',
              '`gene_ridge_smooth_half` 的自相关上升是构造结果（后处理平滑），不得作为空间关系或恢复证据；',
              '本阶段不做同点/邻域基线对照，也不在信息预算之间判胜负。', '',
              '## 3. 边界与未运行', '',
              '- 分位阈值（q75/q80/q90）是操作阈值，取严格大于分位值；逐场报告 area_frac 与 at_threshold_frac，',
              '  不同场的点集大小可以不同。',
              '- 预测场是可见模块/基因的函数；预测与可见模块的相关部分由输入结构决定，不能读作独立生物学证据。',
              '- 相关、自相关与尺度在任一方方差低于合同 sd_floor 时为 NOT_ESTIMABLE，不填 0、不填 1；RMSE 仍保留。',
              '- 数学上可算不等于有足够生物学参照：检出覆盖、单基因主导、有效基因数、边界/孤立点支撑与未确认的物理尺度',
              '  仍会限制解释；这些证据逐项并列，不合成单一分数，也不用统一百分比删除切片/目标/行。',
              '- 不得把"可计算"格子升级为已具生物学评价资格；NOT_TESTABLE/NOT_ESTIMABLE 也不读作场不存在或恢复失败。',
              '- 生物学可用性由主代理结合上述连续证据判断，本清单只提供描述性事实与不可评价状态。',
              '- 未训练/未重训、未做独立 marker 隔离、未读取人工 TLS、未做正式显著性、未做物理尺度换算。',
              '- 已知患者仅 8 个研究内单位；切片数不冒称患者数。', '']
    (WORK / 'REPORT.md').write_text('\n'.join(lines) + '\n')


def cross_stage_consistency(measured_rows, records):
    """The measured reference used for deltas must equal the frozen measured table."""
    reference = {}
    for row in measured_rows:
        reference[(row['section_id'], row['target'])] = row
    compared = 0
    worst = 0.0
    for row in records:
        key = (row['section_id'], row['target'])
        ref = reference.get(key)
        if ref is None:
            raise RuntimeError(f'Prediction row without a measured reference: {key}')
        for field in ['dist_sd', 'dist_frac_zero', 'geo_autocorr_near', 'geo_autocorr_outer',
                      'geo_centroid_shift_pitch', 'geo_patch75_n', 'geo_decay_half_pitch']:
            a, b = numeric(ref.get(field)), numeric(row.get(f'meas_{field}'))
            if a is None or b is None:
                if a is not None or b is not None:
                    raise RuntimeError(f'Estimability mismatch across stages: {key}/{field}')
                continue
            worst = max(worst, abs(a - b))
            compared += 1
    if compared and worst > 1e-12:
        raise RuntimeError(f'Measured descriptors differ between prepare and run stages: {worst}')
    return dict(compared=int(compared), max_abs_difference=worst)


def summarize():
    ctx = setup(verify='none')
    measured = list(csv.DictReader((WORK / 'measured_fields.tsv').open(), delimiter='\t'))
    measured_relations = list(csv.DictReader((WORK / 'measured_module_relations.tsv').open(),
                                             delimiter='\t'))
    records, relations, receipts = [], [], []
    seen = set()
    for split in ctx['splits']:
        out = WORK / 'splits' / split['name']
        if not (out / 'COMPLETE').exists():
            raise RuntimeError(f'Split incomplete: {split["name"]}')
        split_rows = list(csv.DictReader((out / 'prediction_metrics.tsv').open(), delimiter='\t'))
        for r in split_rows:
            key = (r['arm'], r['section_id'], r['target'], r['view'], r['model'])
            if key in seen:
                raise RuntimeError(f'Duplicate arm/section/target/view/model row: {key}')
            seen.add(key)
        records += split_rows
        relations += list(csv.DictReader((out / 'prediction_module_relations.tsv').open(), delimiter='\t'))
        receipts.append(json.loads((out / 'split_receipt.json').read_text()))
    if {r['arm'] for r in records} != {'internal', 'source'}:
        raise RuntimeError('Both arms must be present and kept apart')
    cross_stage = cross_stage_consistency(measured, records)
    bad_columns = check_column_names(list(records[0].keys()) + MACRO_METRICS)
    if bad_columns:
        raise RuntimeError(f'Output column names suggest a pass/fail or significance claim: {bad_columns}')
    base.tsv_write(WORK / 'prediction_metrics.tsv', records)
    base.tsv_write(WORK / 'prediction_module_relations.tsv', relations)
    units = aggregate_units(records)
    base.tsv_write(WORK / 'unit_metrics.tsv', units)
    macro = macro_rows(units, records, 'split') + macro_rows(units, records, 'arm')
    base.tsv_write(WORK / 'macro_summary.tsv', macro)
    base.tsv_write(WORK / 'macro_strata.tsv', strata_rows(units))
    base.tsv_write(WORK / 'within_unit_consistency.tsv',
                   within_unit_rows(measured_relations, relations))
    base.tsv_write(WORK / 'diagnostic_sensitivity.tsv', diagnostic_sensitivity(measured))
    base.tsv_write(WORK / 'measured_summary.tsv', measured_summary(measured))
    base.json_write(WORK / 'cross_stage_consistency.json', cross_stage)
    write_report(measured, records, macro, ctx)
    freeze = {}
    for path in sorted(WORK.rglob('*')):
        if path.is_file() and path.name != 'freeze_manifest.json':
            freeze[str(path.relative_to(WORK))] = dict(bytes=path.stat().st_size, sha256=base.sha256(path))
    base.json_write(WORK / 'freeze_manifest.json',
                    dict(code_sha256=base.sha256(__file__),
                         contract_sha256=base.sha256(WORK / 'contract.json'),
                         protocol_sha256=base.sha256(PROTOCOL),
                         files=freeze, human_tls_labels_read=False, biological_pass_threshold=None,
                         shuffle_reference='descriptive perturbation only, no p-value',
                         evidence_roles=['measurement_eligibility',
                                         'input_echo_relation_not_spatial_increment',
                                         'exploratory_reanalysis_of_developed_predictions']))
    (WORK / 'STATUS').write_text('DESCRIPTIVE_REVIEW_COMPLETE_NO_HUMAN_TLS_NO_PASS_THRESHOLD\n')
    print(json.dumps(dict(measured_rows=len(measured), prediction_rows=len(records),
                          relation_rows=len(relations), unit_rows=len(units),
                          macro_rows=len(macro), states=state_counts(records))), flush=True)


# --------------------------------------------------------------------------- checks
def synthetic_geo(nx=20, ny=20):
    xs, ys = np.meshgrid(np.arange(nx, dtype=float), np.arange(ny, dtype=float))
    xy = np.column_stack([xs.ravel(), ys.ravel()])
    index = lambda i, j: i * ny + j
    edges = []
    for i in range(nx):
        for j in range(ny):
            if i + 1 < nx:
                edges.append((index(i, j), index(i + 1, j)))
            if j + 1 < ny:
                edges.append((index(i, j), index(i, j + 1)))
    edges = np.array(edges, np.int64)
    near = sparse.csr_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])),
                             shape=(len(xy), len(xy)))
    near = near + near.T
    near = sparse.diags(1.0 / np.asarray(near.sum(axis=1)).ravel()) @ near
    outer = sparse.identity(len(xy), format='csr')
    return Geo(xy, 1.0, near.tocsr(), outer.tocsr(),
               [(0.0, 1.5), (1.5, 2.5), (2.5, 3.5), (3.5, 5.0)], 'synthetic')


def brute_components(mask, geo):
    adjacency = defaultdict(list)
    for a, b in geo.edges:
        if mask[a] and mask[b]:
            adjacency[int(a)].append(int(b))
            adjacency[int(b)].append(int(a))
    seen = set()
    count = 0
    for node in np.flatnonzero(mask):
        node = int(node)
        if node in seen:
            continue
        count += 1
        stack = [node]
        seen.add(node)
        while stack:
            for nxt in adjacency[stack.pop()]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
    return count


def check_synthetic(seed, n_perm, sd_floor):
    geo = synthetic_geo()
    n = geo.n
    rng = np.random.default_rng(seed)
    results = []

    def record(name, expected, observed, passed):
        results.append(dict(check=name, expected=str(expected), observed=str(observed),
                            passed=bool(passed)))

    constant = np.full(n, 2.5)
    desc, _, notes = describe_field(constant, geo, seed, [], n_perm, False, sd_floor)
    record('constant_autocorr_not_estimable', 'None', desc['geo_autocorr_near'],
           desc['geo_autocorr_near'] is None)
    record('constant_patch_not_estimable', 'None', desc['geo_patch75_n'], desc['geo_patch75_n'] is None)
    record('constant_patch_reason', 'constant_field', 'geo_patch75_n:constant_field' in notes,
           'geo_patch75_n:constant_field' in notes)
    record('constant_sd_band', 'zero', desc['sd_band'], desc['sd_band'] == 'zero')

    tiny_noise = 2.5 + rng.normal(0, 1e-9, n)
    t_desc, _, t_notes = describe_field(tiny_noise, geo, seed, [], n_perm, False, sd_floor)
    record('float_noise_below_floor_not_estimable', 'NOT_ESTIMABLE_DEGENERATE_VARIANCE',
           t_desc['geo_autocorr_near'], t_desc['geo_autocorr_near'] is None)
    record('float_noise_sd_band', 'below_floor', t_desc['sd_band'],
           t_desc['sd_band'] == 'below_floor')
    record('float_noise_reason', 'degenerate_variance',
           'geo_autocorr_near:degenerate_variance' in t_notes,
           'geo_autocorr_near:degenerate_variance' in t_notes)

    gradient = np.repeat(np.arange(20, dtype=float)[None, :], 20, axis=0).ravel()
    g_desc, g_extra, _ = describe_field(gradient, geo, seed, [], n_perm, True, sd_floor)
    record('gradient_centroid_shift_positive', '>0', g_desc['geo_centroid_shift_pitch'],
           g_desc['geo_centroid_shift_pitch'] > 0)
    record('gradient_autocorr_high', '>0.5', g_desc['geo_autocorr_near'],
           g_desc['geo_autocorr_near'] > 0.5)
    record('gradient_decay_defined', 'not None', g_desc['geo_decay_half_pitch'],
           g_desc['geo_decay_half_pitch'] is not None)

    blobs = np.ones((20, 20))
    blobs[2:12, 2:7] = 3.0
    blobs[2:12, 14:19] = 3.0
    blob_field = blobs.ravel()
    b_desc, b_extra, _ = describe_field(blob_field, geo, seed, [], n_perm, True, sd_floor)
    record('two_blobs_patch_count', 2, b_desc['geo_patch75_n'], b_desc['geo_patch75_n'] == 2)
    record('two_blobs_area_frac', '~0.25', b_desc['geo_patch75_area_frac'],
           abs(b_desc['geo_patch75_area_frac'] - 0.25) < 0.02)
    record('two_blobs_shuffled_more_patches', 'shuffled median > observed',
           (b_desc['geo_patch75_n_shuffled_median'], b_desc['geo_patch75_n']),
           b_desc['geo_patch75_n_shuffled_median'] > b_desc['geo_patch75_n'])

    copy_desc, copy_extra, _ = describe_field(blob_field.copy(), geo, seed, [], n_perm, True, sd_floor)
    rec = recovery_metrics(blob_field, blob_field.copy(), geo, b_extra, copy_extra, [], sd_floor)
    record('copy_rmse_zero', 0.0, rec['rec_rmse'], rec['rec_rmse'] == 0.0)
    record('copy_pearson_one', 1.0, rec['rec_pearson'], rec['rec_pearson'] == 1.0)
    record('copy_jaccard_one', 1.0, rec['rec_jaccard_q75'], rec['rec_jaccard_q75'] == 1.0)
    record('copy_local_detail_one', 1.0, rec['rec_local_pattern_pearson'],
           rec['rec_local_pattern_pearson'] == 1.0)
    record('copy_centroid_zero', 0.0, rec['rec_centroid_displacement_pitch'],
           rec['rec_centroid_displacement_pitch'] == 0.0)
    record('copy_sd_ratio_one', 1.0, rec['rec_sd_ratio'], rec['rec_sd_ratio'] == 1.0)

    affine = 2.0 * blob_field + 0.5
    a_desc, a_extra, _ = describe_field(affine, geo, seed, [], n_perm, False, sd_floor)
    rec = recovery_metrics(blob_field, affine, geo, b_extra, a_extra, [], sd_floor)
    record('affine_pearson_one', 1.0, rec['rec_pearson'], rec['rec_pearson'] == 1.0)
    record('affine_sd_ratio_two', 2.0, rec['rec_sd_ratio'], abs(rec['rec_sd_ratio'] - 2.0) < 1e-9)
    record('affine_slope_half', 0.5, rec['rec_calibration_slope'],
           abs(rec['rec_calibration_slope'] - 0.5) < 1e-9)
    record('affine_rmse_positive', '>0', rec['rec_rmse'], rec['rec_rmse'] > 0)

    permuted = blob_field[rng.permutation(n)]
    p_desc, p_extra, _ = describe_field(permuted, geo, seed, [], n_perm, True, sd_floor)
    rec = recovery_metrics(blob_field, permuted, geo, b_extra, p_extra, [], sd_floor)
    record('permuted_pearson_low', '|r|<0.3', rec['rec_pearson'], abs(rec['rec_pearson']) < 0.3)
    record('permuted_jaccard_within_shuffle_range', 'observed <= shuffled p90',
           (rec['rec_jaccard_q75'], rec['rec_jaccard_q75_shuffled_p90']),
           rec['rec_jaccard_q75'] <= rec['rec_jaccard_q75_shuffled_p90'])

    smooth = np.asarray(neighbor_mean(geo, blob_field))
    s_desc, s_extra, _ = describe_field(smooth, geo, seed, [], n_perm, True, sd_floor)
    rec = recovery_metrics(blob_field, smooth, geo, b_extra, s_extra, [], sd_floor)
    record('smoothing_raises_autocorr_lowers_local_detail', 'autocorr up, local detail < copy',
           (s_desc['geo_autocorr_near'], rec['rec_local_pattern_pearson']),
           s_desc['geo_autocorr_near'] > b_desc['geo_autocorr_near']
           and rec['rec_local_pattern_pearson'] < 1.0)
    record('smoothing_has_no_success_column', 'no pass/success column',
           check_column_names(list(s_desc.keys())), check_column_names(list(s_desc.keys())) == [])

    mask_two = np.zeros(n, bool)
    for i in range(2, 12):
        for j in range(2, 12):
            mask_two[i * 20 + j] = True
    for i in range(14, 19):
        for j in range(14, 19):
            mask_two[i * 20 + j] = True
    stats = patch_stats(mask_two, geo)
    record('patch_two_components', 2, stats['n_patches'], stats['n_patches'] == 2)
    brute = brute_components(mask_two, geo)
    record('patch_matches_bruteforce', brute, stats['n_patches'], brute == stats['n_patches'])

    zero_inflated = np.zeros(n)
    positive = np.flatnonzero(rng.random(n) < 0.15)
    zero_inflated[positive] = 1.0 + rng.random(len(positive))
    zi_desc, zi_extra, _ = describe_field(zero_inflated, geo, seed, [], n_perm, False, sd_floor)
    record('zero_inflated_threshold_is_nonzero_set', '~0.15', zi_desc['geo_patch75_area_frac'],
           abs(zi_desc['geo_patch75_area_frac'] - 0.15) < 0.04)
    record('zero_inflated_mask_equals_nonzero', 'True',
           bool(np.array_equal(zi_extra['masks'][0.75], zero_inflated > 0)),
           bool(np.array_equal(zi_extra['masks'][0.75], zero_inflated > 0)))

    isolated_block = np.ones(n)
    isolated_block[0] = 5.0
    isolated_geo = synthetic_geo()
    isolated_geo.near = sparse.csr_matrix((isolated_geo.n, isolated_geo.n))
    isolated_geo.edges = np.zeros((0, 2), np.int64)
    i_desc, i_extra, _ = describe_field(isolated_block, isolated_geo, seed, [], 2, False, sd_floor)
    record('isolated_spots_keep_own_value', 'detail component all zero',
           float(np.max(np.abs(i_extra['dv']))), float(np.max(np.abs(i_extra['dv']))) == 0.0)
    record('isolated_spots_autocorr_not_estimable', 'None', i_desc['geo_autocorr_near'],
           i_desc['geo_autocorr_near'] is None)

    shuffled_labels = np.arange(n)
    rng.shuffle(shuffled_labels)
    try:
        assert_aligned(shuffled_labels, np.arange(n))
        detected = False
    except ValueError:
        detected = True
    record('shuffled_barcodes_detected', 'raises', detected, detected)
    try:
        assert_visible_target_disjoint(['G1', 'G2'], ['G2', 'G3'])
        rejected = False
    except ValueError:
        rejected = True
    record('target_gene_in_visible_input_rejected', 'raises', rejected, rejected)
    record('denylist_matches_annotation_paths', 'all patterns matched',
           [p for p in DENIED_SUBSTRINGS if not is_denied('/x' + p + '/y')],
           all(is_denied('/x' + p + '/y') for p in DENIED_SUBSTRINGS))
    record('numeric_rejects_nan', 'None', numeric(float('nan')), numeric(float('nan')) is None)
    record('quantile_summary_skips_nan', (1.0, 1.0, 1.0),
           quantile_summary([float('nan'), 1.0]), quantile_summary([float('nan'), 1.0])[0] == 1.0)
    record('fake_success_column_rejected', 'detected', check_column_names(['rec_rmse', 'recovery_pass']),
           check_column_names(['rec_rmse', 'recovery_pass']) == ['recovery_pass'])
    return results


def dense_graph_autocorr(graph, dv, use_stored_weights):
    """Naive row-wise evaluation of Σ_i dv_i·(W dv)_i with the frozen isolated-point rule."""
    indptr, indices = graph.indptr, graph.indices
    stored = graph.data.astype(np.float64)
    total = 0.0
    for i in range(len(dv)):
        lo, hi = indptr[i], indptr[i + 1]
        if lo == hi:
            total += dv[i] * dv[i]
            continue
        weights = stored[lo:hi] if use_stored_weights else np.full(hi - lo, 1.0 / (hi - lo))
        total += dv[i] * float(np.dot(weights, dv[indices[lo:hi]]))
    return total


def compute_dense_checks(y, geo):
    """Independent dense re-implementation of the declared statistics.

    Two weight variants are computed so the residual can be attributed exactly:
      * stored weights: the frozen near/outer graphs keep their float32 weights, as the
        frozen pipeline does; the naive re-implementation over those same weights must
        agree to machine precision (this tests the arithmetic path);
      * uniform float64 weights 1/degree; differences against the stored-weight path are
        caused by float32 rounding of the stored weights and are reported, not silently
        absorbed by a wide tolerance.
    Isolated spots (degree 0) keep their own value for both graphs, exactly as
    base.neighbor_mean does and as the frozen feature construction does for near and outer.
    """
    n = len(y)
    indptr = geo.near.indptr
    indices = geo.near.indices
    stored = geo.near.data.astype(np.float64)
    nm_stored = np.empty(n, np.float64)
    nm_uniform = np.empty(n, np.float64)
    weight_deviation = 0.0
    for i in range(n):
        lo, hi = indptr[i], indptr[i + 1]
        if lo == hi:
            nm_stored[i] = y[i]
            nm_uniform[i] = y[i]
            continue
        neighbours = indices[lo:hi]
        weights = stored[lo:hi]
        nm_stored[i] = float(np.dot(weights, y[neighbours]))
        nm_uniform[i] = float(np.mean(y[neighbours]))
        weight_deviation = max(weight_deviation, float(np.max(np.abs(weights - 1.0 / len(neighbours)))))
    out = dict(weight_deviation_from_uniform=weight_deviation,
               n_edges_frozen=int(len(geo.edges)),
               n_outer_isolated=int(geo.outer_isolated),
               n_pairs_dense=int(len(cKDTree(geo.xy).query_pairs(1.1 * geo.pitch, output_type='ndarray'))))
    dv = y - nm_stored
    den = float(np.sum(dv ** 2))
    for label, graph in [('near', geo.near), ('outer', geo.outer)]:
        stored_total = dense_graph_autocorr(graph, dv, True)
        uniform_total = dense_graph_autocorr(graph, dv, False)
        out[f'autocorr_{label}_stored'] = (stored_total / den) if den > 0 else None
        out[f'autocorr_{label}_uniform'] = (uniform_total / den) if den > 0 else None
    weights = y - y.min()
    centroid = (weights[:, None] * geo.xy).sum(axis=0) / weights.sum()
    out['sd'] = float(np.sqrt(np.mean((y - y.mean()) ** 2)))
    out['centroid_shift'] = float(np.linalg.norm(centroid - geo.xy.mean(axis=0)) / geo.pitch)
    return out


def gene_matrix_provenance(ctx, info, cache):
    """Row->barcode and column->gene-symbol provenance for the frozen visible matrices."""
    target = info['target']
    sid = info['section_id']
    with h5py.File(base.HEST / f"st/{sid}.h5ad") as f:
        genes = base.strings(f['var/_index'][:])
        barcodes = base.strings(f['obs/_index'][:])
        keep = np.flatnonzero(f['obs/in_tissue'][:] == 1)
        x = base.read_matrix(f['X'], len(barcodes), len(genes))[keep]
    assert_aligned(barcodes[keep], cache['barcodes'])
    out = {}
    for view in ctx['c']['views']:
        recomputed, _ = follow.gene_input(x, genes, ctx['vocab'][target][view])
        root = PRED if view == 'small5' else FUP
        stored = np.load(root / f'genes/{sid}_t{info["ti"]}.npy')
        diff = float(np.max(np.abs(recomputed.astype(np.float64) - stored.astype(np.float64))))
        permuted = recomputed[np.random.default_rng(7).permutation(len(recomputed))]
        perm_diff = float(np.max(np.abs(permuted.astype(np.float64) - stored.astype(np.float64))))
        out[view] = dict(max_abs_diff=diff, row_permutation_detected=bool(perm_diff > 1e-3),
                         n_genes=int(stored.shape[1]), n_spots=int(stored.shape[0]),
                         stored_file=str((root / f'genes/{sid}_t{info["ti"]}.npy').relative_to(ROOT)))
    return out


def check_real(ctx, sids, c):
    results = []
    store, relations, timings = [], [], {}
    provenance_sections = set(sids[:2])

    def record(name, expected, observed, passed):
        results.append(dict(check=name, expected=str(expected), observed=str(observed),
                            passed=bool(passed)))

    frozen_metrics = {}
    for row in csv.DictReader((PRED / 'section_metrics.tsv').open(), delimiter='\t'):
        frozen_metrics[(row['split'], row['section_id'], row['target'], row['view'], row['model'])] = row
    for sid in sids:
        row = ctx['byid'][sid]
        cache = base.load_cache(row)
        geo, pitch_recomputed = geo_of(row, cache, c)
        record(f'{sid}_pitch_matches_manifest', row['empirical_pitch_source_pixels'], pitch_recomputed,
               abs(pitch_recomputed - float(row['empirical_pitch_source_pixels'])) < 1e-9)
        record(f'{sid}_isolated_spots_match_manifest', row['near_isolated_spots'], geo.isolated,
               int(row['near_isolated_spots']) == geo.isolated)
        record(f'{sid}_confirmed_unit_matches_identity_table',
               ctx['identity'][sid]['patient_confirmed_within_study'],
               ctx['byid'][sid]['confirmed_unit'],
               bool(ctx['identity'][sid]['patient_confirmed_within_study'] == 'True')
               == ctx['byid'][sid]['confirmed_unit'])
        start = time.time()
        for ti, target in enumerate(c['targets']):
            info = dict(row, ti=ti, target=target)
            sup = support_metrics(ctx, info, cache)
            record(f'{sid}_{target}_readout_matches_frozen', '<=1e-6', sup['sup_readout_max_abs_diff'],
                   sup['sup_readout_max_abs_diff'] <= 1e-6)
            record(f'{sid}_{target}_support_denominators_consistent', 0,
                   sup['sup_genes_absent_from_panel'] + sup['sup_genes_in_panel']
                   - sup['sup_genes_defined'],
                   sup['sup_genes_absent_from_panel'] + sup['sup_genes_in_panel']
                   == sup['sup_genes_defined'])
            record(f'{sid}_{target}_detected_denominator_consistent', 0,
                   sup['sup_genes_in_panel_zero_counts'] + sup['sup_genes_detected']
                   - sup['sup_genes_in_panel'],
                   sup['sup_genes_in_panel_zero_counts'] + sup['sup_genes_detected']
                   == sup['sup_genes_in_panel'])
            if sid in provenance_sections:
                prov = gene_matrix_provenance(ctx, info, cache)
                for view, value in prov.items():
                    record(f'{sid}_{target}_{view}_gene_matrix_matches_frozen', '<=1e-5',
                           value['max_abs_diff'], value['max_abs_diff'] <= 1e-5)
                    record(f'{sid}_{target}_{view}_gene_matrix_row_order_checkable', 'detected',
                           value['row_permutation_detected'], value['row_permutation_detected'])
                    record(f'{sid}_{target}_{view}_gene_matrix_shape', 'n_genes>0',
                           (value['n_spots'], value['n_genes']),
                           value['n_genes'] > 0 and value['n_spots'] == len(cache['y']))
            y, mdesc, mex, _ = measured_descriptors(ctx, info, cache, geo)
            _, mdesc_reshuffled, _, _ = measured_descriptors(ctx, info, cache, geo, do_shuffle=True)
            shared = [k for k in mdesc if k in mdesc_reshuffled]
            differences = [abs(numeric(v) - numeric(mdesc_reshuffled[k])) for k, v in mdesc.items()
                           if k in mdesc_reshuffled and numeric(v) is not None
                           and numeric(mdesc_reshuffled[k]) is not None]
            record(f'{sid}_{target}_measured_descriptor_shuffle_invariance', 'all identical',
                   (len(shared), max(differences) if differences else None),
                   bool(differences) and max(differences) == 0.0)
            dense = compute_dense_checks(y, geo)
            record(f'{sid}_{target}_frozen_edge_set_matches_dense_pairs', dense['n_edges_frozen'],
                   dense['n_pairs_dense'], dense['n_edges_frozen'] == dense['n_pairs_dense'])
            record(f'{sid}_{target}_stored_weights_match_uniform_within_float32', '<=1e-6',
                   dense['weight_deviation_from_uniform'],
                   dense['weight_deviation_from_uniform'] <= 1e-6)
            record(f'{sid}_{target}_outer_isolated_spots', 'recorded', dense['n_outer_isolated'],
                   dense['n_outer_isolated'] == geo.outer_isolated)
            for ring in ['near', 'outer']:
                record(f'{sid}_{target}_autocorr_{ring}_matches_dense_same_weights', '<=1e-12',
                       abs(dense[f'autocorr_{ring}_stored'] - mdesc[f'geo_autocorr_{ring}']),
                       abs(dense[f'autocorr_{ring}_stored'] - mdesc[f'geo_autocorr_{ring}']) <= 1e-12)
                record(f'{sid}_{target}_autocorr_{ring}_uniform_vs_stored_weight_residual', '<=1e-6',
                       abs(dense[f'autocorr_{ring}_uniform'] - mdesc[f'geo_autocorr_{ring}']),
                       abs(dense[f'autocorr_{ring}_uniform'] - mdesc[f'geo_autocorr_{ring}']) <= 1e-6)
            record(f'{sid}_{target}_sd_matches_dense', dense['sd'], mdesc['dist_sd'],
                   abs(dense['sd'] - mdesc['dist_sd']) < 1e-9)
            record(f'{sid}_{target}_centroid_matches_dense', dense['centroid_shift'],
                   mdesc['geo_centroid_shift_pitch'],
                   abs(dense['centroid_shift'] - mdesc['geo_centroid_shift_pitch']) < 1e-9)
            for p_q in PATCH_QS:
                mask = y > float(np.quantile(y, p_q))
                stats = patch_stats(mask, geo)
                brute = brute_components(mask, geo)
                record(f'{sid}_{target}_patch{pkey(p_q)}_matches_bruteforce', brute,
                       None if stats is None else stats['n_patches'],
                       stats is not None and stats['n_patches'] == brute)
            for split in ctx['splits']:
                if sid not in split['test']:
                    continue
                prows, prels = prediction_rows(ctx, split, info, cache, geo, mdesc, mex, y)
                store += prows
                relations += prels
                for p in prows:
                    frozen = frozen_metrics.get((p['split'], sid, target, p['view'], p['model']))
                    if frozen is None:
                        continue
                    for key, value, name in [('rec_rmse', math.sqrt(float(frozen['mse'])), 'rmse'),
                                             ('rec_centered_rmse', math.sqrt(float(frozen['centered_mse'])),
                                              'centered_rmse')]:
                        record(f'{sid}_{target}_{p["view"]}_{p["model"]}_{name}_matches_frozen',
                               value, p[key], abs(value - p[key]) < 1e-5)
                    if frozen['pearson'] != '':
                        frozen_value = float(frozen['pearson'])
                        mine = p['rec_pearson']
                        if mine is None:
                            record(f'{sid}_{target}_{p["view"]}_{p["model"]}_constant_pearson_artifact',
                                   'recompute None, frozen <= 1e-6', frozen_value,
                                   abs(frozen_value) <= 1e-6
                                   and float(p['dist_sd']) < c['sd_floor'])
                        else:
                            record(f'{sid}_{target}_{p["view"]}_{p["model"]}_pearson_matches_frozen',
                                   frozen_value, mine, abs(frozen_value - mine) < 1e-5)
                    if p['model'] == 'train_mean':
                        record(f'{sid}_{target}_{p["view"]}_train_mean_degenerate_state',
                               'NOT_ESTIMABLE_*', p['state'], str(p['state']).startswith('NOT_ESTIMABLE'))
                        record(f'{sid}_{target}_{p["view"]}_train_mean_rmse_kept', 'not None',
                               p['rec_rmse'], p['rec_rmse'] is not None)
                        record(f'{sid}_{target}_{p["view"]}_train_mean_sd_below_floor',
                               f'<{c["sd_floor"]}', p['dist_sd'],
                               float(p['dist_sd']) < c['sd_floor'])
        timings[sid] = time.time() - start
    arms = {(r['arm'], r['split'], r['section_id'], r['target'], r['view'], r['model'])
            for r in store}
    record('no_arm_row_duplication', 'unique', len(arms), len(arms) == len(store))
    groups = defaultdict(set)
    for r in relations:
        groups[(r['model'], r['split'], r['unit_id'], r['target'], r['view'])].add(r['section_id'])
    expected_within = 2 * sum(1 for value in groups.values() if len(value) >= 2)
    real_within = within_unit_rows([], relations)
    record('within_unit_runs_on_real_prediction_relation_rows', expected_within, len(real_within),
           len(real_within) == expected_within)
    record('within_unit_real_profile_columns', "['measured_rel_pearson', 'predicted_rel_pearson']",
           sorted({r['profile_column'] for r in real_within}),
           sorted({r['profile_column'] for r in real_within})
           <= ['measured_rel_pearson', 'predicted_rel_pearson'])
    record('internal_and_source_arms_separated', 'both present',
           sorted({r['arm'] for r in store}), {r['arm'] for r in store} == {'internal', 'source'})
    record('organ_extrapolation_flagged', 'at least one True row',
           sum(1 for r in store if str(r['cross_organ_extrapolation']) == 'True'),
           any(str(r['cross_organ_extrapolation']) == 'True' for r in store))
    record('prediction_rows_role', 'exploratory_reanalysis_of_developed_predictions',
           sorted({r['evidence_role'] for r in store}),
           all(r['evidence_role'] == 'exploratory_reanalysis_of_developed_predictions'
               for r in store))
    record('relation_rows_role', 'input_echo_relation_not_spatial_increment',
           sorted({r['evidence_role'] for r in relations}),
           all(r['evidence_role'] == 'input_echo_relation_not_spatial_increment'
               for r in relations))
    return results, store, relations, timings


def check_prediction_inventory(ctx, c):
    results = []

    def record(name, expected, observed, passed):
        results.append(dict(check=name, expected=str(expected), observed=str(observed),
                            passed=bool(passed)))

    selected = []
    for split in ctx['splits']:
        for sid in split['test']:
            for ti in range(len(c['targets'])):
                for view in c['views']:
                    selected.append(PRED / split['name'] / f'{sid}_t{ti}_{view}.npz')
    record('prediction_packages_selected', 1000, len(selected), len(selected) == 1000)
    record('all_selected_files_exist', 'True', all(p.exists() for p in selected),
           all(p.exists() for p in selected))
    record('train_indices_excluded', 'True',
           all('train_indices' not in p.name for p in selected),
           all('train_indices' not in p.name for p in selected))
    record('prediction_packages_match_freeze_list', 'True',
           all(str(p.relative_to(PRED)) in
               json.loads((PRED / 'freeze_manifest.json').read_text())['files'] for p in selected[:50]),
           all(str(p.relative_to(PRED)) in
               json.loads((PRED / 'freeze_manifest.json').read_text())['files'] for p in selected[:50]))
    return results


def check_guard_writes(ctx):
    results = []

    def record(name, expected, observed, passed):
        results.append(dict(check=name, expected=str(expected), observed=str(observed),
                            passed=bool(passed)))

    probe = ROOT / 'infra/field_biological_review_20261019_write_probe.txt'
    try:
        with probe.open('w') as handle:
            handle.write('probe')
        blocked = False
    except PermissionError:
        blocked = True
    record('write_outside_output_dir_rejected', 'raises', blocked, blocked)
    record('write_probe_not_created', 'absent', probe.exists(), not probe.exists())
    denied = ROOT / 'data/other_sources/zenodo_st_crc_cms/Pathology_SpotAnnotations/probe.txt'
    record('hidden_truth_path_denied', 'True', is_denied(denied), is_denied(denied))
    return results


def aggregation_record(split, arm, sid, unit, extrapolation=False, **over):
    """Minimal record carrying the real prediction-column set, for aggregation contracts."""
    row = dict(section_id=sid, split=split, arm=arm, unit_id=unit, source_key='https://x',
               organ='Bowel', platform='Visium', patient_known='False', confirmed_unit='False',
               target='BOWEL', view='full15', model='gene_ridge', model_class='raw_gene_ridge',
               cross_organ_extrapolation=extrapolation, n_sections=1,
               evidence_role='exploratory_reanalysis_of_developed_predictions', state='DESCRIBED')
    for metric in MACRO_METRICS:
        row[metric] = 0.5 if ('jaccard' in metric or 'pearson' in metric) else 0.1
    row.update(over)
    return row


def check_aggregation_contracts(seed):
    """Counterexamples for the unit/arm/stratum aggregation contracts (review R-01, R-02, R-04)."""
    results = []

    def record(name, expected, observed, passed):
        results.append(dict(check=name, expected=str(expected), observed=str(observed),
                            passed=bool(passed)))

    pairs = [aggregation_record('internal_0', 'internal', 'S1', 'u1', rec_pearson=0.8),
             aggregation_record('internal_0', 'internal', 'S2', 'u1', rec_pearson='')]
    units = aggregate_units(pairs)
    unit = units[0]
    record('partial_unit_metric_is_not_estimable', 'None', unit['rec_pearson'],
           unit['rec_pearson'] is None)
    record('partial_unit_metric_recorded', 'rec_pearson in partial_unit_metrics',
           unit['partial_unit_metrics'], 'rec_pearson' in unit['partial_unit_metrics'])
    record('partial_unit_estimable_section_count', 1, unit['rec_pearson__n'],
           unit['rec_pearson__n'] == 1)
    record('partial_unit_not_reported_as_complete', 'False', unit['complete_for_all_metrics'],
           unit['complete_for_all_metrics'] is False)
    macro = macro_rows(units, pairs, 'split')
    macro_partial = [r for r in macro if r['metric'] == 'rec_pearson']
    record('macro_counts_partial_units', 1,
           macro_partial[0]['n_units_partial'] if macro_partial else None,
           bool(macro_partial) and macro_partial[0]['n_units_partial'] == 1)
    record('macro_excludes_partial_unit_from_mean', 0,
           macro_partial[0]['n_estimable'] if macro_partial else None,
           bool(macro_partial) and macro_partial[0]['n_estimable'] == 0)
    strata = strata_rows(units)
    empty_strata = [r for r in strata if r['metric'] == 'rec_pearson' and r['n_estimable'] == 0]
    record('empty_stratum_row_kept', 'row with reason', len(empty_strata),
           bool(empty_strata) and empty_strata[0]['reason'] == 'no_estimable_units_in_stratum')

    arm = [aggregation_record('source_0', 'source', 'S1', 'u1', extrapolation=True),
           aggregation_record('source_1', 'source', 'S2', 'u2', extrapolation=False)]
    arm_units = aggregate_units(arm)
    arm_macro = [r for r in macro_rows(arm_units, arm, 'arm') if r['arm'] == 'source'
                 and r['metric'] == 'rec_rmse']
    record('arm_level_extrapolation_label_mixed', 'mixed',
           arm_macro[0]['cross_organ_extrapolation'] if arm_macro else None,
           bool(arm_macro) and arm_macro[0]['cross_organ_extrapolation'] == 'mixed')
    split_macro = [r for r in macro_rows(arm_units, arm, 'split') if r['metric'] == 'rec_rmse']
    labels = sorted({str(r['cross_organ_extrapolation']) for r in split_macro})
    record('split_level_extrapolation_label_kept', "['False', 'True']", labels,
           labels == ['False', 'True'])

    measured_like = [dict(unit_id='u1', target='BOWEL', view='full15', section_id='S1',
                          rel_module=f'M{i}', rel_pearson=0.1 * i) for i in range(3)] + \
                    [dict(unit_id='u1', target='BOWEL', view='full15', section_id='S2',
                          rel_module=f'M{i}', rel_pearson=0.12 * i) for i in range(3)]
    measured_rows = within_unit_rows(measured_like, [])
    record('within_unit_reads_measured_columns', 1, len(measured_rows),
           len(measured_rows) == 1)
    predicted_like = []
    for sid, scale in [('S1', 0.1), ('S2', 0.12)]:
        for i in range(3):
            predicted_like.append(dict(unit_id='u1', target='BOWEL', view='full15', section_id=sid,
                                       rel_module=f'M{i}', model='gene_ridge', split='internal_0',
                                       measured_rel_pearson=scale * i,
                                       predicted_rel_pearson=scale * i * 0.5))
    prediction_rows = within_unit_rows([], predicted_like)
    columns = sorted({r['profile_column'] for r in prediction_rows})
    record('within_unit_reads_prediction_columns', "['measured_rel_pearson', 'predicted_rel_pearson']",
           columns, columns == ['measured_rel_pearson', 'predicted_rel_pearson'])
    try:
        within_unit_rows([], [dict(unit_id='u1', target='BOWEL', view='full15', section_id='S1',
                                   rel_module='M0', model='gene_ridge', split='internal_0')])
        closed = False
    except KeyError:
        closed = True
    record('within_unit_fails_closed_on_missing_columns', 'raises KeyError', closed, closed)
    return results


def check_command(sids=None):
    if sids is None:
        sids = ['MISC33', 'TENX154', 'TENX161', 'NCBI642', 'TENX171']
    ctx = setup(verify='sample', sids=sids)
    c = ctx['c']
    out = WORK / 'checks'
    out.mkdir(parents=True, exist_ok=True)
    synthetic = check_synthetic(c['seed'], c['shuffle_permutations'], c['sd_floor'])
    aggregation = check_aggregation_contracts(c['seed'])
    inventory = check_prediction_inventory(ctx, c)
    guard_checks = check_guard_writes(ctx)
    real, store, relations, timings = check_real(ctx, sids, c)
    results = synthetic + aggregation + inventory + guard_checks + real
    failures = [r for r in results if not r['passed']]
    per_section = float(np.mean(list(timings.values()))) if timings else None
    visits = Counter()
    for split in ctx['splits']:
        visits[split['arm']] += len(split['test'])
    n_visits = sum(visits.values())
    extrapolation = dict(measured_section_seconds=timings,
                         mean_seconds_per_section_run_stage=per_section,
                         split_section_visits=dict(visits),
                         projected_run_cpu_seconds=(per_section * n_visits if per_section else None),
                         projected_run_wall_hours=(round(per_section * n_visits / c['workers'] / 3600.0, 2)
                                                   if per_section else None),
                         projected_prepare_cpu_seconds=(round(per_section * 72, 1) if per_section else None),
                         note='run-stage projection uses split test-section visits (72 internal + 53 source); '
                              'prepare reads the same sections once each')
    payload = dict(candidate_hashes=dict(script=base.sha256(__file__),
                                         contract=base.sha256(WORK / 'contract.json')),
                   check_sections=sids, n_checks=len(results), n_failed=len(failures),
                   results=results, failures=failures, timings_seconds=timings,
                   extrapolation=extrapolation, setup_checks=ctx['checks'])
    base.json_write(out / 'checks.json', payload)
    base.tsv_write(out / 'prediction_metrics_sample.tsv', store)
    base.tsv_write(out / 'prediction_module_relations_sample.tsv', relations)
    print(json.dumps(dict(n_checks=len(results), n_failed=len(failures),
                          failures=failures[:8], extrapolation=extrapolation), indent=2), flush=True)
    return 0 if not failures else 1


# --------------------------------------------------------------------------- main
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['prepare', 'run', 'summarize', 'check'])
    parser.add_argument('--sections', nargs='*', default=None)
    args = parser.parse_args()
    if os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', '-1'):
        raise RuntimeError('GPU disabled')
    if args.command == 'prepare':
        prepare()
    elif args.command == 'run':
        run()
    elif args.command == 'summarize':
        summarize()
    else:
        sys.exit(check_command(args.sections))


if __name__ == '__main__':
    main()
