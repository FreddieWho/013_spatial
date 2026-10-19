#!/usr/bin/env python3
"""Proxy-field audit, segment 1 (contract: docs/FIELD_PROXY_AUDIT_CONTRACT_20261019.md, D-277).

Rules fixed before computation (D-277):
  * one readout formula for every map: gene map = log1p(1e4 * raw count / all-gene library of
    that spot); a set map is the arithmetic mean of its gene maps; single-gene,
    leave-one-out and disjoint-subset maps share the same counts, denominator, tissue spots
    and frozen barcode order and only change the averaged gene set
  * measurable gene = in panel, at least 20 tissue spots with a positive count, gene map
    variance > 0; unmeasurable genes stay in the original proxy definition and coverage record
    and are never treated as biological negatives
  * technical control = one within-section regression per measurable gene map on
    [intercept, log1p(non-target library), non-target detection fraction]; raw and residual
    results are always reported separately and residuals never replace the raw maps
  * disjoint subsets need at least 3 measurable genes each (>= 6 in total), are evaluated with
    10 fixed seeds that are all reported, and are compared with a quartile-bin matched
    out-of-set control for the same 10 seeds; insufficient bins make that draw NOT_ESTIMABLE
  * metrics are computed per section first, then aggregated over complete established units
    only (equal unit weight, 0 is a valid value); segment 1 is descriptive only and never fixes
    the members, weights or candidate counts of segment 2

Commands: prepare | audit | summarize | smoke
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
from pathlib import Path

import h5py
import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import field_discovery_20261018 as base                     # noqa: E402
import field_biological_review_20261019 as review          # noqa: E402

WORK = ROOT / 'infra/field_proxy_audit_20261019'
PANEL = ROOT / 'infra/tls_multiscale_gobp_20261002/panel_controls.json'
SUPPORT = ROOT / 'infra/field_biological_review_20261019/measured_fields.tsv'
MANIFEST = ROOT / 'infra/field_source_holdout_20261018/input_manifest.json'
IDENTITY = ROOT / 'infra/field_discovery_review_20261018/identity_review.tsv'
ELIGIBILITY = ROOT / 'infra/field_source_holdout_20261018/source_eligibility.tsv'
CACHE = base.WORK
DRAWS = 10
TABLES = ['section_target_summary.tsv', 'synergy.tsv', 'gene_coherence.tsv',
          'deletion_sensitivity.tsv', 'subset_consistency.tsv', 'subset_draws.tsv',
          'not_testable.tsv']
METRICS = ['n_measurable', 'effective_gene_count', 'top_gene_share', 'proxy_technical_r2',
           'pair_raw.median', 'pair_raw.q10', 'pair_raw.q90', 'pair_residual.median',
           'pair_residual.q10', 'pair_residual.q90', 'technical_r2.median',
           'deletion_corr_vs_measurable_mean', 'deletion_jaccard_q75',
           'subset_corr_median', 'subset_jaccard_q75_median', 'control_corr_median']


def machine_contract():
    return json.loads((WORK / 'contract.json').read_text())


def seed_for(*parts):
    return int.from_bytes(hashlib.sha256('|'.join(str(p) for p in parts).encode()).digest()[:8],
                          'little')


def setup(verify='sample', sids=None):
    c = machine_contract()
    rows = json.loads(MANIFEST.read_text())
    byid = {r['section_id']: r for r in rows}
    panels = json.loads(PANEL.read_text())
    identity = {r['section_id']: r for r in csv.DictReader(IDENTITY.open(), delimiter='\t')}
    eligibility = {r['section_id']: r['status'] for r in
                   csv.DictReader(ELIGIBILITY.open(), delimiter='\t')}
    for row in rows:
        if bool(identity[row['section_id']]['patient_confirmed_within_study'] == 'True') \
                != bool(row['patient_known']):
            raise ValueError(f'Identity flag disagreement: {row["section_id"]}')
        row['confirmed_unit'] = bool(row['patient_known'])
        row['source_eligibility'] = eligibility.get(row['section_id'], 'UNKNOWN')
    frozen_order = json.loads((CACHE / 'contract.json').read_text())['targets']
    if set(frozen_order) != set(c['targets']):
        raise ValueError('Target set differs from the frozen cache order')
    c['frozen_target_order'] = frozen_order
    allowed = {Path(__file__), Path(base.__file__), Path(review.__file__), PANEL, SUPPORT,
               MANIFEST, IDENTITY, ELIGIBILITY, WORK / 'contract.json', CACHE / 'contract.json'}
    freeze = json.loads((CACHE / 'freeze_manifest.json').read_text())
    allowed.add(CACHE / 'freeze_manifest.json')
    allowed.update(CACHE / p for p in freeze['files'])
    allowed.update(base.HEST / f'st/{sid}.h5ad' for sid in byid)
    guard = base.InputGuard(ROOT, allowed, WORK)
    guard.install()
    review.install_denylist_hook()
    base.reserve(c)
    checks = dict(targets={t: len(panels[t]) for t in c['targets']},
                  frozen_target_order=frozen_order,
                  n_sections=len(rows), n_units=len({r['unit_id'] for r in rows}),
                  known_sections=sum(1 for r in rows if r['confirmed_unit']),
                  known_units=len({r['unit_id'] for r in rows if r['confirmed_unit']}),
                  source_not_testable=sum(1 for r in rows
                                          if r['source_eligibility'].startswith('NOT_TESTABLE')),
                  free_bytes=shutil.disk_usage(WORK).free)
    checks['storage_reserve_ok'] = checks['free_bytes'] >= c['min_free_bytes']
    if not checks['storage_reserve_ok']:
        raise RuntimeError(f'Storage reserve violated: {checks["free_bytes"]}')
    if verify in ('all', 'sample'):
        verified = Counter()
        for row in rows:
            if verify == 'sample' and sids is not None and row['section_id'] not in sids:
                continue
            if base.sha256(base.HEST / f'st/{row["section_id"]}.h5ad') != row['expression_sha256']:
                raise ValueError(f'Expression hash changed: {row["section_id"]}')
            with np.load(CACHE / f'cache/{row["section_id"]}.npz', allow_pickle=False) as z:
                if len(z['barcodes']) != row['n_spots']:
                    raise ValueError(f'Cache spot count changed: {row["section_id"]}')
            verified['sections'] += 1
        checks['verified'] = dict(verified)
    return dict(c=c, rows=rows, byid=byid, panels=panels, guard=guard, checks=checks)


# --------------------------------------------------------------------------- statistics
def corr(a, b):
    a = np.asarray(a, np.float64)
    b = np.asarray(b, np.float64)
    a = a - a.mean()
    b = b - b.mean()
    den = math.sqrt(float(a @ a) * float(b @ b))
    return float(a @ b / den) if den > 0 else None


def corr_matrix(matrix):
    x = np.asarray(matrix, np.float64)
    x = x - x.mean(axis=0)
    sd = np.sqrt((x ** 2).sum(axis=0))
    with np.errstate(invalid='ignore', divide='ignore'):
        out = (x.T @ x) / np.outer(sd, sd)
    out[~np.isfinite(out)] = np.nan
    return out


def jaccard_high(a, b, q):
    ma = np.asarray(a) > float(np.quantile(a, q))
    mb = np.asarray(b) > float(np.quantile(b, q))
    union = int(np.count_nonzero(ma | mb))
    return (int(np.count_nonzero(ma & mb)) / union) if union else None


def centroid_shift_pitch(xy, a, b, pitch):
    ca = review.weighted_centroid(xy, np.asarray(a) - np.min(a))
    cb = review.weighted_centroid(xy, np.asarray(b) - np.min(b))
    if ca is None or cb is None:
        return None
    return float(np.linalg.norm(ca - cb) / pitch)


def regression_residual(values, covariates):
    """One multivariate within-section regression: value ~ intercept + covariates."""
    design = np.c_[np.ones(len(covariates)), covariates]
    beta, *_ = np.linalg.lstsq(design, values, rcond=None)
    return values - design @ beta


def variance_explained(raw, residual):
    var = float(np.var(raw))
    return float(1 - np.var(residual) / var) if var > 0 else None


def distribution(values):
    arr = np.asarray([v for v in values if v is not None and np.isfinite(v)], np.float64)
    if len(arr) == 0:
        return dict(n=0, mean=None, median=None, q10=None, q25=None, q75=None, q90=None,
                    frac_positive=None, frac_above_0_5=None)
    return dict(n=int(len(arr)), mean=float(arr.mean()), median=float(np.median(arr)),
                q10=float(np.quantile(arr, 0.10)), q25=float(np.quantile(arr, 0.25)),
                q75=float(np.quantile(arr, 0.75)), q90=float(np.quantile(arr, 0.90)),
                frac_positive=float(np.mean(arr > 0)), frac_above_0_5=float(np.mean(arr > 0.5)))


def expand(prefix, stats):
    return {f'{prefix}.{k}': v for k, v in stats.items()}


# --------------------------------------------------------------------------- section data
def read_section(row, union_genes):
    with h5py.File(base.HEST / f'st/{row["section_id"]}.h5ad') as handle:
        genes = base.strings(handle['var/_index'][:])
        barcodes = base.strings(handle['obs/_index'][:])
        keep = np.flatnonzero(handle['obs/in_tissue'][:] == 1)
        matrix = base.read_matrix(handle['X'], len(barcodes), len(genes))[keep]
    lookup = {}
    for gene in union_genes:
        hits = np.flatnonzero(genes == gene)
        if len(hits) > 1:
            raise ValueError(f'Duplicate gene symbol in the panel: {gene}')
        if len(hits) == 1:
            lookup[gene] = int(hits[0])
    positions = [lookup[g] for g in union_genes if g in lookup]
    present = [g for g in union_genes if g in lookup]
    library = np.asarray(matrix.sum(axis=1)).ravel()
    z_union = np.asarray(base.normalized_subset(matrix, positions, library).todense()) \
        if positions else np.zeros((len(keep), 0))
    return dict(genes=genes, barcodes=barcodes[keep], matrix=matrix, positions=positions,
                present=present, index={g: i for i, g in enumerate(present)},
                z_union=np.asarray(z_union, np.float64), library=library,
                n_genes_total=matrix.shape[1])


def build_bins(union_stats):
    if len(union_stats) < 4:
        return None, None
    detection = np.array([u[1] for u in union_stats])
    expression = np.array([u[2] for u in union_stats])
    return np.quantile(detection, [0.25, 0.5, 0.75]), np.quantile(expression, [0.25, 0.5, 0.75])


def match_cells(needed, pool_by_cell, rng):
    """Sample pool genes for one half only (kept for reference; the control uses the joint rule)."""
    chosen = []
    for cell, count in needed.items():
        candidates = pool_by_cell.get(cell, [])
        if len(candidates) < count:
            return None
        picks = rng.choice(len(candidates), size=count, replace=False)
        chosen.extend(candidates[int(i)] for i in picks)
    return chosen


def match_cells_joint(pool_by_cell, needed_a, needed_b, rng):
    """D-277 no-reuse rule: both halves are sampled in one pass, so a pool gene used by one
    half can never be used by the other half of the same control draw."""
    chosen_a, chosen_b = [], []
    for cell in sorted(set(needed_a) | set(needed_b)):
        na, nb = needed_a.get(cell, 0), needed_b.get(cell, 0)
        candidates = pool_by_cell.get(cell, [])
        if len(candidates) < na + nb:
            return None
        picks = rng.choice(len(candidates), size=na + nb, replace=False)
        chosen_a.extend(candidates[int(i)] for i in picks[:na])
        chosen_b.extend(candidates[int(i)] for i in picks[na:])
    if set(chosen_a) & set(chosen_b):
        return None
    return chosen_a, chosen_b


def audit_section(job):
    ctx = setup(verify='none')
    c, panels = ctx['c'], ctx['panels']
    row = ctx['byid'][job['sid']]
    out = Path(job['out'])
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'COMPLETE').exists() and not job.get('force'):
        return dict(sid=job['sid'], status='already_complete')
    started = time.time()
    union_genes = sorted({g for t in c['targets'] for g in panels[t]})
    data = read_section(row, union_genes)
    cache = base.load_cache(row)
    review.assert_aligned(data['barcodes'], cache['barcodes'])
    xy = np.asarray(cache['xy_source_pixels'], np.float64)
    pitch = float(row['empirical_pitch_source_pixels'])
    if abs(float(np.median(cKDTree(xy).query(xy, k=2)[0][:, 1])) - pitch) > 1e-9:
        raise ValueError(f'Pitch mismatch: {row["section_id"]}')
    identity = dict(section_id=row['section_id'], unit_id=row['unit_id'], organ=row['organ'],
                    platform=row['metadata']['st_technology'], source_key=row['source_key'],
                    source_eligibility=row['source_eligibility'],
                    confirmed_unit=row['confirmed_unit'], n_spots=int(len(cache['barcodes'])),
                    pitch_source_pixels=pitch)
    z_u = data['z_union']
    counts_u = data['matrix'][:, data['positions']].tocsr() if data['positions'] else None
    n_any_u = np.asarray((counts_u > 0).sum(axis=0)).ravel() if counts_u is not None else np.zeros(0)
    det_u = np.asarray((counts_u > 0).mean(axis=0)).ravel() if counts_u is not None else np.zeros(0)
    mean_u = z_u.mean(axis=0) if z_u.size else np.zeros(0)
    measurable_u = (n_any_u >= 20) & (np.std(z_u, axis=0) > 0)
    union_stats = [(gene, float(det_u[i]), float(mean_u[i]))
                   for i, gene in enumerate(data['present']) if measurable_u[i]]
    det_edges, expr_edges = build_bins(union_stats)
    summaries, synergy_rows, gene_rows = [], [], []
    deletion_rows, subset_rows, draw_rows, not_testable = [], [], [], []
    for target in c['targets']:
        member_genes = [g for g in panels[target] if g in data['index']]
        n_panel = len(panels[target])
        base_row = dict(identity, target=target, n_panel=n_panel, n_in_panel=len(member_genes),
                        n_absent=n_panel - len(member_genes))
        if not member_genes:
            not_testable.append(dict(section_id=row['section_id'], target=target,
                                     reason='NO_PANEL_GENES', n_panel=n_panel, n_in_panel=0,
                                     n_measurable=0))
            continue
        cols = [data['index'][g] for g in member_genes]
        z = z_u[:, cols]
        counts = counts_u[:, cols]
        n_any = n_any_u[cols]
        detection = det_u[cols]
        mean_norm = mean_u[cols]
        measurable = (n_any >= 20) & (np.std(z, axis=0) > 0)
        n_measurable = int(np.count_nonzero(measurable))
        total = float(mean_norm.sum())
        share = mean_norm / total if total > 0 else np.zeros_like(mean_norm)
        proxy = z.mean(axis=1)
        keep_cols = np.ones(len(data['genes']), bool)
        keep_cols[[data['positions'][i] for i in cols]] = False
        non_target = data['matrix'][:, np.flatnonzero(keep_cols)]
        non_target_library = np.asarray(non_target.sum(axis=1)).ravel()
        non_target_detection = np.asarray((non_target > 0).sum(axis=1)).ravel() / max(non_target.shape[1], 1)
        covariates = np.c_[np.log1p(non_target_library), non_target_detection]
        proxy_res = regression_residual(proxy, covariates)
        ry = np.asarray(cache['y'][:, c['frozen_target_order'].index(target)], np.float64)
        summary = dict(base_row, n_in_panel_zero_counts=int(np.count_nonzero(n_any == 0)),
                       n_detected_any=int(np.count_nonzero(n_any > 0)),
                       n_measurable=n_measurable,
                       effective_gene_count=(float((share.sum() ** 2) / float((share ** 2).sum()))
                                             if total > 0 else None),
                       top_gene_share=(float(share.max()) if total > 0 else None),
                       top3_gene_share=(float(np.sort(share)[::-1][:3].sum()) if total > 0 else None),
                       det_frac_median=float(np.median(detection)),
                       proxy_mean=float(proxy.mean()), proxy_sd=float(proxy.std()),
                       readout_max_abs_diff=float(np.max(np.abs(proxy - ry))),
                       proxy_technical_r2=variance_explained(proxy, proxy_res),
                       proxy_corr_non_target_depth=corr(proxy, covariates[:, 0]),
                       proxy_corr_non_target_detection=corr(proxy, covariates[:, 1]),
                       n_non_target_genes=int(non_target.shape[1]))
        if n_measurable < 3:
            not_testable.append(dict(section_id=row['section_id'], target=target,
                                     reason='FEWER_THAN_THREE_MEASURABLE', n_panel=n_panel,
                                     n_in_panel=len(member_genes), n_measurable=n_measurable))
            summaries.append(summary)
            continue
        zm = z[:, measurable]
        zm_res = np.column_stack([regression_residual(zm[:, i], covariates)
                                  for i in range(zm.shape[1])])
        genes = [g for g, keep in zip(member_genes, measurable) if keep]
        matrix_raw = corr_matrix(zm)
        matrix_res = corr_matrix(zm_res)
        upper = np.triu_indices(zm.shape[1], 1)
        synth = dict(base_row, n_genes=int(zm.shape[1]), n_pairs=int(len(upper[0])))
        synth.update(expand('pair_raw', distribution(list(matrix_raw[upper]))))
        synth.update(expand('pair_residual', distribution(list(matrix_res[upper]))))
        tech_r2 = [variance_explained(zm[:, i], zm_res[:, i]) for i in range(zm.shape[1])]
        synth.update(expand('technical_r2', distribution(tech_r2)))
        synergy_rows.append(synth)
        summary.update({k: v for k, v in synth.items()
                        if k.startswith('pair_') or k.startswith('technical_r2')})
        diag = np.arange(zm.shape[1])
        coherence_raw = np.array([np.nanmean(matrix_raw[i, diag != i]) for i in diag]) \
            if zm.shape[1] > 1 else np.zeros(zm.shape[1])
        coherence_res = np.array([np.nanmean(matrix_res[i, diag != i]) for i in diag]) \
            if zm.shape[1] > 1 else np.zeros(zm.shape[1])
        for i, gene in enumerate(genes):
            gene_rows.append(dict(section_id=row['section_id'], target=target, gene=gene,
                                  detection_fraction=float(detection[measurable][i]),
                                  spots_detected=int(n_any[measurable][i]),
                                  mean_normalised=float(mean_norm[measurable][i]),
                                  contribution_share=float(share[measurable][i]),
                                  coherence_with_other_genes=float(coherence_raw[i]),
                                  coherence_with_other_genes_residual=float(coherence_res[i]),
                                  technical_r2=tech_r2[i]))
        # ---- leave-one-out on the measurable set
        top = int(np.argmax(share[measurable]))
        remaining = [i for i in range(zm.shape[1]) if i != top]
        reduced = zm[:, remaining].mean(axis=1)
        measurable_mean = zm.mean(axis=1)
        deletion_rows.append(dict(
            base_row, removed_gene=genes[top], removed_share=float(share[measurable][top]),
            removed_spots_detected=int(n_any[measurable][top]),
            n_genes_before=int(zm.shape[1]), n_genes_after=int(zm.shape[1] - 1),
            corr_vs_measurable_mean=corr(measurable_mean, reduced),
            corr_vs_proxy_full_definition=corr(proxy, reduced),
            jaccard_q75_vs_measurable_mean=jaccard_high(measurable_mean, reduced,
                                                        c['operational_quantile']),
            jaccard_q75_vs_proxy_full_definition=jaccard_high(proxy, reduced,
                                                             c['operational_quantile']),
            centroid_shift_pitch_vs_measurable_mean=centroid_shift_pitch(xy, measurable_mean,
                                                                        reduced, pitch),
            mean_ratio=float(reduced.mean() / proxy.mean()) if proxy.mean() else None,
            sd_ratio=float(reduced.std() / proxy.std()) if proxy.std() else None,
            corr_residual=corr(zm_res.mean(axis=1), zm_res[:, remaining].mean(axis=1)),
            state='DESCRIBED'))
        # ---- disjoint subsets: 10 fixed draws + quartile-bin matched control
        if n_measurable < 6:
            not_testable.append(dict(section_id=row['section_id'], target=target,
                                     reason='FEWER_THAN_SIX_MEASURABLE', n_panel=n_panel,
                                     n_in_panel=len(member_genes), n_measurable=n_measurable))
            summaries.append(summary)
            continue
        pool_genes = [g for g in data['present']
                      if g not in set(member_genes) and g in data['index']]
        pool_cols = [data['index'][g] for g in pool_genes
                     if measurable_u[data['index'][g]]]
        pool_by_cell = defaultdict(list)
        if det_edges is not None:
            for j, col in enumerate(pool_cols):
                cell = (int(np.searchsorted(det_edges, det_u[col])),
                        int(np.searchsorted(expr_edges, mean_u[col])))
                pool_by_cell[cell].append(j)
        observed, controls = [], []
        for draw in range(DRAWS):
            rng = np.random.default_rng(seed_for(c['seed'], row['section_id'], target,
                                                 'subset_split', draw))
            order = rng.permutation(zm.shape[1])
            half = zm.shape[1] // 2
            mask = np.zeros(zm.shape[1], bool)
            mask[order[:half]] = True
            map_a, map_b = zm[:, mask].mean(axis=1), zm[:, ~mask].mean(axis=1)
            res_a, res_b = zm_res[:, mask].mean(axis=1), zm_res[:, ~mask].mean(axis=1)
            observed.append(dict(draw=draw, corr=corr(map_a, map_b),
                                 jaccard=jaccard_high(map_a, map_b, c['operational_quantile']),
                                 corr_residual=corr(res_a, res_b),
                                 centroid=centroid_shift_pitch(xy, map_a, map_b, pitch),
                                 n_half_a=int(np.count_nonzero(mask)),
                                 n_half_b=int(np.count_nonzero(~mask))))
            control_corr = control_jac = None
            if pool_cols and det_edges is not None:
                det_m = det_u[cols][measurable]
                mean_m = mean_u[cols][measurable]
                needed = {}
                for i in np.flatnonzero(mask):
                    cell = (int(np.searchsorted(det_edges, det_m[i])),
                            int(np.searchsorted(expr_edges, mean_m[i])))
                    needed[cell] = needed.get(cell, 0) + 1
                needed_b = {}
                for i in np.flatnonzero(~mask):
                    cell = (int(np.searchsorted(det_edges, det_m[i])),
                            int(np.searchsorted(expr_edges, mean_m[i])))
                    needed_b[cell] = needed_b.get(cell, 0) + 1
                control_seed = seed_for(c['seed'], row['section_id'], target, 'control_draw', draw)
                sampled = match_cells_joint(pool_by_cell, needed, needed_b,
                                            np.random.default_rng(control_seed))
                if sampled is not None:
                    chosen_a, chosen_b = sampled
                    ctrl_a = z_u[:, [pool_cols[i] for i in chosen_a]].mean(axis=1)
                    ctrl_b = z_u[:, [pool_cols[i] for i in chosen_b]].mean(axis=1)
                    control_corr = corr(ctrl_a, ctrl_b)
                    control_jac = jaccard_high(ctrl_a, ctrl_b, c['operational_quantile'])
                    shared = len(set(chosen_a) & set(chosen_b))
                else:
                    shared = None
            controls.append(dict(draw=draw, corr=control_corr, jaccard=control_jac,
                                 evaluable=control_corr is not None, shared=shared,
                                 control_seed=control_seed))
        for obs, ctl in zip(observed, controls):
            draw_rows.append(dict(base_row, draw=obs['draw'],
                                  split_seed=seed_for(c['seed'], row['section_id'], target,
                                                      'subset_split', obs['draw']),
                                  n_half_a=obs['n_half_a'], n_half_b=obs['n_half_b'],
                                  corr=obs['corr'], jaccard_q75=obs['jaccard'],
                                  corr_residual=obs['corr_residual'], centroid_shift_pitch=obs['centroid'],
                                  control_corr=ctl['corr'], control_jaccard_q75=ctl['jaccard'],
                                  control_evaluable=ctl['evaluable'],
                                  control_shared_genes=ctl['shared'],
                                  control_seed=ctl['control_seed']))
        obs_corr = distribution([o['corr'] for o in observed])
        obs_jac = distribution([o['jaccard'] for o in observed])
        obs_res = distribution([o['corr_residual'] for o in observed])
        ctl_corr = distribution([c['corr'] for c in controls])
        ctl_jac = distribution([c['jaccard'] for c in controls])
        subset_rows.append(dict(
            base_row, n_draws=DRAWS, n_measurable=n_measurable,
            n_evaluable_control_draws=int(sum(1 for c_ in controls if c_['evaluable'])),
            corr_median=obs_corr['median'], corr_q10=obs_corr['q10'], corr_q90=obs_corr['q90'],
            jaccard_q75_median=obs_jac['median'],
            corr_residual_median=obs_res['median'],
            centroid_shift_pitch_median=distribution([o['centroid'] for o in observed])['median'],
            control_corr_median=ctl_corr['median'], control_corr_q10=ctl_corr['q10'],
            control_corr_q90=ctl_corr['q90'], control_jaccard_q75_median=ctl_jac['median'],
            pool_measurable=len(pool_cols), state='DESCRIBED'))
        summaries.append(summary)
    base.tsv_write(out / 'section_target_summary.tsv', summaries)
    base.tsv_write(out / 'synergy.tsv', synergy_rows)
    base.tsv_write(out / 'gene_coherence.tsv', gene_rows)
    base.tsv_write(out / 'deletion_sensitivity.tsv', deletion_rows)
    base.tsv_write(out / 'subset_consistency.tsv', subset_rows)
    base.tsv_write(out / 'subset_draws.tsv', draw_rows)
    base.tsv_write(out / 'not_testable.tsv', not_testable)
    base.json_write(out / 'section_receipt.json',
                    dict(section_id=row['section_id'], unit_id=row['unit_id'],
                         n_spots=int(len(cache['barcodes'])), seconds=time.time() - started,
                         receipt=ctx['guard'].receipt()))
    (out / 'COMPLETE').write_text('SECTION_AUDIT_DONE\n')
    return dict(sid=row['section_id'], seconds=round(time.time() - started, 1),
                targets=len(summaries))


# --------------------------------------------------------------------------- control recompute
CONTROL_DRAW_COLUMNS = ['control_corr', 'control_jaccard_q75', 'control_evaluable',
                        'control_shared_genes', 'control_seed']
CONTROL_SUMMARY_COLUMNS = ['control_corr_median', 'control_corr_q10', 'control_corr_q90',
                           'control_jaccard_q75_median', 'n_evaluable_control_draws']


def recompute_controls_section(job):
    """Recompute only the matched-control columns under the D-277 no-reuse rule.

    The target's 10 splits, the bins, the technical regression and every other column are
    untouched: the target half masks are regenerated from the same seeds and asserted against
    the stored half sizes before any control value is replaced.
    """
    ctx = setup(verify='none')
    c, panels = ctx['c'], ctx['panels']
    row = ctx['byid'][job['sid']]
    out = Path(job['out'])
    draws_path, subset_path = out / 'subset_draws.tsv', out / 'subset_consistency.tsv'
    if not draws_path.exists() or not subset_path.exists():
        return dict(sid=job['sid'], status='no_subset_rows')
    draws = list(csv.DictReader(draws_path.open(), delimiter='\t'))
    subsets = list(csv.DictReader(subset_path.open(), delimiter='\t'))
    if not draws:
        return dict(sid=job['sid'], status='no_draws')
    union_genes = sorted({g for t_ in c['targets'] for g in panels[t_]})
    data = read_section(row, union_genes)
    cache = base.load_cache(row)
    review.assert_aligned(data['barcodes'], cache['barcodes'])
    z_u = data['z_union']
    counts_u = data['matrix'][:, data['positions']].tocsr()
    n_any_u = np.asarray((counts_u > 0).sum(axis=0)).ravel()
    det_u = np.asarray((counts_u > 0).mean(axis=0)).ravel()
    mean_u = z_u.mean(axis=0)
    measurable_u = (n_any_u >= 20) & (np.std(z_u, axis=0) > 0)
    union_stats = [(g, float(det_u[i]), float(mean_u[i]))
                   for i, g in enumerate(data['present']) if measurable_u[i]]
    det_edges, expr_edges = build_bins(union_stats)
    changed = 0
    for target in c['targets']:
        tdraws = [r for r in draws if r['target'] == target]
        tsubset = [r for r in subsets if r['target'] == target]
        if not tdraws or not tsubset:
            continue
        members = [g for g in panels[target] if g in data['index']]
        cols = [data['index'][g] for g in members]
        n_any = n_any_u[cols]
        measurable = (n_any >= 20) & (np.std(z_u[:, cols], axis=0) > 0)
        zm = z_u[:, cols][:, measurable]
        det_m = det_u[cols][measurable]
        mean_m = mean_u[cols][measurable]
        pool_genes = [g for g in data['present'] if g not in set(members)]
        pool_cols = [data['index'][g] for g in pool_genes if measurable_u[data['index'][g]]]
        pool_by_cell = defaultdict(list)
        if det_edges is not None:
            for j, col in enumerate(pool_cols):
                cell = (int(np.searchsorted(det_edges, det_u[col])),
                        int(np.searchsorted(expr_edges, mean_u[col])))
                pool_by_cell[cell].append(j)
        for r in tdraws:
            draw = int(r['draw'])
            rng = np.random.default_rng(seed_for(c['seed'], row['section_id'], target,
                                                 'subset_split', draw))
            order = rng.permutation(zm.shape[1])
            half = zm.shape[1] // 2
            mask = np.zeros(zm.shape[1], bool)
            mask[order[:half]] = True
            if int(np.count_nonzero(mask)) != int(r['n_half_a']) or \
                    int(np.count_nonzero(~mask)) != int(r['n_half_b']):
                raise ValueError(f'Target split changed while recomputing controls: '
                                 f'{row["section_id"]}/{target}/{draw}')
            needed_a, needed_b = {}, {}
            for i in np.flatnonzero(mask):
                cell = (int(np.searchsorted(det_edges, det_m[i])),
                        int(np.searchsorted(expr_edges, mean_m[i])))
                needed_a[cell] = needed_a.get(cell, 0) + 1
            for i in np.flatnonzero(~mask):
                cell = (int(np.searchsorted(det_edges, det_m[i])),
                        int(np.searchsorted(expr_edges, mean_m[i])))
                needed_b[cell] = needed_b.get(cell, 0) + 1
            control_seed = seed_for(c['seed'], row['section_id'], target, 'control_draw', draw)
            sampled = match_cells_joint(pool_by_cell, needed_a, needed_b,
                                        np.random.default_rng(control_seed)) \
                if (pool_cols and det_edges is not None) else None
            if sampled is None:
                r['control_corr'] = ''
                r['control_jaccard_q75'] = ''
                r['control_evaluable'] = 'False'
                r['control_shared_genes'] = ''
            else:
                chosen_a, chosen_b = sampled
                ctrl_a = z_u[:, [pool_cols[i] for i in chosen_a]].mean(axis=1)
                ctrl_b = z_u[:, [pool_cols[i] for i in chosen_b]].mean(axis=1)
                r['control_corr'] = corr(ctrl_a, ctrl_b)
                r['control_jaccard_q75'] = jaccard_high(ctrl_a, ctrl_b,
                                                        c['operational_quantile'])
                r['control_evaluable'] = 'True'
                r['control_shared_genes'] = 0
            r['control_seed'] = control_seed
            changed += 1
        corr_values = [num_ for num_ in (float(r['control_corr']) for r in tdraws
                                         if r['control_corr'] not in ('', 'None'))] \
            if any(r['control_corr'] not in ('', 'None') for r in tdraws) else []
        corr_stats = distribution(corr_values)
        jac_values = [float(r['control_jaccard_q75']) for r in tdraws
                      if r['control_jaccard_q75'] not in ('', 'None')]
        for s in tsubset:
            s['control_corr_median'] = corr_stats['median']
            s['control_corr_q10'] = corr_stats['q10']
            s['control_corr_q90'] = corr_stats['q90']
            s['control_jaccard_q75_median'] = (distribution(jac_values)['median']
                                               if jac_values else '')
            s['n_evaluable_control_draws'] = sum(1 for r in tdraws
                                                 if r['control_evaluable'] == 'True')
    base.tsv_write(draws_path, draws)
    base.tsv_write(subset_path, subsets)
    return dict(sid=job['sid'], status='controls_recomputed', draws=changed)


def recompute_controls(sids=None, workers=None):
    ctx = setup(verify='none')
    c = ctx['c']
    root = WORK / 'sections'
    jobs = [dict(sid=s, out=str(root / s)) for s in (sids or sorted(ctx['byid']))]
    with multiprocessing.get_context('spawn').Pool(processes=workers or c['workers']) as pool:
        for info in pool.imap_unordered(recompute_controls_section, jobs):
            print(json.dumps(info), flush=True)
    (WORK / 'CONTROLS_RECOMPUTED').write_text('D277_NO_REUSE_CONTROLS\n')


# --------------------------------------------------------------------------- commands
def prepare():
    ctx = setup(verify='all')
    if (WORK / 'PREPARED').exists():
        raise RuntimeError('Prepared outputs exist; preserve this version')
    base.json_write(WORK / 'input_checks.json', ctx['checks'])
    base.json_write(WORK / 'prepare_receipt.json', dict(receipt=ctx['guard'].receipt()))
    (WORK / 'PREPARED').write_text('FROZEN_INPUTS_VERIFIED\n')
    print(json.dumps(ctx['checks'], indent=1), flush=True)


def audit(sids=None, workers=None, out_root=None, force=False):
    ctx = setup(verify='none')
    c = ctx['c']
    root = Path(out_root) if out_root else WORK / 'sections'
    root.mkdir(parents=True, exist_ok=True)
    jobs = [dict(sid=s, out=str(root / s), force=force) for s in (sids or sorted(ctx['byid']))]
    with multiprocessing.get_context('spawn').Pool(processes=workers or c['workers']) as pool:
        for info in pool.imap_unordered(audit_section, jobs):
            print(json.dumps(info), flush=True)
    if not out_root:
        (WORK / 'AUDIT_COMPLETE').write_text('ALL_SELECTED_SECTIONS_DONE\n')


def collect(root):
    collected = {name: [] for name in TABLES}
    for d in sorted(Path(root).iterdir()):
        if not d.is_dir() or not (d / 'COMPLETE').exists():
            continue
        for name in TABLES:
            path = d / name
            if path.exists():
                collected[name] += list(csv.DictReader(path.open(), delimiter='\t'))
    return collected


def summarize(root=None, out=None, arms=True):
    ctx = setup(verify='none')
    c = ctx['c']
    root = Path(root) if root else WORK / 'sections'
    out = Path(out) if out else WORK
    collected = collect(root)
    for name, rows in collected.items():
        base.tsv_write(out / name, rows)
    per_section = [dict(r) for r in collected['section_target_summary.tsv']]
    approved = {(r['section_id'], r['target']): r
                for r in csv.DictReader(SUPPORT.open(), delimiter='\t')}
    for row in per_section:
        a = approved.get((row['section_id'], row['target']), {})
        row['approved_support_source'] = str(SUPPORT.relative_to(ROOT))
        row['approved_genes_in_panel'] = a.get('sup_genes_in_panel')
        row['approved_genes_detected'] = a.get('sup_genes_detected')
        row['approved_effective_gene_count'] = a.get('sup_effective_gene_count')
        row['approved_top_gene_share'] = a.get('sup_top_gene_share')
        row['approved_spots_with_any_gene_frac'] = a.get('sup_spots_with_any_gene_frac')
    deletion = {(r['section_id'], r['target']): r for r in collected['deletion_sensitivity.tsv']}
    subset = {(r['section_id'], r['target']): r for r in collected['subset_consistency.tsv']}
    for row in per_section:
        d = deletion.get((row['section_id'], row['target']))
        s = subset.get((row['section_id'], row['target']))
        row['deletion_corr_vs_measurable_mean'] = d.get('corr_vs_measurable_mean') if d else None
        row['deletion_jaccard_q75'] = d.get('jaccard_q75_vs_measurable_mean') if d else None
        row['subset_corr_median'] = s.get('corr_median') if s else None
        row['subset_jaccard_q75_median'] = s.get('jaccard_q75_median') if s else None
        row['control_corr_median'] = s.get('control_corr_median') if s else None
    base.tsv_write(out / 'section_metrics.tsv', per_section)
    units = aggregate_units(per_section)
    base.tsv_write(out / 'unit_metrics.tsv', units)
    if arms:
        macro = unit_macro(units, c['targets'])
        base.tsv_write(out / 'unit_macro.tsv', macro)
        strata = unit_strata(units, c['targets'])
        base.tsv_write(out / 'unit_strata.tsv', strata)
    write_report(per_section, collected, units, c, out, arms)
    (out / 'STATUS').write_text('PROXY_AUDIT_SEGMENT1_COMPLETE_PENDING_REVIEW\n')
    print(json.dumps(dict(sections=len({r['section_id'] for r in per_section}),
                          section_rows=len(per_section),
                          synergy_rows=len(collected['synergy.tsv']),
                          gene_rows=len(collected['gene_coherence.tsv']),
                          deletion_rows=len(collected['deletion_sensitivity.tsv']),
                          subset_rows=len(collected['subset_consistency.tsv']),
                          draw_rows=len(collected['subset_draws.tsv']),
                          not_testable=len(collected['not_testable.tsv']),
                          unit_rows=len(units))), flush=True)


def aggregate_units(per_section):
    groups = defaultdict(list)
    for r in per_section:
        groups[(r['unit_id'], r['target'])].append(r)
    rows = []
    for (unit, target), group in sorted(groups.items()):
        sections = sorted({g['section_id'] for g in group})
        row = dict(unit_id=unit, target=target, n_sections=len(sections),
                   section_ids=';'.join(sections), organ=group[0]['organ'],
                   platform=group[0]['platform'], source_key=group[0]['source_key'],
                   source_eligibility=group[0]['source_eligibility'],
                   confirmed_unit=group[0]['confirmed_unit'])
        for metric in METRICS:
            values = [review.numeric(g.get(metric)) for g in group]
            present = [v for v in values if v is not None]
            row[f'{metric}__n'] = len(present)
            if not present:
                row[metric] = None
                row[f'{metric}__state'] = 'NOT_ESTIMABLE_UNIT'
            elif len(present) < len(sections):
                row[metric] = None
                row[f'{metric}__state'] = 'INCOMPLETE_UNIT'
            else:
                row[metric] = float(np.mean(present))
                row[f'{metric}__state'] = 'COMPLETE'
        rows.append(row)
    return rows


def unit_macro(units, targets):
    rows = []
    for target in targets:
        for metric in METRICS:
            for scope, selector in [('internal', lambda u: True),
                                    ('source', lambda u: not str(u['source_eligibility'])
                                     .startswith('NOT_TESTABLE'))]:
                group = [u for u in units if u['target'] == target and selector(u)]
                present = [u[metric] for u in group if u[metric] is not None]
                rows.append(dict(target=target, metric=metric, scope=scope, n_units=len(group),
                                 n_complete=len(present),
                                 mean=(float(np.mean(present)) if present else None),
                                 median=(float(np.median(present)) if present else None),
                                 q10=(float(np.quantile(present, 0.10)) if present else None),
                                 q90=(float(np.quantile(present, 0.90)) if present else None)))
    return rows


def unit_strata(units, targets):
    rows = []
    for key in ['source_key', 'organ', 'platform', 'confirmed_unit', 'source_eligibility']:
        for target in targets:
            groups = defaultdict(list)
            for u in units:
                if u['target'] == target:
                    groups[u[key]].append(u)
            for value, group in sorted(groups.items(), key=lambda kv: str(kv[0])):
                for metric in ['pair_raw.median', 'pair_residual.median', 'n_measurable',
                               'top_gene_share', 'subset_corr_median', 'control_corr_median']:
                    present = [u[metric] for u in group if u[metric] is not None]
                    rows.append(dict(stratum=key, value=str(value), target=target, metric=metric,
                                     n_units=len(group), n_complete=len(present),
                                     mean=(float(np.mean(present)) if present else None),
                                     median=(float(np.median(present)) if present else None)))
    return rows


def write_report(per_section, collected, units, c, out, arms):
    def med(rows, field):
        vals = [review.numeric(r.get(field)) for r in rows]
        vals = [v for v in vals if v is not None]
        return float(np.median(vals)) if vals else None

    def show(value, digits=3):
        return 'NA' if value is None else format(value, f'.{digits}f')

    lines = ['# 四个代理场体检（第一段，描述性）', '',
             '只描述代理集合的内部一致性。**不把片内同向写成生物协同、独立病理真值或 HF-01**；'
             '不命名真实场；高值集合是操作阈值（严格大于 q75），不是生物真值。', '',
             f"覆盖：{len({r['section_id'] for r in per_section})} 片 × {len(c['targets'])} 目标；"
             f"单位 {len({r['unit_id'] for r in units})}（跨目标）；内部 72 片/65 单位、"
             '来源 53 片/47 单位；来源不足 19 片仅内部；肺不可作来源；乳腺/肝脏为器官外推。', '',
             '## 1. 代理支持与可测基因', '',
             '| 目标 | 片数 | 面板内中位 | 任意检出中位（已审口径） | 可测基因中位（≥20 点且方差>0） | 已审有效基因数中位 | 已审首位份额中位 | 读出加权首位份额中位 | 原读数复算最大差 |',
             '|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for target in c['targets']:
        rows = [r for r in per_section if r['target'] == target]
        if not rows:
            continue
        worst = max([review.numeric(r['readout_max_abs_diff']) for r in rows
                     if review.numeric(r['readout_max_abs_diff']) is not None])
        lines.append(f"| {target} | {len(rows)} | {show(med(rows, 'n_in_panel'), 0)} | "
                     f"{show(med(rows, 'approved_genes_detected'), 0)} | "
                     f"{show(med(rows, 'n_measurable'), 0)} | "
                     f"{show(med(rows, 'approved_effective_gene_count'), 2)} | "
                     f"{show(med(rows, 'approved_top_gene_share'))} | "
                     f"{show(med(rows, 'top_gene_share'))} | {worst:.1e} |")
    lines += ['', '上述“已审”列直接取自 `infra/field_biological_review_20261019/measured_fields.tsv`'
              '（来源逐行记录在 `approved_support_source`），其有效基因数与首位份额按原始计数总额加权；'
              '“读出加权”列按本表统一读出（log1p 均值）计算，两者定义不同、分列报告，不混用。', '',
              '## 2. 片内基因空间协同（基因对相关分布；原图 vs 单次片内回归残差）', '',
              '| 目标 | 基因对中位 | 原图中位 | 原图 q10–q90 | 残差中位 | 残差 q10–q90 | 技术可解释方差中位 |',
              '|---|---:|---:|---|---:|---|---:|']
    for target in c['targets']:
        rows = [r for r in per_section if r['target'] == target and r.get('pair_raw.n')]
        if not rows:
            continue
        lines.append(f"| {target} | {show(med(rows, 'pair_raw.n'), 0)} | "
                     f"{show(med(rows, 'pair_raw.median'))} | "
                     f"{show(med(rows, 'pair_raw.q10'))}–{show(med(rows, 'pair_raw.q90'))} | "
                     f"{show(med(rows, 'pair_residual.median'))} | "
                     f"{show(med(rows, 'pair_residual.q10'))}–{show(med(rows, 'pair_residual.q90'))} | "
                     f"{show(med(rows, 'technical_r2.median'))} |")
    lines += ['', '技术对照：每张可测基因图**一次**片内回归（截距 + 排除本目标基因后的点级文库 log1p + '
              '非目标检出比例）。该回归也会去掉与细胞密度耦合的生物信号，残差相似或消失均不能单独判定'
              '生物/技术来源；共同的全基因归一化分母仍可能引入依赖。', '',
              '## 3. 删除最大贡献基因（可测集合内）', '',
              '| 目标 | 片数 | 被删份额中位 | 与可测均值相关中位 | q75 Jaccard 中位 | 质心位移中位(pitch) | 剩余基因中位 |',
              '|---|---:|---:|---:|---:|---:|---:|']
    deletion = collected['deletion_sensitivity.tsv']
    for target in c['targets']:
        rows = [r for r in per_section if r['target'] == target
                and r.get('deletion_corr_vs_measurable_mean') is not None]
        drs = [r for r in deletion if r['target'] == target]
        if not rows:
            continue
        lines.append(f"| {target} | {len(drs)} | {show(med(drs, 'removed_share'))} | "
                     f"{show(med(rows, 'deletion_corr_vs_measurable_mean'))} | "
                     f"{show(med(rows, 'deletion_jaccard_q75'))} | "
                     f"{show(med(drs, 'centroid_shift_pitch_vs_measurable_mean'))} | "
                     f"{show(med(drs, 'n_genes_after'), 0)} |")
    lines += ['', '与“原代理全定义图”的相关另记在 `deletion_sensitivity.tsv`'
              '（`corr_vs_proxy_full_definition`）：基因集合不同、公式相同，单独报告。', '',
              '## 4. 互斥子集一致性（10 次固定种子全部报告 + 四分位档匹配的非目标对照）', '',
              '| 目标 | 片数 | 两半相关中位 | q10–q90（片间） | q75 Jaccard 中位 | 匹配对照相关中位 | 可评对照次数中位 |',
              '|---|---:|---:|---|---:|---:|---:|']
    subset_rows = collected['subset_consistency.tsv']
    for target in c['targets']:
        rows = [r for r in per_section if r['target'] == target and r.get('subset_corr_median') is not None]
        srs = [r for r in subset_rows if r['target'] == target]
        if not rows:
            continue
        lines.append(f"| {target} | {len(rows)} | {show(med(rows, 'subset_corr_median'))} | "
                     f"{show(med(srs, 'corr_q10'))}–{show(med(srs, 'corr_q90'))} | "
                     f"{show(med(rows, 'subset_jaccard_q75_median'))} | "
                     f"{show(med(rows, 'control_corr_median'))} | "
                     f"{show(med(srs, 'n_evaluable_control_draws'), 1)} |")
    zero_draws = [r for r in subset_rows if int(r['n_evaluable_control_draws']) == 0]
    lines += ['', '**对照可评覆盖（必须与上表同读，且不得跨目标合并）**', '',
              '| 目标 | 切片×目标 | 有可评对照 | 无可评对照 | 可比切片里观测更高 | 观测减对照中位 |',
              '|---|---:|---:|---:|---:|---:|']
    for target in c['targets']:
        rows = [r for r in subset_rows if r['target'] == target]
        if not rows:
            continue
        comparable = [r for r in rows
                      if review.numeric(r.get('control_corr_median')) is not None]
        diffs = [review.numeric(r['corr_median']) - review.numeric(r['control_corr_median'])
                 for r in comparable]
        lines.append(f"| {target} | {len(rows)} | {len(comparable)} | {len(rows)-len(comparable)} | "
                     f"{(100*np.mean([d > 0 for d in diffs])) if diffs else float('nan'):.0f}% | "
                     f"{(np.median(diffs)) if diffs else float('nan'):+.3f} |")
    lines += ['', f'全体 288 个切片×目标中 {len(zero_draws)} 个没有任何可评对照抽取'
              '（四分位档在本片目标外可测基因中候选不足，逐次记录在 `subset_draws.tsv`）。'
              '**每个目标单独读；不给跨目标合并的差，也不把该差当显著性、生物参照或 HF-01。**'
              '单次对照内两半不重复取同一基因（`control_shared_genes` 全为 0）。', '',
              '每个互斥子集至少 3 个可测基因（合计 ≥6），否则该片该目标为 NOT_ESTIMABLE'
              '（见 `not_testable.tsv`）；对照档位候选不足只让**该次抽取**不可评，'
              '逐次计数保存在 `subset_draws.tsv`，不得改档位或补抽。', '',
              '## 5. 单位与分母（完整单位等权；任一片不可评则该单位该指标不进均值）', '']
    for metric in METRICS:
        sel = [u for u in units if metric in u and u[metric] is not None]
        if not sel:
            continue
        lines.append(f"- {metric}：完整单位 {len(sel)}/{len(units)}，单位均值 "
                     f"{np.mean([u[metric] for u in sel]):.3f}")
    lines += ['', '### 6 个可辨来源（仅 ELIGIBLE；乳腺/肝脏为器官外推）', '',
              '| 来源 | 目标 | 单位数 | 完整 | 两半相关中位 | 匹配对照中位 |', '|---|---|---:|---:|---:|---:|']
    known = [r for r in units if r['source_eligibility'] == 'ELIGIBLE']
    grouped = defaultdict(list)
    for r in known:
        grouped[(r['source_key'], r['target'])].append(r)
    for (source, target), group in sorted(grouped.items()):
        present = [review.numeric(g.get('subset_corr_median')) for g in group]
        present = [v for v in present if v is not None]
        control = [review.numeric(g.get('control_corr_median')) for g in group]
        control = [v for v in control if v is not None]
        organ_flag = '（器官外推）' if source in (
            'https://www.biorxiv.org/content/10.1101/2022.11.09.515839v1.full',
            'https://www.nature.com/articles/s41467-023-43458-x') else ''
        lines.append(f"| {source[:52]}{organ_flag} | {target} | {len(group)} | {len(present)} | "
                     f"{'NA' if not present else format(float(np.median(present)), '.3f')} | "
                     f"{'NA' if not control else format(float(np.median(control)), '.3f')} |")
    lines += ['', '标题型 `source_key`（`title:...`）属来源身份不足的 19 片，**只留在内部描述**，'
              '不得读作第 7 个来源；肺切片全部在其中，故肺来源不可评。', '',
              '完整分母与分层见 `unit_metrics.tsv`（`__n`、`__state`）、`unit_macro.tsv`、'
              '`unit_strata.tsv`；逐片见 `section_metrics.tsv`。', '',
              '## 6. 不可评价记录', '',
              f"- 共 {len(collected['not_testable.tsv'])} 条：" +
              json.dumps(dict(Counter(r['reason'] for r in collected['not_testable.tsv'])),
                         ensure_ascii=False),
              '- 结构化不可评价（不是生物学阴性）：无 panel 基因；可测基因 <3（协同/留一）；<6（两半各需 ≥3）；'
              '常量图；技术残差退化；对照档位候选不足（逐次记录）；条码/坐标不一致（硬失败）。', '',
              '## 7. 边界与第二段', '',
              '- 描述性：无 PASS 阈值、无显著性、无独立病理或生物真值；片内同向 ≠ 生物协同 ≠ HF-01。',
              '- 两半互斥、另配匹配对照，避免同一批基因既构图又自证。',
              '- **第一段不选定第二段的基因成员、权重或候选数**；第二段未授权，本报告不提供拆分结论。', '']
    (out / 'REPORT.md').write_text('\n'.join(lines) + '\n')


def smoke(sids=('MISC33', 'ZEN38')):
    global WORK
    main = WORK
    smoke_dir = main / 'smoke'
    smoke_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(main / 'contract.json', smoke_dir / 'contract.json')
    WORK = smoke_dir
    ctx = setup(verify='sample', sids=sids)
    base.json_write(smoke_dir / 'input_checks.json', ctx['checks'])
    (smoke_dir / 'PREPARED').write_text('SMOKE_ONLY\n')
    started = time.time()
    audit(sids=list(sids), out_root=smoke_dir / 'sections', workers=2)
    summarize(root=smoke_dir / 'sections', out=smoke_dir, arms=False)
    WORK = main
    print(json.dumps(dict(smoke_sections=list(sids), seconds=round(time.time() - started, 1))),
          flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['prepare', 'audit', 'summarize', 'smoke',
                                             'recompute-controls'])
    parser.add_argument('--sections', nargs='*', default=None)
    parser.add_argument('--workers', type=int, default=None)
    parser.add_argument('--force', action='store_true')
    args = parser.parse_args()
    if os.environ.get('CUDA_VISIBLE_DEVICES') not in (None, '', '-1'):
        raise RuntimeError('GPU disabled')
    if args.command == 'prepare':
        prepare()
    elif args.command == 'audit':
        audit(sids=args.sections, workers=args.workers, force=args.force)
    elif args.command == 'summarize':
        summarize()
    elif args.command == 'recompute-controls':
        recompute_controls(sids=args.sections, workers=args.workers)
    else:
        smoke()
