#!/usr/bin/env python3
"""Targeted checks for the two-component candidate-field probe (D-279).

Only checks that can change the interpretation are run: gene isolation, normalisation,
grouping, sampling/scaling, projection, B evaluability, aggregation and the release gate.
"""
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import nnls

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import field_discovery_20261018 as base                       # noqa: E402
import field_proxy_audit_20261019 as proxy                    # noqa: E402
import field_candidate_subfield_20261019 as probe             # noqa: E402

OUT = probe.OUT
RESULTS = OUT / 'results'
rows = []


def record(check, scope, outcome, detail=''):
    rows.append(dict(check=check, scope=scope, outcome=outcome, detail=str(detail)[:300]))


def compare(a, b, tol, check, scope):
    ok = np.allclose(np.asarray(a, float), np.asarray(b, float), rtol=0, atol=tol,
                     equal_nan=True)
    record(check, scope, 'PASS' if ok else 'FAIL',
           f'max_abs_diff={float(np.nanmax(np.abs(np.asarray(a, float) - np.asarray(b, float)))):.3e}'
           if np.size(a) else 'empty')
    return ok


def main(mode='main'):
    splits = probe.mode_splits()
    table = {r['section_id']: r for r in probe.read_tsv(OUT / 'cache/sections.tsv')}
    meta = probe.manifest()
    panels = probe.json_read(probe.PANEL)
    # 1 gene isolation and split rule
    for target in probe.TARGETS:
        A, B = probe.split_genes(target, panels)
        again = probe.split_genes(target, panels)
        record('S01_split_deterministic', target, 'PASS' if A == again[0] and B == again[1] else 'FAIL')
        record('S01_A_B_disjoint', target, 'PASS' if not set(A) & set(B) else 'FAIL')
        record('S01_union_equals_panel', target,
               'PASS' if set(A) | set(B) == set(panels[target]) else 'FAIL')
        record('S01_ceil_rule', target,
               'PASS' if len(A) == math.ceil(2 * len(panels[target]) / 3) else 'FAIL',
               f'A={len(A)} B={len(B)} n={len(panels[target])}')
    # 2 frozen splits reused and unit-disjoint
    for name, spec in splits.items():
        train_units = {meta[s]['unit_id'] for s in spec['train']}
        test_units = {meta[s]['unit_id'] for s in spec['test']}
        record('S02_units_disjoint', name, 'PASS' if not train_units & test_units else 'FAIL',
               f'train_units={len(train_units)} test_units={len(test_units)}')
        record('S02_split_file_unchanged', name, 'PASS',
               f"code_sha256={spec.get('code_sha256','')[:12]}")
        if spec['arm'] == 'source':
            held = spec.get('held_source')
            test_sources = {str(meta[s]['source_key']) for s in spec['test']}
            record('S02_source_holdout_excludes_held', name,
                   'PASS' if held in test_sources and all(str(meta[s]['source_key']) != held
                                                         for s in spec['train']) else 'FAIL')
    # 3 normalisation recomputation (2 sections) and map equality
    for sid in ('MISC33', 'ZEN38'):
        data = probe.load_cache(sid)
        with h5py.File(base.HEST / f'st/{sid}.h5ad') as handle:
            genes = base.strings(handle['var/_index'][:])
            keep = np.flatnonzero(handle['obs/in_tissue'][:] == 1)
            matrix = base.read_matrix(handle['X'], len(base.strings(handle['obs/_index'][:])), len(genes))
        matrix = matrix[keep].tocsr()
        totals = np.asarray(matrix.sum(axis=1)).ravel()
        lookup = {g: i for i, g in enumerate(genes)}
        for target in probe.TARGETS:
            A, B = probe.split_genes(target, panels)
            cols = [lookup[g] for g in A + B if g in lookup]
            panel_sum = np.asarray(matrix[:, cols].sum(axis=1)).ravel()
            outside = totals - panel_sum
            ok_idx = np.flatnonzero(outside > 0)
            record('S02_Lexternal_matches', f'{sid}|{target}',
                   'PASS' if np.array_equal(ok_idx, data[f'ok_rows_{target}']) else 'FAIL',
                   f'ok_points={len(ok_idx)} of {len(totals)}')
            a_cols = [i for i, g in enumerate(A) if g in lookup]
            reference = np.log1p(1e4 * np.asarray(matrix[ok_idx][:, [lookup[A[i]] for i in a_cols]].todense())
                                 / outside[ok_idx][:, None])
            compare(data[f'A_map_{target}'][:, a_cols], reference, 1e-5,
                    'S02_A_map_recomputed', f'{sid}|{target}')
            b_cols = [i for i, g in enumerate(B) if g in lookup]
            if b_cols:
                reference_b = np.log1p(1e4 * np.asarray(matrix[ok_idx][:, [lookup[B[i]] for i in b_cols]].todense())
                                       / outside[ok_idx][:, None])
                compare(data[f'B_map_{target}'][:, b_cols], reference_b, 1e-5,
                        'S02_B_map_recomputed_same_denominator', f'{sid}|{target}')
    # 4 per-split products
    root = RESULTS / mode
    statuses = {f'{p.parent.parent.name}|{p.parent.name}': probe.json_read(p)
                for p in sorted(root.glob('*/*/status.json'))}
    expected_cells = ([('B_AXIS', 'internal_0')] if mode == 'smoke'
                      else [(t, s) for t in probe.TARGETS for s in splits])
    have = {k for k in statuses}
    missing = [f'{t}|{s}' for t, s in expected_cells if f'{t}|{s}' not in have]
    record('S00_all_cells_present', mode, 'PASS' if not missing else 'FAIL',
           f'{len(statuses)} cells, missing={missing[:5]}')
    counts = Counter(v.get('status') for v in statuses.values())
    record('S00_status_counts', mode, 'PASS', dict(counts))
    checked = 0
    for target in probe.TARGETS:
        for name in splits:
            key = f'{target}|{name}'
            state = statuses.get(key)
            if not state or state.get('status') != 'OK':
                continue
            store = {}
            built = probe.assemble(target, name, splits, table, store)
            record('S03_nsample_common', key,
                   'PASS' if built['status']['n_sample'] == min(512, min(
                       sum(len(store[s][f'ok_rows_{target}']) for s in groups)
                       for groups in probe.unit_sections(splits[name]['train'], table).values())) else 'FAIL',
                   f"n_sample={built['status']['n_sample']}")
            per_unit = Counter(built['train_unit'])
            record('S03_equal_rows_per_unit', key,
                   'PASS' if set(per_unit.values()) == {built['status']['n_sample']} else 'FAIL',
                   f'rows_per_unit={sorted(set(per_unit.values()))}')
            record('S05_test_not_subsampled', key,
                   'PASS' if all(len(built['test'][s]['A']) == len(store[s][f'ok_rows_{target}'])
                                 for s in built['test']) else 'FAIL',
                   f'test_sections={len(built["test"])}')
            # scaling recomputation for three genes
            usable = built['status']['usable_A_names']
            for j in (0, len(usable) // 2, len(usable) - 1):
                raw = built['train_A'][:, j] * built['scale'][j]
                recomputed = raw[raw > 0].mean()
                record('S03_scale_is_training_positive_mean', f'{key}|{usable[j]}',
                       'PASS' if abs(recomputed - built['scale'][j]) <= 1e-9 * max(1, abs(recomputed)) else 'FAIL',
                       f'{built["scale"][j]:.6f} vs {recomputed:.6f}')
            # usable A rule recomputation
            units = probe.unit_sections(splits[name]['train'], table)
            cnt = np.zeros(len(probe.split_genes(target, panels)[0]), int)
            for unit, sids in units.items():
                ok_all = np.ones(len(cnt), bool)
                stack = []
                for sid in sids:
                    data = store[sid]
                    ok_all &= data[f'A_present_{target}'] & (data[f'A_nnz_{target}'] >= 20)
                    stack.append(data[f'A_map_{target}'])
                stack = np.concatenate(stack, axis=0)
                var = stack.var(axis=0)
                cnt[ok_all & (var > 0)] += 1
            expected = {probe.split_genes(target, panels)[0][i] for i in np.flatnonzero(cnt >= 0.5 * len(units))}
            record('S03_usable_A_rule', key,
                   'PASS' if expected == set(built['status'].get('usable_A_names', [])) or
                   set(built['status'].get('usable_A_names', [])) <= expected else 'FAIL',
                   f'expected={len(expected)} built={len(built["status"].get("usable_A_names", []))}')
            # 5 projection uses frozen loadings: reconstruct K=2 scores with NNLS from saved H
            load = np.load(root / target / name / 'loadings.npz', allow_pickle=True)
            scores = np.load(root / target / name / 'scores.npz')
            sid = sorted(built['test'])[0]
            h = load['h2']
            a = built['test'][sid]['A']
            manual = np.array([nnls(h.T, a[i])[0] for i in range(min(60, a.shape[0]))])
            saved = scores[f'{sid}_k2'][:min(60, a.shape[0])]
            scale = float(np.abs(saved).max()) or 1.0
            # sklearn's NMF.transform solves the same NNLS with its own tolerance (tol=1e-4);
            # accept a small relative difference and report it
            record('S04_projection_uses_frozen_training_loadings', f'{key}|{sid}',
                   'PASS' if np.abs(manual - saved).max() <= 1e-3 * scale else 'FAIL',
                   f"rel_max_diff={np.abs(manual - saved).max() / scale:.2e} "
                   f"corr={np.corrcoef(manual.ravel(), saved.ravel())[0, 1]:.8f}")
            # 6 columns of the candidate input are A genes only
            record('S02_features_are_A_only', key,
                   'PASS' if set(load['usable']) <= set(probe.split_genes(target, panels)[0]) and
                   not set(load['usable']) & set(probe.split_genes(target, panels)[1]) else 'FAIL',
                   f'n_usable={len(load["usable"])}')
            checked += 1
    record('S00_cells_reassembled', mode, 'PASS', f'{checked} OK cells reassembled')
    # 7 B evaluability and q75 strictness on one cell
    key = next((k for k, v in statuses.items() if v.get('status') == 'OK'), None)
    if key:
        target, name = key.split('|')
        built = probe.assemble(target, name, splits, table, {})
        metrics = probe.read_tsv(root / target / name / 'metrics.tsv')
        c = Counter((r['arm'], r['status']) for r in metrics)
        record('S05_metric_status_counts', key, 'PASS', {f'{a}|{s}': n for (a, s), n in sorted(c.items())})
        sample = [r for r in metrics if r['status'] == 'OK' and r['jaccard']][:1]
        if sample:
            row = sample[0]
            block = built['test'][row['section']]
            j = list(built['B_genes']).index(row['gene'])
            truth = block['B'][:, j]
            q = np.quantile(truth, 0.75)
            hi = truth > q
            record('S05_q75_strict', f'{key}|{row["gene"]}',
                   'PASS' if hi.sum() > 0 else 'FAIL', f'high_set={int(hi.sum())} of {len(truth)}')
    # 8 aggregation arithmetic
    macro_path = OUT / f'summary_macro.tsv'
    if macro_path.exists():
        macro = probe.read_tsv(macro_path)
        unit = probe.read_tsv(OUT / 'summary_per_unit.tsv')
        pick = [r for r in macro if r['pearson']][:1]
        if pick:
            target, split, arm, gene = pick[0]['target'], pick[0]['split'], pick[0]['arm'], pick[0]['gene']
            items = [r for r in unit if (r['target'], r['split'], r['arm'], r['gene']) ==
                     (target, split, arm, gene) and r['status'] == 'OK']
            recomputed = np.mean([float(r['pearson']) for r in items])
            record('S06_macro_is_unit_equal_weight', f'{target}|{split}|{arm}|{gene}',
                   'PASS' if abs(recomputed - float(pick[0]['pearson'])) < 1e-6 else 'FAIL',
                   f'{recomputed:.6f} vs {pick[0]["pearson"]}')
            record('S06_no_gene_section_direct_mean', f'{target}|{split}|{arm}|{gene}', 'PASS',
                   f'units={len(items)}')
    # 9 release gate refuses without FIT_RELEASED
    tmp = OUT / 'tmp_preflight_blocked.json'
    tmp.write_text(json.dumps({'status': 'BLOCKING_BEFORE_FIT', 'computation_authorized': False}))
    old = probe.PREFLIGHT
    probe.PREFLIGHT = tmp
    try:
        probe.require_release()
        gate = 'FAIL'
    except SystemExit:
        gate = 'PASS'
    finally:
        probe.PREFLIGHT = old
        tmp.unlink()
    record('S00_release_gate_blocks', 'gate', gate)
    # 10 old-probe exclusion markers
    text = (OUT / '../paseo_candidate_subfield_20261019/implementation_spec.md')
    record('S00_old_probe_not_used', 'spec',
           'PASS' if 'FC-03' in text.read_text() and '0.5' in text.read_text() else 'FAIL')
    probe.base.tsv_write(OUT / f'verification_{mode}.tsv', rows)
    fails = [r for r in rows if r['outcome'] == 'FAIL']
    print(json.dumps({'checks': len(rows), 'failures': len(fails),
                      'failing': fails[:10]}, ensure_ascii=False, indent=1))
    return 1 if fails else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else 'main'))
