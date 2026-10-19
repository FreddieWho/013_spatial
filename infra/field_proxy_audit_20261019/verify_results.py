#!/usr/bin/env python3
"""Targeted verification for proxy-audit segment 1 (D-277 scope).

Only the checks the task allows: 72 sections x 4 targets coverage, 65 established units,
reuse of the approved support-gene table with its source stated, the frozen readout formula and
barcode order, the single-regression technical control (covariates exclude the target set and
residuals are orthogonal to both covariates), the 10 fixed draws, and the unit denominators.
No re-fitting, no new models, no broad testing.

Usage: python verify_results.py [--work DIR]
"""
import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import field_discovery_20261018 as base                     # noqa: E402
import field_proxy_audit_20261019 as A                      # noqa: E402

APPROVED = ROOT / 'infra/field_biological_review_20261019/measured_fields.tsv'
FROZEN_TARGETS = ['B_AXIS', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
                  'HALLMARK_INTERFERON_GAMMA_RESPONSE',
                  'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']


def rows(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def num(value):
    try:
        return float(value) if value not in (None, '', 'None', 'nan') else None
    except ValueError:
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', default=str(ROOT / 'infra/field_proxy_audit_20261019'))
    args = parser.parse_args()
    W = Path(args.work)
    results = []

    def record(name, expected, observed, passed):
        results.append(dict(check=name, expected=str(expected), observed=str(observed),
                            passed=bool(passed)))

    sections = rows(W / 'section_metrics.tsv')
    units = rows(W / 'unit_metrics.tsv')
    deletion = rows(W / 'deletion_sensitivity.tsv')
    subset = rows(W / 'subset_consistency.tsv')
    draws = rows(W / 'subset_draws.tsv')
    not_testable = rows(W / 'not_testable.tsv')
    ctx = A.setup(verify='none', sids=['MISC33', 'NCBI642'])

    record('section_target_rows', 288, len(sections), len(sections) == 288)
    record('sections_covered', 72, len({r['section_id'] for r in sections}),
           len({r['section_id'] for r in sections}) == 72)
    record('targets_covered', 4, len({r['target'] for r in sections}),
           len({r['target'] for r in sections}) == 4)
    per_target_units = Counter(r['target'] for r in units)
    record('unit_rows_per_target', '65 for each target', dict(per_target_units),
           set(per_target_units.values()) == {65})
    record('units_distinct', 65, len({r['unit_id'] for r in units}),
           len({r['unit_id'] for r in units}) == 65)
    states = Counter()
    for r in units:
        for key, value in r.items():
            if key.endswith('__state'):
                states[value] += 1
    record('unit_states_present', 'COMPLETE/INCOMPLETE/NOT_ESTIMABLE keys',
           dict(states), set(states) <= {'COMPLETE', 'INCOMPLETE_UNIT', 'NOT_ESTIMABLE_UNIT'})

    worst = max((num(r['readout_max_abs_diff']) or 0) for r in sections)
    record('frozen_readout_reproduced', '<=1e-6', worst, worst <= 1e-6)

    approved = {(r['section_id'], r['target']): r for r in rows(APPROVED)}
    joined = [r for r in sections if 'approved_genes_detected' in r]
    record('approved_support_table_joined', 288, len(joined), len(joined) == 288)
    mismatch = 0
    for r in joined:
        a = approved.get((r['section_id'], r['target']))
        if a is None:
            mismatch += 1
            continue
        if num(r['approved_genes_detected']) != num(a['sup_genes_detected']) \
                or num(r['approved_genes_in_panel']) != num(a['sup_genes_in_panel']):
            mismatch += 1
    record('approved_support_values_match', 0, mismatch, mismatch == 0)
    record('support_source_stated', 'one source path recorded',
           sorted({r.get('approved_support_source', '') for r in joined})[:1],
           all(r.get('approved_support_source') for r in joined))

    record('deletion_rows', 288, len(deletion), len(deletion) == 288)
    record('subset_rows', 288, len(subset), len(subset) == 288)
    per_case_draws = Counter((r['section_id'], r['target']) for r in draws)
    record('draws_per_case', '10 for every evaluated case',
           dict(Counter(per_case_draws.values())), set(per_case_draws.values()) == {10})
    halves_ok = all(int(r['n_half_a']) >= 3 and int(r['n_half_b']) >= 3 for r in draws)
    record('each_subset_at_least_three_genes', '>=3 per half', halves_ok, halves_ok)
    seeds = {}
    for r in draws:
        seeds.setdefault((r['section_id'], r['target']), set()).add(r['split_seed'])
    distinct = all(len(v) == 10 for v in seeds.values())
    record('ten_distinct_split_seeds_per_case', '10 distinct', distinct, distinct)
    evaluable = [num(r['n_evaluable_control_draws']) for r in subset]
    record('control_evaluable_recorded', '0..10 recorded', (min(evaluable), max(evaluable)),
           all(v is not None and 0 <= v <= 10 for v in evaluable))

    reasons = Counter(r['reason'] for r in not_testable)
    record('not_testable_reasons_recorded', 'structured reasons only', dict(reasons),
           set(reasons) <= {'NO_PANEL_GENES', 'FEWER_THAN_THREE_MEASURABLE',
                            'FEWER_THAN_SIX_MEASURABLE'})

    # frozen readout formula and barcode order, recomputed for two sample sections
    formula_ok = barcode_ok = True
    for sid in ['MISC33', 'NCBI642']:
        row = ctx['byid'][sid]
        cache = base.load_cache(row)
        union = sorted({g for t in ctx['c']['targets'] for g in ctx['panels'][t]})
        data = A.read_section(row, union)
        if not np.array_equal(data['barcodes'], cache['barcodes']):
            barcode_ok = False
        for target in ctx['c']['targets']:
            frozen = base.target_readout(data['matrix'], data['genes'], ctx['panels'][target])
            members = [g for g in ctx['panels'][target] if g in data['index']]
            proxy = data['z_union'][:, [data['index'][g] for g in members]].mean(axis=1)
            if float(np.max(np.abs(proxy - frozen))) > 1e-6:
                formula_ok = False
    record('barcode_order_matches_cache', 'equal', barcode_ok, barcode_ok)
    record('same_formula_for_every_map', '<=1e-6 vs frozen target_readout', formula_ok, formula_ok)

    # technical control: covariates exclude the target set, and residuals are orthogonal to both
    sid = 'MISC33'
    row = ctx['byid'][sid]
    cache = base.load_cache(row)
    union = sorted({g for t in ctx['c']['targets'] for g in ctx['panels'][t]})
    data = A.read_section(row, union)
    orth = exclusions = 0.0
    for target in ctx['c']['targets']:
        members = [g for g in ctx['panels'][target] if g in data['index']]
        cols = [data['index'][g] for g in members]
        z = data['z_union'][:, cols]
        keep = np.ones(len(data['genes']), bool)
        keep[[data['positions'][i] for i in cols]] = False
        non_target = data['matrix'][:, np.flatnonzero(keep)]
        library = np.asarray(non_target.sum(axis=1)).ravel()
        detection = np.asarray((non_target > 0).sum(axis=1)).ravel() / max(non_target.shape[1], 1)
        # the covariates must not contain the target counts themselves: compare with an
        # independent exclusion built from the union index
        target_mask = np.isin(data['genes'], members)
        rebuilt_library = np.asarray(data['matrix'][:, ~target_mask].sum(axis=1)).ravel()
        exclusions = max(exclusions, float(np.max(np.abs(rebuilt_library - library))))
        cov = np.c_[np.log1p(library), detection]
        for i in range(min(5, z.shape[1])):
            residual = A.regression_residual(z[:, i], cov)
            orth = max(orth, abs(A.corr(residual, cov[:, 0]) or 0.0),
                       abs(A.corr(residual, cov[:, 1]) or 0.0))
    record('residual_orthogonal_to_both_covariates', '<=1e-8', orth, orth <= 1e-8)
    record('technical_covariates_exclude_target_genes', '0.0', exclusions, exclusions == 0.0)

    payload = dict(n_checks=len(results), n_failed=sum(1 for r in results if not r['passed']),
                   results=results)
    (W / 'verify_results.json').write_text(json.dumps(payload, indent=2) + '\n')
    print(json.dumps(payload, indent=2), flush=True)
    return 0 if payload['n_failed'] == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
