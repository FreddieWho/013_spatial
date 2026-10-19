#!/usr/bin/env python3
"""Limited verification for the descriptive spatial-contrast results (D-275 scope).

Checks only what the task allows: aggregation weighting/denominators, barcode and row alignment,
and metric direction/sign recomputation. No re-fitting, no new models, no broad testing.

Usage: python verify_results.py [--work DIR] [--sample N]
"""
import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import field_raw_gene_spatial_20261019 as R          # noqa: E402
import field_discovery_20261018 as base              # noqa: E402
import field_biological_review_20261019 as review    # noqa: E402

MODELS = ['point_ridge', 'point_ridge_smooth_half', 'local_ridge', 'background_ridge',
          'permuted_ridge']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', default=str(ROOT / 'infra/field_raw_gene_spatial_20261019'))
    parser.add_argument('--sample', type=int, default=4)
    args = parser.parse_args()
    W = Path(args.work)
    results = []

    def record(name, expected, observed, passed):
        results.append(dict(check=name, expected=str(expected), observed=str(observed),
                            passed=bool(passed)))

    rows = list(csv.DictReader((W / 'prediction_metrics.tsv').open(), delimiter='\t'))
    paired = list(csv.DictReader((W / 'paired_metrics.tsv').open(), delimiter='\t'))
    units = list(csv.DictReader((W / 'unit_metrics.tsv').open(), delimiter='\t'))
    unit_paired = list(csv.DictReader((W / 'unit_paired.tsv').open(), delimiter='\t'))
    macro = list(csv.DictReader((W / 'macro_summary.tsv').open(), delimiter='\t'))
    claims = list(csv.DictReader((W / 'readout_claim_summary.tsv').open(), delimiter='\t'))
    denom = list(csv.DictReader((W / 'denominators.tsv').open(), delimiter='\t'))

    record('section_rows', 2500, len(rows), len(rows) == 2500)
    record('paired_rows', 2000, len(paired), len(paired) == 2000)
    record('unit_rows', 2240, len(units), len(units) == 2240)
    record('unit_paired_rows', 1792, len(unit_paired), len(unit_paired) == 1792)
    record('readout_claim_rows', 40, len(claims), len(claims) == 40)
    record('denominator_rows', 210, len(denom), len(denom) == 210)

    per = Counter()
    for r in units:
        per[(r['arm'], r['target'], r['model'])] += 1
    internal = {v for (a, t, m), v in per.items() if a == 'internal'}
    source = {v for (a, t, m), v in per.items() if a == 'source'}
    record('unit_group_counts', 'internal 65 / source 47', (sorted(internal), sorted(source)),
           internal == {65} and source == {47})

    ids = defaultdict(set)
    for r in rows:
        ids[r['arm']].add(r['section_id'])
    record('sections_per_arm', 'internal 72 / source 53',
           (len(ids['internal']), len(ids['source'])),
           len(ids['internal']) == 72 and len(ids['source']) == 53)

    confirmed = {r['section_id'] for r in rows if r['confirmed_unit'] == 'True'}
    unknown = {r['section_id'] for r in rows if r['confirmed_unit'] != 'True'}
    record('identity_split', 'confirmed 14 / unknown 58',
           (len(confirmed), len(unknown)), len(confirmed) == 14 and len(unknown) == 58)

    not_testable = {r['section_id'] for r in rows
                    if str(r.get('source_eligibility', '')).startswith('NOT_TESTABLE')}
    record('source_not_testable_sections_kept_internal', 19, len(not_testable),
           len(not_testable) == 19)
    extrapolation = {r['split'] for r in rows if str(r['cross_organ_extrapolation']) == 'True'}
    record('cross_organ_extrapolation_splits', "['source_0', 'source_5']",
           sorted(extrapolation), sorted(extrapolation) == ['source_0', 'source_5'])

    numeric_columns = ['rec_rmse', 'rec_centered_rmse', 'rec_pearson', 'rec_local_pattern_pearson',
                       'rec_jaccard_q75', 'rec_jaccard_q80', 'rec_jaccard_q90',
                       'rec_centroid_displacement_pitch', 'rec_patch75_largest_radius_abs_delta',
                       'rec_sd_ratio']
    bad = 0
    present = 0
    for r in rows:
        for c in numeric_columns:
            v = r.get(c, '')
            if v in ('', 'None'):
                continue
            present += 1
            if not np.isfinite(float(v)):
                bad += 1
    record('finite_numeric_values', 'no non-finite values', (present, bad), bad == 0)

    # barcode / row alignment against the frozen cache for a few stored packages
    ctx = R.setup(verify='none', sids=sorted(ids['internal'])[:1])
    mismatches = 0
    for r in rows[:args.sample * 5]:
        row = ctx['byid'][r['section_id']]
        cache = base.load_cache(row)
        ti = ctx['c']['targets'].index(r['target'])
        with np.load(W / 'splits' / r['split'] / f'{r["section_id"]}_t{ti}.npz',
                     allow_pickle=False) as z:
            if not np.array_equal(z['barcodes'], cache['barcodes']) or \
                    not np.array_equal(z['y'], cache['y'][:, ti]):
                mismatches += 1
    record('barcode_and_readout_alignment', '0 mismatches', mismatches, mismatches == 0)

    # metric direction recomputation on real rows (detail-component Pearson, strict-above Jaccard)
    geo_cache = {}
    recomputed = []
    for r in rows[:args.sample * 5]:
        row = ctx['byid'][r['section_id']]
        ti = ctx['c']['targets'].index(r['target'])
        cache = base.load_cache(row)
        geo = geo_cache.setdefault(r['section_id'], R.geo_for(ctx, row, cache))
        with np.load(W / 'splits' / r['split'] / f'{r["section_id"]}_t{ti}.npz',
                     allow_pickle=False) as z:
            y = z['y'].astype(np.float64)
            p = z[r['model']].astype(np.float64)
        dy = y - base.neighbor_mean(cache['near'], y)
        dp = p - base.neighbor_mean(cache['near'], p)
        manual_local = float(np.corrcoef(dy, dp)[0, 1])
        mask_y = y > np.quantile(y, 0.75)
        mask_p = p > np.quantile(p, 0.75)
        manual_jac = float(np.count_nonzero(mask_y & mask_p) / np.count_nonzero(mask_y | mask_p))
        recomputed.append((abs(manual_local - float(r['rec_local_pattern_pearson'])),
                           abs(manual_jac - float(r['rec_jaccard_q75']))))
    worst_local = max(x[0] for x in recomputed)
    worst_jac = max(x[1] for x in recomputed)
    record('local_detail_pearson_recompute', '<=1e-9', worst_local, worst_local <= 1e-9)
    record('jaccard_q75_strict_above_recompute', '<=1e-9', worst_jac, worst_jac <= 1e-9)

    # complete-unit rule: a two-section unit equals the mean of its sections, partials excluded
    two = [u for u in units if u['n_sections'] == '2']
    ok_mean = 0
    checked = 0
    for u in two[:20]:
        secs = [r for r in rows if r['arm'] == u['arm'] and r['unit_id'] == u['unit_id']
                and r['target'] == u['target'] and r['model'] == u['model']]
        if not secs or u['rec_rmse'] in ('', 'None'):
            continue
        checked += 1
        mean = float(np.mean([float(s['rec_rmse']) for s in secs]))
        if abs(mean - float(u['rec_rmse'])) < 1e-12:
            ok_mean += 1
    record('two_section_unit_equal_weighted_mean', 'all checked units match',
           (checked, ok_mean), checked > 0 and checked == ok_mean)

    partials = [u for u in units if u['partial_unit_metrics']]
    record('partial_units_excluded_from_complete_mean',
           'partial metric value is empty and listed', len(partials),
           all(u[u['partial_unit_metrics'].split(';')[0]] in ('', 'None') for u in partials))

    payload = dict(n_checks=len(results), n_failed=sum(1 for r in results if not r['passed']),
                   results=results)
    (W / 'verify_results.json').write_text(json.dumps(payload, indent=2) + '\n')
    print(json.dumps(payload, indent=2), flush=True)
    return 0 if payload['n_failed'] == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
