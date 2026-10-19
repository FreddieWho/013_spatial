#!/usr/bin/env python3
"""Per-case frozen-baseline agreement audit for the accepted full run (read-only).

The frozen implementation (Git 0e6e5f5) is not modified. This audit re-reads its stored
predictions after `fit` and compares, for every (split, target, section) case, the re-fitted
point / fixed-smoothing predictions with the frozen source-holdout predictions using the
tolerances fixed before the run in the machine contract. Any violation or missing case makes
the audit exit non-zero so the caller stops before evaluate/summarize.

Usage: python audit_frozen_agreement.py [--work DIR] [--pred DIR]
"""
import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_WORK = ROOT / 'infra/field_raw_gene_spatial_20261019'
DEFAULT_PRED = ROOT / 'infra/field_source_holdout_20261018'

MODEL_PAIRS = [('point_ridge', 'gene_ridge',
                'point_prediction_max_abs_diff_vs_frozen',
                'point_section_rmse_abs_diff_vs_frozen'),
               ('point_ridge_smooth_half', 'gene_ridge_smooth_half',
                'smooth_prediction_max_abs_diff_vs_frozen',
                'smooth_section_rmse_abs_diff_vs_frozen')]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work', default=str(DEFAULT_WORK))
    parser.add_argument('--pred', default=str(DEFAULT_PRED))
    args = parser.parse_args()
    work, pred = Path(args.work), Path(args.pred)
    contract = json.loads((work / 'contract.json').read_text())
    tol = contract['tolerances']
    splits = json.loads((work / 'splits.json').read_text())
    rows, failures = [], []
    cases = Counter()
    sections_seen = 0
    for split in splits:
        for ti, target in enumerate(contract['targets']):
            cases[(split['name'], target)] += 1
            diag_path = work / 'splits' / split['name'] / 'fit_diagnostics.tsv'
            if not diag_path.exists():
                failures.append(f'missing fit diagnostics: {split["name"]}/{target}')
                continue
            diag = [r for r in csv.DictReader(diag_path.open(), delimiter='\t')
                    if r['target'] == target]
            if len(diag) != len(contract['fit_models']):
                failures.append(f'incomplete fit diagnostics: {split["name"]}/{target} '
                                f'({len(diag)} rows)')
            for row in diag:
                if row['converged'] != 'True':
                    failures.append(f'non-converged fit: {split["name"]}/{target}/{row["model"]}')
            max_abs = {name: 0.0 for name, *_ in MODEL_PAIRS}
            max_rmse = {name: 0.0 for name, *_ in MODEL_PAIRS}
            n_sections = 0
            for sid in split['test']:
                mine_path = work / 'splits' / split['name'] / f'{sid}_t{ti}.npz'
                frozen_path = pred / split['name'] / f'{sid}_t{ti}_full15.npz'
                if not mine_path.exists() or not frozen_path.exists():
                    failures.append(f'missing package: {split["name"]}/{target}/{sid}')
                    continue
                with np.load(mine_path, allow_pickle=False) as mine, \
                        np.load(frozen_path, allow_pickle=False) as frozen:
                    y = mine['y'].astype(np.float64)
                    if not np.array_equal(mine['y'], frozen['y']):
                        failures.append(f'measured y mismatch: {split["name"]}/{target}/{sid}')
                    for name, frozen_name, tol_max, tol_rmse in MODEL_PAIRS:
                        a = mine[name].astype(np.float64)
                        b = frozen[frozen_name].astype(np.float64)
                        d = float(np.max(np.abs(a - b)))
                        r = abs(float(np.sqrt(np.mean((y - a) ** 2)))
                                - float(np.sqrt(np.mean((y - b) ** 2))))
                        max_abs[name] = max(max_abs[name], d)
                        max_rmse[name] = max(max_rmse[name], r)
                        if d > tol[tol_max]:
                            failures.append(f'{name} max|diff| {d:.3e} > {tol[tol_max]:.1e} '
                                            f'at {split["name"]}/{target}/{sid}')
                        if r > tol[tol_rmse]:
                            failures.append(f'{name} section-RMSE diff {r:.3e} > {tol[tol_rmse]:.1e} '
                                            f'at {split["name"]}/{target}/{sid}')
                n_sections += 1
                sections_seen += 1
            row = dict(split=split['name'], arm=split['arm'], target=target,
                       n_sections=n_sections,
                       point_max_abs_diff=max_abs['point_ridge'],
                       point_section_rmse_diff=max_rmse['point_ridge'],
                       smooth_max_abs_diff=max_abs['point_ridge_smooth_half'],
                       smooth_section_rmse_diff=max_rmse['point_ridge_smooth_half'],
                       point_max_abs_limit=tol['point_prediction_max_abs_diff_vs_frozen'],
                       point_section_rmse_limit=tol['point_section_rmse_abs_diff_vs_frozen'])
            rows.append(row)
    expected_cases = len(splits) * len(contract['targets'])
    if len(cases) != expected_cases:
        failures.append(f'case coverage {len(cases)} != {expected_cases}')
    summary = dict(
        cases_checked=len(cases), cases_expected=expected_cases,
        section_target_packages=sections_seen,
        worst_point_max_abs_diff=max([r['point_max_abs_diff'] for r in rows], default=None),
        worst_point_section_rmse_diff=max([r['point_section_rmse_diff'] for r in rows], default=None),
        worst_smooth_max_abs_diff=max([r['smooth_max_abs_diff'] for r in rows], default=None),
        worst_smooth_section_rmse_diff=max([r['smooth_section_rmse_diff'] for r in rows], default=None),
        limits=dict(point_max_abs=tol['point_prediction_max_abs_diff_vs_frozen'],
                    section_rmse_abs=tol['point_section_rmse_abs_diff_vs_frozen']),
        n_failures=len(failures), failures=failures[:40], status='PASS' if not failures else 'FAIL')
    with (work / 'audit_frozen_agreement.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    (work / 'audit_frozen_agreement.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2), flush=True)
    return 0 if not failures else 1


if __name__ == '__main__':
    sys.exit(main())
