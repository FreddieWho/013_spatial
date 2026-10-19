#!/usr/bin/env python3
"""Bounded end-to-end smoke driver for the descriptive field review (test scaffolding).

It calls the real functions of scripts/field_biological_review_20261019.py
(prepare_section -> run_split -> summarize) with WORK redirected to
infra/field_biological_review_20261019/smoke/, so the previously untested summarize path
is exercised on real data without touching the production output directories and without
running the full 72-section / 1000-package analysis.

Splits are trimmed to the sampled sections; splits without a sampled section get an empty,
explicitly labelled table so that summarize's completeness checks stay honest.

Usage: python scripts/smoke_field_biological_review_20261019.py [--sections S1 S2 ...]
"""
import argparse
import json
import shutil
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import field_biological_review_20261019 as review
import field_discovery_20261018 as base

DEFAULT_SECTIONS = ['ZEN38', 'ZEN42', 'TENX171', 'NCBI642', 'TENX161']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sections', nargs='*', default=DEFAULT_SECTIONS)
    args = parser.parse_args()
    sids = list(args.sections)
    smoke = review.WORK / 'smoke'
    smoke.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(review.WORK / 'contract.json', smoke / 'contract.json')
    review.WORK = smoke
    start = time.time()
    measured, relations = [], []
    for sid in sids:
        result = review.prepare_section(dict(sid=sid))
        measured += result['measured']
        relations += result['relations']
        print(json.dumps(dict(prepared=sid, seconds=round(result['seconds'], 2))), flush=True)
    base.tsv_write(smoke / 'measured_fields.tsv', measured)
    base.tsv_write(smoke / 'measured_module_relations.tsv', relations)
    (smoke / 'PREPARED').write_text('SMOKE_ONLY_TRUNCATED_PREPARE\n')
    ctx = review.setup(verify='none')
    split_summary = []
    for split in ctx['splits']:
        test = [s for s in split['test'] if s in sids]
        trimmed = dict(split, test=test)
        if test:
            info = review.run_split(trimmed)
            split_summary.append(dict(split=split['name'], arm=split['arm'], sections=test,
                                      **{k: v for k, v in info.items() if k != 'split'}))
        else:
            out = smoke / 'splits' / split['name']
            out.mkdir(parents=True, exist_ok=True)
            base.tsv_write(out / 'prediction_metrics.tsv', [])
            base.tsv_write(out / 'prediction_module_relations.tsv', [])
            base.json_write(out / 'split_receipt.json',
                            dict(name=split['name'], arm=split['arm'], n_sections=0, rows=0,
                                 relation_rows=0, trimmed_empty_for_smoke=True))
            (out / 'COMPLETE').write_text('TRIMMED_EMPTY_SMOKE\n')
            split_summary.append(dict(split=split['name'], arm=split['arm'], sections=[],
                                      status='trimmed_empty'))
    review.summarize()
    outputs = {}
    for name in ['measured_fields.tsv', 'measured_module_relations.tsv', 'prediction_metrics.tsv',
                 'prediction_module_relations.tsv', 'unit_metrics.tsv', 'macro_summary.tsv',
                 'macro_strata.tsv', 'within_unit_consistency.tsv', 'diagnostic_sensitivity.tsv',
                 'measured_summary.tsv', 'cross_stage_consistency.json', 'REPORT.md',
                 'freeze_manifest.json', 'STATUS']:
        path = smoke / name
        outputs[name] = dict(exists=path.exists(),
                             bytes=path.stat().st_size if path.exists() else 0,
                             sha256=base.sha256(path) if path.exists() else None)
    summary = dict(sections=sids, split_summary=split_summary, outputs=outputs,
                   seconds=round(time.time() - start, 1),
                   code_sha256=base.sha256(Path(review.__file__)),
                   contract_sha256=base.sha256(review.WORK / 'contract.json'))
    base.json_write(smoke / 'smoke_summary.json', summary)
    print(json.dumps(dict(sections=sids, seconds=summary['seconds'],
                          outputs={k: v['bytes'] for k, v in outputs.items()}), indent=1), flush=True)


if __name__ == '__main__':
    main()
