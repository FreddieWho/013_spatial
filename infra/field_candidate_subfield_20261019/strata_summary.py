#!/usr/bin/env python3
"""Strata summary for the candidate-field probe (D-279): internal, source, organ, platform.

Reads summary_per_unit.tsv (complete units only) and cache/sections.tsv; writes
infra/field_candidate_subfield_20261019/summary_strata.tsv
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'infra/field_candidate_subfield_20261019'


def read_tsv(path):
    import csv
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main(mode='main'):
    units = read_tsv(OUT / f'summary_per_unit.tsv')
    sections = read_tsv(OUT / 'cache/sections.tsv')
    splits = {s.parent.name: json.loads(s.read_text()) for s in
              sorted((ROOT / 'infra/field_source_holdout_20261018').glob('*/split.json'))}
    unit_meta = defaultdict(dict)
    for row in sections:
        unit_meta[row['unit_id']] = dict(organ=row['organ'], platform=row['platform'],
                                         source=row['source_key'], disease=row['disease'],
                                         patient_known=row['patient_known'])
    keys = defaultdict(list)
    for row in units:
        if row['status'] != 'OK':
            continue
        split = splits[row['split']]
        meta = unit_meta.get(row['unit'])
        if not meta:
            continue
        keys[('mode', split['arm'])].append((row['target'], row['arm'], row['gene'], row))
        keys[('source', meta['source'])].append((row['target'], row['arm'], row['gene'], row))
        keys[('organ', meta['organ'])].append((row['target'], row['arm'], row['gene'], row))
        keys[('platform', meta['platform'])].append((row['target'], row['arm'], row['gene'], row))
        keys[('disease', meta['disease'])].append((row['target'], row['arm'], row['gene'], row))
        keys[('patient_known', meta['patient_known'])].append((row['target'], row['arm'], row['gene'], row))
        if split['arm'] == 'source':
            keys[('held_source', split.get('held_source'))].append(
                (row['target'], row['arm'], row['gene'], row))
    out = []
    for (kind, value), items in keys.items():
        grouped = defaultdict(list)
        for target, arm, gene, row in items:
            grouped[(target, arm, gene)].append(row)
        for (target, arm, gene), rows in grouped.items():
            pears = [float(r['pearson']) for r in rows if r['pearson']]
            jacs = [float(r['jaccard']) for r in rows if r['jaccard']]
            out.append(dict(strata_type=kind, strata_value=str(value)[:90], target=target, arm=arm,
                            gene=gene, n_units=len(rows),
                            pearson=(f'{np.mean(pears):.6f}' if pears else ''),
                            jaccard=(f'{np.mean(jacs):.6f}' if jacs else '')))
    out.sort(key=lambda r: (r['strata_type'], r['strata_value'], r['target'], r['arm'], r['gene']))
    path = OUT / f'summary_strata_{mode}.tsv'
    path.write_text('strata_type\tstrata_value\ttarget\tarm\tgene\tn_units\tpearson\tjaccard\n' +
                    '\n'.join('\t'.join(str(r[k]) for k in ('strata_type', 'strata_value', 'target',
                                                            'arm', 'gene', 'n_units', 'pearson',
                                                            'jaccard')) for r in out) + '\n')
    summary = defaultdict(set)
    for row in out:
        summary[row['strata_type']].add(row['strata_value'])
    print(json.dumps({'rows': len(out),
                      'strata': {k: len(v) for k, v in summary.items()}}, ensure_ascii=False))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else 'main'))
