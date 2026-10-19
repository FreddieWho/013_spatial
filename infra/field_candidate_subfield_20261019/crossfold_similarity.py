#!/usr/bin/env python3
"""Cross-fold descriptive similarity of training loadings (D-279).

Only training loadings are used; test maps, test B genes and held-out scores are never read.
A pair is NOT_ESTIMABLE when the two folds share fewer than 6 usable A genes.
Writes infra/field_candidate_subfield_20261019/crossfold_similarity.tsv
"""
import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'infra/field_candidate_subfield_20261019'
RESULTS = OUT / 'results/main'


def read_tsv(path):
    import csv
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def cosine(a, b):
    den = np.linalg.norm(a) * np.linalg.norm(b)
    return float(np.dot(a, b) / den) if den else float('nan')


def main():
    splits = {s.parent.name: json.loads(s.read_text()) for s in
              sorted((ROOT / 'infra/field_source_holdout_20261018').glob('*/split.json'))}
    rows = []
    for target in ['B_AXIS', 'HALLMARK_INTERFERON_GAMMA_RESPONSE',
                   'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
                   'HALLMARK_OXIDATIVE_PHOSPHORYLATION']:
        store = {}
        for name in sorted(splits):
            path = RESULTS / target / name / 'loadings.npz'
            if not path.exists():
                continue
            data = np.load(path, allow_pickle=True)
            store[name] = dict(h1=data['h1'], h2=data['h2'], usable=list(data['usable']))
        for a, b in itertools.combinations(sorted(store), 2):
            shared = [g for g in store[a]['usable'] if g in set(store[b]['usable'])]
            if len(shared) < 6:
                rows.append(dict(target=target, split_a=a, split_b=b, n_shared=len(shared),
                                 cos_comp1='', cos_comp2='', cos_best_pair='',
                                 status='NOT_ESTIMABLE_SHARED_A_LT_6'))
                continue
            ia = {g: i for i, g in enumerate(store[a]['usable'])}
            ib = {g: i for i, g in enumerate(store[b]['usable'])}
            ha = store[a]['h2'][:, [ia[g] for g in shared]]
            hb = store[b]['h2'][:, [ib[g] for g in shared]]
            c11 = cosine(ha[0], hb[0]); c22 = cosine(ha[1], hb[1])
            c12 = cosine(ha[0], hb[1]); c21 = cosine(ha[1], hb[0])
            best = max(0.5 * (c11 + c22), 0.5 * (c12 + c21))
            rows.append(dict(target=target, split_a=a, split_b=b, n_shared=len(shared),
                             cos_comp1=f'{c11:.6f}', cos_comp2=f'{c22:.6f}',
                             cos_best_pair=f'{best:.6f}', status='OK'))
    OUT.joinpath('crossfold_similarity.tsv').write_text(
        'target\tsplit_a\tsplit_b\tn_shared\tcos_comp1\tcos_comp2\tcos_best_pair\tstatus\n' +
        '\n'.join('\t'.join(str(r[k]) for k in ('target', 'split_a', 'split_b', 'n_shared',
                                                'cos_comp1', 'cos_comp2', 'cos_best_pair', 'status'))
                  for r in rows) + '\n')
    ok = [r for r in rows if r['status'] == 'OK']
    print(f'pairs {len(rows)} | estimable {len(ok)} | not estimable {len(rows) - len(ok)}')
    for target in sorted({r['target'] for r in ok}):
        vals = [float(r['cos_best_pair']) for r in ok if r['target'] == target]
        print(f'  {target[:42]:42s} median best-pair cosine {np.median(vals):+.3f}  n={len(vals)}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
