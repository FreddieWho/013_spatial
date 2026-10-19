#!/usr/bin/env python3
"""Read-only feasibility check for the two-component candidate-field probe (D-278).

Reads expression and the frozen manifests only; fits nothing, uses no held-out score.
It answers the questions that can block the implementation:
  * which sections the frozen fold file contains beyond the 72 qualified ones
  * the fixed A/B gene split (name + seed 0) per target and its overlaps
  * per-gene measurability (>= 20 tissue spots with a positive count, variance > 0)
    across the 72 sections, per target, for A and for B genes
  * how many A genes are usable per fold under candidate training-only rules and how many
    B genes are evaluable per test section
"""
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import field_discovery_20261018 as base                      # noqa: E402

MANIFEST = ROOT / 'infra/field_source_holdout_20261018/input_manifest.json'
FOLDS = ROOT / 'infra/tls_hooknet_rna_triage_20261014/split_folds_v1.tsv'
PANEL = ROOT / 'infra/tls_multiscale_gobp_20261002/panel_controls.json'
OUT = ROOT / 'infra/paseo_candidate_subfield_20261019'
TARGETS = ['B_AXIS', 'HALLMARK_INTERFERON_GAMMA_RESPONSE',
           'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
           'HALLMARK_OXIDATIVE_PHOSPHORYLATION']


def split_genes(genes):
    genes = sorted(set(genes))
    order = [genes[i] for i in np.random.default_rng(0).permutation(len(genes))]
    n_a = int(np.ceil(len(genes) * 2 / 3))
    return sorted(order[:n_a]), sorted(order[n_a:])


def main():
    rows = json.loads(MANIFEST.read_text())
    byid = {r['section_id']: r for r in rows}
    folds = {r['section_id']: r for r in csv.DictReader(FOLDS.open(), delimiter='\t')}
    panels = json.loads(PANEL.read_text())
    report = dict(
        fold_file_rows=len(folds), manifest_rows=len(rows),
        fold_file_extra_sections=sorted(set(folds) - set(byid)),
        manifest_missing_from_fold_file=sorted(set(byid) - set(folds)),
        folds=dict(Counter(folds[s]['fold'] for s in byid)),
        units=len({r['unit_id'] for r in rows}),
        source_eligible=sum(1 for r in rows if r['source_key'].startswith('http')),
        targets={})
    union = sorted({g for t in TARGETS for g in panels[t]})
    measurable = {}                      # (sid, gene) -> int nonzero spots
    for i, sid in enumerate(sorted(byid)):
        row = byid[sid]
        with h5py.File(base.HEST / f'st/{sid}.h5ad') as handle:
            genes = base.strings(handle['var/_index'][:])
            barcodes = base.strings(handle['obs/_index'][:])
            keep = np.flatnonzero(handle['obs/in_tissue'][:] == 1)
            matrix = base.read_matrix(handle['X'], len(barcodes), len(genes))[keep]
        lookup = {}
        for gene in union:
            hits = np.flatnonzero(genes == gene)
            if len(hits) == 1:
                lookup[gene] = int(hits[0])
        positions = [lookup[g] for g in union if g in lookup]
        present = [g for g in union if g in lookup]
        counts = matrix[:, positions].tocsr()
        n_any = np.asarray((counts > 0).sum(axis=0)).ravel()
        z = np.asarray(base.normalized_subset(matrix, positions,
                                              np.asarray(matrix.sum(axis=1)).ravel()).todense())
        sd = np.std(np.asarray(z, np.float64), axis=0)
        for j, gene in enumerate(present):
            measurable[(sid, gene)] = int(n_any[j]) if (n_any[j] >= 20 and sd[j] > 0) else int(min(n_any[j], 19))
        print(json.dumps(dict(section=sid, done=i + 1, of=len(byid))), flush=True)
    gene_rows = []
    for target in TARGETS:
        A, B = split_genes(panels[target])
        info = dict(panel=len(set(panels[target])), n_a=len(A), n_b=len(B),
                    overlap=len(set(A) & set(B)))
        for label, genes in [('A', A), ('B', B)]:
            per_gene = {}
            for gene in genes:
                frac = np.mean([measurable.get((sid, gene), 0) >= 20 for sid in byid])
                per_gene[gene] = float(frac)
                gene_rows.append(dict(target=target, side=label, gene=gene,
                                      measurable_section_fraction=float(frac)))
            fracs = np.array(list(per_gene.values()))
            info[f'{label}_measurable_fraction'] = dict(
                n=len(fracs), median=float(np.median(fracs)),
                q10=float(np.quantile(fracs, 0.10)), q90=float(np.quantile(fracs, 0.90)),
                all_sections=int(np.sum(fracs == 1.0)),
                at_least_80pct=int(np.sum(fracs >= 0.8)),
                at_least_half=int(np.sum(fracs >= 0.5)),
                never=int(np.sum(fracs == 0.0)))
        # per-section counts of measurable A and B genes
        a_counts = [sum(1 for g in A if measurable.get((sid, g), 0) >= 20) for sid in byid]
        b_counts = [sum(1 for g in B if measurable.get((sid, g), 0) >= 20) for sid in byid]
        info['per_section_measurable'] = dict(
            A_min=int(min(a_counts)), A_median=float(np.median(a_counts)), A_max=int(max(a_counts)),
            B_min=int(min(b_counts)), B_median=float(np.median(b_counts)), B_max=int(max(b_counts)))
        # per-fold usable A under training-only rules
        info['per_fold_usable_A'] = {}
        for fold in sorted({folds[s]['fold'] for s in byid}):
            train = [s for s in byid if folds[s]['fold'] != fold]
            for rule in [1.0, 0.8, 0.5]:
                usable = sum(1 for g in A
                             if np.mean([measurable.get((s, g), 0) >= 20 for s in train]) >= rule)
                info['per_fold_usable_A'][f'fold{fold}_rule{int(rule*100)}'] = int(usable)
        # B genes evaluable per section (for the verification side)
        info['B_evaluable_section_counts'] = dict(Counter(
            [sum(1 for g in B if measurable.get((sid, g), 0) >= 20) for sid in byid]))
        report['targets'][target] = info
    (OUT / 'feasibility.json').write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    base.tsv_write(OUT / 'feasibility_gene_measurability.tsv', gene_rows)
    print(json.dumps(report, indent=2, ensure_ascii=False)[:4000], flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
