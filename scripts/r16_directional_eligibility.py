"""D-143 local qualification and descriptive direction contrasts; no spatial p."""
from pathlib import Path
import csv
import hashlib
import json
import shutil

import h5py
import numpy as np
from scipy.spatial import cKDTree

from r04.io_contract import _decode_strings
from r16.recovery.corrected import feature_block, fit_linear, tls_components, NotTestable
from r16.recovery.repair_pipeline import axes

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'infra/repair_20260921'
OUT = ROOT / 'infra/directional_eligibility_20260922'


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


def write_table(name, rows):
    if not rows:
        raise ValueError(f'Empty table: {name}')
    with (OUT / name).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def sectors(xy, labels, component):
    """Selection depends only on measured geometry and explicit source labels."""
    center = xy[component].mean(0)
    tumor = xy[labels[:, 1] == 1]
    if not len(tumor):
        return np.nan, [], 'NO_TUM_LABEL'
    distance, nearest = cKDTree(tumor).query(center)
    if not 4 <= distance <= 24:
        return float(distance), [], 'DISTANCE_OUTSIDE_4_24'
    direction = tumor[nearest] - center
    delta = xy - center
    radius = np.linalg.norm(delta, axis=1)
    angle = np.arctan2(delta[:, 1], delta[:, 0]) - np.arctan2(direction[1], direction[0])
    edges = distance * np.array([.25, .5, .75, 1.])
    bins = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        radial = (radius >= lo) & (radius < hi) & (labels[:, 0] == 0)
        groups = []
        for turn in [0., np.pi / 2, -np.pi / 2]:
            offset = (angle - turn + np.pi) % (2 * np.pi) - np.pi
            groups.append(np.flatnonzero(radial & (np.abs(offset) <= np.pi / 6)))
        if all(len(g) >= 3 for g in groups):
            bins.append(groups)
    return float(distance), bins, 'ELIGIBLE' if len(bins) >= 2 else 'INSUFFICIENT_MATCHED_BANDS'


def contrast(values, bins):
    return float(np.mean([values[a].mean() - (values[b].mean() + values[c].mean()) / 2
                          for a, b, c in bins]))


def main():
    assert shutil.disk_usage(ROOT).free >= 1_200_000_000_000
    OUT.mkdir(exist_ok=False)
    hashes = {}

    def track(p):
        hashes[str(p.relative_to(ROOT))] = digest(p)
        return p

    sections = json.loads(track(SOURCE / 'sections.json').read_text())
    definitions = json.loads(track(SOURCE / 'six_programs.json').read_text())
    genes = json.loads(track(SOURCE / 'genes.json').read_text())
    audit_path = track(ROOT / 'infra/structure-registry/gt_source_audit.tsv')
    audit = list(csv.DictReader(audit_path.open(), delimiter='\t'))
    track(ROOT / 'infra/r04/marker_proxy_combined.json')
    coverage, annotation, geometry, effects = [], [], [], []
    registry = [dict(gt_source_id=r['gt_source_id'], structure=r['structure_id'],
                     physical_unit=r['physical_unit_id'], status=r['audit_status'],
                     replayed_in_this_run=r['physical_unit_id'] in
                     {s['section'] for s in sections if s['cohort'] != 'DISCOVERY'})
                for r in audit]
    for s in sections:
        key = s['key']
        source_genes = set(json.loads(track(SOURCE / 'cache' / f'{key}_source_genes.json').read_text()))
        for name, d in definitions.items():
            wanted = set(d['input_genes'] + d['readout_genes'])
            coverage.append(dict(section=s['section'], patient=s['patient'], cohort=s['cohort'],
                                 program=name, source_complete=not (wanted - source_genes),
                                 analysis_complete=not (wanted - (source_genes & set(genes))),
                                 missing_source=';'.join(sorted(wanted - source_genes))))
        if s['cohort'] == 'DISCOVERY':
            continue  # Cache NaNs are not evidence that the original sources lack annotation.
        ident = dict(cohort=s['cohort'], patient=s['patient'], section=s['section'])
        if s['cohort'] == 'ST-CRC':
            sources = [r for r in audit if r['physical_unit_id'] == s['section'] and r['audit_status'] == 'AUDITABLE_GT']
            if not sources:
                annotation.append(dict(**ident, category='NO_AUDITED_SOURCE', n=0, matched_barcodes=0, source=''))
            for r in sources:
                p = track(ROOT / r['path'])
                assert digest(p) == r['sha256']
                rows = list(csv.DictReader(p.open()))
                field = next(k for k in rows[0] if k != 'Barcode')
                lookup = {r['Barcode']: r[field].strip() for r in rows}
                assert len(lookup) == len(rows)
                categories = [lookup.get(b, 'MISSING_BARCODE') or 'EMPTY_UNKNOWN' for b in s['barcodes']]
                for cat in sorted(set(categories)):
                    annotation.append(dict(**ident, category=cat, n=categories.count(cat),
                                           matched_barcodes=sum(b in lookup for b in s['barcodes']), source=str(p.relative_to(ROOT))))
            continue
        xy = np.load(track(SOURCE / 'cache' / f'{key}_coords.npy'))
        labels = np.load(track(SOURCE / 'cache' / f'{key}_labels.npy'))
        alias = s['section'].split('::')[-1]
        p = track(ROOT / 'data/other_sources/zenodo_usz_tls_visium/TLS_VISIUM_USZ/h5ad_preprocessed' / f'{alias}.h5ad')
        with h5py.File(p) as f:
            bc = _decode_strings(f['obs']['_index'][()])
            cats = _decode_strings(f['obs']['ground_truth']['categories'][()])
            codes = f['obs']['ground_truth']['codes'][()]
        assert len(set(bc)) == len(bc)
        lookup = {b: cats[int(c)] if 0 <= c < len(cats) else 'UNASSIGNED' for b, c in zip(bc, codes)}
        raw = [lookup.get(b, 'MISSING_BARCODE') for b in s['barcodes']]
        replay = np.array([[np.nan, np.nan] if c in {'UNASSIGNED', 'MISSING_BARCODE'} else [float(c == 'TLS'), float(c == 'TUM')] for c in raw])
        assert np.array_equal(labels, replay, equal_nan=True)
        for cat in sorted(set(raw)):
            annotation.append(dict(**ident, category=cat, n=raw.count(cat), matched_barcodes=sum(b in lookup for b in s['barcodes']), source=str(p.relative_to(ROOT))))
        selected = []
        for i, component in enumerate(tls_components(xy, labels[:, 0])):
            distance, bins, status = sectors(xy, labels, component)
            geometry.append(dict(**ident, focus=i, n_tls_spots=len(component),
                                 tumor_distance=distance, matched_bands=len(bins), status=status,
                                 band_counts=json.dumps([[len(g) for g in b] for b in bins])))
            if status == 'ELIGIBLE':
                selected.append((i, bins))
        if not selected:
            continue
        path = track(SOURCE / 'cache' / f'{key}_counts.npy')
        assert digest(path) == s['counts_sha256']
        counts = np.load(path, mmap_mode='r')
        for name in ['P-stromal', 'P-mhc2']:
            d = definitions[name]
            assert set(d['input_genes'] + d['readout_genes']) <= source_genes & set(genes)
            for mode in ['raw', 'depth_normalized']:
                b = feature_block(counts, genes, d['input_genes'], d['readout_genes'], axes(), mode)
                X = b['X'][:, :-1]
                residual = b['y'] - X @ fit_linear(X, b['y'])
                sd = float(np.std(residual))
                if sd <= 1e-12:
                    raise NotTestable('CONSTANT_RESIDUAL')
                for i, bins in selected:
                    effects.append(dict(**ident, focus=i, program=name, mode=mode,
                                        residual_sd=sd, directional_delta_sd=contrast(residual / sd, bins)))
        print(key, 'foci', len(selected), 'eligible', flush=True)
    write_table('program_eligibility.tsv', coverage)
    write_table('source_annotations.tsv', annotation)
    write_table('registry_scope.tsv', registry)
    write_table('focus_eligibility.tsv', geometry)
    patients = []
    if effects:
        write_table('focus_effects.tsv', effects)
        # Explicit section then patient aggregation; repeated foci are not independent people.
        for name in ['P-stromal', 'P-mhc2']:
            for mode in ['raw', 'depth_normalized']:
                subset = [e for e in effects if e['program'] == name and e['mode'] == mode]
                for patient in sorted({e['patient'] for e in subset}):
                    rows = [e for e in subset if e['patient'] == patient]
                    means = [np.mean([e['directional_delta_sd'] for e in rows if e['section'] == section])
                             for section in sorted({e['section'] for e in rows})]
                    patients.append(dict(patient=patient, program=name, mode=mode, n_foci=len(rows),
                                         n_sections=len(means), directional_delta_sd=float(np.mean(means))))
        write_table('patient_effects.tsv', patients)
    summary = []
    for name in ['P-stromal', 'P-mhc2']:
        for mode in ['raw', 'depth_normalized']:
            values = [p['directional_delta_sd'] for p in patients if p['program'] == name and p['mode'] == mode]
            summary.append(dict(program=name, mode=mode, n_patients=len(values),
                                median_delta_sd=float(np.median(values)) if values else None,
                                positive_patients=sum(v > 0 for v in values), spatial_p='NOT_CALIBRATED'))
    result = dict(status='COMPLETED', scope='69 source-gene inventories; 22 external annotation replays; USZ directional description',
                  n_foci=len(geometry), n_eligible_foci=sum(r['status'] == 'ELIGIBLE' for r in geometry),
                  summary=summary, hidden_prediction='NOT_RUN', spatial_null='NOT_RUN',
                  general_scanner='NOT_RUN', go_family_extension='NOT_RUN',
                  input_hashes=hashes, code_sha256=digest(Path(__file__)),
                  output_hashes={p.name: digest(p) for p in OUT.glob('*.tsv')},
                  disk_free_bytes=shutil.disk_usage(ROOT).free)
    (OUT / 'receipt.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ['input_hashes', 'output_hashes']}, indent=2))


if __name__ == '__main__':
    main()
