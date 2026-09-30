"""Aggregate complete stage outputs at patient level; no silent partial-stage completion."""
import csv
from itertools import product
import json
from pathlib import Path
import warnings
from collections import defaultdict

import numpy as np

from r16.recovery.corrected import fit_from_moments
from scripts.r16_gobp_stage import SOURCE, OUT, digest, atomic_json


def table(name, rows):
    rows = iter(rows)
    first = next(rows)
    with (OUT / name).open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(first), delimiter='\t')
        w.writeheader()
        w.writerow(first)
        w.writerows(rows)


def mean(a, axis):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmean(a, axis=axis)


def median(a, axis):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        return np.nanmedian(a, axis=axis)


def improvement(base, value, ok):
    return np.divide(100 * (base - value), base, out=np.full_like(base, np.nan), where=ok & (base > 1e-12))


def exact_p(values):
    """All patient rows must be present; do not turn missing patients into a smaller test."""
    signs = np.array(list(product([-1., 1.], repeat=values.shape[0])))
    good = np.isfinite(values).all(0)
    result = np.full(values.shape[1], np.nan)
    v = values[:, good]
    result[good] = np.mean(abs(signs @ v) >= abs(v.sum(0))[None, :] - 1e-12, axis=0)
    return result


def holm(rows, field):
    p = np.array([r[field] if np.isfinite(r[field]) else 1. for r in rows])
    order = np.argsort(p)
    adjusted = np.minimum(1., np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    for i, q in zip(order, adjusted):
        rows[i]['holm_family_p'] = float(q)


def common_hits(summary, metric, threshold, mode):
    sets = []
    for c, n in [('ST-CRC', 7), ('USZ', 8)]:
        sets.append({r['program'] for r in summary if r['mode'] == mode and r['cohort'] == c
                     and r['n_patients'] == n and r[metric] > threshold})
    return sets[0] & sets[1]


def main():
    compute = json.loads((OUT / 'compute_receipt.json').read_text())
    assert compute['status'] == 'COMPLETED' and compute['tasks'] == 44 and compute['programs'] == 6870
    sections = [s for s in json.loads((SOURCE / 'sections.json').read_text()) if s['cohort'] != 'DISCOVERY']
    ids = sorted(json.loads((SOURCE / 'go_representatives.json').read_text()))
    molecular_rows, molecular_summary, mask_rows, mask_summary = [], [], [], []
    direction_rows, direction_summary, distance_rows, distance_summary = [], [], [], []
    receipts = []
    for mode in ['raw', 'depth_normalized']:
        all_data = {}
        for s in sections:
            tag = f"{s['key']}_{mode}"
            rec = json.loads((OUT / 'sections' / f'{tag}.json').read_text())
            path = OUT / 'sections' / f'{tag}.npz'
            assert rec['status'] == 'COMPLETED' and rec['n_programs'] == 6870 and digest(path) == rec['output_sha256']
            with np.load(path) as f:
                assert f['ids'].tolist() == ids
                all_data[s['key']] = {k: f[k] for k in f.files if k != 'ids'}
            receipts.append(rec)
        for cohort in ['ST-CRC', 'USZ']:
            chosen = [s for s in sections if s['cohort'] == cohort]
            patients = sorted({s['patient'] for s in chosen})
            patient_data = []
            for p in patients:
                group = [all_data[s['key']] for s in chosen if s['patient'] == p]
                item = {k: np.mean([g[k] for g in group], axis=0) for k in ['gram', 'rhs', 'yy']}
                item['mask'] = np.mean([mean(g['masked_mse'], axis=1) for g in group], axis=0)
                item['mask_variance'] = np.mean([mean(g['masked_target_variance'], axis=1) for g in group], axis=0)
                item['direction'] = np.mean([mean(g['directional_delta'], axis=1) for g in group if g['directional_delta'].shape[1]], axis=0) if any(g['directional_delta'].shape[1] for g in group) else np.full(len(ids), np.nan)
                item['distance'] = np.mean([mean(g['distance_cv_mse'], axis=1) for g in group], axis=0)
                patient_data.append(item)
            gram = np.array([d['gram'] for d in patient_data])
            rhs = np.array([d['rhs'] for d in patient_data])
            yy = np.array([d['yy'] for d in patient_data])
            mse = np.full((len(patients), len(ids), 3), np.nan)
            for i, p in enumerate(patients):
                others = np.arange(len(patients)) != i
                tg, tr = gram[others].mean(0), rhs[others].mean(0)
                valid = np.isfinite(tg).all((1, 2)) & np.isfinite(gram[i]).all((1, 2))
                for j, n in enumerate([9, 10, 11]):
                    beta = fit_from_moments(tg[valid, :n, :n], tr[valid, :n])
                    errors = yy[i, valid] - 2 * np.sum(beta * rhs[i, valid, :n], axis=1) + np.einsum('bi,bij,bj->b', beta, gram[i, valid, :n, :n], beta)
                    mse[i, valid, j] = np.maximum(errors, 0)
            var = yy - rhs[:, :, 0] ** 2
            mol_deltas = [improvement(mse[:, :, 0], mse[:, :, j], var > 1e-12) for j in [1, 2]]
            for i, p in enumerate(patients):
                for j, g in enumerate(ids):
                    molecular_rows.append(dict(mode=mode, cohort=cohort, patient=p, program=g,
                                               status='COMPUTED' if np.isfinite(mol_deltas[0][i, j]) else 'NOT_TESTABLE',
                                               mse_M0=mse[i, j, 0], mse_M1=mse[i, j, 1], mse_M2=mse[i, j, 2],
                                               delta_M1=mol_deltas[0][i, j], delta_M2=mol_deltas[1][i, j]))
            mm = [median(d, 0) for d in mol_deltas]
            for j, g in enumerate(ids):
                molecular_summary.append(dict(mode=mode, cohort=cohort, program=g, n_patients=int(np.isfinite(mol_deltas[0][:, j]).sum()),
                                              median_delta_M1=mm[0][j], median_delta_M2=mm[1][j]))
            masked = np.array([d['mask'] for d in patient_data])
            mv = np.array([d['mask_variance'] for d in patient_data])
            for i, p in enumerate(patients):
                for j, g in enumerate(ids):
                    mask_rows.append(dict(mode=mode, cohort=cohort, patient=p, program=g,
                                          status='COMPUTED' if np.isfinite(mv[i,j]) and mv[i,j]>1e-12 else 'NOT_TESTABLE_CONSTANT_OR_MISSING_TARGET',
                                          target_var=mv[i, j], **{f'mse_{m}':masked[i, j, k] for k, m in enumerate(['M0','M1','M2','NN','KNN8'])}))
            for mi, model in [(1, 'M1'), (2, 'M2')]:
                for bi, baseline in [(0, 'M0'), (3, 'NN'), (4, 'KNN8')]:
                    delta = improvement(masked[:, :, bi], masked[:, :, mi], mv > 1e-12)
                    pvalues = exact_p(delta)
                    med = median(delta, 0)
                    for j, g in enumerate(ids):
                        mask_summary.append(dict(mode=mode, cohort=cohort, program=g, model=model, baseline=baseline,
                                                 n_patients=int(np.isfinite(delta[:, j]).sum()), median_improvement=med[j],
                                                 positive_patients=int((delta[:, j] > 0).sum()), patient_sign_flip_p=pvalues[j]))
            if cohort == 'USZ':
                direct = np.array([d['direction'] for d in patient_data])
                eligible = np.isfinite(direct).any(1)
                pv = exact_p(direct[eligible])
                med = median(direct, 0)
                for i, p in enumerate(patients):
                    if not eligible[i]:
                        continue
                    for j, g in enumerate(ids):
                        direction_rows.append(dict(mode=mode, patient=p, program=g,
                                                   status='COMPUTED' if np.isfinite(direct[i,j]) else 'NOT_TESTABLE',delta_sd=direct[i, j]))
                for j, g in enumerate(ids):
                    direction_summary.append(dict(mode=mode, program=g, n_patients=int(np.isfinite(direct[:, j]).sum()),
                                                  median_delta_sd=med[j], positive_patients=int((direct[:, j] > 0).sum()),
                                                  negative_patients=int((direct[:, j] < 0).sum()), patient_sign_flip_p=pv[j]))
                dist = np.array([d['distance'] for d in patient_data])
                delta = improvement(dist[:, :, 0], dist[:, :, 1], np.isfinite(dist).all(2))
                pv = exact_p(delta)
                med = median(delta, 0)
                for i, p in enumerate(patients):
                    for j, g in enumerate(ids):
                        distance_rows.append(dict(mode=mode, patient=p, program=g,
                                                  status='COMPUTED' if np.isfinite(delta[i,j]) else 'NOT_TESTABLE_CONSTANT_OR_MISSING_TARGET',
                                                  mse_QC=dist[i, j, 0], mse_QC_distance=dist[i, j, 1], improvement=delta[i, j]))
                for j, g in enumerate(ids):
                    distance_summary.append(dict(mode=mode, program=g, n_patients=int(np.isfinite(delta[:, j]).sum()),
                                                 median_improvement=med[j], positive_patients=int((delta[:, j] > 0).sum()),
                                                 patient_sign_flip_p=pv[j]))
        print('aggregated', mode, flush=True)
    for rows in [mask_summary, direction_summary, distance_summary]:
        holm(rows, 'patient_sign_flip_p')
    for name, rows in [('molecular_patients', molecular_rows), ('molecular_summary', molecular_summary),
                       ('masked_patients', mask_rows), ('masked_summary', mask_summary),
                       ('direction_patients', direction_rows), ('direction_summary', direction_summary),
                       ('distance_patients', distance_rows), ('distance_summary', distance_summary)]:
        table(name + '.tsv', rows)
    mol = {mode:common_hits(molecular_summary, 'median_delta_M1', 20, mode) for mode in ['raw','depth_normalized']}
    masks = {}
    for mode in ['raw','depth_normalized']:
        masks[mode] = {}
        for model in ['M1','M2']:
            rows = [r for r in mask_summary if r['model'] == model and r['baseline'] == 'KNN8']
            hits = common_hits(rows, 'median_improvement', 0, mode)
            masks[mode][model] = hits
    stability = []
    for g in ids:
        stability.append(dict(program=g, molecular_raw_20=g in mol['raw'], molecular_normalized_20=g in mol['depth_normalized'],
                              masked_M1_raw=g in masks['raw']['M1'], masked_M1_normalized=g in masks['depth_normalized']['M1'],
                              masked_M2_raw=g in masks['raw']['M2'], masked_M2_normalized=g in masks['depth_normalized']['M2']))
    table('cross_mode_candidates.tsv', stability)
    distance_counts = {mode:dict(positive_median=int(sum(r['mode']==mode and r['n_patients']==8 and r['median_improvement']>0 for r in distance_summary)),
                                positive_all8=int(sum(r['mode']==mode and r['positive_patients']==8 for r in distance_summary))) for mode in ['raw','depth_normalized']}
    direction_counts = {mode:dict(positive_all3=sum(r['mode']==mode and r['positive_patients']==3 for r in direction_summary),
                                  negative_all3=sum(r['mode']==mode and r['negative_patients']==3 for r in direction_summary)) for mode in ['raw','depth_normalized']}
    grouped = defaultdict(list)
    for r in mask_summary:
        grouped[(r['program'],r['model'])].append(r)
    context = []
    directions = {(r['program'],r['mode']):r for r in direction_summary}
    distances = {(r['program'],r['mode']):r for r in distance_summary}
    for g in ids:
        row = dict(program=g, molecular_robust_gt20=g in mol['raw'] & mol['depth_normalized'])
        for model in ['M1','M2']:
            rows = grouped[(g,model)]
            complete = len(rows)==12 and all(r['n_patients']==(7 if r['cohort']=='ST-CRC' else 8) for r in rows)
            row[model+'_all_conditions_complete'] = complete
            row[model+'_minimum_median_improvement'] = min(r['median_improvement'] for r in rows) if complete else np.nan
            row[model+'_all_patients_all_conditions_positive'] = complete and all(r['positive_patients']==r['n_patients'] for r in rows)
        dr = [directions[(g,m)] for m in ['raw','depth_normalized']]
        row['direction_all3_positive_both_modes'] = all(r['positive_patients']==3 for r in dr)
        row['direction_all3_negative_both_modes'] = all(r['negative_patients']==3 for r in dr)
        di = [distances[(g,m)] for m in ['raw','depth_normalized']]
        row['distance_positive_median_both_modes'] = all(r['n_patients']==8 and r['median_improvement']>0 for r in di)
        context.append(row)
    table('candidate_context.tsv', context)
    all_baselines = {model:dict(
        strictly_positive_minimum_median=sum(r[model+'_minimum_median_improvement']>0 for r in context),
        minimum_median_above_one_percent=sum(r[model+'_minimum_median_improvement']>1 for r in context),
        best_minimum_median_percent=float(np.nanmax([r[model+'_minimum_median_improvement'] for r in context])),
        all_patients_all_conditions_positive=sum(r[model+'_all_patients_all_conditions_positive'] for r in context)) for model in ['M1','M2']}
    result = dict(status='COMPLETED', programs=6870, section_mode_tasks=44,
                  molecular_both_cohorts_gt20={k:len(v) for k,v in mol.items()},
                  molecular_stable_raw_normalized=len(mol['raw'] & mol['depth_normalized']),
                  masked_both_cohorts_positive_vs_KNN8={mode:{model:len(v) for model,v in value.items()} for mode,value in masks.items()},
                  masked_stable_raw_normalized={model:len(masks['raw'][model] & masks['depth_normalized'][model]) for model in ['M1','M2']},
                  distance_counts=distance_counts, direction_counts=direction_counts,
                  masked_both_modes_cohorts_all_three_baselines=all_baselines,
                  molecular45_masked92_overlap=sum(r['molecular_robust_gt20'] and r['M2_minimum_median_improvement']>0 for r in context),
                  stable_distance_positive_median=sum(r['distance_positive_median_both_modes'] for r in context),
                  stable_direction_positive=sum(r['direction_all3_positive_both_modes'] for r in context),
                  stable_direction_negative=sum(r['direction_all3_negative_both_modes'] for r in context),
                  holm_rejections=dict(mask=sum(r['holm_family_p']<=.05 for r in mask_summary),
                                       direction=sum(r['holm_family_p']<=.05 for r in direction_summary),
                                       distance=sum(r['holm_family_p']<=.05 for r in distance_summary)),
                  patient_test_scope='Conditional on independent symmetric patient effects; whole endpoint-family Holm; not pointwise spatial p',
                  actual_spatial_association_p='NOT_CALIBRATED', hidden_structure_AUC='NOT_TESTABLE_NO_POSITIVE_WINDOWS',
                  row_counts={name:len(rows) for name,rows in [('molecular_patients',molecular_rows),('mask_patients',mask_rows),('direction_patients',direction_rows),('distance_patients',distance_rows)]},
                  numerical_output_hashes={p.name:digest(p) for p in OUT.glob('*.tsv')})
    # NumPy comparisons can produce NumPy integer counts in nested result dictionaries.
    result = json.loads(json.dumps(result,default=lambda v:v.item()))
    atomic_json(OUT / 'results.json', result)
    print(json.dumps({k:v for k,v in result.items() if k!='numerical_output_hashes'},indent=2))


if __name__ == '__main__':
    main()
