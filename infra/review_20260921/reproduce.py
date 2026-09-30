"""Read-only implementation audit probes; writes only to this audit directory.

Run from repository root:
LD_LIBRARY_PATH=/opt/anaconda3/lib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MPLCONFIGDIR=/tmp/spatial-review-mpl PYTHONDONTWRITEBYTECODE=1 \
  PYTHONPATH=. python infra/review_20260921/reproduce.py
These probes are diagnostics, not replacement scientific runs.
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
from scipy.interpolate import BSpline
from scipy.linalg import null_space
from scipy.spatial import cKDTree
from scipy.stats import binomtest, spearmanr

from r04.structure_readout import evaluate_shared_specific_readout
from r04.loaders import read_10x_positions
from r16 import axes as A
from r16.section_io import usz_tls_labels
from r16.recovery import mask_tasks as T
from scripts.r16_recovery_spline import ns_basis

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
REC = ROOT / 'infra/r16/recovery_20260918'


def read(path):
    return json.loads(path.read_text())


def main():
    results = {}
    # Basis identity and a strictly linear, noiseless signal.
    x = np.linspace(0, 1, 301)
    design = ns_basis(x)
    prediction = design @ np.linalg.lstsq(design, x, rcond=None)[0]
    knots = np.r_[[0.] * 4, [1/3, 2/3], [1.] * 4]
    results['spline_basis'] = {
        'endpoint_design_rows': design[[0, -1]].tolist(),
        'linear_signal_r2': float(1 - ((x-prediction)**2).sum()/((x-x.mean())**2).sum()),
        'selected_basis_second_derivatives_at_endpoints': [
            BSpline(knots, np.eye(6)[i], 3).derivative(2)([0, 1]).tolist()
            for i in [1, 2, 3]],
    }
    # Same shape, different additive origin, using the production classifier.
    def classify(v):
        rho, p = spearmanr(np.arange(len(v)), v)
        rel = (v[-1]-v[0])/max(abs(v).max(), 1e-9)
        return ('high-to-low' if rho < 0 else 'low-to-high') if p < .05 and abs(rel) >= .15 else ('flat' if abs(rel) < .15 else 'nonmonotonic')
    v = np.linspace(-.02, .02, 8)
    results['shape_offset'] = {'centered': classify(v), 'plus_half': classify(v+.5)}
    geometry = []
    for row in read(ROOT/'infra/r04/role_manifests/external_validation_manifest.json')['rows']:
        barcodes, xy = read_10x_positions(Path(row['coordinate_locator']))
        labels = usz_tls_labels(row['section_id'], barcodes, ROOT)
        distance = cKDTree(xy[labels==1]).query(xy)[0]
        dn = (np.clip(distance/np.quantile(distance,.99),-.1,1)+.1)/1.1
        alias = row['section_id'].split('::')[-1]
        path = ROOT/'data/other_sources/zenodo_usz_tls_visium/TLS_VISIUM_USZ/h5ad_preprocessed'/f'{alias}.h5ad'
        with h5py.File(path) as h:
            cats = [c.decode() if isinstance(c,bytes) else c for c in h['obs/ground_truth/categories'][:]]
            codes = h['obs/ground_truth/codes'][:]
            idx = [b.decode() if isinstance(b,bytes) else b for b in h['obs/_index'][:]]
        lut = {b:(cats[codes[i]] if codes[i]>=0 else 'UNKNOWN') for i,b in enumerate(idx)}
        keep = np.array([lut.get(b)=='TUM' for b in barcodes])
        xg = dn[keep]; Xg = ns_basis(xg)
        bg = np.linalg.lstsq(Xg,xg,rcond=None)[0]
        pg = Xg@bg
        qs = np.quantile(xg,[0,1/3,2/3,1])
        tg = np.r_[[qs[0]]*4,qs[1:3],[qs[-1]]*4]
        spl = BSpline(tg,np.eye(6),3)
        N = null_space(spl.derivative(2)(qs[[0,-1]]))
        natural = spl(xg)@N
        pnat = natural@np.linalg.lstsq(natural,xg,rcond=None)[0]
        bins = np.linspace(0,.75,9); bc=[]; bv=[]
        for lo,hi in zip(bins[:-1],bins[1:]):
            if ((xg>=lo)&(xg<hi)).sum()>=10:
                c=(lo+hi)/2;bc.append(c);bv.append(np.r_[1,spl(c)[1:4]]@bg)
        geometry.append({'section':alias,'n_tum':int(keep.sum()),
                         'injected_truth':'strictly_linear_increasing',
                         'legacy_shape':classify(np.array(bv)),
                         'legacy_r2':float(1-((xg-pg)**2).sum()/((xg-xg.mean())**2).sum()),
                         'natural_r2':float(1-((xg-pnat)**2).sum()/((xg-xg.mean())**2).sum()),
                         'finite_bin_centers':bc})
    (OUT/'spline_geometry_injection.json').write_text(json.dumps(geometry,indent=2))
    # T5 probe passes although the actual hidden-window predictor changes.
    coords = np.array([(i, j) for i in range(10) for j in range(10)], dtype=float)
    win = {'center': (4.5, 4.5), 'half': 1.5}
    hidden = T.window_mask(coords, win)
    rng = np.random.default_rng(21)
    inp = rng.normal(size=100)
    target = 2*inp + .01*rng.normal(size=100)
    beta = np.linalg.lstsq(np.c_[np.ones((~hidden).sum()), inp[~hidden]], target[~hidden], rcond=None)[0]
    changed = inp.copy(); changed[hidden] = 999
    f1, _ = T.visible_features(coords, inp, win)
    f2, _ = T.visible_features(coords, changed, win)
    p1 = np.c_[np.ones(hidden.sum()), inp[hidden]] @ beta
    p2 = np.c_[np.ones(hidden.sum()), changed[hidden]] @ beta
    results['t5_probe_scope'] = {
        'production_neighbor_probe_passes': bool(np.array_equal(f1, f2, equal_nan=True)),
        'actual_M1_prediction_max_change': float(abs(p1-p2).max()),
    }
    # L009 experimental units and repeated measurements.
    f = pd.read_csv(REC/'masked_structure_fine.tsv', sep='\t')
    f = f[f.mhc2_hit.isin([True, False]) & f.b_hit.isin([True, False])].copy()
    f['delta'] = f.mhc2_hit.astype(float)-f.b_hit.astype(float)
    pats = f.groupby('patient').delta.mean()
    results['l009'] = {
        'rows': len(f), 'patients': f.patient.nunique(),
        'foci': f.groupby(['patient', 'focus']).ngroups,
        'pooled_delta': float(f.delta.mean()), 'equal_patient_delta': float(pats.mean()),
        'patient_deltas': pats.to_dict(),
        'patient_sign_test_p_diagnostic_only': binomtest(int((pats>0).sum()), int((pats!=0).sum())).pvalue,
    }
    # Rank1 specific-increment is not a test of absolute predictive information.
    tr = np.linspace(-2, 2, 200)[:, None]
    te = np.linspace(-1.9, 1.9, 180)[:, None]
    ytr = {'TLS': (tr[:, 0] > 0).astype(float), 'TUMOR_STROMA_BOUNDARY': (tr[:, 0] > .3).astype(float)}
    yte = {'TLS': (te[:, 0] > 0).astype(float), 'TUMOR_STROMA_BOUNDARY': (te[:, 0] > .3).astype(float)}
    ro = evaluate_shared_specific_readout(tr, te, ytr, yte)
    results['r04_rank1_counterexample'] = {
        'shared_absolute_auc': {k: v['auc'] for k, v in ro['shared_only']['metrics'].items()},
        'specific_auc_increment_used_by_final_gate': ro['specific_increment']['auc_delta_by_target'],
    }
    # Frozen definitions and baseline target overlap.
    defs = read(REC/'programs/program_definitions.json')
    proxy = read(ROOT/'infra/r04/marker_proxy_combined.json')
    axes = {k: v['voted_genes'] for k, v in proxy['classes'].items() if k not in A.RETIRED_AXES}
    axes['Plasma'] = A.PLASMA_GENES
    results['t4_overlap'] = {
        pn: {side: {ax: sorted(set(gs)&set(d[side])) for ax, gs in axes.items() if set(gs)&set(d[side])}
             for side in ['input_genes', 'readout_genes']}
        for pn, d in defs.items()}
    t1 = read(REC/'audit/pc_matching_reanalysis.json')
    results['t1_artifact'] = {'corrected': t1['corrected'], 'reproduced': [g for g in t1['groups'] if g['grade_new_rule']=='EXPLORATORY_REPRODUCED']}
    # GO duplicates are deduplicated AFTER random splitting.
    go = read(REC/'programs_go/definitions.json')
    mapping = read(REC/'programs_go/id_map.json')
    mismatch = []
    for rep, members in mapping.items():
        for gid in members:
            if set(go[gid]['input_genes']) != set(go[rep]['input_genes']):
                mismatch.append({'id': gid, 'representative': rep,
                                 'declared_input': go[gid]['input_genes'],
                                 'executed_input': go[rep]['input_genes']})
    axgenes = set().union(*(set(g) for g in axes.values()))
    results['go_definitions'] = {
        'ids': len(go), 'representatives': len(mapping),
        'duplicate_extra_ids': len(go)-len(mapping),
        'duplicate_ids_with_different_declared_split': len(mismatch),
        'examples': mismatch[:3],
        'representatives_target_overlaps_axes': sum(bool(set(go[r]['readout_genes']) & axgenes) for r in mapping),
        'unfrozen_timestamp_representatives': sum(not go[r]['frozen_for_external_at'] for r in mapping),
    }
    inc = pd.read_csv(REC/'programs_go/incremental_prediction_full.tsv', sep='\t')
    gm = inc.groupby(['program','cohort']).dM1_pct.median().unstack()
    hits = gm[(gm['ST-CRC']>20) & (gm['USZ']>20)]
    summary = pd.read_csv(REC/'programs_go/go_summary_full.tsv', sep='\t')
    results['go_counts'] = {'increment_rows': len(inc), 'unique_representatives': inc.program.nunique(),
                           'both_cohorts_over_20_unique_reps': len(hits),
                           'both_cohorts_over_20_expanded_ids': sum(len(mapping[r]) for r in hits.index),
                           'summary_rows': len(summary),
                           'finite_dM1': int(np.isfinite(inc.dM1_pct).sum())}
    # Small authoritative artifacts and source versions, not new execution claims.
    gate = read(ROOT/'infra/r04/r04_final_gate_20260911.json')
    results['r04_gate_hashes'] = {
        r['path']: hashlib.sha256((ROOT/r['path']).read_bytes()).hexdigest()==r['sha256']
        for r in gate['evidence']}
    results['r04_unnamed_gate'] = gate['unnamed_field_reproducibility_gate']
    manifests = {}
    for role in ['training', 'internal_validation', 'external_validation']:
        rows = read(ROOT/f'infra/r04/role_manifests/{role}_manifest.json')['rows']
        manifests[role] = {'sections':len(rows),'patients':len({r['patient_id'] for r in rows}),
                           'lineages': sorted({r['lineage'] for r in rows})}
    results['manifest_units'] = manifests
    paths = list((ROOT/'r16').rglob('*.py')) + list((ROOT/'r04').rglob('*.py')) + list((ROOT/'scripts').glob('*.py'))
    paths += [ROOT/'STATUS.md', ROOT/'docs/plan.md', ROOT/'docs/decisions.md', ROOT/'docs/roadmap.md']
    (OUT/'source_hashes.json').write_text(json.dumps({str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)}, indent=2))
    (OUT/'probes.json').write_text(json.dumps(results, indent=2, ensure_ascii=False))
    print('Saved probes.json, source_hashes.json, spline_geometry_injection.json')


if __name__ == '__main__':
    main()
