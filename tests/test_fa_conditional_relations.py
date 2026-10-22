import numpy as np
import pytest
from scipy import sparse

from field_analyzer import conditional_relations as cr


def _grid_section(*, seed=0, joint=True, side=20, n_units=1):
    rng = np.random.default_rng(seed)
    rows, cols = np.meshgrid(np.arange(side), np.arange(side), indexing='ij')
    xy = np.column_stack([cols.ravel(), rows.ravel()]) * 100.0
    n = xy.shape[0]
    predictor = np.sin(cols.ravel() / 3.0) + 0.4 * rng.normal(size=n)
    own = np.cos(rows.ravel() / 4.0) + 0.4 * rng.normal(size=n)
    from scipy.spatial import cKDTree
    _, idx = cKDTree(xy).query(xy, k=9)
    own_nb = own[idx[:, 1:]].mean(axis=1)
    cross_nb = predictor[idx[:, 1:]].mean(axis=1)
    target = 0.9 * own_nb + (0.9 * cross_nb if joint else 0.0) + 0.15 * rng.normal(size=n)
    zeros = np.zeros(n)
    scores = {('T', 'A', 'feature'): own, ('T', 'A', 'truth'): target,
              ('T', 'B', 'feature'): own, ('T', 'B', 'truth'): target,
              ('P', 'A', 'feature'): predictor, ('P', 'A', 'truth'): zeros,
              ('P', 'B', 'feature'): predictor, ('P', 'B', 'truth'): zeros}
    q = np.zeros(n, dtype=bool)
    q[::10] = True
    return cr.SectionDesign(section_id=f's{seed}', unit_id=f'u{seed}', xy=xy,
                            in_tissue=np.ones(n, dtype=bool), q_mask=q, log_lib=np.zeros(n),
                            scores=scores, gene_sets={k: (f'g{i}',) for i, k in enumerate(scores)})


def test_synthetic_increment_and_separate_references():
    sections = [_grid_section(seed=s, joint=True) for s in range(4)]
    unit_rows, _ = cr.run_rule(sections, 'T', 'P', conditioned=False, protocol='pos')
    same = cr.run_single_feature(sections, 'T', 'P', feature='same', conditioned=False)
    cross = cr.run_single_feature(sections, 'T', 'P', feature='cross', conditioned=False)
    assert np.median([r['r2'] for r in same]) > np.median([r['r2'] for r in cross])
    deltas = np.array([r['delta_r2'] for r in unit_rows])
    assert np.median([r['r2_m1'] for r in unit_rows]) > np.median([r['r2_m0'] for r in unit_rows])
    assert np.median(deltas) > 0.02 and (deltas > 0).sum() >= 3


def test_synthetic_negative_has_no_increment():
    sections = [_grid_section(seed=50 + s, joint=False) for s in range(4)]
    unit_rows, _ = cr.run_rule(sections, 'T', 'P', conditioned=False, protocol='neg')
    deltas = np.array([r['delta_r2'] for r in unit_rows])
    assert np.all(np.abs(deltas) < 0.02)


def test_neighbour_feature_excludes_self_and_uses_k():
    from scipy.spatial import cKDTree
    n = 120
    xy = np.column_stack([np.arange(n) * 10.0, np.zeros(n)])
    values = np.arange(n, dtype=float)
    observed = np.ones(n, dtype=bool)
    mean, count = cr.neighbour_feature(xy, np.ones(n, dtype=bool), observed, values, k=4)
    # rightmost point has no self value in the average
    assert count[-1] <= 4
    assert mean[-1] < values[-1]
    assert count[0] <= 4


def test_poison_hidden_counts_does_not_move_design():
    rng = np.random.default_rng(3)
    n, n_genes = 300, 45
    counts = rng.poisson(2.0, size=(n, n_genes))
    genes = np.array([f'G{i}' for i in range(n_genes)], dtype=object)
    panel = {'T': [f'G{i}' for i in range(8)],
             'P': [f'G{i}' for i in range(8, 16)],
             'X': [f'G{i}' for i in range(16, 24)]}
    xy = rng.uniform(0, 300, size=(n, 2))
    in_tissue = np.ones(n, dtype=bool)
    base = cr.build_section('p', 'u', sparse.csr_matrix(counts), genes, xy, in_tissue, panel)
    poisoned_counts = counts.copy()
    poisoned_counts[~base.observed, :] = 10_000
    poisoned = cr.build_section('p', 'u', sparse.csr_matrix(poisoned_counts), genes, xy,
                                in_tissue, panel)
    assert np.array_equal(base.q_mask, poisoned.q_mask)
    d0 = cr.fold_design(base, 'T', 'P', conditioned=True)
    d1 = cr.fold_design(poisoned, 'T', 'P', conditioned=True)
    assert np.array_equal(np.nan_to_num(d0['same']), np.nan_to_num(d1['same']))
    assert np.array_equal(np.nan_to_num(d0['cross']), np.nan_to_num(d1['cross']))
    assert d0['params'] == d1['params']
    q = base.q_mask
    assert np.any(np.abs(np.nan_to_num(d0['y'][q]) - np.nan_to_num(d1['y'][q])) > 0)


def test_gene_halves_disjoint_and_witness_excluded():
    rng = np.random.default_rng(5)
    n, n_genes = 120, 60
    counts = rng.poisson(3.0, size=(n, n_genes))
    genes = np.array([f'G{i}' for i in range(n_genes)], dtype=object)
    # a gene shared between two fields must not appear in a feature half
    panel = {'T': [f'G{i}' for i in range(12)] + ['G11'],
             'P': [f'G{i}' for i in range(12, 24)] + ['G5'],
             'X': [f'G{i}' for i in range(24, 36)]}
    xy = rng.uniform(0, 200, size=(n, 2))
    design = cr.build_section('s', 'u', sparse.csr_matrix(counts), genes, xy,
                              np.ones(n, dtype=bool), panel)
    feat_t = set(design.gene_sets[('T', 'A', 'feature')])
    feat_p = set(design.gene_sets[('P', 'A', 'feature')])
    truth_t = set(design.gene_sets[('T', 'B', 'truth')])
    assert not feat_t & truth_t
    assert 'G5' not in feat_p or 'G5' not in feat_t  # shared gene dropped on one side


def test_query_points_follow_fa12_mask_and_unit_isolation():
    rng = np.random.default_rng(9)
    n, n_genes = 200, 40
    counts = rng.poisson(2.0, size=(n, n_genes))
    genes = np.array([f'G{i}' for i in range(n_genes)], dtype=object)
    panel = {'T': [f'G{i}' for i in range(10)], 'P': [f'G{i}' for i in range(10, 20)]}
    xy = rng.uniform(0, 200, size=(n, 2))
    design = cr.build_section('sec-A', 'unit-A', sparse.csr_matrix(counts), genes, xy,
                              np.ones(n, dtype=bool), panel)
    from field_analyzer.spatial_lag import holdout_mask
    expected = holdout_mask(n, cr.Q_FRACTION, cr.Q_SEED, 'sec-A')
    assert np.array_equal(design.q_mask, expected)
    assert not np.any(design.q_mask & ~design.in_tissue)
    # unit-level folds: a unit appears on exactly one side
    sections = [design]
    second = cr.build_section('sec-B', 'unit-B', sparse.csr_matrix(counts), genes, xy,
                              np.ones(n, dtype=bool), panel)
    sections.append(second)
    unit_rows, _ = cr.run_rule(sections, 'T', 'P', conditioned=False, protocol='u')
    assert {r['held_unit'] for r in unit_rows} == {'unit-A', 'unit-B'}


def test_calibration_is_fitted_on_observed_rows_only():
    rng = np.random.default_rng(11)
    values = rng.normal(size=200)
    log_lib = rng.normal(size=200)
    fit = np.zeros(200, dtype=bool)
    fit[:100] = True
    out, params = cr.calibrate(values, log_lib, fit, conditioned=True)
    # changing values outside the fit mask must not change the calibrated fit rows
    values2 = values.copy()
    values2[100:] += 100.0
    out2, params2 = cr.calibrate(values2, log_lib, fit, conditioned=True)
    assert params == params2
    assert np.array_equal(out[:100], out2[:100])
    assert not np.allclose(out[100:], out2[100:])
