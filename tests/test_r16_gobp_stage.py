import numpy as np
import pytest

from r16.recovery.corrected import feature_block, masked_predict
from r16.recovery.gobp_stage import FeatureEngine, mask_predictions, geometry, sign_flip_p


def fixture_data():
    rng = np.random.default_rng(420)
    counts = rng.poisson(.7, (144, 35)).astype(float)
    counts *= rng.integers(1, 12, (144, 1))
    xy = np.array([(i, j) for i in range(12) for j in range(12)], float)
    genes = [f'g{i}' for i in range(35)]
    axes = {str(i): genes[4 + i * 3:7 + i * 3] for i in range(6)}
    axes['0'] += ['g0', 'g2']
    defs = {'test': {'input_genes': genes[:2], 'readout_genes': genes[2:4]}}
    return counts, xy, genes, axes, defs


@pytest.mark.parametrize('mode', ['raw', 'depth_normalized'])
def test_sparse_engine_and_real_mask_match_direct_oracle(mode):
    counts, xy, genes, axes, defs = fixture_data()
    engine = FeatureEngine(counts, genes, axes, defs)
    ids, x, y, bad = next(engine.blocks(mode))
    direct = feature_block(counts, genes, genes[:2], genes[2:4], axes, mode)
    assert ids == ['test'] and not bad
    np.testing.assert_allclose(x[:, 0], direct['X'], atol=1e-12)
    np.testing.assert_allclose(y[:, 0], direct['y'], atol=1e-12)
    hidden = (xy[:, 0] >= 4) & (xy[:, 0] <= 7) & (xy[:, 1] >= 4) & (xy[:, 1] <= 7)
    geo = geometry(xy, hidden)
    pred = mask_predictions(x, y, geo)
    oracle = masked_predict(counts, xy, genes, genes[:2], genes[2:4], axes, hidden, mode)
    for i, name in enumerate(['M0', 'M1', 'M2', 'NN', 'KNN8']):
        np.testing.assert_allclose(pred[:, 0, i], oracle[name], atol=2e-8, rtol=2e-8)
    changed = counts.copy()
    changed[hidden] = 1e6
    _, xa, ya, _ = next(FeatureEngine(changed, genes, axes, defs).blocks(mode))
    altered = mask_predictions(xa, ya, geo)
    np.testing.assert_array_equal(pred, altered)


def test_small_patient_test_cannot_reject_even_perfect_direction():
    assert sign_flip_p(np.array([1., 2., 3.])) == .25
    assert sign_flip_p(np.array([-1., -2., -3.])) == .25
    assert sign_flip_p(np.zeros(3)) == 1.


def test_missing_source_and_empty_axis_are_explicit():
    counts, _, genes, axes, defs = fixture_data()
    engine = FeatureEngine(counts, genes, axes, defs, unavailable={'g0'})
    ids, _, _, bad = next(engine.blocks('raw'))
    assert not ids and bad == [('test', 'MISSING_PROGRAM_GENES')]
    axes['0'] = ['g0', 'g2']
    ids, _, _, bad = next(FeatureEngine(counts, genes, axes, defs).blocks('raw'))
    assert not ids and bad == [('test', 'NOT_SEPARABLE')]
