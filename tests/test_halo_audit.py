"""Conclusion-affecting regression tests; no real-data fitting."""
import importlib.util
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / f'{name}.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_rank_deficient_residual_matches_least_squares():
    m = module('gobp_halo_screen')
    rng = np.random.default_rng(43)
    xy = rng.normal(size=(90, 2))
    X = m.design(np.zeros(90), xy)
    y = rng.normal(size=(90, 3))
    expected = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    actual = m.residualize(y, m.residual_operator(X))
    np.testing.assert_allclose(actual, expected, atol=1e-12)


def test_redundant_covariates_do_not_change_residual():
    m = module('gobp_halo_screen')
    rng = np.random.default_rng(4)
    X = np.column_stack([np.ones(80), rng.normal(size=(80, 2))])
    y = rng.normal(size=80)
    X2 = np.column_stack([X, X[:, 0], X[:, 1], np.zeros(80)])
    np.testing.assert_allclose(m.residualize(y, m.residual_operator(X)),
                               m.residualize(y, m.residual_operator(X2)), atol=1e-12)


def test_sigma_recovers_analytic_gaussian_kernel():
    m = module('gobp_halo_null')
    sigma = 300.0
    corr = {(lo, hi): np.exp(-((lo+hi)/2)**2/(4*sigma**2))
            for lo, hi in m.LAG_BANDS}
    np.testing.assert_allclose(m.sigma_from_bands(corr), sigma, rtol=1e-12)


def test_unfit_correlation_is_unknown():
    m = module('gobp_halo_null')
    assert np.isnan(m.sigma_from_bands({(100, 200): .3}))
    assert np.isnan(m.sigma_from_bands({(100, 200): -.3, (200, 400): -.2}))


def test_visium_six_neighbors_have_same_pitch():
    m = module('tls_pool_expand')
    # Doubled-column Visium metadata: row and column have equal parity.
    xy = m.hex_xy([2, 2, 2, 1, 1, 3, 3], [2, 0, 4, 1, 3, 1, 3], 100)
    np.testing.assert_allclose(np.linalg.norm(xy[1:] - xy[0], axis=1), 100)
