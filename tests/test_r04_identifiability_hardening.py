from __future__ import annotations

import numpy as np
import pytest

from r04.diagnostics import (
    canonical_permutation_cutoffs,
    mnsf_spatial_effect_matrix,
    platform_convergence_summary,
)
from r04.models.mnsf import (
    _inducing_projection,
    _rate_initial_components,
)
from r04.synthetic import make_identifiability_calibration_sections
from scripts.r04_k_calibration import gene_crossfit_splits
from scripts.r04_k_calibration import required_recovered_directions
from scripts.r04_k_calibration import _fit_intercept_nb


def test_inducing_value_projection_is_identity_at_inducing_points() -> None:
    inducing = np.array([
        [-1.0, -0.5],
        [-0.2, 0.7],
        [0.6, -0.4],
        [1.0, 0.8],
    ])
    basis, _, _ = _inducing_projection(inducing, inducing, 1.7, 1e-6)
    np.testing.assert_allclose(basis, np.eye(len(inducing)), atol=2e-5)


def test_rate_initialisation_has_k_invariant_total_mass() -> None:
    rng = np.random.default_rng(11)
    library = rng.lognormal(size=30)
    counts = rng.poisson(library[:, None] * np.linspace(1.0, 4.0, 12)).astype(float)
    totals = []
    for factors in (1, 2, 4, 6):
        loading, amplitude, baseline, dispersion = _rate_initial_components(
            counts, library, factors, seed=19
        )
        assert loading.shape == (12, factors)
        assert np.all(loading > 0)
        np.testing.assert_allclose(loading.sum(axis=0), 1.0, atol=1e-8)
        assert np.all(amplitude > 0) and np.all(baseline > 0) and np.all(dispersion > 0)
        totals.append(float(amplitude.sum() + baseline.sum()))
    np.testing.assert_allclose(totals, np.repeat(totals[0], len(totals)), rtol=1e-10)


def test_correct_specification_has_two_independent_gene_effects() -> None:
    _, _, metadata = make_identifiability_calibration_sections(
        mode="two_independent",
        domain="disconnected",
        generator="mnsf_correct",
        sections=2,
        spots_per_section=30,
        genes=24,
        seed=23,
    )
    loading = np.asarray(metadata["true_loading"], dtype=float)
    assert np.linalg.matrix_rank(loading - loading.mean(axis=0, keepdims=True)) == 2
    cosine = np.dot(loading[:, 0], loading[:, 1]) / (
        np.linalg.norm(loading[:, 0]) * np.linalg.norm(loading[:, 1])
    )
    assert cosine < 0.5


def test_crescent_control_places_both_fields_on_observed_support() -> None:
    _, _, metadata = make_identifiability_calibration_sections(
        mode="two_independent",
        domain="crescent",
        generator="mnsf_correct",
        sections=3,
        spots_per_section=60,
        genes=24,
        seed=20260807,
    )
    assert min(metadata["field_standard_deviation"]) > 0.15


def test_mnsf_effect_matrix_preserves_continuous_two_direction_truth() -> None:
    _, fields, metadata = make_identifiability_calibration_sections(
        mode="two_independent", genes=24, seed=29
    )
    effect = mnsf_spatial_effect_matrix(
        fields,
        np.asarray(metadata["true_loading"]),
        np.asarray(metadata["true_amplitude"]),
    )
    assert effect.shape == (len(fields), 24)
    assert np.linalg.matrix_rank(effect - effect.mean(axis=0, keepdims=True)) == 2


def test_direction_matched_null_and_platform_reject_rising_loss() -> None:
    rng = np.random.default_rng(31)
    truth = rng.normal(size=(80, 2))
    estimate = truth + rng.normal(scale=0.02, size=truth.shape)
    groups = np.repeat(np.arange(4), 20)
    cutoff = canonical_permutation_cutoffs(
        truth, estimate, groups=groups, seed=37, draws=40
    )
    assert cutoff.shape == (2,)
    assert np.all(cutoff >= 0) and np.all(cutoff <= 1)

    falling_then_flat = list(np.linspace(10.0, 5.0, 100)) + [5.0] * 150
    rising = list(np.linspace(5.0, 6.0, 150))
    assert platform_convergence_summary(
        falling_then_flat, window=50, stable_windows=2, relative_improvement=1e-4
    )["converged"]
    assert not platform_convergence_summary(
        rising, window=50, stable_windows=2, relative_improvement=0.05
    )["converged"]
    stochastic_plateau = [5.0] * 50 + [5.00002] * 50 + [4.99999] * 50
    assert platform_convergence_summary(
        stochastic_plateau, window=50, stable_windows=2, relative_improvement=1e-4
    )["converged"]


def test_gene_crossfit_never_uses_evaluation_genes_for_adaptation() -> None:
    splits = gene_crossfit_splits(24, folds=3, seed=41)
    assert len(splits) == 3
    observed_evaluation = []
    for adaptation, evaluation in splits:
        assert len(adaptation) > 0 and len(evaluation) > 0
        assert set(adaptation).isdisjoint(evaluation)
        observed_evaluation.extend(evaluation.tolist())
    assert sorted(observed_evaluation) == list(range(24))


def test_gene_crossfit_rejects_degenerate_requests() -> None:
    with pytest.raises(ValueError):
        gene_crossfit_splits(3, folds=1, seed=1)
    with pytest.raises(ValueError):
        gene_crossfit_splits(3, folds=4, seed=1)


def test_k_recovery_requires_only_representable_truth_directions() -> None:
    assert required_recovered_directions(1, 2) == 1
    assert required_recovered_directions(2, 2) == 2
    assert required_recovered_directions(4, 2) == 2


def test_intercept_nb_baseline_is_positive_and_gene_aligned() -> None:
    sections, _, _ = make_identifiability_calibration_sections(
        mode="no_field", sections=2, spots_per_section=20, genes=8, seed=43
    )
    rate, dispersion = _fit_intercept_nb(sections)
    assert rate.shape == dispersion.shape == (8,)
    assert np.all(rate > 0) and np.all(dispersion > 0)
