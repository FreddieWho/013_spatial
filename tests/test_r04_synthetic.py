from __future__ import annotations

import numpy as np
import pytest

from r04.diagnostics import (
    canonical_subspace_correlations,
    effective_factor_count,
    platform_converged,
)
from r04.synthetic import (
    irregular_domain_mask,
    make_identifiability_calibration_sections,
    make_irregular_continuous_field_sections,
)


@pytest.mark.parametrize("domain", ["crescent", "branch", "disconnected"])
def test_irregular_calibration_domains_are_nonempty_and_nonrectangular(domain: str) -> None:
    sections, truth = make_irregular_continuous_field_sections(
        domain=domain, sections=2, spots_per_section=30, genes=10, seed=11
    )
    assert len(sections) == 2
    assert truth.shape == (60, 2)
    assert np.isfinite(truth).all()
    assert np.unique(sections[0].coords[:, 0]).size > 5
    assert irregular_domain_mask(np.array([[0.0, 0.0], [0.9, 0.9]]), domain).dtype == bool


@pytest.mark.parametrize("mode,k_true", [
    ("no_field", 0),
    ("one_disconnected", 1),
    ("two_independent", 2),
    ("two_collinear", 2),
])
def test_identifiability_calibration_declares_continuous_truth(mode: str, k_true: int) -> None:
    sections, truth, metadata = make_identifiability_calibration_sections(
        mode=mode, domain="disconnected", sections=2, spots_per_section=30, genes=10, seed=19
    )
    assert len(sections) == 2
    assert truth.shape == (60, k_true)
    assert metadata["k_true"] == k_true
    assert metadata["field_definition"].startswith("smooth continuous")
    assert all(section.library_size is not None for section in sections)


def test_subspace_diagnostic_is_invariant_to_factor_rotation() -> None:
    rng = np.random.default_rng(22)
    truth = rng.normal(size=(80, 2))
    rotation = np.array([[0.0, -1.0], [1.0, 0.0]])
    estimate = truth @ rotation
    correlations = canonical_subspace_correlations(truth, estimate)
    assert correlations.tolist() == pytest.approx([1.0, 1.0], abs=1e-6)


def test_effective_factor_count_uses_smallest_one_se_supported_k() -> None:
    result = effective_factor_count(
        {1: 0.805, 2: 0.82, 4: 0.821},
        {1: 0.01, 2: 0.02, 4: 0.02},
        {1: 0.30, 2: 0.25, 4: 0.01},
        null_cutoff=0.05,
        stability_by_k={1: 1.0, 2: 1.0, 4: 1.0},
    )
    assert result["selected_k"] == 1
    assert result["status"] == "K_SELECTED"


def test_platform_convergence_requires_a_recent_plateau() -> None:
    assert platform_converged([10.0, 9.999, 9.9985, 9.9982, 9.9981], window=5)
    assert not platform_converged([10.0, 9.0, 8.0, 7.0, 6.0], window=5)
