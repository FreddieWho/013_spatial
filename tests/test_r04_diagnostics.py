from __future__ import annotations

import pytest

from r04.diagnostics import deterministic_objective_platform_summary


def test_deterministic_objective_platform_uses_final_two_intervals() -> None:
    result = deterministic_objective_platform_summary([
        {"step": 6000, "objective_mean": 100.0},
        {"step": 6300, "objective_mean": 99.95},
        {"step": 6600, "objective_mean": 99.90},
        {"step": 6900, "objective_mean": 99.85},
    ])

    assert result["converged"] is True
    assert result["steps"] == [6300, 6600, 6900]
    assert result["relative_improvements"] == pytest.approx([
        0.0005002501250625,
        0.0005005005005005,
    ])


def test_deterministic_objective_platform_allows_small_rising_objective() -> None:
    result = deterministic_objective_platform_summary([
        {"step": 6000, "objective_mean": 100.0},
        {"step": 6300, "objective_mean": 99.95},
        {"step": 6600, "objective_mean": 99.96},
    ])

    assert result["converged"] is True


def test_deterministic_objective_platform_rejects_large_rising_objective() -> None:
    result = deterministic_objective_platform_summary([
        {"step": 6000, "objective_mean": 100.0},
        {"step": 6300, "objective_mean": 99.95},
        {"step": 6600, "objective_mean": 100.2},
    ])

    assert result["converged"] is False


def test_deterministic_objective_platform_requires_three_evaluations() -> None:
    result = deterministic_objective_platform_summary([
        {"step": 7200, "objective_mean": 99.0},
        {"step": 7500, "objective_mean": 98.99},
    ])

    assert result["converged"] is False
    assert result["required_evaluations"] == 3
