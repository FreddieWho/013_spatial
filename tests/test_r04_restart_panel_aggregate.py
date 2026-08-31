from __future__ import annotations

import copy

import numpy as np
import pytest

from scripts.r04_restart_panel_aggregate import (
    aggregate_restart_scores,
    pairwise_subspace_stability,
    validated_reused_subspace,
)


def _cell(*, fold: int, k: int, offset: float) -> dict:
    patients = {"p1": -10.0 + offset, "p2": -20.0 + offset}
    return {
        "schema": "r04.real_k_cell.v1",
        "status": "FIT_AND_SCORED",
        "fold": fold,
        "k_model": k,
        "input_hash": f"fold-{fold}",
        "validation_patients": ["p1", "p2"],
        "patient_scores": patients,
        "parameters": {
            "gene_folds": 2,
            "inducing_points": 16,
            "inference_steps": 2400,
            "lengthscale": 3.0,
            "nonspatial_rank": 1,
            "optimization_schedule": "joint",
            "optimization_seed": 1000 + fold + k,
            "shared_steps": 0,
            "split_seed": 20260817,
            "steps": 7200,
        },
        "fit_platform": {"converged": False},
        "inference_platform": {"converged": True},
        "fit_diagnostics": {
            "deterministic_objective_platform": {"converged": True},
            "factor_collapse_warning": False,
            "field_collapse_warning": False,
            "loading_collapse_warning": False,
        },
    }


def _entries() -> list[dict]:
    entries = []
    for restart in range(5):
        for fold, delta in ((0, -2.0 - restart), (4, 4.0 + restart)):
            k0 = _cell(fold=fold, k=0, offset=0.0)
            k3 = _cell(fold=fold, k=3, offset=delta)
            if restart == 0:
                k0["status"] = "FIT_AND_SCORED_DUAL_GATE_DIAGNOSTIC"
                k3["status"] = "FIT_AND_SCORED_DUAL_GATE_DIAGNOSTIC"
            entries.append(
                {
                    "restart_index": restart,
                    "fold": fold,
                    "k0": k0,
                    "k3": k3,
                }
            )
    return entries


def test_restart_score_aggregate_preserves_fold_directions() -> None:
    result = aggregate_restart_scores(_entries())

    assert result["restart_indices"] == [0, 1, 2, 3, 4]
    assert result["folds"]["0"]["mean_delta_by_restart"] == [
        -2.0,
        -3.0,
        -4.0,
        -5.0,
        -6.0,
    ]
    assert result["folds"]["4"]["mean_delta_by_restart"] == [
        4.0,
        5.0,
        6.0,
        7.0,
        8.0,
    ]
    assert result["folds"]["0"]["direction_preserved_across_restarts"]
    assert result["folds"]["4"]["direction_preserved_across_restarts"]
    assert result["platform_counts"] == {
        "cells": 20,
        "legacy_fit_converged": 0,
        "dense_objective_converged": 20,
        "inference_converged": 20,
    }


def test_restart_score_aggregate_rejects_patient_mismatch() -> None:
    entries = _entries()
    entries[0] = copy.deepcopy(entries[0])
    entries[0]["k3"]["patient_scores"].pop("p2")
    entries[0]["k3"]["validation_patients"].remove("p2")

    with pytest.raises(ValueError, match="patient score keys differ"):
        aggregate_restart_scores(entries)


def test_restart_score_aggregate_rejects_protocol_drift() -> None:
    entries = _entries()
    entries[-1] = copy.deepcopy(entries[-1])
    entries[-1]["k3"]["parameters"]["inference_steps"] = 1200

    with pytest.raises(ValueError, match="protocol differs"):
        aggregate_restart_scores(entries)


def test_restart_score_aggregate_surfaces_explicit_protocol_exception() -> None:
    entries = _entries()
    entries[0] = copy.deepcopy(entries[0])
    entries[0]["k0"]["parameters"]["steps"] = 7800

    result = aggregate_restart_scores(
        entries,
        allowed_protocol_exceptions=[
            {
                "restart_index": 0,
                "fold": 0,
                "k_model": 0,
                "field": "steps",
                "reference": 7200,
                "observed": 7800,
                "reason": "documented continuation",
            }
        ],
    )

    assert result["protocol"]["steps"] == 7200
    assert result["protocol_exceptions"][0]["observed"] == 7800


def test_restart_score_aggregate_replays_current_dense_platform_rule() -> None:
    entries = _entries()
    entries[0] = copy.deepcopy(entries[0])
    entries[0]["k0"]["fit_objective_platform"] = {
        "converged": False,
        "steps": [7201, 7500, 7800],
        "objective_means": [100.0, 100.01, 100.0],
        "relative_improvements": [-0.0001, 0.0001],
    }

    result = aggregate_restart_scores(entries)

    assert result["platform_counts"]["dense_objective_converged"] == 20
    assert result["pairs"][0]["k0_platform"]["dense_objective_source"].startswith(
        "top_level_"
    )


def test_pairwise_subspace_stability_ignores_rotation_and_factor_order() -> None:
    rng = np.random.default_rng(7)
    base = rng.normal(size=(80, 3))
    rotated = base @ np.asarray(
        [[0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0]]
    )
    result = pairwise_subspace_stability({0: base, 1: rotated})

    assert result["pair_count"] == 1
    assert result["pairs"][0]["canonical_correlations"] == pytest.approx(
        [1.0, 1.0, 1.0]
    )
    assert result["minimum_canonical_correlation"] == pytest.approx(1.0)
    assert [item["mean"] for item in result["canonical_correlation_by_direction"]] == pytest.approx(
        [1.0, 1.0, 1.0]
    )


def test_reused_subspace_requires_matching_source_hashes() -> None:
    provenance = [{"path": "cell.json", "sha256": "a" * 64}]
    training = {"path": "training.json", "sha256": "b" * 64}
    previous = {
        "schema": "r04.restart_stability_panel.v1",
        "input_provenance": provenance,
        "training_manifest": training,
        "subspace_stability": {"0": {"value": 1}, "4": {"value": 2}},
    }

    assert validated_reused_subspace(
        previous, input_provenance=provenance, training_manifest=training
    )["4"]["value"] == 2
    with pytest.raises(ValueError, match="provenance differs"):
        validated_reused_subspace(
            previous,
            input_provenance=[{"path": "cell.json", "sha256": "c" * 64}],
            training_manifest=training,
        )
