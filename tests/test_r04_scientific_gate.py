from __future__ import annotations

import json
from pathlib import Path

from scripts.r04_scientific_gate import (
    _calibration_result_reasons,
    _input_validation_reasons,
    _nested_provenance_reasons,
    _recomputed_aggregate_reasons,
)


def test_scientific_gate_fails_closed_on_unbound_success_statuses() -> None:
    reasons = _input_validation_reasons(
        {"schema": "legacy", "status": "PASS_OBJECTIVE_ORACLE"},
        [],
        {"schema": "legacy", "final_k_status": "K_SELECTED"},
        {"schema": "legacy", "status": "K0_SELECTED"},
    )
    assert "objective_schema_invalid" in reasons
    assert "objective_gradient_oracle_missing" in reasons
    assert "calibration_matrix_incomplete_or_duplicated" in reasons
    assert "k_calibration_schema_invalid" in reasons
    assert "k0_gate_schema_invalid" in reasons


def test_scientific_gate_rejects_numeric_oracle_outside_tolerance() -> None:
    objective = {
        "schema": "r04.objective_oracle.v1",
        "fit_infer_dense_scaling": "per_spot_gene_sum",
        "dense_objective": 10.0,
        "full_batch_objective": 10.0,
        "mean_gene_batch_objective": 10.1,
        "objective_abs_error": 0.1,
        "mean_batch_gradient_max_abs_error": 0.2,
        "gradient_finite": True,
        "batch_count": 20,
    }
    reasons = _input_validation_reasons(objective, [], None, None)
    assert "objective_or_gradient_oracle_tolerance_failed" in reasons


def test_scientific_gate_rejects_nan_canonical_correlation() -> None:
    calibration = {
        "schema": "r04.identifiability_calibration.v2",
        "model": "mnsf",
        "seed": 20260807,
        "truth": {
            "mode": "two_independent",
            "domain": "disconnected",
            "generator": "mnsf_correct",
            "generator_version": "independent_gene_effects_v2",
            "k_eff_true": 2,
        },
        "platform": {"converged": True},
        "subspace": {
            "canonical_correlations": [0.9, float("nan")],
            "direction_matched_null_95": [0.2, 0.2],
        },
    }
    assert _calibration_result_reasons(calibration) == [
        "planted_effect_subspace_not_recovered_above_null"
    ]


def test_scientific_gate_rejects_short_calibration_contract() -> None:
    calibration = {
        "schema": "r04.identifiability_calibration.v2",
        "status": "CALIBRATION_COMPLETE_NOT_VALIDATED",
        "model": "mnsf",
        "seed": 20260807,
        "steps": 1,
        "n_factors_model": 2,
        "config": {
            "inducing_points": 16,
            "nonspatial_rank": 0,
            "lengthscale": 3.0,
            "learning_rate": 0.05,
            "null_draws": 200,
            "platform_window": 50,
        },
        "truth": {
            "mode": "two_independent",
            "domain": "disconnected",
            "generator": "mnsf_correct",
            "generator_version": "independent_gene_effects_v2",
        },
    }
    reasons = _input_validation_reasons({}, [calibration], None, None)
    assert "calibration_execution_contract_invalid" in reasons


def test_scientific_gate_recomputes_nested_source_hashes(tmp_path) -> None:
    source = tmp_path / "source.json"
    source.write_text('{"schema":"r04.k_calibration.v2"}\n', encoding="utf-8")
    aggregate = {
        "input_provenance": [{
            "path": str(source),
            "bytes": source.stat().st_size,
            "sha256": "0" * 64,
        }]
    }
    assert _nested_provenance_reasons(
        aggregate, label="k_calibration"
    ) == ["k_calibration_source_hash_mismatch"]


def test_scientific_gate_recomputes_nested_aggregate_semantics(
    tmp_path, monkeypatch
) -> None:
    source = tmp_path / "source.json"
    source.write_text("{}\n", encoding="utf-8")
    aggregate = {
        "status": "forged",
        "input_provenance": [{"path": str(source)}],
    }
    monkeypatch.setattr(
        "scripts.r04_scientific_gate.aggregate_k_shards",
        lambda shards: {"status": "recomputed"},
    )
    assert _recomputed_aggregate_reasons(
        aggregate, label="k_calibration"
    ) == ["k_calibration_recomputed_aggregate_mismatch"]


def test_current_delivered_k_aggregates_replay_from_nested_sources() -> None:
    root = Path(__file__).resolve().parents[1]
    for label, relative_path in (
        ("k_calibration", "infra/r04/k_calibration_aggregate_v4_20260813.json"),
        ("k0_gate", "infra/r04/k0_gate_aggregate_v4_20260813.json"),
    ):
        aggregate = json.loads((root / relative_path).read_text(encoding="utf-8"))
        assert _recomputed_aggregate_reasons(aggregate, label=label) == []
