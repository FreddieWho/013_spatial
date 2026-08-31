from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from r04.continuation_platform_audit import audit_continuation_manifest, replace_final_live_evaluation
from r04.diagnostics import deterministic_objective_platform_summary

PROTOCOL = {"checkpoint_step": 8400, "mc_draws": 4, "repeat_seed_offset": 7000}


def _manifest(tmp_path: Path) -> dict[str, object]:
    cells = []
    for restart in (2, 3, 4):
        for factors in (0, 3):
            for fold in (0, 4):
                seed = 1000 + restart * 10 + factors + fold
                checkpoint_dir = f"checkpoints/r{restart}/k{factors}/f{fold}"
                source = {
                    "schema": "r04.late_lr_cell.v1", "status": "FIT_ONLY_DIAGNOSTIC",
                    "restart_index": restart, "k_model": factors, "fold": fold,
                    "input_hash": f"input-{restart}-{fold}", "training_patients": ["P1", "P2"],
                    "manifest_hash": "manifest-1", "optimization_seed": seed,
                    "checkpoint_dir": checkpoint_dir, "protocol": PROTOCOL,
                    "config_hash": f"config-{restart}-{factors}-{fold}",
                    "fit_platform": {"converged": False},
                    "fit_diagnostics": {
                        "steps": 8400, "start_step": 7200, "evaluation_mc_draws": 4,
                        "environment_hash": "environment-1",
                        "evaluation_trace": [
                            {"step": 7600, "objective_mean": 10.0,
                             "objective_kind": "dense_full_panel_fixed_stateless_mc", "mc_draws": 4},
                            {"step": 8000, "objective_mean": 9.995,
                             "objective_kind": "dense_full_panel_fixed_stateless_mc", "mc_draws": 4},
                            {"step": 8400, "objective_mean": 9.99,
                             "objective_kind": "dense_full_panel_fixed_stateless_mc", "mc_draws": 4},
                        ],
                    },
                }
                source_path = tmp_path / f"source-r{restart}-k{factors}-f{fold}.json"
                source_path.write_text(json.dumps(source), encoding="utf-8")
                checkpoint_root = tmp_path / checkpoint_dir
                checkpoint_root.mkdir(parents=True)
                for name, content in {
                    "checkpoint": 'model_checkpoint_path: "ckpt-8400"\n',
                    "checkpoint_metadata.json": json.dumps({
                        "checkpoint_schema": "r04.mnsf_checkpoint.v2",
                        "input_hash": source["input_hash"],
                        "config_hash": source["config_hash"],
                        "environment_hash": "environment-1",
                    }),
                    "ckpt-8400.index": "index",
                    "ckpt-8400.data-00000-of-00001": "data",
                }.items():
                    (checkpoint_root / name).write_text(content, encoding="utf-8")
                checkpoint = str(Path(checkpoint_dir) / "ckpt-8400")
                repeat = {
                    "schema": "r04.checkpoint_audit.v1", "status": "DIAGNOSTIC_ONLY",
                    "factors": factors, "fold": fold, "mc_draws": 4,
                    "manifest_hash": "manifest-1", "gene_count": 2,
                    "folds": 5, "group_seed": 20260807, "inducing_points": 16,
                    "lengthscale": 3.0, "ridge": 1e-5,
                    "checkpoints": [
                        {"checkpoint": checkpoint, "global_step": 8400, "mc_draws": 4,
                         "seed": seed + 7000, "objective_mean": 9.99,
                         "objective_sd": 0.0, "objective_draws": [9.99] * 4,
                         "n_genes": 2, "n_spots": 10,
                         "nuisance_gauge": {"rank": 1, "h_min": [0.1], "h_mean": [0.2],
                                            "h_max": [0.3], "baseline_v_projection": [0.4]}},
                        {"checkpoint": checkpoint, "global_step": 8400, "mc_draws": 4,
                         "seed": seed + 7000, "objective_mean": 9.99,
                         "objective_sd": 0.0, "objective_draws": [9.99] * 4,
                         "n_genes": 2, "n_spots": 10,
                         "nuisance_gauge": {"rank": 1, "h_min": [0.1], "h_mean": [0.2],
                                            "h_max": [0.3], "baseline_v_projection": [0.4]}},
                    ],
                }
                repeat_path = tmp_path / f"audit-r{restart}-k{factors}-f{fold}.json"
                repeat_path.write_text(json.dumps(repeat), encoding="utf-8")
                cells.append({
                    "cell_id": f"r{restart}_k{factors}_f{fold}",
                    "restart_index": restart, "k_model": factors, "fold": fold,
                    "source_cell": str(source_path.relative_to(tmp_path)),
                    "checkpoint_repeat_audit": str(repeat_path.relative_to(tmp_path)),
                    "input_hash": source["input_hash"], "training_patients": source["training_patients"],
                    "manifest_hash": "manifest-1", "optimization_seed": seed, "protocol": PROTOCOL,
                })
    return {
        "checkpoint_audit_protocol": {
            "folds": 5, "group_seed": 20260807, "inducing_points": 16,
            "lengthscale": 3.0, "ridge": 1e-5, "mc_draws": 4,
            "repeat_seed_offset": 7000,
        },
        "cells": cells,
    }


def test_real_checkpoint_audit_schema_and_frozen_replacement(tmp_path: Path) -> None:
    result = audit_continuation_manifest(_manifest(tmp_path), tmp_path)
    assert result["status"] == "ALL_PASS"
    assert result["cells"][0]["exact_repeat"] is True
    assert result["cells"][0]["dense"]["objective_means"][-1] == 9.99
    trace = replace_final_live_evaluation(
        [{"step": 8400, "objective_mean": 1.0}],
        {"schema": "r04.checkpoint_audit.v1", "checkpoints": [
            {"global_step": 8400, "objective_mean": 2.0},
            {"global_step": 8400, "objective_mean": 2.0},
        ]},
    )
    assert trace[-1]["objective_mean"] == 2.0


def test_dense_threshold_boundary_is_off_platform(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path)
    source_path = tmp_path / manifest["cells"][0]["source_cell"]
    source = json.loads(source_path.read_text())
    source["fit_diagnostics"]["evaluation_trace"] = [
        {"step": step, "objective_mean": objective,
         "objective_kind": "dense_full_panel_fixed_stateless_mc", "mc_draws": 4}
        for step, objective in (
            (7600, 1_000_000.0),
            (8000, 998_993.474),
            (8400, 997_987.0),
        )
    ]
    source_path.write_text(json.dumps(source), encoding="utf-8")
    repeat_path = tmp_path / manifest["cells"][0]["checkpoint_repeat_audit"]
    repeat = json.loads(repeat_path.read_text(encoding="utf-8"))
    for checkpoint in repeat["checkpoints"]:
        checkpoint["objective_mean"] = 997_987.0
        checkpoint["objective_sd"] = 0.0
        checkpoint["objective_draws"] = [997_987.0] * 4
    repeat_path.write_text(json.dumps(repeat), encoding="utf-8")
    result = audit_continuation_manifest(manifest, tmp_path)
    assert result["status"] == "STOPPED_OFF_PLATFORM"
    assert result["inference_authorized"] is False


@pytest.mark.parametrize("mutation,expected", [
    (lambda m: m["cells"].pop(), "exactly_12_entries_required"),
    (lambda m: m["cells"].append(copy.deepcopy(m["cells"][0])), "duplicate_cell"),
    (lambda m: m["cells"][0].update(input_hash="wrong"), "input_hash_mismatch"),
])
def test_panel_invalid_or_duplicate_fails_closed(tmp_path: Path, mutation, expected: str) -> None:
    manifest = _manifest(tmp_path)
    mutation(manifest)
    result = audit_continuation_manifest(manifest, tmp_path)
    assert result["status"] == "INVALID_AUDIT"
    assert any(expected in problems for problems in result["failing_cells"].values())


def test_repeat_mismatch_is_off_platform(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path)
    path = tmp_path / manifest["cells"][0]["checkpoint_repeat_audit"]
    audit = json.loads(path.read_text())
    audit["checkpoints"][1]["objective_draws"] = [9.98] * 4
    audit["checkpoints"][1]["objective_mean"] = 9.98
    audit["checkpoints"][1]["objective_sd"] = 0.0
    path.write_text(json.dumps(audit), encoding="utf-8")
    result = audit_continuation_manifest(manifest, tmp_path)
    assert result["status"] == "STOPPED_OFF_PLATFORM"
    assert "repeat_objective_draws_mismatch" in result["failing_cells"]["r2_k0_f0"]


def test_checkpoint_bundle_is_recorded(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path)
    checkpoint_dir = tmp_path / "checkpoints/r2/k0/f0"
    (checkpoint_dir / "ckpt-8400.index").write_text("changed", encoding="utf-8")
    result = audit_continuation_manifest(manifest, tmp_path)

    assert result["status"] == "ALL_PASS"
    assert result["cells"][0]["checkpoint_bundle"]["files"]


def test_hash_and_pair_provenance_mismatch_are_invalid(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path)
    source_ref = manifest["cells"][0]["source_cell"]
    source_path = tmp_path / source_ref
    source_path.write_text(source_path.read_text() + " ", encoding="utf-8")
    manifest["cells"][0]["source_cell"] = {"path": source_ref, "sha256": "0" * 64}
    manifest["cells"][2]["input_hash"] = "different"
    result = audit_continuation_manifest(manifest, tmp_path)
    assert result["status"] == "INVALID_AUDIT"
    assert "invalid_artifact:ValueError" in result["artifact_errors"]["r2_k0_f0"]
    assert "input_hash_mismatch" in result["provenance_errors"]["r2_k3_f0"]


def test_pair_provenance_comes_from_source_artifacts(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path)
    entry = manifest["cells"][2]
    source_path = tmp_path / entry["source_cell"]
    source = json.loads(source_path.read_text(encoding="utf-8"))
    source["input_hash"] = "different"
    source_path.write_text(json.dumps(source), encoding="utf-8")
    entry["input_hash"] = "different"
    result = audit_continuation_manifest(manifest, tmp_path)

    assert result["status"] == "INVALID_AUDIT"
    assert "k_pair_provenance_mismatch" in result["provenance_errors"]["r2_k0_f0"]
    assert "k_pair_provenance_mismatch" in result["provenance_errors"]["r2_k3_f0"]
    assert result["passed_cell_count"] == 10


@pytest.mark.parametrize("mutation", [
    lambda audit: audit["checkpoints"][0].update(objective_draws="not-a-list"),
    lambda audit: audit["checkpoints"][0].update(objective_mean=123.0),
])
def test_malformed_repeat_payload_is_invalid(tmp_path: Path, mutation) -> None:
    manifest = _manifest(tmp_path)
    path = tmp_path / manifest["cells"][0]["checkpoint_repeat_audit"]
    audit = json.loads(path.read_text(encoding="utf-8"))
    mutation(audit)
    path.write_text(json.dumps(audit), encoding="utf-8")

    result = audit_continuation_manifest(manifest, tmp_path)
    assert result["status"] == "INVALID_AUDIT"
    assert result["inference_authorized"] is False


def test_repeat_ridge_drift_is_invalid(tmp_path: Path) -> None:
    manifest = _manifest(tmp_path)
    path = tmp_path / manifest["cells"][0]["checkpoint_repeat_audit"]
    audit = json.loads(path.read_text(encoding="utf-8"))
    audit["ridge"] = 1e-4
    path.write_text(json.dumps(audit), encoding="utf-8")

    result = audit_continuation_manifest(manifest, tmp_path)
    assert result["status"] == "INVALID_AUDIT"
    assert "checkpoint_audit_ridge_mismatch" in result["provenance_errors"]["r2_k0_f0"]


def test_diagnostic_threshold_function_keeps_boundary_fail() -> None:
    result = deterministic_objective_platform_summary([
        {"step": 1, "objective_mean": 1_000_000.0},
        {"step": 2, "objective_mean": 998_993.474},
        {"step": 3, "objective_mean": 997_987.0},
    ])
    assert result["converged"] is False
