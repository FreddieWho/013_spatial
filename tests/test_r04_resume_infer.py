from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from r04.continuation_platform_audit import checkpoint_bundle_fingerprint
from scripts.r04_resume_infer import (
    EXPECTED_FROZEN_INFERENCE_PROTOCOL,
    _checkpoint_environment_hash,
    _objective_trace,
    _require_continuation_audit,
    _require_frozen_inference_endpoint,
)


def test_objective_trace_accepts_continuation_cell_audit() -> None:
    trace = _objective_trace({
        "fit_diagnostics": {
            "evaluation_trace": [
                {"step": 6600, "objective_mean": 10.0},
                {"step": 6900, "objective_mean": 9.99},
            ]
        }
    })

    assert trace[-1]["step"] == 6900


def test_checkpoint_environment_hash_is_read_from_checkpoint_metadata(tmp_path: Path) -> None:
    checkpoint_dir = tmp_path / "checkpoints"
    checkpoint_dir.mkdir()
    (checkpoint_dir / "checkpoint_metadata.json").write_text(
        json.dumps({"environment_hash": "source-hash"}), encoding="utf-8"
    )
    assert _checkpoint_environment_hash(checkpoint_dir) == "source-hash"
    (checkpoint_dir / "checkpoint_metadata.json").write_text(
        json.dumps({"environment_hash": ""}), encoding="utf-8"
    )
    with pytest.raises(RuntimeError, match="no environment hash"):
        _checkpoint_environment_hash(checkpoint_dir)


def test_objective_trace_accepts_checkpoint_audit() -> None:
    trace = _objective_trace({
        "checkpoints": [
            {"global_step": 6600, "objective_mean": 10.0},
            {"global_step": 7200, "objective_mean": 9.99},
        ]
    })

    assert trace == [
        {"step": 6600, "objective_mean": 10.0},
        {"step": 7200, "objective_mean": 9.99},
    ]


def _checkpoint_dir(tmp_path: Path) -> Path:
    root = tmp_path / "checkpoints"
    root.mkdir()
    for name, content in {
        "checkpoint": 'model_checkpoint_path: "ckpt-8400"\n',
        "checkpoint_metadata.json": "{}",
        "ckpt-8400.index": "index",
        "ckpt-8400.data-00000-of-00001": "data",
    }.items():
        (root / name).write_text(content, encoding="utf-8")
    return root


def test_v1_continuation_source_requires_panel_audit(tmp_path: Path) -> None:
    source = {
        "schema": "r04.late_lr_cell.v1",
        "restart_index": 2,
        "fit_diagnostics": {"start_step": 7200, "steps": 8400},
    }
    with pytest.raises(RuntimeError, match="require a panel audit"):
        _require_continuation_audit(source, None, tmp_path / "cell.json", _checkpoint_dir(tmp_path))


def test_stopped_panel_cannot_authorize_inference(tmp_path: Path) -> None:
    source_path = tmp_path / "cell.json"
    source = {
        "schema": "r04.late_lr_cell.v1",
        "restart_index": 2,
        "fit_diagnostics": {"start_step": 7200, "steps": 8400},
    }
    source_path.write_text(json.dumps(source), encoding="utf-8")
    audit_path = tmp_path / "panel.json"
    audit_path.write_text(json.dumps({
        "schema": "r04.continuation_platform_audit.v2",
        "status": "STOPPED_OFF_PLATFORM",
        "paired_gate": False,
        "inference_authorized": False,
        "cells": [],
    }), encoding="utf-8")

    with pytest.raises(RuntimeError, match="not ALL_PASS"):
        _require_continuation_audit(
            source, audit_path, source_path, _checkpoint_dir(tmp_path)
        )


def test_continuation_gate_binds_source_and_checkpoint_hashes(tmp_path: Path) -> None:
    checkpoint_dir = _checkpoint_dir(tmp_path)
    source_path = tmp_path / "cell.json"
    source = {
        "schema": "r04.late_lr_cell.v1",
        "restart_index": 2,
        "k_model": 3,
        "fold": 0,
        "fit_diagnostics": {"start_step": 7200, "steps": 8400},
    }
    source_path.write_text(json.dumps(source), encoding="utf-8")
    source_sha256 = hashlib.sha256(source_path.read_bytes()).hexdigest()
    repeat_path = tmp_path / "repeat.json"
    repeat_path.write_text(json.dumps({"repeat": "artifact"}), encoding="utf-8")
    repeat_sha256 = hashlib.sha256(repeat_path.read_bytes()).hexdigest()
    audit_path = tmp_path / "panel.json"
    audit_path.write_text(json.dumps({
        "schema": "r04.continuation_platform_audit.v2",
        "status": "ALL_PASS",
        "paired_gate": True,
        "inference_authorized": True,
        "cells": [{
            "identity": {"restart_index": 2, "k_model": 3, "fold": 0},
            "source_sha256": source_sha256,
            "checkpoint_repeat_audit_path": str(repeat_path),
            "checkpoint_repeat_audit_sha256": repeat_sha256,
            "checkpoint_bundle": checkpoint_bundle_fingerprint(checkpoint_dir),
            "dense": {"converged": True},
            "exact_repeat": True,
        }],
    }), encoding="utf-8")

    assert _require_continuation_audit(source, audit_path, source_path, checkpoint_dir)
    (checkpoint_dir / "ckpt-8400.index").write_text("changed", encoding="utf-8")
    with pytest.raises(RuntimeError, match="checkpoint hash mismatch"):
        _require_continuation_audit(source, audit_path, source_path, checkpoint_dir)


def test_continuation_gate_rechecks_repeat_artifact(tmp_path: Path) -> None:
    checkpoint_dir = _checkpoint_dir(tmp_path)
    source_path = tmp_path / "cell.json"
    source = {
        "schema": "r04.late_lr_cell.v1",
        "restart_index": 2,
        "k_model": 3,
        "fold": 0,
        "fit_diagnostics": {"start_step": 7200, "steps": 8400},
    }
    source_path.write_text(json.dumps(source), encoding="utf-8")
    repeat_path = tmp_path / "repeat.json"
    repeat_path.write_text("repeat", encoding="utf-8")
    audit_path = tmp_path / "panel.json"
    audit_path.write_text(json.dumps({
        "schema": "r04.continuation_platform_audit.v2",
        "status": "ALL_PASS",
        "paired_gate": True,
        "inference_authorized": True,
        "cells": [{
            "identity": {"restart_index": 2, "k_model": 3, "fold": 0},
            "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
            "checkpoint_repeat_audit_path": str(repeat_path),
            "checkpoint_repeat_audit_sha256": hashlib.sha256(repeat_path.read_bytes()).hexdigest(),
            "checkpoint_bundle": checkpoint_bundle_fingerprint(checkpoint_dir),
            "dense": {"converged": True},
            "exact_repeat": True,
        }],
    }), encoding="utf-8")

    repeat_path.write_text("tampered", encoding="utf-8")
    with pytest.raises(RuntimeError, match="repeat audit hash"):
        _require_continuation_audit(source, audit_path, source_path, checkpoint_dir)
    repeat_path.unlink()
    with pytest.raises(RuntimeError, match="repeat audit is missing"):
        _require_continuation_audit(source, audit_path, source_path, checkpoint_dir)


def test_audited_endpoint_rejects_additional_fit_steps_before_scoring() -> None:
    row = {"checkpoint_bundle": {"step": 8400}}

    _require_frozen_inference_endpoint(row, fit_steps=8400, expected_start_step=8400)
    with pytest.raises(RuntimeError, match="without more fit steps"):
        _require_frozen_inference_endpoint(row, fit_steps=9000, expected_start_step=8400)


def test_audited_endpoint_rejects_inference_protocol_drift() -> None:
    row = {"checkpoint_bundle": {"step": 8400}}
    protocol = dict(EXPECTED_FROZEN_INFERENCE_PROTOCOL)
    protocol["lengthscale"] = 4.0

    with pytest.raises(RuntimeError, match="protocol does not match"):
        _require_frozen_inference_endpoint(
            row,
            fit_steps=8400,
            expected_start_step=8400,
            inference_protocol=protocol,
        )


def test_audited_endpoint_rejects_inference_learning_rate_drift() -> None:
    row = {"checkpoint_bundle": {"step": 8400}}
    protocol = dict(EXPECTED_FROZEN_INFERENCE_PROTOCOL)
    protocol["learning_rate"] = 0.02
    with pytest.raises(RuntimeError, match="protocol does not match"):
        _require_frozen_inference_endpoint(
            row,
            fit_steps=8400,
            expected_start_step=8400,
            inference_protocol=protocol,
        )


def test_v2_declared_continuation_also_requires_panel_audit(tmp_path: Path) -> None:
    source = {
        "schema": "r04.late_lr_cell.v2",
        "continuation_panel": True,
        "fit_diagnostics": {"start_step": 8400, "steps": 9000},
    }

    with pytest.raises(RuntimeError, match="require a panel audit"):
        _require_continuation_audit(source, None, tmp_path / "cell.json", None)


def test_practical_acceptance_allows_only_the_named_stopped_panel_cell(
    tmp_path: Path,
) -> None:
    checkpoint_dir = _checkpoint_dir(tmp_path)
    source = {
        "schema": "r04.late_lr_cell.v1",
        "restart_index": 3,
        "k_model": 3,
        "fold": 0,
        "fit_diagnostics": {"start_step": 7200, "steps": 8400},
    }
    source_path = tmp_path / "cell.json"
    source_path.write_text(json.dumps(source), encoding="utf-8")
    repeat_path = tmp_path / "repeat.json"
    repeat_path.write_text(json.dumps({"repeat": "artifact"}), encoding="utf-8")

    cells = []
    for restart in (2, 3, 4):
        for factors in (0, 3):
            for fold in (0, 4):
                cell_id = f"r{restart}_k{factors}_f{fold}"
                exceptional = cell_id == "r3_k3_f0"
                row = {
                    "cell_id": cell_id,
                    "identity": {
                        "restart_index": restart,
                        "k_model": factors,
                        "fold": fold,
                    },
                    "exact_repeat": True,
                    "dense": {
                        "converged": not exceptional,
                        "relative_improvement_threshold": 0.001,
                        "relative_improvements": (
                            [0.0010065261463332008, -0.000307587704458387]
                            if exceptional else [0.0002, -0.0001]
                        ),
                        "objective_means": (
                            [5307.040283203125, 5301.6986083984375, 5303.329345703125]
                            if exceptional else [10.0, 9.998, 9.999]
                        ),
                        "steps": [7600, 8000, 8400],
                    },
                }
                if exceptional:
                    row.update({
                        "source_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
                        "checkpoint_repeat_audit_path": str(repeat_path),
                        "checkpoint_repeat_audit_sha256": hashlib.sha256(repeat_path.read_bytes()).hexdigest(),
                        "checkpoint_bundle": checkpoint_bundle_fingerprint(checkpoint_dir),
                    })
                cells.append(row)
    audit = {
        "schema": "r04.continuation_platform_audit.v2",
        "status": "STOPPED_OFF_PLATFORM",
        "cell_count": 12,
        "passed_cell_count": 11,
        "paired_gate": False,
        "inference_authorized": False,
        "artifact_errors": {},
        "provenance_errors": {},
        "failing_cells": {"r3_k3_f0": ["dense_objective_off_platform"]},
        "panel_source_hash": "panel-source",
        "cells": cells,
    }
    audit_path = tmp_path / "panel.json"
    audit_path.write_text(json.dumps(audit), encoding="utf-8")
    acceptance = {
        "schema": "r04.practical_inference_acceptance.v1",
        "status": "PRACTICAL_PLATFORM_ACCEPTED_WITH_ONE_MARGINAL_EXCEPTION",
        "decision_id": "D-084",
        "strict_audit_sha256": hashlib.sha256(audit_path.read_bytes()).hexdigest(),
        "strict_audit_status": "STOPPED_OFF_PLATFORM",
        "panel_source_hash": "panel-source",
        "accepted_exception": {
            "cell_id": "r3_k3_f0",
            "identity": {"restart_index": 3, "k_model": 3, "fold": 0},
            "reason": "MARGINAL_FIRST_WINDOW_OVERRUN",
            "threshold": 0.001,
            "first_window": 0.0010065261463332008,
            "net_7600_8400": 0.0006992480369416418,
            "final_window_abs": 0.000307587704458387,
            "exact_repeat": True,
        },
        "scope": {
            "checkpoint_step": 8400,
            "allowed": ["held_out_inference", "held_out_nb_scoring"],
            "fit_updates_allowed": False,
        },
        "selected_k": None,
    }
    acceptance_path = tmp_path / "acceptance.json"
    acceptance_path.write_text(json.dumps(acceptance), encoding="utf-8")

    row = _require_continuation_audit(
        source,
        audit_path,
        source_path,
        checkpoint_dir,
        practical_acceptance_path=acceptance_path,
    )
    assert row is not None
    assert row["inference_authority"] == "PRACTICAL_POSTHOC_ACCEPTANCE"
    assert row["strict_audit_status"] == "STOPPED_OFF_PLATFORM"
