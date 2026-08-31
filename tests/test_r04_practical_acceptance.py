from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from r04.practical_acceptance import validate_practical_inference_acceptance


def _strict_audit() -> dict[str, object]:
    cells = []
    for restart in (2, 3, 4):
        for factors in (0, 3):
            for fold in (0, 4):
                cell_id = f"r{restart}_k{factors}_f{fold}"
                exceptional = cell_id == "r3_k3_f0"
                cells.append({
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
                })
    return {
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


def _acceptance(audit_path: Path, audit: dict[str, object]) -> dict[str, object]:
    return {
        "schema": "r04.practical_inference_acceptance.v1",
        "status": "PRACTICAL_PLATFORM_ACCEPTED_WITH_ONE_MARGINAL_EXCEPTION",
        "decision_id": "D-084",
        "strict_audit_sha256": hashlib.sha256(audit_path.read_bytes()).hexdigest(),
        "strict_audit_status": "STOPPED_OFF_PLATFORM",
        "panel_source_hash": audit["panel_source_hash"],
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


def _write_audit(tmp_path: Path) -> tuple[Path, dict[str, object]]:
    audit = _strict_audit()
    path = tmp_path / "panel.json"
    path.write_text(json.dumps(audit), encoding="utf-8")
    return path, audit


def test_current_marginal_exception_authorizes_only_frozen_inference(tmp_path: Path) -> None:
    audit_path, audit = _write_audit(tmp_path)
    result = validate_practical_inference_acceptance(
        _acceptance(audit_path, audit), audit, audit_path
    )
    assert result["accepted_cell_id"] == "r3_k3_f0"
    assert result["checkpoint_step"] == 8400
    assert result["inference_authority"] == "PRACTICAL_POSTHOC_ACCEPTANCE"


@pytest.mark.parametrize("mutation", [
    lambda acceptance: acceptance.update(selected_k=3),
    lambda acceptance: acceptance["scope"].update(fit_updates_allowed=True),
    lambda acceptance: acceptance["accepted_exception"].update(cell_id="r2_k3_f0"),
    lambda acceptance: acceptance["accepted_exception"].update(first_window=0.0011),
])
def test_acceptance_scope_or_exception_tamper_fails_closed(
    tmp_path: Path, mutation
) -> None:
    audit_path, audit = _write_audit(tmp_path)
    acceptance = _acceptance(audit_path, audit)
    mutation(acceptance)
    with pytest.raises(ValueError):
        validate_practical_inference_acceptance(acceptance, audit, audit_path)


def test_acceptance_is_bound_to_exact_strict_audit_bytes(tmp_path: Path) -> None:
    audit_path, audit = _write_audit(tmp_path)
    acceptance = _acceptance(audit_path, audit)
    audit_path.write_text(audit_path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="strict audit hash"):
        validate_practical_inference_acceptance(acceptance, audit, audit_path)


def test_strict_audit_remains_stopped_and_unique_failure_is_required(tmp_path: Path) -> None:
    audit_path, audit = _write_audit(tmp_path)
    acceptance = _acceptance(audit_path, audit)
    audit["status"] = "ALL_PASS"
    audit_path.write_text(json.dumps(audit), encoding="utf-8")
    acceptance["strict_audit_sha256"] = hashlib.sha256(audit_path.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="strict audit state"):
        validate_practical_inference_acceptance(acceptance, audit, audit_path)
