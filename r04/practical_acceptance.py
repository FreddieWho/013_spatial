"""Fail-closed validation for a post-hoc frozen-inference decision.

This module never changes the strict continuation audit.  It only validates a
separate, user-approved decision that permits held-out inference from the
already frozen checkpoints.
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Any


EXPECTED_CELLS = {
    (restart, factors, fold)
    for restart in (2, 3, 4)
    for factors in (0, 3)
    for fold in (0, 4)
}
EXPECTED_EXCEPTION_ID = "r3_k3_f0"
EXPECTED_EXCEPTION_IDENTITY = {"restart_index": 3, "k_model": 3, "fold": 0}
EXPECTED_ALLOWED_ACTIONS = ["held_out_inference", "held_out_nb_scoring"]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _finite_number(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
    )


def _same_number(actual: Any, expected: float) -> bool:
    return _finite_number(actual) and math.isclose(
        float(actual), expected, rel_tol=1e-12, abs_tol=1e-12
    )


def validate_practical_inference_acceptance(
    acceptance: dict[str, Any],
    strict_audit: dict[str, Any],
    strict_audit_path: str | Path,
) -> dict[str, Any]:
    """Validate the one accepted marginal exception and return its scope."""
    audit_path = Path(strict_audit_path)
    if acceptance.get("schema") != "r04.practical_inference_acceptance.v1":
        raise ValueError("practical acceptance schema mismatch")
    if acceptance.get("status") != "PRACTICAL_PLATFORM_ACCEPTED_WITH_ONE_MARGINAL_EXCEPTION":
        raise ValueError("practical acceptance status mismatch")
    if acceptance.get("decision_id") != "D-084":
        raise ValueError("practical acceptance decision mismatch")
    if acceptance.get("strict_audit_sha256") != _sha256(audit_path):
        raise ValueError("strict audit hash mismatch")
    if acceptance.get("strict_audit_status") != "STOPPED_OFF_PLATFORM":
        raise ValueError("practical acceptance rewrites the strict audit status")
    if acceptance.get("selected_k") is not None:
        raise ValueError("practical acceptance cannot select K")

    required_state = {
        "schema": "r04.continuation_platform_audit.v2",
        "status": "STOPPED_OFF_PLATFORM",
        "cell_count": 12,
        "passed_cell_count": 11,
        "paired_gate": False,
        "inference_authorized": False,
        "artifact_errors": {},
        "provenance_errors": {},
        "failing_cells": {EXPECTED_EXCEPTION_ID: ["dense_objective_off_platform"]},
    }
    if any(strict_audit.get(key) != value for key, value in required_state.items()):
        raise ValueError("strict audit state does not match the accepted 11/12 panel")
    if acceptance.get("panel_source_hash") != strict_audit.get("panel_source_hash"):
        raise ValueError("panel source hash mismatch")

    cells = strict_audit.get("cells")
    if not isinstance(cells, list) or len(cells) != 12:
        raise ValueError("strict audit cells are incomplete")
    identities: set[tuple[int, int, int]] = set()
    exception_rows: list[dict[str, Any]] = []
    for row in cells:
        if not isinstance(row, dict) or row.get("exact_repeat") is not True:
            raise ValueError("strict audit cell is not exactly repeatable")
        identity = row.get("identity")
        if not isinstance(identity, dict):
            raise ValueError("strict audit cell identity is invalid")
        try:
            identity_tuple = (
                int(identity["restart_index"]),
                int(identity["k_model"]),
                int(identity["fold"]),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("strict audit cell identity is invalid") from exc
        identities.add(identity_tuple)
        dense = row.get("dense")
        if not isinstance(dense, dict):
            raise ValueError("strict audit dense result is invalid")
        if dense.get("converged") is not True:
            exception_rows.append(row)
    if identities != EXPECTED_CELLS or len(identities) != len(cells):
        raise ValueError("strict audit panel identities do not match the frozen panel")
    if len(exception_rows) != 1 or exception_rows[0].get("cell_id") != EXPECTED_EXCEPTION_ID:
        raise ValueError("strict audit does not contain the unique accepted exception")

    row = exception_rows[0]
    if row.get("identity") != EXPECTED_EXCEPTION_IDENTITY:
        raise ValueError("accepted exception identity mismatch")
    dense = row["dense"]
    improvements = dense.get("relative_improvements")
    objectives = dense.get("objective_means")
    steps = dense.get("steps")
    threshold = dense.get("relative_improvement_threshold")
    if (
        not isinstance(improvements, list)
        or len(improvements) != 2
        or not all(_finite_number(value) for value in improvements)
        or not isinstance(objectives, list)
        or len(objectives) != 3
        or not all(_finite_number(value) for value in objectives)
        or steps != [7600, 8000, 8400]
        or not _finite_number(threshold)
    ):
        raise ValueError("accepted exception diagnostics are invalid")
    threshold_value = float(threshold)
    first_window = float(improvements[0])
    final_window_abs = abs(float(improvements[1]))
    net_change = (float(objectives[0]) - float(objectives[2])) / abs(float(objectives[0]))
    if not (
        threshold_value > 0
        and threshold_value < first_window <= threshold_value * 1.01
        and abs(net_change) <= threshold_value
        and final_window_abs <= threshold_value
    ):
        raise ValueError("strict audit exception is not the accepted marginal case")

    declared = acceptance.get("accepted_exception")
    if not isinstance(declared, dict):
        raise ValueError("accepted exception declaration is missing")
    expected_declaration = {
        "cell_id": EXPECTED_EXCEPTION_ID,
        "identity": EXPECTED_EXCEPTION_IDENTITY,
        "reason": "MARGINAL_FIRST_WINDOW_OVERRUN",
        "threshold": threshold_value,
        "first_window": first_window,
        "net_7600_8400": net_change,
        "final_window_abs": final_window_abs,
        "exact_repeat": True,
    }
    for key, expected in expected_declaration.items():
        actual = declared.get(key)
        if isinstance(expected, float):
            if not _same_number(actual, expected):
                raise ValueError(f"accepted exception {key} mismatch")
        elif actual != expected:
            raise ValueError(f"accepted exception {key} mismatch")

    scope = acceptance.get("scope")
    if not isinstance(scope, dict) or scope != {
        "checkpoint_step": 8400,
        "allowed": EXPECTED_ALLOWED_ACTIONS,
        "fit_updates_allowed": False,
    }:
        raise ValueError("practical acceptance scope mismatch")
    return {
        "accepted_cell_id": EXPECTED_EXCEPTION_ID,
        "accepted_identity": EXPECTED_EXCEPTION_IDENTITY,
        "checkpoint_step": 8400,
        "inference_authority": "PRACTICAL_POSTHOC_ACCEPTANCE",
        "strict_audit_status": "STOPPED_OFF_PLATFORM",
    }
