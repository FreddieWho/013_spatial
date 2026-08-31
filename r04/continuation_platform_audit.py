"""Pure, read-only audit and aggregation for the R-04 continuation panel."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable

from .diagnostics import deterministic_objective_platform_summary

EXPECTED_CELLS = {(r, k, f) for r in (2, 3, 4) for k in (0, 3) for f in (0, 4)}
EXPECTED_REPEAT_PROTOCOL = {
    "folds": 5,
    "group_seed": 20260807,
    "inducing_points": 16,
    "lengthscale": 3.0,
    "ridge": 1e-5,
}
EXPECTED_MANIFEST_PROTOCOL = {
    **EXPECTED_REPEAT_PROTOCOL,
    "mc_draws": 4,
    "repeat_seed_offset": 7000,
}


def _canonical(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def checkpoint_bundle_fingerprint(
    checkpoint_dir: str | Path, *, step: int = 8400,
) -> dict[str, Any]:
    """Hash the state file, metadata, index, and data shards for one checkpoint."""
    root = Path(checkpoint_dir)
    required = [
        root / "checkpoint",
        root / "checkpoint_metadata.json",
        root / f"ckpt-{step}.index",
    ]
    shards = sorted(root.glob(f"ckpt-{step}.data-*"))
    if any(not path.is_file() for path in required) or not shards:
        raise ValueError(f"checkpoint bundle is incomplete: {root}/ckpt-{step}")
    state = (root / "checkpoint").read_text(encoding="utf-8").splitlines()
    if not state or state[0].strip() != f'model_checkpoint_path: "ckpt-{step}"':
        raise ValueError(f"checkpoint state does not select ckpt-{step}: {root}")
    files = [
        {"name": path.name, "sha256": _sha256(path)}
        for path in [*required, *shards]
    ]
    bundle_sha256 = hashlib.sha256(_canonical(files).encode("utf-8")).hexdigest()
    return {"step": step, "sha256": bundle_sha256, "files": files}


def _read_ref(ref: Any, base_dir: Path) -> tuple[dict[str, Any], str, str]:
    if isinstance(ref, str):
        path, declared = base_dir / ref, None
    elif isinstance(ref, dict) and ("path" in ref or "file" in ref):
        path = base_dir / str(ref.get("path", ref.get("file")))
        declared = ref.get("sha256")
    else:
        raise ValueError("audit references must point to JSON files")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    actual = _sha256(path)
    if declared is not None and declared != actual:
        raise ValueError(f"declared hash mismatch: {path}")
    return value, actual, str(path)


def _source_nested(source: dict[str, Any], key: str) -> Any:
    diagnostics = source.get("fit_diagnostics")
    if not isinstance(diagnostics, dict) or key not in diagnostics:
        raise ValueError(f"source cell is missing fit_diagnostics.{key}")
    return diagnostics[key]


def _validated_source_trace(source: dict[str, Any], mc_draws: int) -> list[dict[str, Any]]:
    trace = _source_nested(source, "evaluation_trace")
    if not isinstance(trace, list) or len(trace) < 3 or any(not isinstance(item, dict) for item in trace):
        raise ValueError("source evaluation trace is invalid")
    recent = trace[-3:]
    if [item.get("step") for item in recent] != [7600, 8000, 8400]:
        raise ValueError("source evaluation trace has unexpected final steps")
    if any(
        item.get("objective_kind") != "dense_full_panel_fixed_stateless_mc"
        or item.get("mc_draws") != mc_draws
        for item in recent
    ):
        raise ValueError("source evaluation trace protocol mismatch")
    return trace


def _checkpoint_entries(audit: dict[str, Any]) -> list[dict[str, Any]]:
    if audit.get("schema") != "r04.checkpoint_audit.v1":
        raise ValueError("checkpoint audit schema mismatch")
    checkpoints = audit.get("checkpoints")
    if not isinstance(checkpoints, list) or len(checkpoints) != 2:
        raise ValueError("checkpoint audit must contain exactly two checkpoints")
    if any(not isinstance(item, dict) for item in checkpoints):
        raise ValueError("checkpoint audit contains an invalid checkpoint")
    return checkpoints


def _is_finite_number(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def _validate_checkpoint_payload(audit: dict[str, Any]) -> None:
    mc_draws = audit.get("mc_draws")
    if not isinstance(mc_draws, int) or isinstance(mc_draws, bool) or mc_draws < 1:
        raise ValueError("checkpoint audit mc_draws is invalid")
    for item in _checkpoint_entries(audit):
        draws = item.get("objective_draws")
        if (
            not isinstance(draws, list)
            or len(draws) != mc_draws
            or any(not _is_finite_number(value) for value in draws)
        ):
            raise ValueError("checkpoint objective_draws are invalid")
        objective_mean = item.get("objective_mean")
        objective_sd = item.get("objective_sd")
        if not _is_finite_number(objective_mean) or not _is_finite_number(objective_sd):
            raise ValueError("checkpoint objective summary is invalid")
        expected_mean = sum(float(value) for value in draws) / len(draws)
        expected_sd = (
            math.sqrt(
                sum((float(value) - expected_mean) ** 2 for value in draws)
                / (len(draws) - 1)
            )
            if len(draws) > 1
            else 0.0
        )
        if not math.isclose(float(objective_mean), expected_mean, rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError("checkpoint objective_mean is inconsistent with draws")
        if not math.isclose(float(objective_sd), expected_sd, rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError("checkpoint objective_sd is inconsistent with draws")
        for field in ("global_step", "mc_draws", "seed", "n_genes", "n_spots"):
            value = item.get(field)
            if not isinstance(value, int) or isinstance(value, bool):
                raise ValueError(f"checkpoint {field} is invalid")
        if item["mc_draws"] != mc_draws or item["n_genes"] < 1 or item["n_spots"] < 1:
            raise ValueError("checkpoint dimensions or draw count are invalid")
        if not isinstance(item.get("checkpoint"), str) or not item["checkpoint"]:
            raise ValueError("checkpoint path is invalid")
        gauge = item.get("nuisance_gauge")
        if not isinstance(gauge, dict):
            raise ValueError("checkpoint nuisance_gauge is invalid")
        rank = gauge.get("rank")
        if not isinstance(rank, int) or isinstance(rank, bool) or rank < 0:
            raise ValueError("checkpoint nuisance rank is invalid")
        for field in ("h_min", "h_mean", "h_max", "baseline_v_projection"):
            values = gauge.get(field)
            if (
                not isinstance(values, list)
                or len(values) != rank
                or any(not _is_finite_number(value) for value in values)
            ):
                raise ValueError(f"checkpoint nuisance {field} is invalid")


def _frozen_checkpoint(audit: dict[str, Any]) -> dict[str, Any]:
    return _checkpoint_entries(audit)[0]


def _same_path(left: str, right: str, base_dir: Path) -> bool:
    return (base_dir / left).resolve() == (base_dir / right).resolve()


def _repeat_audit(
    audit: dict[str, Any], *, factors: int, fold: int, source_checkpoint_dir: str,
    source_mc_draws: int, optimization_seed: int, source_manifest_hash: str,
    base_dir: Path,
) -> tuple[list[str], list[str]]:
    _validate_checkpoint_payload(audit)
    checkpoints = _checkpoint_entries(audit)
    provenance_errors: list[str] = []
    repeat_errors: list[str] = []
    if audit.get("status") != "DIAGNOSTIC_ONLY":
        provenance_errors.append("checkpoint_audit_status_mismatch")
    if audit.get("factors") != factors or audit.get("fold") != fold:
        provenance_errors.append("checkpoint_audit_identity_mismatch")
    if audit.get("mc_draws") != source_mc_draws:
        provenance_errors.append("checkpoint_audit_mc_draws_mismatch")
    if audit.get("manifest_hash") != source_manifest_hash:
        provenance_errors.append("checkpoint_audit_manifest_hash_mismatch")
    for field, expected_value in EXPECTED_REPEAT_PROTOCOL.items():
        if audit.get(field) != expected_value:
            provenance_errors.append(f"checkpoint_audit_{field}_mismatch")
    expected = str(Path(source_checkpoint_dir) / "ckpt-8400")
    for item in checkpoints:
        if item.get("global_step") != 8400:
            provenance_errors.append("checkpoint_step_mismatch")
        if not _same_path(str(item.get("checkpoint", "")), expected, base_dir):
            provenance_errors.append("checkpoint_path_mismatch")
        if item.get("mc_draws") != source_mc_draws:
            provenance_errors.append("repeat_mc_draws_mismatch")
        if item.get("seed") != optimization_seed + 7000:
            provenance_errors.append("repeat_seed_mismatch")
    if checkpoints[0].get("checkpoint") != checkpoints[1].get("checkpoint"):
        provenance_errors.append("repeat_checkpoint_mismatch")
    if checkpoints[0].get("objective_draws") != checkpoints[1].get("objective_draws"):
        repeat_errors.append("repeat_objective_draws_mismatch")
    for field in ("objective_mean", "objective_sd", "n_genes", "n_spots", "nuisance_gauge"):
        if checkpoints[0].get(field) != checkpoints[1].get(field):
            repeat_errors.append(f"repeat_{field}_mismatch")
    if audit.get("gene_count") != checkpoints[0].get("n_genes"):
        provenance_errors.append("checkpoint_audit_gene_count_mismatch")
    return sorted(set(provenance_errors)), sorted(set(repeat_errors))


def replace_final_live_evaluation(
    trace_or_audit: Iterable[dict[str, Any]] | dict[str, Any],
    checkpoint_audit: dict[str, Any] | None = None,
) -> list[dict[str, Any]]:
    """Replace a pre-update final point with the frozen entry from outer audit."""
    if isinstance(trace_or_audit, dict):
        trace: Iterable[dict[str, Any]] = []
        audit = trace_or_audit
    else:
        trace, audit = trace_or_audit, checkpoint_audit
    if not isinstance(audit, dict):
        raise ValueError("checkpoint audit is required")
    frozen = _frozen_checkpoint(audit)
    step = int(frozen["global_step"])
    result = [dict(item) for item in trace if int(item.get("step", -1)) != step]
    result.append({"step": step, "objective_mean": float(frozen["objective_mean"]), "source": "frozen_checkpoint"})
    return result


def _panel_hash(rows: list[dict[str, Any]]) -> str:
    values = [{"identity": row["identity"], "source_sha256": row["source_sha256"]} for row in rows]
    payload = _canonical(sorted(values, key=_canonical)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def audit_continuation_manifest(manifest: dict[str, Any], base_dir: str | Path = ".") -> dict[str, Any]:
    """Validate exactly twelve declared cells without modifying any artifact."""
    root = Path(base_dir)
    entries = manifest.get("cells")
    if not isinstance(entries, list):
        raise ValueError("manifest must contain a cells list")
    artifact_errors: dict[str, list[str]] = {}
    provenance_errors: dict[str, list[str]] = {}
    off_platform: dict[str, list[str]] = {}
    rows: list[dict[str, Any]] = []
    seen: set[tuple[int, int, int]] = set()
    if manifest.get("checkpoint_audit_protocol") != EXPECTED_MANIFEST_PROTOCOL:
        artifact_errors["manifest"] = ["checkpoint_audit_protocol_mismatch"]
    for raw in entries:
        if not isinstance(raw, dict):
            artifact_errors[str(len(artifact_errors))] = ["entry_not_object"]
            continue
        key = (raw.get("restart_index", -1), raw.get("k_model", -1), raw.get("fold", -1))
        cell_id = str(raw.get("cell_id", f"r{key[0]}_k{key[1]}_f{key[2]}"))
        problems: list[str] = []
        if key not in EXPECTED_CELLS:
            problems.append("invalid_identity")
        if key in seen:
            problems.append("duplicate_cell")
        seen.add(key)
        try:
            source, source_hash, _ = _read_ref(raw["source_cell"], root)
            repeat, repeat_hash, repeat_path = _read_ref(raw["checkpoint_repeat_audit"], root)
            if source.get("schema") not in {"r04.late_lr_cell.v1", "r04.late_lr_cell.v2"}:
                problems.append("source_schema_mismatch")
            if source.get("status") != "FIT_ONLY_DIAGNOSTIC":
                problems.append("source_status_mismatch")
            if (source.get("restart_index"), source.get("k_model"), source.get("fold")) != key:
                problems.append("source_identity_mismatch")
            if "input_hash" in raw and raw.get("input_hash") != source.get("input_hash"):
                problems.append("input_hash_mismatch")
            if "training_patients" in raw and raw.get("training_patients") != source.get("training_patients"):
                problems.append("patients_mismatch")
            if "manifest_hash" in raw and raw.get("manifest_hash") != source.get("manifest_hash"):
                problems.append("manifest_hash_mismatch")
            if "optimization_seed" in raw and raw.get("optimization_seed") != source.get("optimization_seed"):
                problems.append("seed_mismatch")
            if _source_nested(source, "steps") != 8400:
                problems.append("fit_steps_mismatch")
            if _source_nested(source, "start_step") != 7200:
                problems.append("fit_start_step_mismatch")
            if "protocol" in raw and raw.get("protocol") != source.get("protocol"):
                problems.append("protocol_mismatch")
            source_mc_draws = int(_source_nested(source, "evaluation_mc_draws"))
            source_trace = _validated_source_trace(source, source_mc_draws)
            checkpoint_root = root / str(source.get("checkpoint_dir", ""))
            checkpoint_metadata = json.loads(
                (checkpoint_root / "checkpoint_metadata.json").read_text(encoding="utf-8")
            )
            if not isinstance(checkpoint_metadata, dict):
                raise ValueError("checkpoint metadata is not an object")
            expected_checkpoint_metadata = {
                "checkpoint_schema": "r04.mnsf_checkpoint.v2",
                "input_hash": source.get("input_hash"),
                "config_hash": source.get("config_hash"),
                "environment_hash": _source_nested(source, "environment_hash"),
            }
            if any(
                checkpoint_metadata.get(field) != expected_value
                for field, expected_value in expected_checkpoint_metadata.items()
            ):
                problems.append("checkpoint_metadata_mismatch")
            repeat_problems, repeat_errors = _repeat_audit(
                repeat, factors=int(key[1]), fold=int(key[2]),
                source_checkpoint_dir=str(source.get("checkpoint_dir", "")),
                source_mc_draws=source_mc_draws,
                optimization_seed=int(source.get("optimization_seed", -1)),
                source_manifest_hash=str(source.get("manifest_hash", "")), base_dir=root,
            )
            problems.extend(repeat_problems)
            checkpoint_bundle = checkpoint_bundle_fingerprint(
                checkpoint_root, step=8400
            )
            dense = deterministic_objective_platform_summary(
                replace_final_live_evaluation(source_trace, repeat)
            )
            rows.append({
                "cell_id": cell_id,
                "identity": {"restart_index": key[0], "k_model": key[1], "fold": key[2]},
                "source_sha256": source_hash,
                "checkpoint_repeat_audit_sha256": repeat_hash,
                "checkpoint_repeat_audit_path": str(Path(repeat_path).resolve()),
                "checkpoint_bundle": checkpoint_bundle,
                "input_hash": source.get("input_hash"),
                "training_patients": source.get("training_patients"),
                "manifest_hash": source.get("manifest_hash"),
                "optimization_seed": source.get("optimization_seed"),
                "dense": dense,
                "legacy": source.get("fit_platform"), "exact_repeat": not repeat_errors,
            })
            if repeat_errors:
                off_platform[cell_id] = repeat_errors
            if not dense.get("converged", False):
                off_platform.setdefault(cell_id, []).append("dense_objective_off_platform")
            if problems:
                provenance_errors[cell_id] = sorted(set(problems))
        except (KeyError, TypeError, ValueError, OSError, json.JSONDecodeError) as exc:
            artifact_errors[cell_id] = [f"invalid_artifact:{type(exc).__name__}"]
    missing = EXPECTED_CELLS - seen
    if len(entries) != 12:
        artifact_errors.setdefault("manifest", []).append("exactly_12_entries_required")
    if missing:
        artifact_errors.setdefault("manifest", []).append("missing_cells")
    for restart in (2, 3, 4):
        for fold in (0, 4):
            pair = [
                row for row in rows
                if row["identity"]["restart_index"] == restart
                and row["identity"]["fold"] == fold
            ]
            if len(pair) == 2 and (
                len({row.get("input_hash") for row in pair}) != 1
                or len({_canonical(row.get("training_patients")) for row in pair}) != 1
                or len({row.get("manifest_hash") for row in pair}) != 1
            ):
                for row in pair:
                    provenance_errors.setdefault(row["cell_id"], []).append(
                        "k_pair_provenance_mismatch"
                    )
    panel_source_hash = _panel_hash(rows) if rows else hashlib.sha256(b"[]").hexdigest()
    all_pass = (
        len(entries) == 12 and not missing and not artifact_errors and not provenance_errors
        and not off_platform and len(rows) == 12
        and all(row["dense"].get("converged") is True and row["exact_repeat"] is True for row in rows)
    )
    status = "ALL_PASS" if all_pass else ("INVALID_AUDIT" if artifact_errors or provenance_errors else "STOPPED_OFF_PLATFORM")
    failing: dict[str, list[str]] = {}
    for bucket in (artifact_errors, provenance_errors, off_platform):
        for cell_id, errors in bucket.items():
            failing.setdefault(cell_id, []).extend(errors)
    failing = {cell_id: sorted(set(errors)) for cell_id, errors in failing.items()}
    failed_cell_ids = set(artifact_errors) | set(provenance_errors) | set(off_platform)
    passed_cell_count = sum(row["cell_id"] not in failed_cell_ids for row in rows)
    return {
        "schema": "r04.continuation_platform_audit.v2", "status": status,
        "cell_count": len(entries), "passed_cell_count": passed_cell_count,
        "panel_source_hash": panel_source_hash, "cells": rows,
        "artifact_errors": artifact_errors, "provenance_errors": provenance_errors,
        "off_platform_cells": off_platform, "failing_cells": failing,
        "paired_gate": all_pass, "inference_authorized": all_pass,
    }
