#!/usr/bin/env python3
"""Run one restart's four frozen held-out inference cells sequentially."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any

from r04.runtime import atomic_json, resource_status
from r04.continuation_platform_audit import checkpoint_bundle_fingerprint
from scripts.r04_resume_infer import EXPECTED_FROZEN_INFERENCE_PROTOCOL


def _cell_specs(
    restart_index: int, output_root: Path, source_root: Path
) -> list[dict[str, Any]]:
    return [
        {
            "cell_id": f"r{restart_index}_k{factors}_f{fold}",
            "factors": factors,
            "fold": fold,
            "fit_audit": source_root / f"k{factors}_fold{fold}" / "cell.json",
            "checkpoint_dir": (
                source_root
                / f"k{factors}_fold{fold}"
                / "checkpoints"
                / f"k{factors}"
                / f"fold{fold}"
            ),
            "output_dir": output_root / f"k{factors}_fold{fold}",
        }
        for fold in (0, 4)
        for factors in (0, 3)
    ]


def _resume_command(
    *,
    python: str,
    repo: Path,
    restart_index: int,
    spec: dict[str, Any],
    panel_audit: Path,
    acceptance: Path,
) -> list[str]:
    return [
        python,
        str(repo / "scripts/r04_resume_infer.py"),
        "--manifest-json", str(repo / "infra/r04/role_manifests/training_manifest.json"),
        "--gene-list", str(repo / "infra/r04/gene_universe.txt"),
        "--source-checkpoint-dir", str(spec["checkpoint_dir"]),
        "--fit-audit-json", str(spec["fit_audit"]),
        "--output-dir", str(spec["output_dir"]),
        "--factors", str(spec["factors"]),
        "--fold", str(spec["fold"]),
        "--folds", "5",
        "--group-seed", "20260807",
        "--restart-index", str(restart_index),
        "--gene-folds", "2",
        "--split-seed", "20260817",
        "--fit-steps", "8400",
        "--expected-start-step", "8400",
        "--inference-steps", "2400",
        "--checkpoint-steps", "600",
        "--inducing-points", "16",
        "--lengthscale", "3.0",
        "--learning-rate", "0.01",
        "--gene-batch-size", "512",
        "--diagnostic-interval", "100",
        "--continuation-panel-audit", str(panel_audit),
        "--practical-acceptance", str(acceptance),
        "--allow-source-environment-hash",
    ]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _expected_provenance(
    spec: dict[str, Any],
    *,
    restart_index: int,
    panel_audit: Path,
    acceptance: Path,
) -> dict[str, Any]:
    return {
        "identity": {
            "restart_index": restart_index,
            "k_model": spec["factors"],
            "fold": spec["fold"],
        },
        "source_sha256": _sha256(Path(spec["fit_audit"])),
        "checkpoint_bundle_sha256": checkpoint_bundle_fingerprint(
            spec["checkpoint_dir"], step=8400
        )["sha256"],
        "strict_audit_sha256": _sha256(panel_audit),
        "practical_acceptance_sha256": _sha256(acceptance),
        "protocol": dict(EXPECTED_FROZEN_INFERENCE_PROTOCOL),
    }


def _completed_output(path: Path, expected_provenance: dict[str, Any]) -> bool:
    summary = path / "r04_resume_infer.json"
    cell_path = path / "cell.json"
    if not summary.is_file() or not cell_path.is_file():
        return False
    try:
        result = json.loads(summary.read_text(encoding="utf-8"))
        cell = json.loads(cell_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    platform = cell.get("inference_platform")
    return (
        result.get("status") == "FROZEN_HELDOUT_DIAGNOSTIC_NOT_K_SELECTION"
        and result.get("inference_authority") == "PRACTICAL_POSTHOC_ACCEPTANCE"
        and result.get("strict_audit_status") == "STOPPED_OFF_PLATFORM"
        and result.get("selected_k") is None
        and result.get("frozen_inference_provenance") == expected_provenance
        and isinstance(platform, dict)
        and platform.get("split_count") == 2
    )


def _child_environment(repo: Path, *, use_gpu: bool) -> dict[str, str]:
    """Build the child environment without changing the inference protocol."""
    env = dict(os.environ)
    if not use_gpu:
        env["CUDA_VISIBLE_DEVICES"] = "-1"
    env["PYTHONPATH"] = str(repo)
    return env


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--restart-index", type=int, choices=(2, 3, 4), required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--panel-audit", type=Path, required=True)
    parser.add_argument("--practical-acceptance", type=Path, required=True)
    parser.add_argument(
        "--use-gpu",
        action="store_true",
        help="leave GPU visibility unchanged; default is the historical CPU-only mode",
    )
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[1]
    args.output_root.mkdir(parents=True, exist_ok=True)
    progress_path = args.output_root / "panel_progress.json"
    specs = _cell_specs(args.restart_index, args.output_root, args.source_root)
    progress: dict[str, Any] = {
        "schema": "r04.frozen_inference_panel_progress.v1",
        "status": "RUNNING",
        "restart_index": args.restart_index,
        "strict_audit_status": "STOPPED_OFF_PLATFORM",
        "inference_authority": "PRACTICAL_POSTHOC_ACCEPTANCE",
        "selected_k": None,
        "completed_cells": [],
        "failed_cells": {},
    }
    atomic_json(progress_path, progress)
    env = _child_environment(repo, use_gpu=args.use_gpu)

    for spec in specs:
        cell_id = str(spec["cell_id"])
        for required in (spec["fit_audit"], spec["checkpoint_dir"], args.panel_audit, args.practical_acceptance):
            if not Path(required).exists():
                progress["failed_cells"][cell_id] = f"MISSING_INPUT:{required}"
                atomic_json(progress_path, progress)
                break
        else:
            expected_provenance = _expected_provenance(
                spec,
                restart_index=args.restart_index,
                panel_audit=args.panel_audit,
                acceptance=args.practical_acceptance,
            )
            if _completed_output(spec["output_dir"], expected_provenance):
                progress["completed_cells"].append(cell_id)
                atomic_json(progress_path, progress)
                continue
            final_files = [
                spec["output_dir"] / "r04_resume_infer.json",
                spec["output_dir"] / "cell.json",
            ]
            if any(path.exists() for path in final_files):
                progress["failed_cells"][cell_id] = "STALE_OUTPUT_PROVENANCE_MISMATCH"
                atomic_json(progress_path, progress)
                continue
            if resource_status(args.output_root) == "BLOCKED_STORAGE":
                progress["failed_cells"][cell_id] = "BLOCKED_STORAGE"
                progress["status"] = "INFERENCE_INCOMPLETE"
                atomic_json(progress_path, progress)
                return 2
            command = _resume_command(
                python=sys.executable,
                repo=repo,
                restart_index=args.restart_index,
                spec=spec,
                panel_audit=args.panel_audit,
                acceptance=args.practical_acceptance,
            )
            completed = subprocess.run(command, cwd=repo, env=env, check=False)
            if completed.returncode == 0 and _completed_output(
                spec["output_dir"], expected_provenance
            ):
                progress["completed_cells"].append(cell_id)
            else:
                progress["failed_cells"][cell_id] = f"RETURN_CODE:{completed.returncode}"
            atomic_json(progress_path, progress)

    progress["status"] = (
        "ALL_CELLS_COMPLETE"
        if len(progress["completed_cells"]) == 4 and not progress["failed_cells"]
        else "INFERENCE_INCOMPLETE"
    )
    atomic_json(progress_path, progress)
    return 0 if progress["status"] == "ALL_CELLS_COMPLETE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
