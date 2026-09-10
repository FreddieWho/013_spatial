#!/usr/bin/env python3
"""Assemble the current R-04 K/readout/resource evidence handoff.

This does not recompute scores or alter historical artifacts.  It validates
the upstream schemas, copies their machine-readable summaries, and records
the hashes and the current conditional K=2 bridge status in one small handoff.
Validation GT remains represented only by the training-role readout status;
internal/external validation GT is never loaded here.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Mapping

from r04.runtime import atomic_json


def _read_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source(path: Path, *, project_root: Path, value: Mapping[str, object]) -> dict[str, object]:
    return {
        "path": str(path.resolve().relative_to(project_root.resolve())),
        "sha256": _sha256(path),
        "schema": value.get("schema"),
        "status": value.get("status"),
    }


def _inner_summary(value: object) -> dict[str, object]:
    if not isinstance(value, dict):
        return {"status": "NOT_AVAILABLE"}
    return {
        "status": value.get("status"),
        "valid_inner_folds": value.get("valid_inner_folds"),
        "total_inner_folds": value.get("total_inner_folds"),
        "specific_minus_shared_mean_mse_delta": value.get(
            "specific_minus_shared_mean_mse_delta"
        ),
        "mean_auc_delta_by_target": value.get("mean_auc_delta_by_target"),
    }


def build_final_artifact(
    aggregate: Mapping[str, object],
    structure_readout: Mapping[str, object],
    k2_stage: Mapping[str, object],
    k2_manifest: Mapping[str, object],
    *,
    source_artifacts: list[dict[str, object]],
    created_at: str = "2026-09-01",
) -> dict[str, object]:
    """Build the current evidence handoff without inventing a K result."""
    if aggregate.get("schema") != "r04.k_semantics_and_downstream_robustness.v1":
        raise ValueError("unsupported aggregate semantics schema")
    if structure_readout.get("schema") != "r04.structure_readout_panel.v1":
        raise ValueError("unsupported structure readout schema")
    if k2_stage.get("schema") != "r04.k2_bridge_stage_a.v1":
        raise ValueError("unsupported K=2 stage schema")
    if k2_manifest.get("schema") != "r04.k2_bridge_stage_a_manifest.v1":
        raise ValueError("unsupported K=2 manifest schema")
    if k2_manifest.get("status") != "FROZEN_PRECOMPUTE_MANIFEST":
        raise ValueError("K=2 manifest is not frozen")
    if k2_stage.get("status") != "BLOCKED_CPU_RUNTIME":
        raise ValueError("this handoff expects the current CPU runtime block")

    evidence = aggregate.get("available_evidence")
    if not isinstance(evidence, dict):
        raise ValueError("aggregate lacks available evidence")
    score_summary = evidence.get("score_summary")
    subspace_summary = evidence.get("subspace_summary_by_fold")
    entries = structure_readout.get("entries")
    if not isinstance(score_summary, dict) or not isinstance(subspace_summary, dict):
        raise ValueError("aggregate lacks score/subspace summaries")
    if not isinstance(entries, list) or len(entries) != 10:
        raise ValueError("structure readout must contain ten K=3 training entries")

    readout_summary = []
    for entry in entries:
        if not isinstance(entry, dict):
            raise ValueError("structure readout entry is not an object")
        coverage = entry.get("coverage_training")
        paired = coverage.get("paired") if isinstance(coverage, dict) else None
        paired_patients = paired.get("patients", []) if isinstance(paired, dict) else []
        readout_summary.append({
            "restart_index": entry.get("restart_index"),
            "fold": entry.get("fold"),
            "paired_patient_count": len(paired_patients),
            "paired_patients": list(paired_patients),
            "full_k3": _inner_summary(
                entry.get("full_k3", {}).get("inner_crossfit")
                if isinstance(entry.get("full_k3"), dict)
                else None
            ),
            "stable_rank2_sensitivity": _inner_summary(
                entry.get("stable_rank2_sensitivity", {}).get("inner_crossfit")
                if isinstance(entry.get("stable_rank2_sensitivity"), dict)
                else None
            ),
        })

    source = copy.deepcopy(source_artifacts)
    return {
        "schema": "r04.k_semantics_and_downstream_robustness.v2",
        "status": "K_SEMANTICS_CORRECTED_K2_BRIDGE_BLOCKED_CPU_RUNTIME",
        "created_at": created_at,
        "source_artifacts": source,
        "k_interpretation": "latent_spatial_effect_dimension_not_structure_count",
        "k_model": copy.deepcopy(aggregate.get("k_model")),
        "k_eff": copy.deepcopy(aggregate.get("k_eff")),
        "selected_k": None,
        "score_summary": copy.deepcopy(score_summary),
        "subspace_stability_by_fold": copy.deepcopy(subspace_summary),
        "third_direction_interpretation": (
            "descriptive instability signal; not a biological structure count and "
            "not sufficient to assign K_eff"
        ),
        "downstream_k_robust": structure_readout.get("downstream_k_robust"),
        "k2_bridge_required": structure_readout.get("k2_bridge_required"),
        "structure_readout": {
            "status": structure_readout.get("status"),
            "scope": structure_readout.get("scope"),
            "downstream_k_robust": structure_readout.get("downstream_k_robust"),
            "k2_bridge_required": structure_readout.get("k2_bridge_required"),
            "reason": structure_readout.get("k2_bridge_reason"),
            "training_role_entry_count": len(entries),
            "per_entry_summary": readout_summary,
            "internal_external_validation_gt": structure_readout.get(
                "internal_external_validation_gt"
            ),
            "factor_naming": "forbidden",
        },
        "k2_bridge": {
            "required": structure_readout.get("k2_bridge_required"),
            "manifest_status": k2_manifest.get("status"),
            "manifest_entry_count": len(k2_manifest.get("entries", [])),
            "stage_status": k2_stage.get("status"),
            "blocked_reason": k2_stage.get("blocked_reason"),
            "attempted_cell": k2_stage.get("attempted_cell"),
            "scientific_result": "NOT_AVAILABLE_NO_K2_CELL_COMPLETED",
            "gpu_required_for_bounded_completion": k2_stage.get(
                "gpu_required_for_bounded_completion"
            ),
        },
        "spatial_null_status": "NOT_RUN",
        "composition_split_status": "NOT_RUN",
        "structure_naming_status": "FORBIDDEN_NOT_RUN",
        "independent_lineage_reproduction_status": "NOT_RUN",
        "claim_boundary": [
            "K_model is representation capacity, not a count of biological structures",
            "K_eff and selected_k remain unresolved",
            "the readout is training-role descriptive evidence only and does not establish structure identity",
            "the K=2 bridge has no scientific result because its first CPU cell stopped before a checkpoint",
            "spatial nulls, composition controls and independent lineage replication remain pending",
            "internal and external validation GT remains sealed",
        ],
    }


# Canonical K closure (D-107). Recorded D-106 paired means, recomputed
# fail-closed from the evidence cells. Tolerances bind the record; they are
# not a discovery rule.
K_CLOSURE_FROZEN_STEPS = 8400
K_CLOSURE_INPUT_HASH_PREFIX = "044d28c5"
K_CLOSURE_OBJECTIVE_HASH_PREFIX = "b5a3fb68"
K_CLOSURE_MANIFEST_HASH_PREFIX = "5a0dc1de"
K_CLOSURE_RECORDED_MEANS = {
    ("k02", 0): -16.61,
    ("k02", 4): 71.96,
    ("k3", 0): -28.96,
    ("k3", 4): 72.73,
}
K_CLOSURE_RECORD_TOL = 0.005


def _require_finite_scores(cell: Mapping[str, object], label: str) -> dict[str, float]:
    import math
    scores = cell.get("patient_scores")
    if not isinstance(scores, dict) or len(scores) != 6:
        raise ValueError(f"{label} must carry six patient scores")
    out = {}
    for patient, value in scores.items():
        number = float(value)
        if not math.isfinite(number):
            raise ValueError(f"{label} has non-finite score for {patient}")
        out[str(patient)] = number
    return out


def _check_evidence_cell(cell: Mapping[str, object], *, fold: int, k: int,
                         label: str, input_prefix: str = K_CLOSURE_INPUT_HASH_PREFIX,
                         objective_prefix: str = K_CLOSURE_OBJECTIVE_HASH_PREFIX
                         ) -> dict[str, float]:
    if cell.get("status") != "FIT_AND_SCORED":
        raise ValueError(f"{label} is not FIT_AND_SCORED")
    if cell.get("k_model") != k or cell.get("fold") != fold:
        raise ValueError(f"{label} k/fold mismatch")
    fit = cell.get("fit_diagnostics")
    if not isinstance(fit, dict) or fit.get("steps") != K_CLOSURE_FROZEN_STEPS:
        raise ValueError(f"{label} fit steps != {K_CLOSURE_FROZEN_STEPS}")
    platform = cell.get("inference_platform")
    if not isinstance(platform, dict) or platform.get("converged") is not True:
        raise ValueError(f"{label} inference not converged")
    if list(platform.get("split_converged", [])) != [True, True]:
        raise ValueError(f"{label} inference splits not both converged")
    if not str(cell.get("input_hash", "")).startswith(input_prefix):
        raise ValueError(f"{label} input hash anchor mismatch")
    objective_hash = ""
    if isinstance(fit, dict):
        objective_hash = str(fit.get("objective_input_hash", ""))
    if not objective_hash.startswith(objective_prefix):
        raise ValueError(f"{label} objective-input hash anchor mismatch")
    return _require_finite_scores(cell, label), str(cell.get("input_hash", "")), objective_hash


def _check_evidence_top(top: Mapping[str, object], label: str) -> None:
    if top.get("status") != "WIRING_ONLY_NOT_SCIENTIFIC":
        raise ValueError(f"{label} top status mismatch")
    if not str(top.get("manifest_hash", "")).startswith(K_CLOSURE_MANIFEST_HASH_PREFIX):
        raise ValueError(f"{label} manifest hash anchor mismatch")
    if top.get("selected_k") is not None:
        raise ValueError(f"{label} selected_k must stay null")


def _check_checkpoint_audit(audit: Mapping[str, object], *, fold: int,
                            factors: int, label: str,
                            bind_in: set[str] | None = None,
                            bind_obj: set[str] | None = None) -> None:
    if audit.get("fold") != fold or audit.get("factors") != factors:
        raise ValueError(f"{label} fold/factors mismatch")
    checkpoints = audit.get("checkpoints")
    if not isinstance(checkpoints, list) or len(checkpoints) != 1:
        raise ValueError(f"{label} must carry exactly one checkpoint")
    entry = checkpoints[0]
    if not isinstance(entry, dict):
        raise ValueError(f"{label} checkpoint entry malformed")
    if entry.get("global_step") != K_CLOSURE_FROZEN_STEPS:
        raise ValueError(f"{label} checkpoint step mismatch")
    if entry.get("objective_input_hash_status") != "VERIFIED":
        raise ValueError(f"{label} checkpoint hash not VERIFIED")
    if bind_in is not None:
        if entry.get("input_hash", "") not in bind_in:
            raise ValueError(f"{label} checkpoint not bound to fold cells")
    elif not str(entry.get("input_hash", "")).startswith(K_CLOSURE_INPUT_HASH_PREFIX):
        raise ValueError(f"{label} checkpoint input hash mismatch")
    if bind_obj is not None:
        if entry.get("objective_input_hash", "") not in bind_obj:
            raise ValueError(f"{label} checkpoint not bound to fold cells")
    elif not str(entry.get("objective_input_hash", "")).startswith(
            K_CLOSURE_OBJECTIVE_HASH_PREFIX):
        raise ValueError(f"{label} checkpoint objective hash mismatch")


def build_k_closure(
    runs: Mapping[str, Mapping[int, Mapping[str, object]]],
    audits: list[Mapping[str, object]],
    *,
    source_artifacts: list[dict[str, object]],
    created_at: str = "2026-09-11",
) -> dict[str, object]:
    """Recompute the D-106 K closure fail-closed from evidence artifacts.

    ``runs`` maps ``{"k02": {fold: {"top", "cells": {k: cell}}}, ...}``.
    Raises on any protocol/hash/status deviation. The recorded D-106 means
    are binding targets, not fitted parameters.
    """
    import math
    scores: dict[tuple[str, int], dict[str, float]] = {}
    # Fold 0 carries the preflight-verified anchors. Fold 4 is a different
    # split (different input hash by construction); its binding is
    # cross-consistency: all four fold-4 cells share one input/objective
    # hash, and both fold-4 checkpoint audits (VERIFIED recomputations)
    # pin the same hashes. Recorded here after verification, not assumed.
    fold4_in: set[str] = set()
    fold4_obj: set[str] = set()
    for tag in ("k02", "k3"):
        for fold in (0, 4):
            run = runs[tag][fold]
            _check_evidence_top(run["top"], f"{tag}/fold{fold} top")
            high_k = 2 if tag == "k02" else 3
            k0, k0_in, k0_obj = _check_evidence_cell(
                run["cells"][0], fold=fold, k=0, label=f"{tag}/fold{fold}/k0",
                **({} if fold == 0 else {"input_prefix": "", "objective_prefix": ""}))
            khi, khi_in, khi_obj = _check_evidence_cell(
                run["cells"][high_k], fold=fold, k=high_k,
                label=f"{tag}/fold{fold}/k{high_k}",
                **({} if fold == 0 else {"input_prefix": "", "objective_prefix": ""}))
            if sorted(k0) != sorted(khi):
                raise ValueError(f"{tag}/fold{fold} patient sets differ")
            if fold == 4:
                fold4_in.update([k0_in, khi_in])
                fold4_obj.update([k0_obj, khi_obj])
            scores[(tag, fold)] = {p: khi[p] - k0[p] for p in sorted(k0)}
    if len(audits) != 4:
        raise ValueError("closure needs exactly four checkpoint audits")
    if len(fold4_in) != 1 or len(fold4_obj) != 1:
        raise ValueError("fold-4 cells disagree on input/objective hash")
    for audit, (fold, factors) in zip(audits, [(0, 0), (0, 2), (4, 0), (4, 2)]):
        _check_checkpoint_audit(audit, fold=fold, factors=factors,
                                label=f"audit k{factors} fold{fold}",
                                bind_in=fold4_in if fold == 4 else None,
                                bind_obj=fold4_obj if fold == 4 else None)
    means = {key: sum(d.values()) / len(d) for key, d in scores.items()}
    for key, recorded in K_CLOSURE_RECORDED_MEANS.items():
        if abs(means[key] - recorded) > K_CLOSURE_RECORD_TOL:
            raise ValueError(f"{key} mean {means[key]:.4f} drifts from D-106 {recorded}")
    for fold in (0, 4):
        m02 = means[("k02", fold)]
        m3 = means[("k3", fold)]
        if (m02 > 0) != (m3 > 0):
            raise ValueError(f"fold{fold} K2/K3 mean direction disagree")
    agree4 = sum(1 for p in scores[("k02", 4)]
                 if (scores[("k02", 4)][p] > 0) == (scores[("k3", 4)][p] > 0))
    if agree4 < 5:
        raise ValueError("fold-4 per-patient K2/K3 direction agreement < 5/6")
    k2_minus_k3_fold0 = means[("k02", 0)] - means[("k3", 0)]
    if not k2_minus_k3_fold0 > 0:
        raise ValueError("fold-0 K2 must remain less negative than K3 per D-107")
    return {
        "schema": "r04.k_closure_canonical.v1",
        "status": "K_CLOSED_WORKING_K3_RANK2_READOUT",
        "created_at": created_at,
        "source_artifacts": copy.deepcopy(source_artifacts),
        "working_k_model": 3,
        "primary_readout_rank": 2,
        "global_k_eff": "NOT_IDENTIFIABLE",
        "selected_k": None,
        "k_search_closed_for_r04": True,
        "third_direction": "NO_STABLE_POSITIVE_INCREMENT",
        "paired_means": {f"{tag}_minus_k0_fold{fold}": means[(tag, fold)]
                           for tag, fold in means},
        "paired_deltas_by_patient": {
            f"{tag}_minus_k0_fold{fold}": deltas
            for (tag, fold), deltas in scores.items()
        },
        "fold4_direction_agreement": f"{agree4}/6",
        "fold0_k2_minus_k3_mean": k2_minus_k3_fold0,
        "claim_boundary": [
            "K=3 is an overcomplete working container, not a Keff=2 proof",
            "the third direction shows no stable positive held-out increment",
            "fold heterogeneity lives within effective rank <= 2",
            "K=2 evidence is single-restart Torch; cross-seed third-direction "
            "instability comes from the TF five-seed panel",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["legacy-handoff", "closure"],
                        default="legacy-handoff")
    parser.add_argument("--aggregate", type=Path, required=False)
    parser.add_argument("--structure-readout", type=Path, required=False)
    parser.add_argument("--k2-stage", type=Path, required=False)
    parser.add_argument("--k2-manifest", type=Path, required=False)
    parser.add_argument("--run-dirs", nargs="*", type=Path, default=[],
                        help="closure mode: <fold0_k02> <fold4_k02> <fold0_k3> <fold4_k3>")
    parser.add_argument("--audits", nargs="*", type=Path, default=[],
                        help="closure mode: four checkpoint-audit JSONs")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    project_root = Path(".").resolve()
    if args.mode == "closure":
        if len(args.run_dirs) != 4 or len(args.audits) != 4:
            raise SystemExit("closure mode needs 4 --run-dirs and 4 --audits")
        run_keys = [("k02", 0), ("k02", 4), ("k3", 0), ("k3", 4)]
        runs: dict[str, dict[int, dict[str, object]]] = {"k02": {}, "k3": {}}
        evidence_paths: list[Path] = []
        evidence_values: list[Mapping[str, object]] = []
        for (tag, fold), raw in zip(run_keys, args.run_dirs):
            root = raw.resolve()
            top = _read_json(root / "r04_real_k_search.json")
            cells = {}
            for k in ((0, 2) if tag == "k02" else (0, 3)):
                cell_path = root / f"cells/k{k}/fold{fold}/cell.json"
                cells[k] = _read_json(cell_path)
                evidence_paths.append(cell_path)
                evidence_values.append(cells[k])
            evidence_paths.append(root / "r04_real_k_search.json")
            evidence_values.append(top)
            runs[tag][fold] = {"top": top, "cells": cells}
        audits = []
        for raw in args.audits:
            path = raw.resolve()
            value = _read_json(path)
            audits.append(value)
            evidence_paths.append(path)
            evidence_values.append(value)
        result = build_k_closure(
            runs, audits,
            source_artifacts=[
                _source(p, project_root=project_root, value=v)
                for p, v in zip(evidence_paths, evidence_values)
            ],
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        atomic_json(args.output, result)
        return 0
    for name in ("aggregate", "structure_readout", "k2_stage", "k2_manifest"):
        if getattr(args, name) is None:
            raise SystemExit(f"legacy-handoff mode needs --{name.replace('_', '-')}")
    paths = [
        args.aggregate.resolve(),
        args.structure_readout.resolve(),
        args.k2_stage.resolve(),
        args.k2_manifest.resolve(),
    ]
    values = [_read_json(path) for path in paths]
    result = build_final_artifact(
        values[0], values[1], values[2], values[3],
        source_artifacts=[
            _source(path, project_root=project_root, value=value)
            for path, value in zip(paths, values)
        ],
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(args.output, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
