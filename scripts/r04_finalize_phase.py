#!/usr/bin/env python3
"""R-04 repaired gate (D-140): separate specificity from a common residual field.

The registered composition artifact contains specific-versus-shared AUC deltas,
not an absolute common-residual readout. Candidate-specific verdicts are retained,
but these inputs cannot establish the absence of every residual field. The
scientific status therefore remains R04_BLOCKED_IDENTIFIABILITY while the
operational stop and closed K search are preserved; no next node is authorized.

A specific candidate needs finite, ordered uncertainty bounds, an adjusted
increment at least EFFECT_FLOOR, and replication in two distinct patients across
both data folds. Shared patients cannot replicate with themselves. These
spot-bootstrap intervals are exploratory and do not become patient-level
confirmatory intervals through this gate. Unnamed loading-subspace recurrence
alone is not evidence of a residual biological field. No model training occurs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from r04.runtime import atomic_json

FOLDS = (0, 4)
RANKS = ("rank1", "rank2")
EFFECT_FLOOR = 0.02
OUTER_FILES = {
    0: "infra/r04/explore_outer_validation_fold0_20260905.json",
    4: "infra/r04/explore_outer_validation_fold4_20260905.json",
}


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _read_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _check_registered(path: Path, record: dict, label: str) -> dict:
    if not path.is_file():
        raise ValueError(f"missing evidence artifact: {label}")
    if _sha256(path) != record.get("sha256"):
        raise ValueError(f"hash mismatch: {label}")
    return _read_json(path)


def _short(patient: str) -> str:
    return patient.split("::")[-1][:8]


def _candidate_verdict(hits: list[tuple[str, float]]) -> str:
    """hits: (patient, delta) with delta>=EFFECT_FLOOR and q025>0.

    SURVIVES needs two disjoint patients with the same sign; tested but
    unreplicated candidates do not survive; untested ones are unjudgeable."""
    if not hits:
        return "DOES_NOT_SURVIVE"
    by_sign: dict[bool, set[str]] = {True: set(), False: set()}
    for patient, delta in hits:
        by_sign[delta > 0].add(patient)
    if any(len(patients) >= 2 for patients in by_sign.values()):
        return "SURVIVES"
    return "DOES_NOT_SURVIVE"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--k-closure", type=Path,
                        default=Path("infra/r04/k_closure_canonical_20260911.json"))
    parser.add_argument("--composition-audit", type=Path,
                        default=Path("infra/r04/composition_final_audit_20260911.json"))
    parser.add_argument("--unnamed-audit", type=Path,
                        default=Path("infra/r04/unnamed_field_audit_20260911.json"))
    parser.add_argument("--output", type=Path,
                        default=Path("infra/repair_20260921/r04_final_gate.json"))
    args = parser.parse_args()
    root = args.project_root.resolve()
    evidence: list[dict[str, object]] = []

    def register(relative: str, label: str) -> dict:
        path = root / relative
        if not path.is_file():
            raise ValueError(f"missing evidence artifact: {label}")
        record = {"path": relative, "sha256": _sha256(path)}
        evidence.append({**record, "role": label})
        return _read_json(path)

    closure = register(str(args.k_closure), "k_closure_canonical")
    if closure.get("status") != "K_CLOSED_WORKING_K3_RANK2_READOUT":
        raise ValueError("K closure is not in the closed state")
    comp = register(str(args.composition_audit), "composition_final_audit_nested")
    if comp.get("status") != "NESTED_AUDIT_COMPLETE":
        raise ValueError("composition audit incomplete")
    unnamed = register(str(args.unnamed_audit), "unnamed_field_audit")
    outer = {}
    for fold in FOLDS:
        outer[fold] = register(OUTER_FILES[fold], f"outer_validation_fold{fold}")

    comp_folds = {f["fold"]: f for f in comp["folds"]}
    if sorted(comp_folds) != [0, 4]:
        raise ValueError("composition audit must cover folds 0 and 4")
    train_patients: dict[int, set[str]] = {}
    for fold in FOLDS:
        patients = set()
        for rank_rec in comp_folds[fold]["ranks"]:
            for person in rank_rec.get("paired_patients", []):
                patients.add(str(person))
        train_patients[fold] = patients
    patient_overlap = sorted(train_patients[0] & train_patients[4])
    for fold in FOLDS:
        heldout: set[str] = set()
        for entry in outer[fold].get("entries", []):
            if isinstance(entry, dict):
                heldout.update(str(p) for p in entry.get("heldout_patients", []))
        if train_patients[fold] & heldout:
            raise ValueError(f"fold{fold} training/heldout patients overlap")

    candidates: list[dict[str, object]] = []
    tested_anywhere = False
    invalid_uncertainty = False
    for structure in ("TLS", "TUMOR_STROMA_BOUNDARY"):
        for rank in RANKS:
            hits: list[tuple[str, float]] = []
            hit_folds = set()
            candidate_invalid = False
            computed = 0
            detail: dict[str, object] = {}
            for fold in FOLDS:
                rank_rec = next(r for r in comp_folds[fold]["ranks"] if r["rank"] == rank)
                folds_hit = []
                for inner in rank_rec["inner_folds"]:
                    if inner.get("status") != "COMPUTED":
                        continue
                    adj = inner.get("adjusted", {})
                    ci = inner.get("adjusted_bootstrap_ci", {}).get("q025_q975", {}).get(structure)
                    delta = adj.get("auc_delta", {}).get(structure)
                    valid = (isinstance(delta, (int, float)) and math.isfinite(delta)
                             and isinstance(ci, (list, tuple)) and len(ci) == 2
                             and all(isinstance(v, (int, float)) and math.isfinite(v) for v in ci)
                             and ci[0] <= ci[1])
                    if not valid:
                        invalid_uncertainty = True
                        candidate_invalid = True
                        continue
                    computed += 1
                    hit = delta >= EFFECT_FLOOR and ci[0] > 0
                    entry = {
                        "test_patients": [_short(p) for p in inner["test_patients"]],
                        "adjusted_auc_delta": delta,
                        "bootstrap_q025_q975": ci,
                        "stable_residual": bool(hit),
                    }
                    folds_hit.append(entry)
                    if hit:
                        hit_folds.add(fold)
                        for person in inner["test_patients"]:
                            hits.append((str(person), delta))
                detail[str(fold)] = folds_hit
            if computed:
                tested_anywhere = True
            verdict = (_candidate_verdict(hits) if computed and not candidate_invalid else "NOT_IDENTIFIABLE")
            if verdict == "SURVIVES" and len(hit_folds) < 2:
                verdict = "DOES_NOT_SURVIVE"
            candidates.append({
                "candidate": f"{structure}_{rank}_residual",
                "verdict": verdict,
                "n_computed_inner_folds": computed,
                "detail": detail,
            })

    unnamed_verdicts = {c["candidate"]: c["verdict"] for c in unnamed["candidates"]}
    rank3_reproduced = any(
        c["verdict"] == "SURVIVES"
        for c in candidates if c["candidate"].endswith("_rank2_residual"))
    isolated_0544 = {
        "fold4_TSB_rank3_0.544": ("REPRODUCED" if rank3_reproduced
                                  else "ISOLATED_EXPLORATORY_SIGNAL"),
        "note": ("rank1/rank2, fold-0, and nested-adjusted do not reproduce it; "
                 "third direction not pursued (D-109)"),
    }
    tls_sections = {}
    for fold in FOLDS:
        tls_sections[str(fold)] = outer[fold].get("tls_note", outer[fold].get("status"))
    surviving = [c["candidate"] for c in candidates if c["verdict"] == "SURVIVES"]
    not_identifiable = [c["candidate"] for c in candidates
                        if c["verdict"] == "NOT_IDENTIFIABLE"]
    uncertainty_ok = bool(tested_anywhere) and not invalid_uncertainty and all(
        inner.get("status") != "COMPUTED" or "adjusted_bootstrap_ci" in inner
        for fold in FOLDS
        for rank_rec in comp_folds[fold]["ranks"]
        for inner in rank_rec["inner_folds"])
    # Specific-versus-shared contrast cannot establish absence of a shared field.
    # These registered artifacts contain no absolute common-residual readout.
    specificity_status = ("SURVIVES" if surviving and uncertainty_ok else
                          "NOT_IDENTIFIABLE" if not_identifiable or not uncertainty_ok else
                          "NO_REPLICATED_SPECIFIC_INCREMENT")
    status = "R04_BLOCKED_IDENTIFIABILITY"
    next_nodes = []
    blocked_nodes = ["R-05:NOT_TRIGGERED_SHARED_FIELD_UNTESTED", "R-06", "R-07"]
    gate = {
        "schema": "r04.final_gate.v2",
        "phase": "R04",
        "status": status,
        "created_at": "2026-09-21",
        "operational_status": "CLOSED_NO_NEW_COMPUTE_AUTHORIZED",
        "specific_residual_status": specificity_status,
        "shared_residual_field_status": "NOT_TESTED",
        "claim_boundary": "A null specific-versus-shared AUC contrast does not reject a common residual field; D-140.",
        "working_k_model": 3,
        "primary_readout_rank": 2,
        "global_k_eff": "NOT_IDENTIFIABLE",
        "selected_k": None,
        "k_search_closed_for_r04": True,
        "surviving_candidates": surviving,
        "composition_gate": {
            "method": "nested patient-level crossfit residualizer (D-108)",
            "effect_floor": EFFECT_FLOOR,
            "effect_floor_provenance": ("project null-equivalence band; "
                                          "D-110 synthesis"),
            "training_patient_overlap_folds_0_4": [_short(p) for p in patient_overlap],
            "overlap_note": ("a29 trains in both folds: it cannot replicate "
                               "with itself; replication needs disjoint patients"),
            "candidate_verdicts": [
                {"candidate": c["candidate"], "verdict": c["verdict"]}
                for c in candidates],
        },
        "unnamed_field_reproducibility_gate": {
            "rank1_subspace": unnamed_verdicts.get("rank1_subspace"),
            "rank2_subspace": unnamed_verdicts.get("rank2_subspace"),
            "rank3": "DESCRIPTIVE_ONLY",
            "note": ("a surviving loading subspace alone is not a surviving "
                     "residual field"),
        },
        "anchored_readout_summary": {
            "outer_validation": "null finding both folds (TSB); TLS one section "
                                "per fold, descriptive only",
            "isolated_signal": isolated_0544,
        },
        "uncertainty_gate": "PASS" if uncertainty_ok else "FAIL",
        "coordinate_null": "NOT_TRIGGERED",
        "independent_lineage": "NOT_TRIGGERED_NO_SURVIVING_CANDIDATE",
        "next_nodes_authorized": next_nodes,
        "blocked_or_not_triggered_nodes": blocked_nodes,
        "evidence": evidence,
    }
    out = root / args.output
    out.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(out, gate)
    print(gate["status"], "| surviving:", surviving, "| uncertain:", not_identifiable)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
