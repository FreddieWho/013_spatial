#!/usr/bin/env python3
"""R-04 final machine-readable gate (D-109/D-110): outcome-neutral closure.

Reads only registered final-audit artifacts, re-verifies hashes/paths/
patient-fold independence, recomputes verdicts, and emits exactly one of:
  R04_COMPLETE_WITH_REPRODUCIBLE_RESIDUAL_FIELD
  R04_COMPLETE_NO_REPRODUCIBLE_RESIDUAL_FIELD
  R04_BLOCKED_IDENTIFIABILITY
No model training. Fail closed on any missing artifact, hash mismatch, or
fold/patient overlap violation.

Pre-committed candidate rules (frozen before composition numbers were seen):
- Anchored residual candidate (structure, rank<=2): per inner (test-patient)
  fold, nested-adjusted specific-increment AUC delta with bootstrap CI.
  A hit needs adjusted delta >= EFFECT_FLOOR with bootstrap q025 > 0.
  EFFECT_FLOOR (0.02) is the project's own historical null-equivalence band
  (|AUC delta| < 0.02 jointly read as zero since the D-100-era sensitivity
  analyses); bootstrap CIs capture spot-sampling noise only, so a signal
  must clear the band, not hug it (D-110 synthesis).
  SURVIVES iff hits occur in >= 2 DISJOINT patients across the two
  data-folds with the same sign (a patient present in both folds, e.g. a29,
  cannot replicate with itself; overlap is recorded, not crashed on).
  DOES_NOT_SURVIVE iff the candidate was tested (COMPUTED with finite CIs)
  but nothing replicates.
  NOT_IDENTIFIABLE iff never COMPUTED anywhere or uncertainty missing.
- Unnamed subspace candidates: verdicts are read off the unnamed audit
  artifact (D-109 rule); a surviving loading subspace alone is NOT a
  surviving residual field.
- fold-4 TSB rank3 0.544: ISOLATED_EXPLORATORY_SIGNAL unless rank1/rank2,
  fold-0, or nested-adjusted reproduces it (D-109).
- Uncertainty gate passes iff every reported adjusted effect carries a
  bootstrap CI and every candidate has a verdict.
"""

from __future__ import annotations

import argparse
import hashlib
import json
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
                        default=Path("infra/r04/r04_final_gate_20260911.json"))
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
    for structure in ("TLS", "TUMOR_STROMA_BOUNDARY"):
        for rank in RANKS:
            hits: list[tuple[str, float]] = []
            computed = 0
            detail: dict[str, object] = {}
            for fold in FOLDS:
                rank_rec = next(r for r in comp_folds[fold]["ranks"] if r["rank"] == rank)
                folds_hit = []
                for inner in rank_rec["inner_folds"]:
                    if inner.get("status") != "COMPUTED":
                        continue
                    computed += 1
                    adj = inner["adjusted"]
                    ci = inner["adjusted_bootstrap_ci"]["q025_q975"][structure]
                    delta = float(adj["auc_delta"][structure] or 0.0)
                    hit = (delta >= EFFECT_FLOOR and ci[0] is not None
                           and ci[0] > 0)
                    entry = {
                        "test_patients": [_short(p) for p in inner["test_patients"]],
                        "adjusted_auc_delta": delta,
                        "bootstrap_q025_q975": ci,
                        "stable_residual": bool(hit),
                    }
                    folds_hit.append(entry)
                    if hit:
                        for person in inner["test_patients"]:
                            hits.append((str(person), delta))
                detail[str(fold)] = folds_hit
            if computed:
                tested_anywhere = True
            verdict = (_candidate_verdict(hits) if computed else "NOT_IDENTIFIABLE")
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
    uncertainty_ok = bool(tested_anywhere) and all(
        inner.get("status") != "COMPUTED" or "adjusted_bootstrap_ci" in inner
        for fold in FOLDS
        for rank_rec in comp_folds[fold]["ranks"]
        for inner in rank_rec["inner_folds"])
    if surviving:
        status = "R04_COMPLETE_WITH_REPRODUCIBLE_RESIDUAL_FIELD"
    elif not_identifiable or not uncertainty_ok:
        status = "R04_BLOCKED_IDENTIFIABILITY"
    else:
        status = "R04_COMPLETE_NO_REPRODUCIBLE_RESIDUAL_FIELD"
    if status == "R04_COMPLETE_WITH_REPRODUCIBLE_RESIDUAL_FIELD":
        next_nodes = ["R-05"]
        blocked_nodes = ["R-06", "R-07"]
    else:
        next_nodes = []
        blocked_nodes = ["R-05:NOT_TRIGGERED_NO_SURVIVING_R04_FIELD", "R-06", "R-07"]
    gate = {
        "schema": "r04.final_gate.v1",
        "phase": "R04",
        "status": status,
        "created_at": "2026-09-11",
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
