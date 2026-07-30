from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

from scripts.r02_build_registry import build_registry
from scripts.r02_validate_gate import evaluate_gate


ROOT = Path(__file__).resolve().parents[1]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


@pytest.fixture()
def built(tmp_path: Path) -> Path:
    output = tmp_path / "structure-registry"
    build_registry(ROOT, output)
    return output


def test_builds_metadata_only_fail_closed_control_plane(built: Path) -> None:
    ontology = read_tsv(built / "structure_ontology.tsv")
    candidates = read_tsv(built / "tls_candidate_summary.tsv")
    instances = read_tsv(built / "structure_instances.tsv")
    audits = read_tsv(built / "gt_source_audit.tsv")
    leakage = read_tsv(built / "leakage_selection_audit.tsv")
    input_policy = read_tsv(built / "input_policy.tsv")
    splits = read_tsv(built / "outer_splits.tsv")

    assert {row["structure_id"] for row in ontology} == {
        "TLS",
        "BLOOD_VESSEL",
        "NECROSIS",
        "TUMOR_STROMA_BOUNDARY",
    }
    assert len(candidates) == 57
    assert len({row["source_tls_id"] for row in candidates}) == 57
    assert "P10_8" in {row["source_tls_id"] for row in candidates}
    assert all(row["source_record_id"].startswith("Table S4!") for row in candidates)
    assert {row["gt_source_id"] for row in candidates} == {"ATLAS_TABLE_S4_TLS_IDS"}
    breakdown: dict[str, int] = {}
    for row in candidates:
        breakdown[row["logical_unit_id"]] = breakdown.get(row["logical_unit_id"], 0) + 1
    assert breakdown == {
        "HTAN_VANDERBILT_CRC": 44,
        "GEO::GSE226997": 4,
        "GEO::GSE274103": 8,
        "GEO::GSE274557": 1,
    }
    assert len(instances) == 57
    assert {row["confirmation_status"] for row in instances} == {
        "NONCONFIRMATORY_SOURCE_REPORTED_ID"
    }
    assert not any(row["gt_geometry_locator"] for row in instances)
    assert sum(row["source_class"] == "STOMICS_TLS_ANNOTATION_CANDIDATE" for row in audits) == 19
    assert {"ATLAS_TABLE_S4_TLS_IDS", "UPSTREAM_CRC_TLS_ID_SUMMARY"} <= {
        row["gt_source_id"] for row in audits
    }
    assert not any(row["audit_status"] == "AUDITABLE_GT" for row in audits)
    assert {row["risk_type"] for row in leakage} == {
        "OUTCOME_IN_SAMPLE_TITLE",
        "PREPROCESSING_SELECTION_BIAS",
    }
    assert (
        "identity, site, treatment, stage, sample and file metadata",
        "FORBIDDEN_AS_PREDICTIVE_INPUT",
    ) in {(row["channel_class"], row["policy"]) for row in input_policy}
    assert len(splits) == 121
    assert len({row["identity_envelope_id"] for row in splits}) == 61
    assert len({row["block_group_id"] for row in splits if row["block_group_id"]}) == 32
    assert all(
        not row["block_id"]
        for row in splits
        if row["identity_granularity"] == "patient_linked_physical_specimen"
    )

    gate = evaluate_gate(ROOT, built)
    assert gate["status"] == "HARD_BLOCKED_NO_AUDITABLE_GT"
    assert gate["confirmatory_instance_count"] == 0
    assert gate["second_structure_frozen"] is False
    assert gate["tls_candidate_count"] == 57
    assert gate["stomics_candidate_source_count"] == 19
    assert gate["outer_patient_group_count"] == 61
    assert gate["explicit_block_group_count"] == 32


def test_false_green_instance_self_report_is_rejected(built: Path) -> None:
    path = built / "structure_instances.tsv"
    rows = read_tsv(path)
    rows[0]["confirmation_status"] = "CONFIRMATORY"
    rows[0]["gt_source_id"] = "forged"
    rows[0]["gt_geometry_locator"] = "row=1"
    _write_tsv(path, rows)

    gate = evaluate_gate(ROOT, built)
    assert gate["status"] == "HARD_BLOCKED_GT_INTEGRITY"
    assert any("unverified_gt_source" in item for item in gate["integrity_errors"])
    assert gate["confirmatory_instance_count"] == 0


@pytest.mark.parametrize(
    ("field", "value", "needle"),
    [
        ("geometry_locator", "", "missing_geometry"),
        ("provenance_status", "UNKNOWN", "unknown_provenance"),
        ("overlap_status", "UNKNOWN", "unknown_overlap"),
        ("same_assay_status", "SAME_ASSAY", "circular_same_assay_gt"),
    ],
)
def test_claimed_gt_source_fails_closed(
    built: Path, field: str, value: str, needle: str
) -> None:
    path = built / "gt_source_audit.tsv"
    rows = read_tsv(path)
    rows[0].update(
        {
            "physical_link_status": "VERIFIED",
            "geometry_locator": "Table S2!I4:J4",
            "provenance_status": "KNOWN",
            "overlap_status": "KNOWN_NONOVERLAP",
            "same_assay_status": "INDEPENDENT_ASSAY",
            "audit_status": "AUDITABLE_GT",
        }
    )
    rows[0][field] = value
    _write_tsv(path, rows)

    gate = evaluate_gate(ROOT, built)
    assert gate["status"] == "HARD_BLOCKED_GT_INTEGRITY"
    assert any(needle in item for item in gate["integrity_errors"])


def test_fully_self_reported_gt_has_no_supported_verifier(built: Path) -> None:
    audit_path = built / "gt_source_audit.tsv"
    audits = read_tsv(audit_path)
    physical_id = read_tsv(built / "outer_splits.tsv")[0]["physical_unit_id"]
    source = audits[0]
    source.update(
        {
            "physical_unit_id": physical_id,
            "physical_link_status": "VERIFIED",
            "geometry_locator": source["path"] + "#sheet=Table S2;range=I4:J4",
            "provenance_status": "KNOWN",
            "provenance_locator": source["path"] + "#sheet=Table S2",
            "overlap_status": "KNOWN_NONOVERLAP",
            "overlap_locator": source["path"] + "#declared-nonoverlap",
            "same_assay_status": "INDEPENDENT_ASSAY",
            "audit_status": "AUDITABLE_GT",
        }
    )
    _write_tsv(audit_path, audits)
    policy_path = built / "r02_policy.json"
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    policy["supported_confirmatory_gt_verifiers"] = [source["source_class"]]
    policy_path.write_text(
        json.dumps(policy, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    gate = evaluate_gate(ROOT, built)
    assert gate["status"] == "HARD_BLOCKED_GT_INTEGRITY"
    assert any("unsupported_gt_verifier" in item for item in gate["integrity_errors"])
    assert any(
        "policy_declares_unimplemented_gt_verifier" in item
        for item in gate["integrity_errors"]
    )
    assert gate["auditable_gt_source_count"] == 0


def test_specimen_locator_cannot_be_promoted_to_block(built: Path) -> None:
    path = built / "outer_splits.tsv"
    rows = read_tsv(path)
    specimen = next(
        row for row in rows if row["identity_granularity"] == "patient_linked_physical_specimen"
    )
    specimen["block_id"] = specimen["physical_specimen_id"]
    _write_tsv(path, rows)

    gate = evaluate_gate(ROOT, built)
    assert gate["status"] == "HARD_BLOCKED_SPLIT_INTEGRITY"
    assert any("specimen_as_block" in item for item in gate["integrity_errors"])


def test_identity_envelope_cannot_cross_fold(built: Path) -> None:
    path = built / "outer_splits.tsv"
    rows = read_tsv(path)
    patient_counts = Counter(row["patient_id"] for row in rows if row["block_id"])
    patient = next(key for key, count in patient_counts.items() if count > 1)
    target = next(row for row in rows if row["patient_id"] == patient)
    target["outer_fold"] = "forged-fold"
    _write_tsv(path, rows)

    gate = evaluate_gate(ROOT, built)
    assert gate["status"] == "HARD_BLOCKED_SPLIT_INTEGRITY"
    assert any("cross_fold_identity_envelope" in item for item in gate["integrity_errors"])


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("patient_id", "FORGED_PATIENT"),
        ("block_id", "FORGED_BLOCK"),
        ("relation_type", "CONFIRMED_STRUCTURE_CORRESPONDENCE"),
        ("evidence_status", "CONFIRMED"),
    ],
)
def test_cross_section_link_self_report_is_rejected(
    built: Path, field: str, value: str
) -> None:
    path = built / "cross_section_links.tsv"
    rows = read_tsv(path)
    rows[0][field] = value
    _write_tsv(path, rows)

    gate = evaluate_gate(ROOT, built)
    assert gate["status"] == "HARD_BLOCKED_SPLIT_INTEGRITY"
    assert any("cross_section_link_mismatch" in item for item in gate["integrity_errors"])


def test_table_s4_source_id_cannot_be_replaced_by_count_ordinal(built: Path) -> None:
    path = built / "tls_candidate_summary.tsv"
    rows = read_tsv(path)
    rows[0]["source_tls_id"] = "SYNTHETIC_ORDINAL_1"
    _write_tsv(path, rows)

    gate = evaluate_gate(ROOT, built)
    assert gate["status"] == "HARD_BLOCKED_GT_INTEGRITY"
    assert "table_s4_candidate_id_set_mismatch" in gate["integrity_errors"]


def test_build_is_byte_deterministic_and_stored_gate_recomputes(tmp_path: Path) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    build_registry(ROOT, first)
    build_registry(ROOT, second)

    first_files = sorted(path.name for path in first.iterdir() if path.is_file())
    second_files = sorted(path.name for path in second.iterdir() if path.is_file())
    assert first_files == second_files
    assert all((first / name).read_bytes() == (second / name).read_bytes() for name in first_files)

    stored = json.loads((first / "r02_gate.json").read_text(encoding="utf-8"))
    assert stored == evaluate_gate(ROOT, first)


def test_build_cli_runs_without_pythonpath_override(tmp_path: Path) -> None:
    output = tmp_path / "cli-registry"
    completed = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/r02_build_registry.py"),
            "--root",
            str(ROOT),
            "--output",
            str(output),
        ],
        cwd=ROOT,
        env={"PATH": str(Path(sys.executable).parent)},
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    assert json.loads((output / "r02_gate.json").read_text(encoding="utf-8"))[
        "status"
    ] == "HARD_BLOCKED_NO_AUDITABLE_GT"


def _write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    assert rows
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
