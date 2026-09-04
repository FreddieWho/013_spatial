from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from r04.structure_gt import labels_for_spots, load_training_gt_catalog


def _write_inputs(tmp_path: Path, *, duplicate: bool = False) -> tuple[Path, Path]:
    root = tmp_path / "project"
    registry_dir = root / "infra" / "structure-registry"
    gt_dir = root / "data" / "gt"
    registry_dir.mkdir(parents=True)
    gt_dir.mkdir(parents=True)
    (gt_dir / "a.csv").write_text(
        "Barcode,pathology_annotation\n"
        "a1,lymphoid_follicle\n"
        "a2,carcinoma\n"
        "a3,normal_mucosa\n",
        encoding="utf-8",
    )
    (gt_dir / "b.csv").write_text(
        "Barcode,pathology_annotation\n"
        "b1,carcinoma_border\n"
        "b2,carcinoma\n",
        encoding="utf-8",
    )
    a_locator = "data/gt/a.csv#pathology_annotation=lymphoid_follicle"
    b_locator = "data/gt/b.csv#pathology_annotation=carcinoma_border,carcinoma_edge"
    rows = [
        ["TLS", "section-a", a_locator, "CONFIRMATORY", "training-fold distance fields"],
        ["TUMOR_STROMA_BOUNDARY", "section-b", b_locator, "CONFIRMATORY", "training-fold localization benchmarks"],
        ["TLS", "section-c", a_locator, "NONCONFIRMATORY_CONTEXT_PRENEOPLASTIC", "discovery and context reference only"],
    ]
    if duplicate:
        rows.append(["TLS", "section-a", a_locator, "CONFIRMATORY", "training-fold distance fields"])
    (registry_dir / "structure_instances.tsv").write_text(
        "structure_id\tphysical_unit_id\tgt_geometry_locator\tconfirmation_status\tallowed_use\n"
        + "\n".join("\t".join(row) for row in rows) + "\n",
        encoding="utf-8",
    )
    manifest = root / "infra" / "r04" / "role_manifests" / "training_manifest.json"
    manifest.parent.mkdir(parents=True)
    manifest.write_text(json.dumps({
        "schema": "r04.role_manifest.v1",
        "status": "READY",
        "role": "training",
        "rows": [{"section_id": value} for value in ("section-a", "section-b", "section-c")],
    }), encoding="utf-8")
    return registry_dir / "structure_instances.tsv", manifest


def test_training_gt_keeps_missing_spots_unknown_and_aligns_exact_barcodes(tmp_path: Path) -> None:
    registry, manifest = _write_inputs(tmp_path)
    catalog = load_training_gt_catalog(registry, manifest)

    labels = labels_for_spots(
        catalog,
        "TLS",
        ["section-a", "section-a", "section-a", "section-b"],
        ["a1", "a2", "missing", "b1"],
    )
    assert labels[:2].tolist() == [1.0, 0.0]
    assert np.isnan(labels[2:]).all()
    boundary = labels_for_spots(
        catalog,
        "TUMOR_STROMA_BOUNDARY",
        ["section-b", "section-b", "section-b"],
        ["b1", "b2", "missing"],
    )
    assert boundary[:2].tolist() == [1.0, 0.0]
    assert np.isnan(boundary[2])
    assert catalog["structures"]["TLS"]["section_count"] == 1


def test_training_gt_rejects_duplicate_source_barcodes(tmp_path: Path) -> None:
    registry, manifest = _write_inputs(tmp_path)
    gt = tmp_path / "project" / "data" / "gt" / "a.csv"
    gt.write_text(
        "Barcode,pathology_annotation\na1,lymphoid_follicle\na1,carcinoma\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="duplicate training GT barcode"):
        load_training_gt_catalog(registry, manifest)


def test_training_gt_requires_training_role_manifest(tmp_path: Path) -> None:
    registry, manifest = _write_inputs(tmp_path)
    value = json.loads(manifest.read_text())
    value["role"] = "internal_validation"
    manifest.write_text(json.dumps(value), encoding="utf-8")
    with pytest.raises(ValueError, match="READY training manifest"):
        load_training_gt_catalog(registry, manifest)
