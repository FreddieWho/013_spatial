from __future__ import annotations

import json

import numpy as np
from scipy import sparse

from r04.diagnostics import continuous_factor_diagnostics
from r04.loaders import load_section_from_row
from r04.synthetic import make_overlapping_sections
from scripts.r04_fit import _align_sections, _manifest_roles
from scripts.r04_partition_manifests import partition_manifest


def test_role_partition_is_pure_and_hash_linked() -> None:
    rows = [
        {"section_id": "T1", "patient_id": "P1", "primary_role": "training"},
        {"section_id": "V1", "patient_id": "P2", "primary_role": "external_validation"},
    ]
    parent = {"status": "READY", "input_manifest_hash": "parent-hash", "rows": rows}
    partitions = partition_manifest(parent)
    assert set(partitions) == {"training", "external_validation"}
    assert partitions["training"]["parent_manifest_hash"] == "parent-hash"
    assert {row["primary_role"] for row in partitions["training"]["rows"]} == {"training"}
    assert _manifest_roles(parent) == ("external_validation", "training")


def test_alignment_can_report_missing_genes_without_zero_fill() -> None:
    sections, _ = make_overlapping_sections(sections=1, spots_per_section=4, genes=4)
    section = sections[0]
    reduced = type(section)(
        section.section_id,
        section.patient_id,
        section.block_id,
        section.lineage,
        section.barcode,
        section.coords,
        sparse.csr_matrix(section.counts.toarray()[:, :3]),
        tuple(section.gene_id[:3]),
        section.library_size,
        section.metadata,
    )
    aligned, coverage = _align_sections([reduced], tuple(section.gene_id), allow_missing=True)
    assert aligned[0].gene_id == tuple(section.gene_id[:3])
    assert aligned[0].counts.shape[1] == 3
    assert coverage[0]["missing_genes"] == [section.gene_id[3]]


def test_factor_diagnostics_expose_duplicate_directions() -> None:
    loading = np.column_stack([np.arange(5), np.arange(5) * 1.0, np.array([1, -1, 1, -1, 1])])
    fields = np.column_stack([np.arange(6), np.arange(6) * 1.0, -np.arange(6)])
    diagnostics = continuous_factor_diagnostics(loading, fields, signed=True)
    assert diagnostics["max_offdiagonal_loading_cosine"] > 0.99
    assert diagnostics["loading_collapse_warning"] is True
    assert diagnostics["effective_rank_participation"] < 3.0


def test_cached_panel_loader_preserves_gene_order_and_coordinates(tmp_path) -> None:
    matrix_path = tmp_path / "panel.npz"
    metadata_path = tmp_path / "panel.json"
    coords_path = tmp_path / "coords.csv"
    sparse.save_npz(matrix_path, sparse.csr_matrix(np.array([[1, 0], [0, 2]], dtype=np.int32)))
    metadata_path.write_text(json.dumps({
        "barcodes": ["BC1", "BC2"],
        "gene_id": ["G2", "G1"],
        "coords": [[4, 5], [6, 7]],
        "library_size": [10, 20],
    }), encoding="utf-8")
    coords_path.write_text("barcode,row,col\nBC1,4,5\nBC2,6,7\n", encoding="utf-8")
    section = load_section_from_row({
        "section_id": "C1",
        "patient_id": "P1",
        "block_id": "B1",
        "lineage": "cached",
        "matrix_locator": str(matrix_path),
        "matrix_kind": "r04_panel_npz",
        "matrix_metadata_locator": str(metadata_path),
        "coordinate_locator": str(coords_path),
    })
    assert section.gene_id == ("G2", "G1")
    np.testing.assert_array_equal(section.coords, np.array([[4, 5], [6, 7]], dtype=float))
    np.testing.assert_array_equal(section.library_size, np.array([10, 20], dtype=float))
