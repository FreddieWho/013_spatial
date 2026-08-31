from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse

from r04.anchor_eval import evaluate_anchor_correspondence
from r04.candidates import aggregate_candidate_stability, evaluate_candidate_gate, match_model_factors
from r04.composition import crossfit_composition_adjustment, ilr_transform
from r04.io_contract import read_h5ad_counts, reject_target_metadata, section_content_hash, validate_count_matrix
from r04.gene_selection import select_training_genes
from r04.loaders import read_10x_positions
from r04.reference import build_reference_panel, read_allowlisted_h5ad
from r04.runtime import config_fingerprint
from r04.serialization import read_candidate_registry, write_candidate_registry
from scripts.r04_build_input_manifest import build_input_rows
from scripts.r04_fit import _model_checkpoint_dir
from r04.spatial import matern32_kernel, normalize_coordinates
from r04.types import CandidateField, FieldFit, R04ContractError, SectionData


def _fit(factor: str, section: str, loading: np.ndarray, field: np.ndarray, input_hash: str = "input") -> FieldFit:
    return FieldFit("test", factor, section, loading, field, np.ones(len(field)) * 0.1, 2.0, float(np.var(field)), input_hash)


def test_fit_checkpoint_paths_are_model_isolated() -> None:
    base = Path("run/checkpoints")
    assert _model_checkpoint_dir(base, "mnsf") == "run/checkpoints/mnsf"
    assert _model_checkpoint_dir(base, "signed") == "run/checkpoints/signed"
    assert _model_checkpoint_dir(None, "mnsf") is None
    with pytest.raises(ValueError, match="unsupported checkpoint model"):
        _model_checkpoint_dir(base, "other")


def test_section_contract_keeps_barcode_opaque_and_counts_integer() -> None:
    data = SectionData(
        "S1", "P1", "B1", "L1", ("AA-1", "AA-2"), np.array([[0, 0], [1, 0]], float),
        sparse.csr_matrix([[1, 0], [0, 2]], dtype=np.int32), ("G1", "G2"),
    )
    assert data.barcode == ("AA-1", "AA-2")
    with pytest.raises(R04ContractError):
        SectionData("S1", "P1", None, "L1", ("AA",), np.array([[0, 0]], float), np.array([[0.5]]), ("G1",))
    with pytest.raises(R04ContractError):
        validate_count_matrix(np.array([[0.5]]))


def test_geometry_is_continuous_not_a_cluster_label() -> None:
    coords = np.array([[0, 0], [1, 0], [0, 1]], float)
    normalized, scale = normalize_coordinates(coords)
    assert scale > 0
    kernel = matern32_kernel(normalized, normalized, 2.0)
    assert kernel.shape == (3, 3)
    assert np.allclose(np.diag(kernel), 1.0)


def test_composition_adjustment_removes_composition_only_signal() -> None:
    composition = np.array([[0.1 + i / 100, 0.9 - i / 100] for i in range(20)])
    field = np.log(composition[:, 0] / composition[:, 1])
    groups = np.repeat(["P1", "P2", "P3", "P4"], 5)
    result = crossfit_composition_adjustment(field, composition, groups=groups)
    assert np.var(result.adjusted) < np.var(field) * 0.05
    assert result.uncertainty_available is False
    assert len(np.unique(result.fold_ids)) == 4
    assert ilr_transform(composition).shape == (20, 1)
    draws = np.stack([composition, np.clip(composition + 0.005, 1e-3, 1.0)], axis=0)
    with_draws = crossfit_composition_adjustment(field, composition, groups=groups, composition_draws=draws)
    assert with_draws.uncertainty_available is True
    assert with_draws.adjusted_draws is not None and with_draws.adjusted_draws.shape == (2, 20)


def test_candidate_union_preserves_overlapping_fields() -> None:
    field = np.linspace(-1, 1, 10)
    mnsf = [_fit("M1", "S1", np.array([1.0, 0.0]), field), _fit("M2", "S1", np.array([0.0, 1.0]), -field)]
    signed = [_fit("Z1", "S1", np.array([1.0, 0.0]), field), _fit("Z2", "S1", np.array([0.0, 1.0]), -field)]
    candidates = match_model_factors(mnsf, signed)
    assert [candidate.tier for candidate in candidates] == ["BOTH", "BOTH"]
    assert all(candidate.member_factor_ids for candidate in candidates)
    assert all(candidate.restart_fraction == 0.0 and candidate.patient_fraction == 0.0 for candidate in candidates)
    stable = aggregate_candidate_stability(candidates[0], restart_total=5, restart_hits=4, patient_total=10, patient_hits=8)
    assert stable.restart_fraction == 0.8 and stable.patient_fraction == 0.8
    assert stable.diagnostics["restart_total"] == 5


def test_candidate_matching_aggregates_factor_across_sections() -> None:
    field = np.linspace(-1, 1, 10)
    mnsf = [_fit("M1", section, np.array([1.0, 0.0]), field) for section in ("S1", "S2")]
    signed = [_fit("Z1", section, np.array([0.8, 0.6]), field) for section in ("S1", "S2")]
    candidates = match_model_factors(mnsf, signed)
    assert len(candidates) == 1
    assert candidates[0].tier == "BOTH"
    assert candidates[0].loading_cosine == pytest.approx(0.8)
    assert candidates[0].diagnostics["matched_sections"] == ["S1", "S2"]


def test_candidate_matching_rejects_different_input_hashes() -> None:
    field = np.linspace(-1, 1, 10)
    mnsf = [_fit("M1", "S1", np.array([1.0, 0.0]), field, "input-a")]
    signed = [_fit("Z1", "S1", np.array([1.0, 0.0]), field, "input-b")]
    candidates = match_model_factors(mnsf, signed)
    assert [candidate.tier for candidate in candidates] == ["MNSF_ONLY", "SIGNED_ONLY"]


def test_section_content_hash_binds_counts_and_coordinates() -> None:
    base = SectionData(
        "S1", "P1", "B1", "L1", ("AA-1", "AA-2"), np.array([[0, 0], [1, 0]], float),
        sparse.csr_matrix([[1, 0], [0, 2]], dtype=np.int32), ("G1", "G2"),
    )
    changed_counts = replace(base, counts=sparse.csr_matrix([[2, 0], [0, 2]], dtype=np.int32))
    changed_coords = replace(base, coords=np.array([[0, 0], [2, 0]], float))
    assert section_content_hash(base) != section_content_hash(changed_counts)
    assert section_content_hash(base) != section_content_hash(changed_coords)


def test_training_gene_selection_is_patient_grouped_and_deterministic() -> None:
    sections = []
    for patient, offset in (("P1", 0), ("P2", 1), ("P3", 2)):
        sections.append(SectionData(
            f"S{patient}", patient, f"B{patient}", "L1",
            tuple(f"c{i}" for i in range(20)),
            np.column_stack([np.arange(20), np.zeros(20)]),
            sparse.csr_matrix(np.column_stack([
                np.ones(20, dtype=np.int32),
                np.ones(20, dtype=np.int32) * (offset + 1),
                np.arange(20, dtype=np.int32) % 3,
                np.arange(20, dtype=np.int32) % 2,
            ])),
            ("G1", "G2", "G3", "MT-ND1"),
        ))
    first = select_training_genes(sections, n_genes=2)
    second = select_training_genes(sections, n_genes=2)
    assert first.gene_id == second.gene_id
    assert "MT-ND1" not in first.gene_id


def test_10x_headerless_positions_are_read_without_using_first_row_as_header(tmp_path) -> None:
    path = tmp_path / "tissue_positions_list.csv"
    path.write_text("bc1,1,4,5,100,200\nbc2,0,6,7,110,210\n", encoding="utf-8")
    barcodes, coords = read_10x_positions(path)
    assert barcodes == ("bc1",)
    assert coords.tolist() == [[4.0, 5.0]]


def test_config_fingerprint_changes_with_model_settings() -> None:
    assert config_fingerprint({"factors": 2}) != config_fingerprint({"factors": 3})


def test_candidate_registry_roundtrip_rechecks_hash(tmp_path) -> None:
    candidate = CandidateField("C1", "MNSF_ONLY", ("M1",), 1, 1, 0, 0, input_hash="input", field_hash="field")
    path = tmp_path / "registry.json"
    registry_hash = write_candidate_registry(path, [candidate])
    loaded_hash, loaded = read_candidate_registry(path)
    assert loaded_hash == registry_hash
    assert loaded[0].member_factor_ids == ("M1",)
    path.write_text(path.read_text(encoding="utf-8").replace("\"field_hash\": \"field\"", "\"field_hash\": \"tampered\""), encoding="utf-8")
    with pytest.raises(ValueError, match="hash"):
        read_candidate_registry(path)


def test_gate_is_fail_closed_for_boundary_and_missing_replication() -> None:
    candidate = CandidateField("C1", "MNSF_ONLY", ("M1",), 1, 1, 1, 1)
    result = evaluate_candidate_gate(
        candidate, heldout_delta=0.2, heldout_ci_low=0.1,
        spatial_variance_probability=0.99, null_fdr=0.04, lengthscale=0.75,
        max_patient_weight=0.2, independent_lineages=2, groups_per_lineage=3,
        single_model=True,
    )
    assert result.status == "SCIENTIFIC_FAIL_CANDIDATE_GATE"
    assert "lengthscale_at_boundary" in result.reasons
    assert "insufficient_independent_lineages" in result.reasons
    complete = evaluate_candidate_gate(
        candidate, heldout_delta=0.2, heldout_ci_low=0.1,
        spatial_variance_probability=0.99, null_fdr=0.04, lengthscale=2.0,
        max_patient_weight=0.2, independent_lineages=3, groups_per_lineage=3,
        single_model=True, uncertainty_complete=False,
    )
    assert "uncertainty_does_not_cover_loading_dispersion_and_lengthscale" in complete.reasons


def test_gate_can_pass_only_with_explicit_supported_uncertainty() -> None:
    candidate = CandidateField("C1", "BOTH", ("M1", "Z1"), 0.9, 0.8, 0.8, 0.8, input_hash="input", field_hash="field")
    result = evaluate_candidate_gate(
        candidate, heldout_delta=0.2, heldout_ci_low=0.1,
        spatial_variance_probability=0.99, null_fdr=0.04, lengthscale=2.0,
        max_patient_weight=0.2, independent_lineages=2, groups_per_lineage=3,
        uncertainty_complete=True, lengthscale_supported=True,
    )
    assert result.status == "PASS_R04_STABLE_FIELD"


def test_gate_cli_blocks_missing_metrics(tmp_path, monkeypatch) -> None:
    from scripts.r04_gate import main

    candidate = CandidateField("C1", "MNSF_ONLY", ("M1",), 1, 1, 0.8, 0.8, input_hash="input", field_hash="field")
    registry = tmp_path / "registry.json"
    registry_hash = write_candidate_registry(registry, [candidate])
    metrics = tmp_path / "metrics.json"
    metrics.write_text(
        '{"candidate_registry_hash": "' + registry_hash + '", "input_hash": "input", "field_hash": "field", '
        '"uncertainty_complete": true, "lengthscale_supported": true}',
        encoding="utf-8",
    )
    output = tmp_path / "gate.json"
    monkeypatch.setattr(sys, "argv", ["r04_gate", "--registry", str(registry), "--candidate-id", "C1", "--metrics", str(metrics), "--output", str(output)])
    assert main() == 2
    assert "required_gate_fields_missing" in output.read_text(encoding="utf-8")


def test_anchor_hash_mismatch_blocks_before_correspondence() -> None:
    result = evaluate_anchor_correspondence(
        candidate_registry_hash="a", expected_registry_hash="b",
        field_values=np.arange(6), anchor_values=np.arange(6), groups=["P1", "P2", "P3"],
    )
    assert result.status == "BLOCKED_CANDIDATE_HASH_MISMATCH"


def test_target_metadata_is_not_a_model_input() -> None:
    with pytest.raises(R04ContractError):
        reject_target_metadata({"ground_truth": "TLS"})


def test_h5ad_reference_reads_only_allowlisted_counts_and_annotations(tmp_path) -> None:
    import h5py

    path = tmp_path / "reference.h5ad"
    with h5py.File(path, "w") as handle:
        handle.create_dataset("X", data=np.arange(12, dtype=np.int32).reshape(4, 3))
        obs = handle.create_group("obs")
        obs.create_dataset("_index", data=np.asarray([b"c1", b"c2", b"c3", b"c4"]))
        obs.create_dataset("donor", data=np.asarray([b"p1", b"p1", b"p2", b"p2"]))
        obs.create_dataset("cell_type", data=np.asarray([b"CD4", b"CD8", b"NK", b"B"]))
        obs.create_dataset("response", data=np.asarray([b"R", b"NR", b"R", b"NR"]))
        var = handle.create_group("var")
        var.create_dataset("_index", data=np.asarray([b"g1", b"g2", b"g3"]))
        obsm = handle.create_group("obsm")
        obsm.create_dataset("spatial", data=np.asarray([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]))
    matrix, genes, donors, broad, source = read_allowlisted_h5ad(
        path, count_key="X", donor_column="donor", cell_type_column="cell_type", source_id="TEST",
        per_donor_type=1, max_total=4,
    )
    assert matrix.shape == (4, 3)
    assert genes == ("g1", "g2", "g3")
    assert set(broad) == {"T_CD4", "T_CD8", "NK", "B_plasma"}
    assert source == "TEST"
    counts, count_barcodes, count_genes, coords = read_h5ad_counts(path)
    assert counts.shape == (4, 3) and count_genes == genes and count_barcodes[0] == "c1" and coords.shape == (4, 2)
    with pytest.raises(ValueError):
        read_allowlisted_h5ad(path, count_key="layers/response", donor_column="donor", cell_type_column="cell_type", source_id="TEST")


def test_dense_h5ad_rejects_fractional_counts_before_integer_cast(tmp_path) -> None:
    import h5py

    path = tmp_path / "fractional.h5ad"
    with h5py.File(path, "w") as handle:
        handle.create_dataset("X", data=np.asarray([[0.0, 1.25], [2.0, 0.0]], dtype=np.float32))
        obs = handle.create_group("obs")
        obs.create_dataset("_index", data=np.asarray([b"c1", b"c2"]))
        var = handle.create_group("var")
        var.create_dataset("_index", data=np.asarray([b"g1", b"g2"]))
        var.create_dataset("gene_ids", data=np.asarray([b"ENSG000001", b"ENSG000002"]))
        obsm = handle.create_group("obsm")
        obsm.create_dataset("spatial", data=np.asarray([[0.0, 0.0], [1.0, 0.0]]))

    with pytest.raises(R04ContractError, match="integer raw counts"):
        read_h5ad_counts(path)


def test_h5ad_prefers_stable_var_gene_ids_when_available(tmp_path) -> None:
    import h5py

    path = tmp_path / "gene_ids.h5ad"
    with h5py.File(path, "w") as handle:
        handle.create_dataset("X", data=np.asarray([[1, 0]], dtype=np.int32))
        obs = handle.create_group("obs")
        obs.create_dataset("_index", data=np.asarray([b"c1"]))
        var = handle.create_group("var")
        var.create_dataset("_index", data=np.asarray([b"GENE1", b"GENE2"]))
        var.create_dataset("gene_ids", data=np.asarray([b"ENSG000001", b"ENSG000002"]))
        obsm = handle.create_group("obsm")
        obsm.create_dataset("spatial", data=np.asarray([[0.0, 0.0]]))
    _, _, genes, _ = read_h5ad_counts(path)
    assert genes == ("ENSG000001", "ENSG000002")


def test_reference_panel_uses_common_genes_and_donor_caps() -> None:
    panel = build_reference_panel(
        [sparse.csr_matrix([[1, 2, 3], [4, 5, 6]]), sparse.csr_matrix([[7, 8], [9, 10]])],
        [("g1", "g2", "g3"), ("g2", "g3")],
        [("p1", "p1"), ("p2", "p2")],
        [("CD4", "CD4"), ("CD8", "CD8")],
        ["A", "B"], per_donor_type=1, max_total=10,
    )
    assert panel.gene_id == ("g2", "g3")
    assert panel.counts.toarray().tolist() == [[2, 3], [7, 8]]


def test_input_manifest_requires_explicit_locator_and_frozen_split(tmp_path) -> None:
    split = tmp_path / "splits.tsv"
    split.write_text(
        "physical_unit_id\tlogical_unit_id\tpatient_id\tblock_id\tleakage_group_id\tprimary_role\touter_fold\tblock_level_eligible\tsplit_status\n"
        "U1\tL1\tP1\tB1\tG1\tdiscovery\tF1\tyes\tFROZEN_R02\n"
        "U2\tL1\tP2\tB2\tG2\ttraining\tF2\tno\tFROZEN_R02\n",
        encoding="utf-8",
    )
    locator = tmp_path / "locators.tsv"
    locator.write_text(
        "physical_unit_id\tmatrix_locator\tcoordinate_locator\tmatrix_kind\tread_authorization\n"
        f"U1\t{tmp_path / 'counts.h5'}\t{tmp_path / 'coords.csv'}\t10x_h5\tMOLECULES_AND_COORDINATES\n",
        encoding="utf-8",
    )
    rows = build_input_rows(split, locator, tmp_path)
    assert len(rows) == 1 and rows[0]["section_id"] == "U1"
    locator.write_text(locator.read_text(encoding="utf-8").replace("U1\t", "U3\t"), encoding="utf-8")
    with pytest.raises(ValueError, match="no explicit R04 locator"):
        build_input_rows(split, locator, tmp_path)


def test_input_manifest_can_fail_closed_on_missing_assets(tmp_path) -> None:
    split = tmp_path / "splits.tsv"
    split.write_text(
        "physical_unit_id\tlogical_unit_id\tpatient_id\tblock_id\tleakage_group_id\tprimary_role\touter_fold\tblock_level_eligible\tsplit_status\n"
        "U1\tL1\tP1\tB1\tG1\ttraining\tF1\tyes\tFROZEN_R02\n",
        encoding="utf-8",
    )
    locator = tmp_path / "locators.tsv"
    locator.write_text(
        "physical_unit_id\tmatrix_locator\tcoordinate_locator\tmatrix_kind\tread_authorization\n"
        f"U1\t{tmp_path / 'missing.h5'}\t{tmp_path / 'coords.csv'}\t10x_h5\tMOLECULES_AND_COORDINATES\n",
        encoding="utf-8",
    )
    with pytest.raises(FileNotFoundError, match="matrix locator does not exist"):
        build_input_rows(split, locator, tmp_path, require_existing=True)
