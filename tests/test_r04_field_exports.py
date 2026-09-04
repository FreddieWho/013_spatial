from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from r04.field_exports import (
    SCHEMA,
    read_heldout_field_export,
    reconstructed_linear_effect,
    write_heldout_field_export,
)
from r04.diagnostics import mnsf_spatial_effect_matrix
from r04.types import FieldFit, SectionData


def _fixtures() -> tuple[list[SectionData], list[list[FieldFit]]]:
    sections: list[SectionData] = []
    groups: list[list[FieldFit]] = []
    for section_index, n_spots in enumerate((3, 2)):
        section_id = f"section-{section_index}"
        section = SectionData(
            section_id=section_id,
            patient_id=f"patient-{section_index}",
            block_id=f"block-{section_index}",
            lineage="training",
            barcode=tuple(f"{section_id}-barcode-{i}" for i in range(n_spots)),
            coords=np.arange(n_spots * 2, dtype=float).reshape(n_spots, 2),
            counts=np.zeros((n_spots, 6), dtype=int),
            gene_id=tuple(f"gene-{i}" for i in range(6)),
        )
        fields = [
            np.arange(n_spots, dtype=float) - (n_spots - 1) / 2,
            np.asarray([1.0, -1.0, 0.0][:n_spots], dtype=float),
        ]
        groups.append([
            FieldFit(
                model_id="mNSF_NB_VI",
                factor_id=f"MNSF_{factor + 1:02d}",
                section_id=section_id,
                loading=np.ones(3),
                field_mean=value,
                field_sd=np.full(n_spots, 0.1),
                lengthscale=3.0,
                spatial_variance=float(np.var(value)),
                input_hash="input",
            )
            for factor, value in enumerate(fields)
        ])
        sections.append(section)
    return sections, groups


def _write_fixture(path: Path) -> None:
    sections, groups = _fixtures()
    write_heldout_field_export(
        path,
        sections=sections,
        field_groups=groups,
        evaluation_loading=np.asarray([[1.0, 0.0], [0.0, 2.0], [1.0, 1.0]]),
        factor_amplitude=np.asarray([1.0, 2.0]),
        frozen_gene_ids=tuple(f"gene-{i}" for i in range(6)),
        adaptation_gene_ids=("gene-0", "gene-2", "gene-4"),
        evaluation_gene_ids=("gene-1", "gene-3", "gene-5"),
        adaptation_indices=np.asarray([0, 2, 4]),
        evaluation_indices=np.asarray([1, 3, 5]),
        provenance={"source": "synthetic", "validation_gt": "NOT_READ"},
    )


def test_field_export_round_trip_preserves_spot_alignment_and_effect(tmp_path: Path) -> None:
    path = tmp_path / "heldout.npz"
    _write_fixture(path)

    exported = read_heldout_field_export(path)
    assert exported["metadata"]["schema"] == SCHEMA
    assert exported["metadata"]["effect_representation"] == (
        "centered_mnsf_spatial_rate_contribution_primary"
    )
    assert exported["barcodes"].tolist() == [
        "section-0-barcode-0",
        "section-0-barcode-1",
        "section-0-barcode-2",
        "section-1-barcode-0",
        "section-1-barcode-1",
    ]
    effect = reconstructed_linear_effect(exported)
    field = np.asarray(exported["field_mean"], dtype=float)
    loading = np.asarray(exported["evaluation_loading"], dtype=float)
    assert effect == pytest.approx(field @ loading.T)
    rate = np.asarray(exported["spatial_rate_effect_centered"], dtype=float)
    assert rate.shape == (5, 3)
    for section_id in np.unique(exported["section_ids"]):
        rows = exported["section_ids"] == section_id
        assert rate[rows].mean(axis=0) == pytest.approx(0.0, abs=1e-6)


def test_field_export_fails_closed_on_tampered_array(tmp_path: Path) -> None:
    path = tmp_path / "heldout.npz"
    _write_fixture(path)
    original = path.read_bytes()
    path.write_bytes(original[:-1] + bytes([original[-1] ^ 1]))

    with pytest.raises(ValueError, match="hash mismatch"):
        read_heldout_field_export(path)


def test_field_export_requires_a_partition_of_frozen_genes(tmp_path: Path) -> None:
    sections, groups = _fixtures()
    with pytest.raises(ValueError, match="partition the frozen universe"):
        write_heldout_field_export(
            tmp_path / "heldout.npz",
            sections=sections,
            field_groups=groups,
            evaluation_loading=np.ones((3, 2)),
            factor_amplitude=np.ones(2),
            frozen_gene_ids=tuple(f"gene-{i}" for i in range(6)),
            adaptation_gene_ids=("gene-0", "gene-2", "gene-4"),
            evaluation_gene_ids=("gene-1", "gene-3", "gene-4"),
            adaptation_indices=np.asarray([0, 2, 4]),
            evaluation_indices=np.asarray([1, 3, 5]),
            provenance={},
        )


def test_linear_effect_is_invariant_to_orthogonal_basis_rotation() -> None:
    rng = np.random.default_rng(23)
    field = rng.normal(size=(20, 3))
    loading = rng.normal(size=(11, 3))
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))

    original = field @ loading.T
    rotated = (field @ q) @ (loading @ q).T
    assert rotated == pytest.approx(original)


def test_rate_effect_is_invariant_to_factor_permutation() -> None:
    rng = np.random.default_rng(31)
    field = rng.normal(size=(25, 3))
    loading = np.abs(rng.normal(size=(9, 3)))
    amplitude = np.asarray([0.4, 1.1, 2.0])
    permutation = np.asarray([2, 0, 1])

    original = mnsf_spatial_effect_matrix(field, loading, amplitude)
    permuted = mnsf_spatial_effect_matrix(
        field[:, permutation], loading[:, permutation], amplitude[permutation]
    )
    assert permuted == pytest.approx(original)
