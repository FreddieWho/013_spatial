from __future__ import annotations

import numpy as np
import pytest

from r04.diagnostics import mnsf_spatial_effect_matrix
from r04.structure_readout import effect_coordinates, evaluate_shared_specific_readout


def _export(field: np.ndarray, loading: np.ndarray, *, rate: bool = False) -> dict[str, object]:
    n, k = field.shape
    genes = np.asarray([f"gene-{i}" for i in range(loading.shape[0])])
    result: dict[str, object] = {
        "field_mean": field,
        "evaluation_loading": loading,
        "factor_amplitude": np.ones(loading.shape[1]),
        "evaluation_gene_ids": genes,
        "section_ids": np.asarray(["section"] * n),
        "barcodes": np.asarray([f"barcode-{i}" for i in range(n)]),
    }
    if rate:
        values = mnsf_spatial_effect_matrix(field, loading, np.ones(k))
        values -= values.mean(axis=0, keepdims=True)
        result["spatial_rate_effect_centered"] = values
    return result


def test_effect_coordinates_use_common_rotation_stable_gene_basis() -> None:
    rng = np.random.default_rng(7)
    train_field = rng.normal(size=(80, 3))
    test_field = rng.normal(size=(35, 3))
    loading = rng.normal(size=(12, 3))
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))

    train_z, test_z, metadata = effect_coordinates(
        _export(train_field, loading), _export(test_field, loading), max_rank=3
    )
    rotated_train_z, rotated_test_z, rotated_metadata = effect_coordinates(
        _export(train_field @ q, loading @ q),
        _export(test_field @ q, loading @ q),
        max_rank=3,
    )

    assert metadata["representation"] == (
        "rotation_stable_centered_linear_spot_by_gene_effect"
    )
    assert metadata["factor_naming"] == "forbidden"
    assert rotated_metadata["singular_values"] == pytest.approx(
        metadata["singular_values"]
    )
    assert rotated_train_z @ rotated_train_z.T == pytest.approx(train_z @ train_z.T)
    assert rotated_test_z @ rotated_test_z.T == pytest.approx(test_z @ test_z.T)


def test_shared_plus_specific_residual_recovers_two_distinct_targets() -> None:
    rng = np.random.default_rng(11)
    train_field = rng.normal(size=(180, 3))
    heldout_field = rng.normal(size=(90, 3))
    loading = rng.normal(size=(18, 3))
    train_z, heldout_z, representation = effect_coordinates(
        _export(train_field, loading), _export(heldout_field, loading)
    )
    train_labels = {
        "TLS": train_z[:, 0] + 1.2 * train_z[:, 1] + rng.normal(scale=0.1, size=len(train_z)),
        "TUMOR_STROMA_BOUNDARY": train_z[:, 0] - 0.9 * train_z[:, 2] + rng.normal(scale=0.1, size=len(train_z)),
    }
    heldout_labels = {
        "TLS": heldout_z[:, 0] + 1.2 * heldout_z[:, 1] + rng.normal(scale=0.1, size=len(heldout_z)),
        "TUMOR_STROMA_BOUNDARY": heldout_z[:, 0] - 0.9 * heldout_z[:, 2] + rng.normal(scale=0.1, size=len(heldout_z)),
    }

    result = evaluate_shared_specific_readout(
        train_z, heldout_z, train_labels, heldout_labels,
        representation=representation,
    )

    assert result["shared_plus_specific_residual"]["mean_mse"] < (
        result["shared_only"]["mean_mse"]
    )
    assert result["specific_increment"]["mean_mse_delta_specific_minus_shared"] < 0


def test_primary_rate_effect_path_avoids_dense_gene_svd() -> None:
    rng = np.random.default_rng(17)
    loading = np.abs(rng.normal(size=(20, 3)))
    train = _export(rng.normal(size=(100, 3)), loading, rate=True)
    heldout = _export(rng.normal(size=(40, 3)), loading, rate=True)
    train_z, heldout_z, metadata = effect_coordinates(train, heldout, max_rank=2)
    assert train_z.shape == (100, 2)
    assert heldout_z.shape == (40, 2)
    assert metadata["representation"] == "centered_mnsf_spatial_rate_contribution"
    assert metadata["rank"] == 2


def test_shared_only_is_sufficient_for_fully_shared_targets() -> None:
    rng = np.random.default_rng(19)
    train_z = rng.normal(size=(160, 2))
    heldout_z = rng.normal(size=(80, 2))
    common_train = train_z[:, 0] - 0.4 * train_z[:, 1]
    common_test = heldout_z[:, 0] - 0.4 * heldout_z[:, 1]
    train_labels = {"A": common_train, "B": common_train}
    heldout_labels = {"A": common_test, "B": common_test}

    result = evaluate_shared_specific_readout(train_z, heldout_z, train_labels, heldout_labels)

    assert result["shared_plus_specific_residual"]["mean_mse"] < 1e-10
    assert result["specific_increment"]["mean_mse_delta_specific_minus_shared"] < 1e-10


def test_readout_requires_two_aligned_structure_targets() -> None:
    with pytest.raises(ValueError, match="exactly two"):
        evaluate_shared_specific_readout(
            np.ones((3, 2)), np.ones((2, 2)), {"only": [1, 2, 3]}, {"only": [1, 2]}
        )
