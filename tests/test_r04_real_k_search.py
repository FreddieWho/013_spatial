from __future__ import annotations

import numpy as np
import pytest

from r04.runtime import resource_status
from scripts.r04_real_k_search import (
    _aggregate,
    _combine_inference_platforms,
    _is_wiring_only,
    _mark_wiring_only,
    _optimization_seed,
    _parse_fold_values,
    _parse_gene_batch_size,
)


def test_crossfit_platform_requires_every_gene_split_to_converge() -> None:
    combined = _combine_inference_platforms([
        {"converged": False, "net_relative_improvement": 0.002},
        {"converged": True, "net_relative_improvement": 0.0002},
    ])
    assert combined["converged"] is False
    assert combined["split_converged"] == [False, True]
    assert combined["split_count"] == 2


def _cell(k: int, fold: int) -> dict[str, object]:
    return {
        "schema": "r04.real_k_cell.v1",
        "status": "FIT_AND_SCORED",
        "k_model": k,
        "fold": fold,
        "patient_scores": {
            f"p{fold * 2}": float(-10 + k * 0.5),
            f"p{fold * 2 + 1}": float(-11 + k * 0.5),
        },
    }


def test_real_k_aggregate_uses_k0_and_grouped_patient_bootstrap() -> None:
    result = _aggregate(
        [_cell(k, fold) for k in (0, 1, 2) for fold in range(2)],
        [0, 1, 2],
        seed=17,
    )
    assert result["selected_k"] == 2
    assert result["delta_to_k0_by_k"]["2"]["ci_low"] > 0
    assert result["null_status"] == "NOT_RUN"
    assert result["next_k_values"] == [3, 4]


def test_real_k_aggregate_keeps_k0_when_increment_ci_crosses_zero() -> None:
    cells = []
    for k in (0, 1):
        for fold in range(2):
            cells.append({
                "schema": "r04.real_k_cell.v1",
                "status": "FIT_AND_SCORED",
                "k_model": k,
                "fold": fold,
                "patient_scores": {f"p{fold * 2}": -10.0, f"p{fold * 2 + 1}": -10.0},
            })
    result = _aggregate(cells, [0, 1], seed=19)
    assert result["selected_k"] == 0
    assert result["status"] == "K0_SELECTED_REAL"


def test_real_k_aggregate_does_not_select_unconverged_cells() -> None:
    cells = []
    for k in (0, 1):
        for fold in range(2):
            cell = _cell(k, fold)
            cell["fit_platform"] = {"converged": False}
            cell["inference_platform"] = {"converged": False}
            cells.append(cell)
    result = _aggregate(cells, [0, 1], seed=23)
    assert result["selected_k"] is None
    assert result["status"] == "K_NOT_IDENTIFIABLE_NOT_CONVERGED"


def test_resource_default_threshold_is_1_2_tb(monkeypatch) -> None:
    monkeypatch.delenv("R04_STORAGE_RESERVE_BYTES", raising=False)
    monkeypatch.setattr("r04.runtime.available_bytes", lambda path: 1_200_000_000_000)
    assert resource_status() == "OK"
    monkeypatch.setattr("r04.runtime.available_bytes", lambda path: 1_199_999_999_999)
    assert resource_status() == "BLOCKED_STORAGE"


def test_resource_threshold_override_is_explicit_and_validated(monkeypatch) -> None:
    monkeypatch.setenv("R04_STORAGE_RESERVE_BYTES", "100")
    monkeypatch.setattr("r04.runtime.available_bytes", lambda path: 100)
    assert resource_status() == "OK"
    monkeypatch.setattr("r04.runtime.available_bytes", lambda path: 99)
    assert resource_status() == "BLOCKED_STORAGE"
    monkeypatch.setenv("R04_STORAGE_RESERVE_BYTES", "not-an-integer")
    with pytest.raises(ValueError, match="must be an integer"):
        resource_status()


def test_parse_fold_values_selects_explicit_available_folds() -> None:
    assert _parse_fold_values("4,0,4", [0, 1, 2, 3, 4]) == [0, 4]
    with pytest.raises(ValueError, match="not present"):
        _parse_fold_values("0,5", [0, 1, 2, 3, 4])


def test_optimization_seed_preserves_canonical_and_separates_restart() -> None:
    from r04.runtime import derived_seed

    canonical = _optimization_seed(20260807, 3, 0, 0)
    restart = _optimization_seed(20260807, 3, 0, 1)
    assert canonical == derived_seed(20260807, "real-k", 3, 0)
    assert restart == _optimization_seed(20260807, 3, 0, 1)
    assert restart != canonical


def test_optimization_seed_rejects_negative_restart() -> None:
    with pytest.raises(ValueError, match="restart index"):
        _optimization_seed(20260807, 3, 0, -1)


def test_fold_subset_is_marked_wiring_only_in_final_result() -> None:
    result = {"selected_k": 3, "status": "K_SELECTED_REAL"}
    assert _is_wiring_only(
        fold_values="0,4",
        max_folds=None,
        max_genes=None,
        max_spots_per_section=None,
    )
    marked = _mark_wiring_only(result, True)
    assert marked["wiring_only"] is True
    assert marked["screen_selected_k"] == 3
    assert marked["selected_k"] is None
    assert marked["status"] == "WIRING_ONLY_NOT_SCIENTIFIC"


def test_gene_batch_size_accepts_full_and_positive_integer() -> None:
    assert _parse_gene_batch_size("full") is None
    assert _parse_gene_batch_size("512") == 512
    with pytest.raises(ValueError, match="positive"):
        _parse_gene_batch_size("0")
