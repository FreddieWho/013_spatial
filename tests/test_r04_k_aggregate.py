from __future__ import annotations

import copy
from functools import lru_cache

import numpy as np
import pytest

from scripts.r04_k_aggregate import aggregate_k_shards
from scripts.r04_k_contract import GENERATOR_CONTRACT
from r04.io_contract import sections_content_hash
from r04.synthetic import make_identifiability_calibration_sections


@lru_cache(maxsize=None)
def _replayed(seed: int) -> tuple[dict, str]:
    sections, _, truth = make_identifiability_calibration_sections(
        mode="two_independent",
        domain="disconnected",
        generator="mnsf_correct",
        sections=3,
        spots_per_section=50,
        genes=24,
        seed=seed,
    )
    return truth, sections_content_hash(sections)


def _recovery(support: bool, required: int) -> list[dict]:
    evaluation_folds = (list(range(12)), list(range(12, 24)))
    universe = set(range(24))
    return [
        {
            "adaptation_genes": sorted(universe - set(evaluation)),
            "evaluation_genes": evaluation,
            "canonical_correlations": [0.9 if support else 0.1] * required,
            "direction_matched_null_95": [0.2] * required,
            "inference_platform": {"converged": support},
        }
        for evaluation in evaluation_folds
    ]


def _shard(k: int, seed: int, scores: list[float], support: bool = True) -> dict:
    assert len(scores) == 2
    truth, input_data_hash = _replayed(seed)
    return {
        "schema": "r04.k_calibration.v2",
        "status": "K_CALIBRATION_COMPLETE",
        "formal_restart_authorized": False,
        "truth": copy.deepcopy(truth),
        "input_data_hash": input_data_hash,
        "generator_contract": dict(GENERATOR_CONTRACT),
        "leakage_control": "heldout_section_fields_adapted_on_disjoint_gene_folds",
        "steps": 600,
        "inference_steps": 250,
        "section_folds": 2,
        "gene_folds": 2,
        "null_draws": 50,
        "k_values": [k],
        "seeds": [seed],
        "holdouts": [0, 1],
        "runs": [
            {
                "holdout": holdout,
                "data_seed": seed,
                "optimization_seed": seed + k * 100 + holdout,
                "split_seed": seed + holdout * 1009,
                "k_model": k,
                "scores": scores,
                "baseline_scores": [-11.0, -11.1],
                "fit_platform": {"converged": support},
                "recovery_pass": support,
                "required_recovered_directions": min(k, 2),
                "recovery": _recovery(support, min(k, 2)),
            }
            for holdout in (0, 1)
        ],
    }


def test_aggregate_selects_smallest_supported_one_se_k() -> None:
    shards = []
    for seed in (11, 13, 17):
        shards.extend([
            _shard(1, seed, [-10.0, -10.1]),
            _shard(2, seed, [-8.0, -8.1]),
            _shard(3, seed, [-8.2, -8.0]),
            _shard(4, seed, [-9.0, -9.1]),
        ])
    result = aggregate_k_shards(shards)
    assert result["one_se_candidate"]["selected_k"] == 2
    assert result["final_k_status"] == "K_SELECTED"


def test_aggregate_keeps_boundary_winner_unidentified() -> None:
    shards = []
    for seed in (11, 13, 17):
        shards.extend([
            _shard(1, seed, [-3.0, -3.1]),
            _shard(2, seed, [-2.0, -2.1]),
        ])
    result = aggregate_k_shards(shards)
    assert result["one_se_candidate"]["search_boundary"]
    assert result["final_k_status"] == "K_SEARCH_BOUNDARY"


def test_aggregate_rejects_duplicate_calibration_cell() -> None:
    shards = [_shard(k, seed, [-8.0, -8.1]) for seed in (11, 13, 17) for k in (1, 2)]
    shards.append(copy.deepcopy(shards[0]))
    with pytest.raises(ValueError, match="duplicate"):
        aggregate_k_shards(shards)


def test_aggregate_rejects_k_dependent_gene_fold() -> None:
    shards = [_shard(k, seed, [-8.0, -8.1]) for seed in (11, 13, 17) for k in (1, 2)]
    target = next(item for item in shards if item["seeds"] == [11] and item["k_values"] == [2])
    target["runs"][0]["recovery"].reverse()
    with pytest.raises(ValueError, match="folds change across K"):
        aggregate_k_shards(shards)


def test_aggregate_rejects_data_change_across_k() -> None:
    shards = [_shard(k, seed, [-8.0, -8.1]) for seed in (11, 13, 17) for k in (1, 2)]
    target = next(item for item in shards if item["seeds"] == [11] and item["k_values"] == [2])
    target["input_data_hash"] = "f" * 64
    with pytest.raises(ValueError, match="input-data hash does not match"):
        aggregate_k_shards(shards)


def test_aggregate_rejects_truth_not_replayed_from_seed() -> None:
    shards = [_shard(k, seed, [-8.0, -8.1]) for seed in (11, 13, 17) for k in (1, 2)]
    shards[0]["truth"]["true_amplitude"][0] += 1.0
    with pytest.raises(ValueError, match="truth metadata does not match"):
        aggregate_k_shards(shards)


def test_aggregate_rejects_nonfinite_heldout_score() -> None:
    shards = [_shard(k, seed, [-8.0, -8.1]) for seed in (11, 13, 17) for k in (1, 2)]
    shards[0]["runs"][0]["scores"][0] = float("nan")
    with pytest.raises(ValueError, match="scores must be finite"):
        aggregate_k_shards(shards)


def test_aggregate_standard_error_uses_independent_data_seeds() -> None:
    shards = []
    expected_seed_means = []
    for seed, score in zip((11, 13, 17), (-7.0, -8.0, -9.0)):
        shards.extend([
            _shard(1, seed, [-10.0, -10.1]),
            _shard(2, seed, [score, score - 0.2]),
        ])
        expected_seed_means.append(score - 0.1)
    result = aggregate_k_shards(shards)
    expected = float(np.std(expected_seed_means, ddof=1) / np.sqrt(3))
    assert result["independent_resampling_unit"] == "data_seed"
    assert result["heldout_score_se_by_k"]["2"] == pytest.approx(expected)
