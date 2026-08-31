from __future__ import annotations

import copy
from functools import lru_cache

import pytest

from r04.io_contract import sections_content_hash
from r04.synthetic import make_identifiability_calibration_sections
from scripts.r04_k0_gate import evaluate_no_field_shards
from scripts.r04_k_contract import GENERATOR_CONTRACT


@lru_cache(maxsize=None)
def _replayed(seed: int) -> tuple[dict, str]:
    sections, _, truth = make_identifiability_calibration_sections(
        mode="no_field",
        domain="disconnected",
        generator="mnsf_correct",
        sections=3,
        spots_per_section=50,
        genes=24,
        seed=seed,
    )
    return truth, sections_content_hash(sections)


def _recovery() -> list[dict]:
    evaluation_folds = (list(range(12)), list(range(12, 24)))
    universe = set(range(24))
    return [
        {
            "adaptation_genes": sorted(universe - set(evaluation)),
            "evaluation_genes": evaluation,
            "canonical_correlations": [],
            "direction_matched_null_95": [],
            "inference_platform": {"converged": True},
        }
        for evaluation in evaluation_folds
    ]


def _shard(seed: int, k: int, score: float, baseline: float) -> dict:
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
        "seeds": [seed],
        "k_values": [k],
        "holdouts": [0, 1],
        "runs": [
            {
                "holdout": holdout,
                "data_seed": seed,
                "optimization_seed": seed + k * 100 + holdout,
                "split_seed": seed + holdout * 1009,
                "k_model": k,
                "scores": [score, score],
                "baseline_scores": [baseline, baseline],
                "fit_platform": {"converged": True},
                "recovery_pass": True,
                "required_recovered_directions": 0,
                "recovery": _recovery(),
            }
            for holdout in (0, 1)
        ],
    }


def test_no_field_gate_selects_zero_when_spatial_models_do_not_improve() -> None:
    shards = []
    for seed in (11, 13, 17):
        shards.extend([_shard(seed, 1, -10.2, -10.0), _shard(seed, 2, -10.4, -10.0)])
    result = evaluate_no_field_shards(shards)
    assert result["status"] == "K0_SELECTED"
    assert result["one_se_candidate"]["selected_k"] == 0


def test_no_field_gate_rejects_k_dependent_baseline() -> None:
    shards = []
    for seed in (11, 13, 17):
        shards.extend([_shard(seed, 1, -10.2, -10.0), _shard(seed, 2, -10.4, -10.0)])
    target = next(item for item in shards if item["seeds"] == [11] and item["k_values"] == [2])
    target["runs"][0]["baseline_scores"][0] = -9.0
    with pytest.raises(ValueError, match="baseline scores change across K"):
        evaluate_no_field_shards(shards)
