#!/usr/bin/env python3
"""Aggregate leakage-free K calibration shards across data/optimisation seeds."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from r04.runtime import atomic_json
from scripts.r04_k_contract import LEAKAGE_CONTROL, validate_k_shards


def aggregate_k_shards(shards: list[dict[str, object]]) -> dict[str, object]:
    validated = validate_k_shards(
        shards, expected_mode="two_independent", expected_k_eff_true=2
    )
    cells = validated["cells"]
    runs = list(cells.values())
    k_values = validated["k_values"]
    data_seeds = validated["data_seeds"]
    scores_by_k = {
        k: [float(score) for run in runs if int(run["k_model"]) == k for score in run["scores"]]
        for k in k_values
    }
    seed_means_by_k = {
        k: {
            seed: float(np.mean([
                float(score)
                for holdout in (0, 1)
                for score in cells[(k, seed, holdout)]["scores"]
            ]))
            for seed in data_seeds
        }
        for k in k_values
    }
    means = {
        k: float(np.mean(list(values.values())))
        for k, values in seed_means_by_k.items()
    }
    ses = {
        k: float(np.std(list(values.values()), ddof=1) / np.sqrt(len(values)))
        if len(values) > 1 else 0.0
        for k, values in seed_means_by_k.items()
    }
    support_fraction: dict[int, float] = {}
    for k in k_values:
        seed_support = []
        for seed in data_seeds:
            shard_runs = [cells[(k, seed, holdout)] for holdout in (0, 1)]
            seed_support.append(all(
                bool(run["fit_platform"]["converged"])
                and bool(run["recovery_pass"])
                and all(bool(item["inference_platform"]["converged"]) for item in run["recovery"])
                for run in shard_runs
            ))
        support_fraction[k] = float(np.mean(seed_support))
    best_k = max(k_values, key=means.get)
    paired_differences = {
        k: [
            seed_means_by_k[k][seed] - seed_means_by_k[best_k][seed]
            for seed in data_seeds
        ]
        for k in k_values
    }
    predictive_eligible = [
        k for k in k_values if means[k] >= means[best_k] - ses[best_k]
    ]
    eligible = [k for k in predictive_eligible if support_fraction[k] >= 2 / 3]
    selected = min(eligible) if eligible else None
    boundary = selected == max(k_values) if selected is not None else best_k == max(k_values)
    if selected is None:
        final_status = "K_NOT_IDENTIFIABLE"
    elif boundary:
        final_status = "K_SEARCH_BOUNDARY"
    elif len(data_seeds) >= 3:
        final_status = "K_SELECTED"
    else:
        final_status = "K_NOT_IDENTIFIABLE"
    return {
        "schema": "r04.k_calibration_aggregate.v1",
        "status": "K_CALIBRATION_AGGREGATED",
        "calibration_cell_count": len(cells),
        "data_seeds": data_seeds,
        "data_optimisation_seeds": data_seeds,
        "k_values": k_values,
        "leakage_control": LEAKAGE_CONTROL,
        "execution_contract": validated["execution_contract"],
        "generator_contract": validated["generator_contract"],
        "input_data_hash_by_seed": {
            str(seed): value for seed, value in validated["input_hash_by_seed"].items()
        },
        "generator_replay_verified": validated["generator_replay_verified"],
        "scores_by_k": {str(k): values for k, values in scores_by_k.items()},
        "independent_resampling_unit": "data_seed",
        "heldout_score_mean_by_seed_and_k": {
            str(k): {str(seed): value for seed, value in values.items()}
            for k, values in seed_means_by_k.items()
        },
        "paired_difference_to_best_by_k": {
            str(k): values for k, values in paired_differences.items()
        },
        "heldout_score_mean_by_k": {str(k): value for k, value in means.items()},
        "heldout_score_se_by_k": {str(k): value for k, value in ses.items()},
        "support_fraction_by_k": {str(k): value for k, value in support_fraction.items()},
        "one_se_candidate": {
            "best_predictive_k": best_k,
            "predictive_eligible_k": predictive_eligible,
            "eligible_k": eligible,
            "selected_k": selected,
            "search_boundary": boundary,
        },
        "k_eff_calibration_truth": int(next(iter(validated["truth_by_seed"].values()))["k_eff_true"]),
        "final_k_status": final_status,
        "formal_restart_authorized": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    shards = [json.loads(path.read_text(encoding="utf-8")) for path in args.input]
    result = aggregate_k_shards(shards)
    result["input_provenance"] = [
        {
            "path": str(path),
            "schema": shard.get("schema"),
            "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
        }
        for path, shard, payload in (
            (path, shard, path.read_bytes()) for path, shard in zip(args.input, shards)
        )
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(args.output, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
