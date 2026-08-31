#!/usr/bin/env python3
"""Select K_eff=0 on no-field controls when spatial models add no held-out value."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from r04.runtime import atomic_json
from scripts.r04_k_contract import LEAKAGE_CONTROL, validate_k_shards


def evaluate_no_field_shards(shards: list[dict[str, object]]) -> dict[str, object]:
    validated = validate_k_shards(
        shards, expected_mode="no_field", expected_k_eff_true=0
    )
    cells = validated["cells"]
    seeds = validated["data_seeds"]
    k_values = validated["k_values"]
    by_k = {
        k: [
            float(score)
            for seed in seeds for holdout in (0, 1)
            for score in cells[(k, seed, holdout)]["scores"]
        ]
        for k in k_values
    }
    baseline = [
        float(score)
        for seed in seeds for holdout in (0, 1)
        for score in cells[(min(k_values), seed, holdout)]["baseline_scores"]
    ]
    all_scores = {0: baseline, **by_k}
    seed_means: dict[int, dict[int, float]] = {
        0: {
            seed: float(np.mean([
                float(score)
                for holdout in (0, 1)
                for score in cells[(min(k_values), seed, holdout)]["baseline_scores"]
            ]))
            for seed in seeds
        }
    }
    seed_means.update({
        k: {
            seed: float(np.mean([
                float(score)
                for holdout in (0, 1)
                for score in cells[(k, seed, holdout)]["scores"]
            ]))
            for seed in seeds
        }
        for k in k_values
    })
    means = {
        k: float(np.mean(list(values.values()))) for k, values in seed_means.items()
    }
    ses = {
        k: float(np.std(list(values.values()), ddof=1) / np.sqrt(len(values)))
        if len(values) > 1 else 0.0
        for k, values in seed_means.items()
    }
    support = {0: 1.0}
    for k in k_values:
        seed_pass = []
        for seed in seeds:
            runs = [cells[(k, seed, holdout)] for holdout in (0, 1)]
            seed_pass.append(all(
                bool(run["fit_platform"]["converged"])
                and bool(run["recovery_pass"])
                and all(bool(item["inference_platform"]["converged"]) for item in run["recovery"])
                for run in runs
            ))
        support[k] = float(np.mean(seed_pass))
    best_k = max(means, key=means.get)
    paired_differences = {
        k: [seed_means[k][seed] - seed_means[best_k][seed] for seed in seeds]
        for k in sorted(seed_means)
    }
    predictive_eligible = [
        k for k in sorted(means) if means[k] >= means[best_k] - ses[best_k]
    ]
    eligible = [k for k in predictive_eligible if support[k] >= 2 / 3]
    selected = min(eligible) if eligible else None
    status = "K0_SELECTED" if selected == 0 and len(seeds) >= 3 else "K0_NOT_CONFIRMED"
    return {
        "schema": "r04.k0_gate.v1",
        "status": status,
        "calibration_cell_count": len(cells),
        "seeds": seeds,
        "candidate_k": [0, *k_values],
        "independent_resampling_unit": "data_seed",
        "heldout_score_mean_by_seed_and_k": {
            str(k): {str(seed): value for seed, value in values.items()}
            for k, values in seed_means.items()
        },
        "paired_difference_to_best_by_k": {
            str(k): values for k, values in paired_differences.items()
        },
        "leakage_control": LEAKAGE_CONTROL,
        "execution_contract": validated["execution_contract"],
        "generator_contract": validated["generator_contract"],
        "input_data_hash_by_seed": {
            str(seed): value for seed, value in validated["input_hash_by_seed"].items()
        },
        "generator_replay_verified": validated["generator_replay_verified"],
        "heldout_score_mean_by_k": {str(k): value for k, value in means.items()},
        "heldout_score_se_by_k": {str(k): value for k, value in ses.items()},
        "support_fraction_by_k": {str(k): value for k, value in support.items()},
        "one_se_candidate": {
            "best_predictive_k": best_k,
            "predictive_eligible_k": predictive_eligible,
            "eligible_k": eligible,
            "selected_k": selected,
        },
        "leakage_control": "heldout_section_fields_adapted_on_disjoint_gene_folds",
        "formal_restart_authorized": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    shards = [json.loads(path.read_text(encoding="utf-8")) for path in args.input]
    result = evaluate_no_field_shards(shards)
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
    return 0 if result["status"] == "K0_SELECTED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
