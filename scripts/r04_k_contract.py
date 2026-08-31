"""Validation shared by the positive-K and K=0 calibration aggregators."""

from __future__ import annotations

import json
from itertools import product

import numpy as np

from r04.io_contract import sections_content_hash
from r04.synthetic import make_identifiability_calibration_sections


LEAKAGE_CONTROL = "heldout_section_fields_adapted_on_disjoint_gene_folds"
EXECUTION_CONTRACT = {
    "steps": 600,
    "inference_steps": 250,
    "section_folds": 2,
    "gene_folds": 2,
    "null_draws": 50,
}
GENERATOR_CONTRACT = {
    "generator": "mnsf_correct",
    "domain": "disconnected",
    "sections": 3,
    "spots_per_section": 50,
    "genes": 24,
}


def _canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _truth_values_equivalent(observed: object, expected: object) -> bool:
    """Compare discrete metadata exactly and floating truth within roundoff."""
    if isinstance(observed, dict) and isinstance(expected, dict):
        return observed.keys() == expected.keys() and all(
            _truth_values_equivalent(observed[key], expected[key])
            for key in observed
        )
    if isinstance(observed, list) and isinstance(expected, list):
        return len(observed) == len(expected) and all(
            _truth_values_equivalent(left, right)
            for left, right in zip(observed, expected)
        )
    if isinstance(observed, bool) or isinstance(expected, bool):
        return observed is expected
    if isinstance(observed, float) or isinstance(expected, float):
        try:
            left = float(observed)
            right = float(expected)
        except (TypeError, ValueError):
            return False
        return bool(
            np.isfinite(left)
            and np.isfinite(right)
            and np.isclose(left, right, rtol=1e-12, atol=1e-12)
        )
    return observed == expected


def validate_k_shards(
    shards: list[dict[str, object]],
    *,
    expected_mode: str,
    expected_k_eff_true: int,
) -> dict[str, object]:
    """Fail closed unless every K is evaluated on the same data and folds."""
    if not shards:
        raise ValueError("at least one K shard is required")

    cells: dict[tuple[int, int, int], dict[str, object]] = {}
    truth_by_seed: dict[int, str] = {}
    truth_payload_by_seed: dict[int, dict[str, object]] = {}
    input_hash_by_seed: dict[int, str] = {}
    regenerated_by_seed: dict[int, tuple[str, dict[str, object]]] = {}
    k_values: set[int] = set()
    data_seeds: set[int] = set()

    for shard in shards:
        if shard.get("schema") != "r04.k_calibration.v2":
            raise ValueError("all inputs must be leakage-free v2 K shards")
        if (
            shard.get("status") != "K_CALIBRATION_COMPLETE"
            or shard.get("formal_restart_authorized") is not False
        ):
            raise ValueError("K shard completion status or authorization is invalid")
        if shard.get("leakage_control") != LEAKAGE_CONTROL:
            raise ValueError("K shard leakage-control contract is missing")
        if shard.get("generator_contract") != GENERATOR_CONTRACT:
            raise ValueError("K shard generator contract does not match calibration")
        for key, expected in EXECUTION_CONTRACT.items():
            if shard.get(key) != expected:
                raise ValueError(f"K shard {key} must equal {expected}")

        shard_k = shard.get("k_values")
        shard_seeds = shard.get("seeds")
        if not isinstance(shard_k, list) or len(shard_k) != 1:
            raise ValueError("each K shard must contain exactly one K value")
        if not isinstance(shard_seeds, list) or len(shard_seeds) != 1:
            raise ValueError("each K shard must contain exactly one data seed")
        k = int(shard_k[0])
        seed = int(shard_seeds[0])
        if k < 1:
            raise ValueError("spatial K values must be positive")
        k_values.add(k)
        data_seeds.add(seed)

        truth = shard.get("truth")
        if not isinstance(truth, dict):
            raise ValueError("K shard truth metadata is missing")
        if truth.get("mode") != expected_mode or truth.get("k_eff_true") != expected_k_eff_true:
            raise ValueError("K shard truth does not match the requested calibration")
        if truth.get("generator") != "mnsf_correct" or truth.get("domain") != "disconnected":
            raise ValueError("K shard truth generator or domain is inconsistent")
        truth_key = _canonical(truth)
        if seed not in regenerated_by_seed:
            regenerated_sections, _, regenerated_truth = (
                make_identifiability_calibration_sections(
                    mode=expected_mode,
                    domain=GENERATOR_CONTRACT["domain"],
                    generator=GENERATOR_CONTRACT["generator"],
                    sections=GENERATOR_CONTRACT["sections"],
                    spots_per_section=GENERATOR_CONTRACT["spots_per_section"],
                    genes=GENERATOR_CONTRACT["genes"],
                    seed=seed,
                )
            )
            regenerated_by_seed[seed] = (
                sections_content_hash(regenerated_sections),
                regenerated_truth,
            )
        expected_input_hash, expected_truth = regenerated_by_seed[seed]
        if not _truth_values_equivalent(truth, expected_truth):
            raise ValueError("truth metadata does not match the replayed generator")
        if seed in truth_by_seed and truth_by_seed[seed] != truth_key:
            raise ValueError("truth metadata changes across K for the same data seed")
        truth_by_seed[seed] = truth_key
        truth_payload_by_seed[seed] = truth

        input_hash = shard.get("input_data_hash")
        if not isinstance(input_hash, str) or len(input_hash) != 64:
            raise ValueError("K shard input-data hash is missing or malformed")
        if input_hash != expected_input_hash:
            raise ValueError("input-data hash does not match the replayed generator")
        if seed in input_hash_by_seed and input_hash_by_seed[seed] != input_hash:
            raise ValueError("input data change across K for the same data seed")
        input_hash_by_seed[seed] = input_hash

        runs = shard.get("runs")
        if not isinstance(runs, list) or not runs:
            raise ValueError("K shard contains no runs")
        run_holdouts: list[int] = []
        for run in runs:
            holdout = int(run["holdout"])
            run_holdouts.append(holdout)
            cell = (k, seed, holdout)
            if cell in cells:
                raise ValueError("duplicate K-by-seed-by-holdout calibration cell")
            if int(run.get("k_model", -1)) != k:
                raise ValueError("run K does not match shard K")
            if int(run.get("data_seed", -1)) != seed:
                raise ValueError("run data seed does not match shard seed")
            if int(run.get("optimization_seed", -1)) != seed + k * 100 + holdout:
                raise ValueError("run optimisation seed violates the frozen contract")
            expected_split_seed = seed + holdout * 1009
            if int(run.get("split_seed", -1)) != expected_split_seed:
                raise ValueError("gene split seed must be independent of K")

            scores = run.get("scores")
            baselines = run.get("baseline_scores")
            recovery = run.get("recovery")
            if not all(isinstance(value, list) and len(value) == 2 for value in (scores, baselines, recovery)):
                raise ValueError("each run must contain exactly two gene cross-fit folds")
            try:
                score_values = np.asarray(scores, dtype=float)
                baseline_values = np.asarray(baselines, dtype=float)
            except (TypeError, ValueError):
                raise ValueError("K shard scores must be numeric") from None
            if not np.isfinite(score_values).all() or not np.isfinite(baseline_values).all():
                raise ValueError("K shard scores must be finite")
            required = min(k, expected_k_eff_true)
            if int(run.get("required_recovered_directions", -1)) != required:
                raise ValueError("required recovered directions do not match K and truth")
            recovered_splits: list[bool] = []
            for split in recovery:
                adaptation = [int(value) for value in split.get("adaptation_genes", [])]
                evaluation = [int(value) for value in split.get("evaluation_genes", [])]
                if set(adaptation) & set(evaluation):
                    raise ValueError("adaptation and evaluation genes overlap")
                if set(adaptation) | set(evaluation) != set(range(GENERATOR_CONTRACT["genes"])):
                    raise ValueError("gene cross-fit does not cover the frozen gene universe")
                try:
                    correlations = np.asarray(
                        split.get("canonical_correlations", []), dtype=float
                    )
                    cutoffs = np.asarray(
                        split.get("direction_matched_null_95", []), dtype=float
                    )
                except (TypeError, ValueError):
                    raise ValueError("recovery diagnostics must be numeric") from None
                if (
                    len(correlations) < required
                    or len(cutoffs) < required
                    or not np.isfinite(correlations[:required]).all()
                    or not np.isfinite(cutoffs[:required]).all()
                ):
                    raise ValueError("recovery diagnostics are incomplete or non-finite")
                recovered_splits.append(bool(np.all(correlations[:required] > cutoffs[:required])))
                inference_platform = split.get("inference_platform")
                if (
                    not isinstance(inference_platform, dict)
                    or not isinstance(inference_platform.get("converged"), bool)
                ):
                    raise ValueError("inference platform status is invalid")
            expected_recovery_pass = all(recovered_splits)
            if (
                not isinstance(run.get("recovery_pass"), bool)
                or run.get("recovery_pass") != expected_recovery_pass
            ):
                raise ValueError("recovery_pass does not match numerical diagnostics")
            fit_platform = run.get("fit_platform")
            if (
                not isinstance(fit_platform, dict)
                or not isinstance(fit_platform.get("converged"), bool)
            ):
                raise ValueError("fit platform status is invalid")
            cells[cell] = run

        declared_holdouts = sorted(int(value) for value in shard.get("holdouts", []))
        if sorted(run_holdouts) != declared_holdouts:
            raise ValueError("shard holdout declaration does not match its runs")

    expected_cells = set(product(sorted(k_values), sorted(data_seeds), (0, 1)))
    if set(cells) != expected_cells:
        raise ValueError("K shards do not form a complete K-by-seed-by-holdout grid")

    for seed, holdout in product(sorted(data_seeds), (0, 1)):
        reference_splits = None
        reference_baseline = None
        for k in sorted(k_values):
            run = cells[(k, seed, holdout)]
            split_signature = [
                (
                    tuple(int(value) for value in split["adaptation_genes"]),
                    tuple(int(value) for value in split["evaluation_genes"]),
                )
                for split in run["recovery"]
            ]
            baseline = np.asarray(run["baseline_scores"], dtype=float)
            if reference_splits is None:
                reference_splits = split_signature
                reference_baseline = baseline
            elif split_signature != reference_splits:
                raise ValueError("evaluation gene folds change across K")
            elif not np.allclose(baseline, reference_baseline, rtol=0.0, atol=1e-12):
                raise ValueError("K=0 baseline scores change across K")

    return {
        "cells": cells,
        "k_values": sorted(k_values),
        "data_seeds": sorted(data_seeds),
        "truth_by_seed": truth_payload_by_seed,
        "input_hash_by_seed": input_hash_by_seed,
        "generator_replay_verified": True,
        "execution_contract": dict(EXECUTION_CONTRACT),
        "generator_contract": dict(GENERATOR_CONTRACT),
    }
