#!/usr/bin/env python3
"""Check dense and gene-minibatch objective equivalence before calibration."""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np

from r04.objectives import gene_batch_sum_objective
from r04.runtime import atomic_json


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    import tensorflow as tf
    import tensorflow_probability as tfp

    counts = tf.constant(np.arange(1, 1 + 8 * 6, dtype=np.float32).reshape(8, 6))
    log_mu = tf.Variable(np.linspace(-0.2, 1.0, 8 * 6, dtype=np.float32).reshape(8, 6))
    theta = tf.constant(np.full((8, 6), 3.0, dtype=np.float32))

    def objective(value: object, indices: np.ndarray) -> object:
        selected_counts = tf.gather(counts, indices, axis=1)
        selected_mu = tf.gather(value, indices, axis=1)
        selected_theta = tf.gather(theta, indices, axis=1)
        dist = tfp.distributions.NegativeBinomial(
            total_count=selected_theta,
            logits=selected_mu - tf.math.log(selected_theta),
        )
        return gene_batch_sum_objective(-dist.log_prob(selected_counts), 6, len(indices))

    with tf.GradientTape() as tape:
        dense = objective(log_mu, np.arange(6, dtype=np.int32))
    dense_gradient = tape.gradient(dense, log_mu)
    full_batch = objective(log_mu, np.arange(6, dtype=np.int32))
    batch_values = []
    batch_gradients = []
    for index in combinations(range(6), 3):
        indices = np.asarray(index, dtype=np.int32)
        with tf.GradientTape() as tape:
            batch_objective = objective(log_mu, indices)
        batch_values.append(float(batch_objective.numpy()))
        batch_gradients.append(tape.gradient(batch_objective, log_mu).numpy())
    mean_batch_gradient = np.mean(batch_gradients, axis=0)
    objective_error = float(abs(float(dense.numpy()) - float(np.mean(batch_values))))
    gradient_error = float(np.max(np.abs(dense_gradient.numpy() - mean_batch_gradient)))
    result = {
        "schema": "r04.objective_oracle.v1",
        "dense_objective": float(dense.numpy()),
        "full_batch_objective": float(full_batch.numpy()),
        "mean_gene_batch_objective": float(np.mean(batch_values)),
        "objective_abs_error": objective_error,
        "mean_batch_gradient_max_abs_error": gradient_error,
        "gradient_finite": bool(np.isfinite(dense_gradient.numpy()).all()),
        "gradient_max_abs": float(np.max(np.abs(dense_gradient.numpy()))),
        "batch_count": len(batch_values),
        "fit_infer_dense_scaling": "per_spot_gene_sum",
        "status": "PASS_OBJECTIVE_ORACLE" if objective_error < 1e-5 and gradient_error < 1e-5 else "SCIENTIFIC_FAIL_OBJECTIVE_ORACLE",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(args.output, result)
    return 0 if result["status"] == "PASS_OBJECTIVE_ORACLE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
