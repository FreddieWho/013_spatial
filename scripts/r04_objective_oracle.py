#!/usr/bin/env python3
"""Check dense and gene-minibatch objective equivalence before calibration."""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import torch

from r04.objectives import gene_batch_sum_objective
from r04.runtime import atomic_json


def _nb_log_prob(y: torch.Tensor, mu: torch.Tensor, theta: torch.Tensor) -> torch.Tensor:
    """NB log prob with mean mu and dispersion theta (total_count)."""
    # Use double precision for the oracle to keep the exhaustive mean error <1e-5
    y_f = y.double()
    theta_f = theta.double().clamp(min=1e-8)
    mu_f = mu.double().clamp(min=1e-8)
    term = torch.lgamma(y_f + theta_f) - torch.lgamma(theta_f) - torch.lgamma(y_f + 1)
    term = term + theta_f * (torch.log(theta_f) - torch.log(theta_f + mu_f))
    term = term + y_f * (torch.log(mu_f) - torch.log(theta_f + mu_f))
    return term


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    counts = torch.tensor(np.arange(1, 1 + 8 * 6, dtype=np.float64).reshape(8, 6), dtype=torch.float64)
    log_mu = torch.nn.Parameter(torch.tensor(np.linspace(-0.2, 1.0, 8 * 6, dtype=np.float64).reshape(8, 6), dtype=torch.float64))
    theta = torch.tensor(np.full((8, 6), 3.0, dtype=np.float64), dtype=torch.float64)

    def objective(value: torch.Tensor, indices: np.ndarray) -> torch.Tensor:
        idx = torch.as_tensor(indices, dtype=torch.long)
        selected_counts = counts[:, idx]
        selected_log_mu = value[:, idx]
        selected_theta = theta[:, idx]
        selected_mu = torch.exp(selected_log_mu)
        lp = _nb_log_prob(selected_counts, selected_mu, selected_theta)
        return gene_batch_sum_objective(-lp, 6, len(indices))

    dense = objective(log_mu, np.arange(6, dtype=np.int32))
    dense_grad_t = torch.autograd.grad(dense, log_mu, retain_graph=False)[0]
    dense_gradient = dense_grad_t.detach().cpu().numpy()
    # Recompute full_batch for output parity (original computed separately)
    full_batch = objective(log_mu, np.arange(6, dtype=np.int32))
    batch_values = []
    batch_gradients = []
    for index in combinations(range(6), 3):
        indices = np.asarray(index, dtype=np.int32)
        batch_objective = objective(log_mu, indices)
        batch_values.append(float(batch_objective.detach().cpu().item()))
        g = torch.autograd.grad(batch_objective, log_mu, retain_graph=False)[0]
        batch_gradients.append(g.detach().cpu().numpy())
    mean_batch_gradient = np.mean(batch_gradients, axis=0)
    objective_error = float(abs(float(dense.detach().cpu().item()) - float(np.mean(batch_values))))
    gradient_error = float(np.max(np.abs(dense_gradient - mean_batch_gradient)))
    result = {
        "schema": "r04.objective_oracle.v1",
        "dense_objective": float(dense.detach().cpu().item()),
        "full_batch_objective": float(full_batch.detach().cpu().item()),
        "mean_gene_batch_objective": float(np.mean(batch_values)),
        "objective_abs_error": objective_error,
        "mean_batch_gradient_max_abs_error": gradient_error,
        "gradient_finite": bool(np.isfinite(dense_gradient).all()),
        "gradient_max_abs": float(np.max(np.abs(dense_gradient))),
        "batch_count": len(batch_values),
        "fit_infer_dense_scaling": "per_spot_gene_sum",
        "status": "PASS_OBJECTIVE_ORACLE" if objective_error < 1e-5 and gradient_error < 1e-5 else "SCIENTIFIC_FAIL_OBJECTIVE_ORACLE",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(args.output, result)
    return 0 if result["status"] == "PASS_OBJECTIVE_ORACLE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
