"""Objective helpers shared by dense and gene-minibatch R-04 fits."""

from __future__ import annotations

from typing import Any

import torch

from .diagnostics import gene_batch_sum_scale


def gene_batch_sum_objective(
    per_spot_gene_loss: Any,
    total_genes: int,
    batch_genes: int,
) -> Any:
    """Estimate a dense per-spot sum from a sampled gene loss matrix.

    ``per_spot_gene_loss`` has shape ``(spots, sampled_genes)`` and contains
    positive losses (or negative log probabilities).  The operation sums over
    sampled genes, averages over spots, then applies the finite-population
    Horvitz--Thompson factor.  It is differentiable when the input is a
    torch.Tensor.
    """
    if isinstance(per_spot_gene_loss, torch.Tensor):
        if per_spot_gene_loss.dim() != 2:
            raise ValueError("per_spot_gene_loss must be a two-dimensional tensor")
        # sum over genes, mean over spots, then HT scale
        value = per_spot_gene_loss.sum(dim=1).mean()
        return gene_batch_sum_scale(value, total_genes, batch_genes)
    # Fallback for numpy arrays (used in diagnostics/tests)
    import numpy as np

    arr = np.asarray(per_spot_gene_loss)
    if arr.ndim != 2:
        raise ValueError("per_spot_gene_loss must be a two-dimensional tensor")
    value = arr.sum(axis=1).mean()
    return gene_batch_sum_scale(value, total_genes, batch_genes)
