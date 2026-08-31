"""Objective helpers shared by dense and gene-minibatch R-04 fits."""

from __future__ import annotations

from typing import Any

from .diagnostics import gene_batch_sum_scale


def gene_batch_sum_objective(
    per_spot_gene_loss: Any,
    total_genes: int,
    batch_genes: int,
    *,
    tf_module: Any | None = None,
) -> Any:
    """Estimate a dense per-spot sum from a sampled gene loss matrix.

    ``per_spot_gene_loss`` has shape ``(spots, sampled_genes)`` and contains
    positive losses (or negative log probabilities).  The operation sums over
    sampled genes, averages over spots, then applies the finite-population
    Horvitz--Thompson factor.  It is differentiable when the input is a
    TensorFlow tensor.
    """
    if len(getattr(per_spot_gene_loss, "shape", ())) != 2:
        raise ValueError("per_spot_gene_loss must be a two-dimensional tensor")
    tf = tf_module
    if tf is None:
        import tensorflow as tf

    return gene_batch_sum_scale(
        tf.reduce_mean(tf.reduce_sum(per_spot_gene_loss, axis=1)),
        total_genes,
        batch_genes,
    )
