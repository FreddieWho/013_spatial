"""Deterministic, read-only diagnostics for frozen R-04 mNSF checkpoints."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Sequence

import numpy as np

from .models.mnsf import _diagonal_gp_kl, _inducing, _inducing_projection
from .spatial import normalize_coordinates
from .types import SectionData


_VALUE_SUFFIX = "/.ATTRIBUTES/VARIABLE_VALUE"


def _softplus(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    return np.logaddexp(0.0, values)


def _softmax(values: np.ndarray, axis: int) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    shifted = values - np.max(values, axis=axis, keepdims=True)
    weights = np.exp(shifted)
    return weights / np.sum(weights, axis=axis, keepdims=True)


def nuisance_gauge_summary(
    raw_b: np.ndarray,
    raw_v: np.ndarray | None,
    raw_h: np.ndarray | None,
) -> dict[str, object]:
    """Summarize the positive additive nuisance decomposition.

    For rank one, ``b + h @ v.T`` is unchanged by ``h -> h-c`` and
    ``b -> b+c*v`` while positivity permits the shift.  The summary reports
    the coordinates that can drift along this gauge direction; it does not
    modify the checkpoint.
    """
    if raw_v is None or raw_h is None or raw_v.size == 0:
        return {
            "rank": 0,
            "h_min": [],
            "h_mean": [],
            "h_max": [],
            "baseline_v_projection": [],
        }
    baseline = _softplus(raw_b) + 1e-5
    loading = _softmax(raw_v, axis=0)
    component = _softplus(raw_h) + 1e-5
    projection = np.sum(loading * baseline[:, None], axis=0) / np.maximum(
        np.sum(loading * loading, axis=0), 1e-12
    )
    return {
        "rank": int(component.shape[1]),
        "h_min": np.min(component, axis=0).tolist(),
        "h_mean": np.mean(component, axis=0).tolist(),
        "h_max": np.max(component, axis=0).tolist(),
        "baseline_v_projection": projection.tolist(),
    }


def _checkpoint_tensor(reader: object, name: str) -> np.ndarray | None:
    key = f"{name}{_VALUE_SUFFIX}"
    try:
        return np.asarray(reader.get_tensor(key))
    except (KeyError, ValueError):
        return None


def _indexed_checkpoint_tensors(reader: object, prefix: str) -> list[np.ndarray]:
    names = [
        name
        for name, _ in reader.get_variable_to_shape_map().items()
        if name.startswith(prefix) and name.endswith(_VALUE_SUFFIX)
    ]
    indexed: list[tuple[int, str]] = []
    pattern = re.compile(rf"^{re.escape(prefix)}_(\d+){re.escape(_VALUE_SUFFIX)}$")
    for name in names:
        match = pattern.match(name)
        if match:
            indexed.append((int(match.group(1)), name))
    return [np.asarray(reader.get_tensor(name)) for _, name in sorted(indexed)]


def evaluate_checkpoint(
    sections: Sequence[SectionData],
    checkpoint_path: str | Path,
    *,
    factors: int,
    inducing_points: int,
    lengthscale: float,
    ridge: float,
    mc_draws: int = 8,
    seed: int = 20260820,
) -> dict[str, object]:
    """Evaluate a frozen checkpoint on the exact training objective.

    The gene panel is always dense.  The only remaining stochasticity is the
    fixed number of variational GP draws, using stateless seeds.  For K=0 the
    result is fully deterministic.  This function never writes or restores a
    checkpoint and never changes model parameters.
    """
    if not sections:
        raise ValueError("at least one section is required")
    if mc_draws < 1:
        raise ValueError("mc_draws must be positive")
    import tensorflow as tf
    import tensorflow_probability as tfp

    reader = tf.train.load_checkpoint(str(checkpoint_path))
    raw_w = _checkpoint_tensor(reader, "raw_w")
    raw_a = _checkpoint_tensor(reader, "raw_a")
    raw_b = _checkpoint_tensor(reader, "raw_b")
    raw_theta = _checkpoint_tensor(reader, "raw_theta")
    raw_v = _checkpoint_tensor(reader, "raw_v")
    raw_h = _checkpoint_tensor(reader, "raw_h")
    if raw_w is None or raw_a is None or raw_b is None or raw_theta is None:
        raise ValueError(f"checkpoint is missing mNSF parameter tensors: {checkpoint_path}")

    counts = np.concatenate([
        section.counts.toarray() if hasattr(section.counts, "toarray")
        else np.asarray(section.counts)
        for section in sections
    ], axis=0).astype(np.float32)
    library = (
        np.concatenate([
            np.asarray(section.library_size, dtype=np.float32)
            for section in sections
        ])
        if all(section.library_size is not None for section in sections)
        else counts.sum(axis=1).astype(np.float32)
    )
    library = np.maximum(library, 1.0)
    basis_list: list[np.ndarray] = []
    prior_inverse_list: list[np.ndarray] = []
    logdet_list: list[float] = []
    section_ranges: list[tuple[int, int]] = []
    offset = 0
    for section in sections:
        coords, _ = normalize_coordinates(section.coords)
        number = min(inducing_points, len(coords))
        inducing = _inducing(coords, number)
        basis, prior_inverse, logdet = _inducing_projection(
            coords, inducing, lengthscale, ridge
        )
        basis_list.append(basis.astype(np.float32))
        prior_inverse_list.append(prior_inverse.astype(np.float32))
        logdet_list.append(float(logdet))
        section_ranges.append((offset, offset + len(coords)))
        offset += len(coords)

    q_loc = _indexed_checkpoint_tensors(reader, "q_loc")
    q_scale = _indexed_checkpoint_tensors(reader, "q_scale")
    if len(q_loc) != len(sections) or len(q_scale) != len(sections):
        raise ValueError("checkpoint GP tensor count does not match sections")
    k = int(factors)
    if raw_w.shape != (counts.shape[1], k):
        raise ValueError(
            f"checkpoint K/genes mismatch: raw_w={raw_w.shape}, expected={(counts.shape[1], k)}"
        )
    basis = [tf.constant(value) for value in basis_list]
    y = tf.constant(counts)
    log_l = tf.constant(np.log(library), dtype=tf.float32)
    raw_w_tf = tf.constant(raw_w, dtype=tf.float32)
    raw_a_tf = tf.constant(raw_a, dtype=tf.float32)
    raw_b_tf = tf.constant(raw_b, dtype=tf.float32)
    raw_theta_tf = tf.constant(raw_theta, dtype=tf.float32)
    raw_v_tf = tf.constant(raw_v, dtype=tf.float32) if raw_v is not None else None
    raw_h_tf = tf.constant(raw_h, dtype=tf.float32) if raw_h is not None else None
    loc_tf = [tf.constant(value, dtype=tf.float32) for value in q_loc]
    scale_tf = [tf.nn.softplus(tf.constant(value, dtype=tf.float32)) + 1e-4 for value in q_scale]
    loading = tf.nn.softmax(raw_w_tf, axis=0) if k else tf.zeros((counts.shape[1], 0))
    amplitude = tf.nn.softplus(raw_a_tf) + 1e-5
    baseline = tf.nn.softplus(raw_b_tf) + 1e-5
    theta = tf.nn.softplus(raw_theta_tf) + 1e-3
    losses: list[float] = []
    for draw in range(mc_draws):
        kl = tf.constant(0.0, dtype=tf.float32)
        samples = []
        for index, loc in enumerate(loc_tf):
            epsilon = tf.random.stateless_normal(
                tf.shape(loc), seed=[int(seed), int(draw * 1000 + index + 1)]
            )
            samples.append(loc + scale_tf[index] * epsilon)
            kl += _diagonal_gp_kl(
                tf,
                loc,
                scale_tf[index],
                tf.constant(prior_inverse_list[index]),
                tf.constant(logdet_list[index], dtype=tf.float32),
            )
        field_parts = [
            basis[index] @ samples[index]
            for index in range(len(samples))
        ]
        if field_parts:
            field = tf.concat([
                part - tf.reduce_mean(part, axis=0, keepdims=True)
                for part in field_parts
            ], axis=0)
        else:
            field = tf.zeros((counts.shape[0], 0), dtype=tf.float32)
        components = []
        if k:
            components.append(
                tf.math.log(loading[None, :, :] + 1e-8)
                + tf.math.log(amplitude[None, None, :])
                + field[:, None, :]
            )
            spatial = tf.reduce_logsumexp(components[0], axis=2)
            components = [spatial]
        if raw_v_tf is not None and raw_h_tf is not None and raw_h_tf.shape[1] > 0:
            nuisance_loading = tf.nn.softmax(raw_v_tf, axis=0)
            nuisance = (tf.nn.softplus(raw_h_tf) + 1e-5) @ tf.transpose(nuisance_loading)
            components.append(tf.math.log(nuisance + 1e-8))
        components.append(tf.broadcast_to(tf.math.log(baseline[None, :]), tf.shape(y)))
        log_mu = log_l[:, None] + tf.reduce_logsumexp(tf.stack(components, axis=2), axis=2)
        distribution = tfp.distributions.NegativeBinomial(
            total_count=theta[None, :],
            logits=log_mu - tf.math.log(theta[None, :]),
        )
        likelihood = tf.reduce_mean(tf.reduce_sum(-distribution.log_prob(y), axis=1))
        losses.append(float((likelihood + 1e-4 * kl / float(counts.shape[0])).numpy()))
    return {
        "checkpoint": str(checkpoint_path),
        "global_step": int(reader.get_tensor(f"global_step{_VALUE_SUFFIX}")),
        "mc_draws": mc_draws,
        "seed": seed,
        "objective_mean": float(np.mean(losses)),
        "objective_sd": float(np.std(losses, ddof=1)) if len(losses) > 1 else 0.0,
        "objective_draws": losses,
        "n_spots": int(counts.shape[0]),
        "n_genes": int(counts.shape[1]),
        "nuisance_gauge": nuisance_gauge_summary(raw_b, raw_v, raw_h),
    }
