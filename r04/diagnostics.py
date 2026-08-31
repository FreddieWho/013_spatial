"""Diagnostics for continuous-field factor identifiability.

These metrics describe shared loading directions and continuous field
directions.  They never turn spots into labels and are intentionally reported
as warnings rather than silently repairing a collapsed factorisation.
"""

from __future__ import annotations

import numpy as np


def gene_batch_mean_scale(value, total_genes: int, batch_genes: int):
    """Apply the Horvitz--Thompson scale to a per-gene batch mean.

    This helper is retained for diagnostics that explicitly estimate a sum of
    per-gene contributions.  It must not be applied to a scalar mean over
    both spots and genes: that would introduce an extra factor of
    ``batch_genes``.  Model likelihoods use :func:`gene_batch_sum_scale`.
    """
    if total_genes < 1 or batch_genes < 1 or batch_genes > total_genes:
        raise ValueError("invalid total or batch gene count")
    return value * (float(total_genes) / float(batch_genes))


def gene_batch_sum_scale(value, total_genes: int, batch_genes: int):
    """Scale a per-spot sum over a uniform gene batch.

    ``value`` must already be averaged over spots but summed over the sampled
    genes.  The result estimates the dense per-spot gene sum.  Keeping this
    operation separate from ``gene_batch_mean_scale`` prevents an easy but
    consequential confusion between a mean over ``(spot, gene)`` and a sum
    over genes.
    """
    if total_genes < 1 or batch_genes < 1 or batch_genes > total_genes:
        raise ValueError("invalid total or batch gene count")
    return value * (float(total_genes) / float(batch_genes))


def _normalise_columns(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    norms = np.linalg.norm(values, axis=0, keepdims=True)
    return values / np.maximum(norms, 1e-12)


def pairwise_cosine(values: np.ndarray, *, absolute: bool = False) -> np.ndarray:
    """Return a symmetric factor-by-factor cosine matrix."""
    cosine = _normalise_columns(values).T @ _normalise_columns(values)
    cosine = np.clip(cosine, -1.0, 1.0)
    if absolute:
        cosine = np.abs(cosine)
    np.fill_diagonal(cosine, 1.0)
    return cosine


def _center_columns(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or values.shape[0] < 1 or values.shape[1] < 1:
        raise ValueError("values must be a non-empty two-dimensional matrix")
    return values - values.mean(axis=0, keepdims=True)


def canonical_subspace_correlations(
    truth: np.ndarray,
    estimate: np.ndarray,
) -> np.ndarray:
    """Return canonical correlations between two continuous field spaces.

    Factor order, sign and rotations inside an equally supported subspace are
    nuisance symmetries.  This metric therefore tests the span of the fields,
    not whether an arbitrary factor number happened to match a planted column.
    """
    left = _center_columns(truth)
    right = _center_columns(estimate)
    if left.shape[0] != right.shape[0]:
        raise ValueError("truth and estimate must have the same number of rows")
    def _basis(values: np.ndarray) -> np.ndarray:
        u, singular, _ = np.linalg.svd(values, full_matrices=False)
        if not len(singular):
            return u[:, :0]
        tolerance = max(values.shape) * np.finfo(float).eps * singular[0]
        return u[:, singular > tolerance]

    q_left = _basis(left)
    q_right = _basis(right)
    if q_left.shape[1] == 0 or q_right.shape[1] == 0:
        return np.empty(0, dtype=float)
    singular = np.linalg.svd(q_left.T @ q_right, compute_uv=False)
    return np.clip(singular, 0.0, 1.0)


def principal_angles(truth: np.ndarray, estimate: np.ndarray) -> np.ndarray:
    """Return principal angles in radians for two continuous field spaces."""
    return np.arccos(np.clip(canonical_subspace_correlations(truth, estimate), 0.0, 1.0))


def canonical_permutation_cutoffs(
    truth: np.ndarray,
    estimate: np.ndarray,
    *,
    groups: np.ndarray | None = None,
    seed: int,
    draws: int = 200,
    quantile: float = 0.95,
) -> np.ndarray:
    """Return direction-matched canonical-correlation permutation cutoffs.

    Rows are permuted within ``groups`` so a section-specific support cannot be
    exchanged with another section.  The first observed canonical direction is
    compared with the first null direction, the second with the second, and so
    on; a maximum first-direction null is not reused for every direction.
    """
    truth = np.asarray(truth, dtype=float)
    estimate = np.asarray(estimate, dtype=float)
    observed = canonical_subspace_correlations(truth, estimate)
    if draws < 1 or not 0.0 < quantile < 1.0:
        raise ValueError("draws and quantile must define a non-empty permutation null")
    if not len(observed):
        return np.empty(0, dtype=float)
    if groups is None:
        groups_array = np.zeros(len(estimate), dtype=int)
    else:
        groups_array = np.asarray(groups)
        if groups_array.ndim != 1 or len(groups_array) != len(estimate):
            raise ValueError("groups must have one entry per row")
    rng = np.random.default_rng(seed)
    null = np.empty((draws, len(observed)), dtype=float)
    group_indices = [np.flatnonzero(groups_array == value) for value in np.unique(groups_array)]
    for draw in range(draws):
        permutation = np.arange(len(estimate))
        for indices in group_indices:
            permutation[indices] = rng.permutation(indices)
        correlation = canonical_subspace_correlations(truth, estimate[permutation])
        if len(correlation) != len(observed):
            raise RuntimeError("permutation changed the estimable subspace rank")
        null[draw] = correlation
    return np.quantile(null, quantile, axis=0)


def field_effect_matrix(field_matrix: np.ndarray, loading: np.ndarray) -> np.ndarray:
    """Map continuous field values and gene effects to a spot-by-gene effect."""
    fields = _center_columns(field_matrix)
    loadings = np.asarray(loading, dtype=float)
    if loadings.ndim != 2 or fields.shape[1] != loadings.shape[1]:
        raise ValueError("field and loading factor dimensions differ")
    return fields @ loadings.T


def mnsf_spatial_effect_matrix(
    field_matrix: np.ndarray,
    loading: np.ndarray,
    amplitude: np.ndarray,
) -> np.ndarray:
    """Return the spot-by-gene spatial rate contribution of an mNSF model."""
    fields = np.asarray(field_matrix, dtype=float)
    loadings = np.asarray(loading, dtype=float)
    amplitudes = np.asarray(amplitude, dtype=float)
    if fields.ndim != 2 or loadings.ndim != 2 or amplitudes.ndim != 1:
        raise ValueError("fields, loading and amplitude must be matrices, matrix and vector")
    if fields.shape[1] != loadings.shape[1] or len(amplitudes) != fields.shape[1]:
        raise ValueError("mNSF factor dimensions differ")
    if np.any(loadings < 0) or np.any(amplitudes <= 0):
        raise ValueError("mNSF loading and amplitude must be non-negative")
    return np.exp(np.clip(fields, -30.0, 30.0)) @ (loadings * amplitudes).T


def effective_factor_count(
    heldout_score_by_k: dict[int, float],
    score_se_by_k: dict[int, float],
    added_signal_by_k: dict[int, float],
    null_cutoff: float,
    stability_by_k: dict[int, float],
    *,
    minimum_stability: float = 2 / 3,
) -> dict[str, object]:
    """Select the smallest supported K using a one-standard-error rule.

    Scores are assumed to be larger-is-better held-out predictive scores.
    ``added_signal_by_k`` is the weakest newly added singular/effect direction
    for that K, compared with a coordinate-permutation null.  This function
    reports a decision; it does not manufacture a factor when the evidence is
    absent.
    """
    if not heldout_score_by_k:
        raise ValueError("at least one K score is required")
    keys = sorted(heldout_score_by_k)
    if set(keys) != set(score_se_by_k) or set(keys) != set(added_signal_by_k) or set(keys) != set(stability_by_k):
        raise ValueError("K score, uncertainty, signal and stability keys must agree")
    if any(k < 1 for k in keys):
        raise ValueError("K must be positive")
    best_k = max(keys, key=lambda k: heldout_score_by_k[k])
    tolerance = max(float(score_se_by_k[best_k]), 0.0)
    eligible = [
        k for k in keys
        if heldout_score_by_k[k] >= heldout_score_by_k[best_k] - tolerance
        and added_signal_by_k[k] > null_cutoff
        and stability_by_k[k] >= minimum_stability
    ]
    selected = min(eligible) if eligible else None
    return {
        "selected_k": selected,
        "best_predictive_k": best_k,
        "one_se_tolerance": tolerance,
        "eligible_k": eligible,
        "status": "K_SELECTED" if selected is not None else "K_NOT_IDENTIFIABLE",
    }


def platform_converged(
    losses: list[float] | np.ndarray,
    *,
    window: int = 5,
    relative_improvement: float = 0.001,
) -> bool:
    """Return whether the recent loss window has reached a plateau."""
    values = np.asarray(losses, dtype=float)
    if window < 2 or relative_improvement < 0:
        raise ValueError("invalid convergence parameters")
    if values.size < window or not np.isfinite(values[-window:]).all():
        return False
    recent = values[-window:]
    denominator = max(abs(float(recent[0])), 1e-12)
    improvement = float((recent[0] - recent[-1]) / denominator)
    return 0.0 <= improvement <= relative_improvement


def platform_convergence_summary(
    losses: list[float] | np.ndarray,
    *,
    window: int = 50,
    stable_windows: int = 2,
    relative_improvement: float = 0.001,
) -> dict[str, object]:
    """Assess a plateau from consecutive non-rising block means.

    ``stable_windows=2`` requires three consecutive blocks: both recent block
    transitions must improve by no more than the threshold.  A rising loss is
    never labelled as convergence.
    """
    values = np.asarray(losses, dtype=float)
    if window < 2 or stable_windows < 1 or relative_improvement < 0:
        raise ValueError("invalid convergence parameters")
    required = window * (stable_windows + 1)
    if values.size < required or not np.isfinite(values[-required:]).all():
        return {
            "converged": False,
            "required_steps": required,
            "observed_steps": int(values.size),
            "block_means": [],
            "relative_improvements": [],
        }
    recent = values[-required:]
    block_means = np.asarray([
        recent[index * window:(index + 1) * window].mean()
        for index in range(stable_windows + 1)
    ])
    improvements = np.asarray([
        (block_means[index] - block_means[index + 1])
        / max(abs(float(block_means[index])), 1e-12)
        for index in range(stable_windows)
    ])
    net_improvement = float(
        (block_means[0] - block_means[-1]) / max(abs(float(block_means[0])), 1e-12)
    )
    converged = bool(
        np.all(np.abs(improvements) <= relative_improvement)
        and net_improvement >= 0.0
    )
    return {
        "converged": converged,
        "required_steps": required,
        "observed_steps": int(values.size),
        "block_means": block_means.tolist(),
        "relative_improvements": improvements.tolist(),
        "net_relative_improvement": net_improvement,
    }


def deterministic_objective_platform_summary(
    trace: list[dict[str, object]],
    *,
    stable_windows: int = 2,
    relative_improvement: float = 0.001,
) -> dict[str, object]:
    """Assess consecutive fixed-draw dense-objective evaluations.

    Unlike the minibatch loss gate, each trace point is a full-panel objective
    evaluated with the same stateless Monte Carlo definition.  A signed
    relative improvement is retained for reporting (negative means the
    objective rose), while the platform criterion accepts either direction as
    long as the absolute relative change is no more than
    ``relative_improvement`` in each final ``stable_windows`` interval.
    """
    if stable_windows < 1 or relative_improvement < 0:
        raise ValueError("invalid deterministic objective convergence parameters")
    points: list[tuple[int, float]] = []
    for item in trace:
        try:
            step = int(item["step"])
            objective = float(item["objective_mean"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("invalid deterministic objective trace") from exc
        if not np.isfinite(objective):
            raise ValueError("deterministic objective trace contains non-finite values")
        points.append((step, objective))
    required = stable_windows + 1
    if len(points) < required:
        return {
            "converged": False,
            "criterion": "dense_full_panel_fixed_stateless_mc",
            "required_evaluations": required,
            "observed_evaluations": len(points),
            "steps": [step for step, _ in points],
            "objective_means": [objective for _, objective in points],
            "relative_improvements": [],
            "relative_improvement_threshold": relative_improvement,
        }
    recent = points[-required:]
    improvements = [
        (recent[index][1] - recent[index + 1][1])
        / max(abs(recent[index][1]), 1e-12)
        for index in range(stable_windows)
    ]
    converged = bool(
        all(abs(value) <= relative_improvement for value in improvements)
    )
    return {
        "converged": converged,
        "criterion": "dense_full_panel_fixed_stateless_mc",
        "required_evaluations": required,
        "observed_evaluations": len(points),
        "steps": [step for step, _ in recent],
        "objective_means": [objective for _, objective in recent],
        "relative_improvements": improvements,
        "relative_improvement_threshold": relative_improvement,
    }


def continuous_factor_diagnostics(
    loading: np.ndarray,
    field_matrix: np.ndarray | None = None,
    *,
    signed: bool = False,
    high_cosine: float = 0.98,
    inactive_energy: float = 0.01,
) -> dict[str, object]:
    """Summarise factor duplication, field similarity and effective rank."""
    loading = np.asarray(loading, dtype=float)
    if loading.ndim != 2 or loading.shape[1] < 1:
        raise ValueError("loading must be a two-dimensional gene-by-factor matrix")
    loading_cosine = pairwise_cosine(loading, absolute=signed)
    off_diagonal = loading_cosine[~np.eye(loading.shape[1], dtype=bool)]
    singular = np.linalg.svd(loading, compute_uv=False)
    singular_sq = singular * singular
    total = float(singular_sq.sum())
    probabilities = singular_sq / max(total, 1e-12)
    entropy_rank = float(np.exp(-np.sum(probabilities * np.log(np.maximum(probabilities, 1e-12)))))
    participation_rank = float(total * total / max(float(np.sum(singular_sq * singular_sq)), 1e-12))
    factor_energy = np.sum(loading * loading, axis=0)
    factor_energy = factor_energy / max(float(factor_energy.sum()), 1e-12)

    result: dict[str, object] = {
        "n_factors": int(loading.shape[1]),
        "max_offdiagonal_loading_cosine": float(np.max(off_diagonal)) if len(off_diagonal) else 0.0,
        "loading_cosine": loading_cosine.tolist(),
        "effective_rank_entropy": entropy_rank,
        "effective_rank_participation": participation_rank,
        "factor_energy": factor_energy.tolist(),
        "inactive_factor_ids": [int(i) for i, value in enumerate(factor_energy) if value < inactive_energy],
        "loading_collapse_warning": bool(len(off_diagonal) and np.max(off_diagonal) >= high_cosine),
    }
    if field_matrix is not None:
        fields = np.asarray(field_matrix, dtype=float)
        if fields.ndim != 2 or fields.shape[1] != loading.shape[1]:
            raise ValueError("field_matrix must have one column per loading factor")
        centered = fields - fields.mean(axis=0, keepdims=True)
        field_cosine = pairwise_cosine(centered, absolute=True)
        result.update({
            "max_offdiagonal_field_cosine": float(np.max(field_cosine[~np.eye(fields.shape[1], dtype=bool)]))
            if fields.shape[1] > 1 else 0.0,
            "field_cosine": field_cosine.tolist(),
            "field_collapse_warning": bool(
                fields.shape[1] > 1 and np.max(field_cosine[~np.eye(fields.shape[1], dtype=bool)]) >= high_cosine
            ),
        })
    result["factor_collapse_warning"] = bool(
        result["loading_collapse_warning"] or result.get("field_collapse_warning", False)
        or result["inactive_factor_ids"]
    )
    return result
