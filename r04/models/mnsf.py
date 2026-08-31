"""A compact variational multi-section NSF-style count model.

This implementation is intentionally explicit rather than hiding the model
behind a clustering API.  It uses shared non-negative gene loadings and
section-specific variational GP basis weights with a negative-binomial count
likelihood.  Large runs must use the pinned TensorFlow/TFP environment.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import time
from typing import Sequence

import numpy as np
from scipy.special import logsumexp

from ..runtime import atomic_json, config_fingerprint, seed_everything
from ..io_contract import sections_content_hash
from ..spatial import matern32_kernel, normalize_coordinates
from ..types import FieldFit, SectionData
from ..diagnostics import (
    continuous_factor_diagnostics,
    deterministic_objective_platform_summary,
)
from ..objectives import gene_batch_sum_objective


@dataclass(frozen=True)
class MNSFConfig:
    factors: int = 4
    inducing_points: int = 64
    nonspatial_rank: int = 1
    lengthscale: float = 2.0
    steps: int = 300
    learning_rate: float = 0.01
    learning_rate_final: float | None = None
    learning_rate_decay_steps: int | None = None
    posterior_draws: int = 200
    seed: int = 20260807
    ridge: float = 1e-5
    checkpoint_dir: str | None = None
    resume_checkpoint_dir: str | None = None
    checkpoint_steps: int = 500
    config_hash: str = ""
    environment_hash: str = ""
    gene_batch_size: int | None = None
    diagnostic_interval: int = 100
    evaluation_interval: int = 0
    evaluation_mc_draws: int = 4
    optimization_schedule: str = "joint"
    shared_steps: int = 0
    initialization: str = "rate_k_invariant"
    gp_parameterization: str = "inducing_values_v2"
    # ``auto`` keeps the historical eager path on CPU and uses the compiled
    # TensorFlow step when a visible GPU is present.  The mode is an execution
    # detail: it must not change the objective, seed schedule, or protocol.
    execution_mode: str = "auto"


class MissingTensorFlowBackend(RuntimeError):
    """Raised instead of silently substituting clustering or graph smoothing."""


def _diagonal_gp_kl(tf, loc, scale, prior_inverse, logdet_prior):
    """KL[q(u)||N(0,K)] for a diagonal q covariance, per factor jointly."""
    factors = tf.cast(tf.shape(loc)[1], tf.float32)
    dimension = tf.cast(tf.size(loc), tf.float32)
    quadratic = tf.reduce_sum((prior_inverse @ loc) * loc)
    trace = tf.reduce_sum(
        (scale * scale) * tf.linalg.diag_part(prior_inverse)[:, None]
    )
    logdet_q = tf.reduce_sum(tf.math.log(scale * scale))
    return 0.5 * (quadratic + trace - dimension + factors * logdet_prior - logdet_q)


def _require_tf():
    try:
        import tensorflow as tf
        import tensorflow_probability as tfp
    except ImportError as exc:
        raise MissingTensorFlowBackend(
            "TensorFlow and tensorflow_probability are required for the R04 NB field model"
        ) from exc
    # Let more than one explicitly isolated worker share a GPU without making
    # TensorFlow reserve all device memory at import time.  This is allocator
    # configuration only; it does not alter model arithmetic.
    for device in tf.config.list_physical_devices("GPU"):
        try:
            tf.config.experimental.set_memory_growth(device, True)
        except RuntimeError:
            # The runtime may already be initialized by an embedding process.
            # In that case TensorFlow has already fixed the allocator policy.
            pass
    return tf, tfp


def _resolve_execution_mode(tf, requested: str) -> str:
    """Resolve an explicit or automatic execution mode without changing math."""
    if requested not in {"auto", "eager", "compiled"}:
        raise ValueError("execution_mode must be 'auto', 'eager', or 'compiled'")
    if requested == "auto":
        return (
            "compiled"
            if tf.config.list_logical_devices("GPU")
            else "eager"
        )
    return requested


def _inducing(coords: np.ndarray, number: int) -> np.ndarray:
    order = np.lexsort((coords[:, 1], coords[:, 0]))
    indices = np.linspace(0, len(order) - 1, number).round().astype(int)
    return coords[order[indices]]


def _softplus_inverse(value: np.ndarray) -> np.ndarray:
    value = np.maximum(np.asarray(value, dtype=float), 1e-8)
    return value + np.log(-np.expm1(-value))


def _rate_initial_components(
    counts: np.ndarray,
    library: np.ndarray,
    factors: int,
    *,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build a positive, data-scale, K-invariant initial rate decomposition."""
    counts = np.asarray(counts, dtype=float)
    library = np.asarray(library, dtype=float)
    if counts.ndim != 2 or library.shape != (len(counts),) or factors < 0:
        raise ValueError("counts, library and factors have incompatible shapes")
    if np.any(counts < 0) or np.any(library <= 0):
        raise ValueError("counts and library must be non-negative and positive")
    rate = counts / library[:, None]
    gene_mean = np.maximum(rate.mean(axis=0), 1e-4)
    if factors == 0:
        # The null keeps the same count model and nuisance channel but has no
        # coordinate-dependent component.  Start at the observed gene rate so
        # K=0 is not penalised by a missing spatial mass at step zero.
        return (
            np.zeros((counts.shape[1], 0), dtype=float),
            np.empty((0,), dtype=float),
            gene_mean,
            np.full(counts.shape[1], 20.0, dtype=float),
        )
    baseline = np.maximum(0.25 * gene_mean, 1e-4)
    spatial_mass = max(float(gene_mean.sum() - baseline.sum()), 1e-4)
    amplitude = np.full(factors, spatial_mass / factors, dtype=float)
    rng = np.random.default_rng(seed)
    loading = rng.lognormal(mean=0.0, sigma=0.5, size=(counts.shape[1], factors))
    loading /= loading.sum(axis=0, keepdims=True)
    dispersion = np.full(counts.shape[1], 20.0, dtype=float)
    return loading, amplitude, baseline, dispersion


def _inducing_projection(
    coords: np.ndarray,
    inducing: np.ndarray,
    lengthscale: float,
    ridge: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return Kxz Kzz^-1, Kzz^-1 and log|Kzz| for inducing values."""
    if ridge <= 0:
        raise ValueError("ridge must be positive")
    kzz = matern32_kernel(inducing, inducing, lengthscale)
    kzz = kzz + np.eye(len(inducing), dtype=float) * ridge
    prior_inverse = np.linalg.inv(kzz)
    sign, logdet = np.linalg.slogdet(kzz)
    if sign <= 0 or not np.isfinite(logdet):
        raise ValueError("inducing covariance is not positive definite")
    cross = matern32_kernel(coords, inducing, lengthscale)
    return cross @ prior_inverse, prior_inverse, float(logdet)


def predict_log_mean_from_fields(
    sections: Sequence[SectionData],
    field_groups: Sequence[Sequence[FieldFit]],
    loading: np.ndarray,
    factor_amplitude: np.ndarray,
    gene_baseline: np.ndarray,
    *,
    nonspatial_component: np.ndarray | None = None,
    nonspatial_loading: np.ndarray | None = None,
) -> np.ndarray:
    """Evaluate the fitted continuous-field NB mean at observed coordinates."""
    if len(sections) != len(field_groups):
        raise ValueError("sections and field groups must have equal length")
    loading = np.asarray(loading, dtype=float)
    amplitude = np.asarray(factor_amplitude, dtype=float)
    baseline = np.asarray(gene_baseline, dtype=float)
    if loading.ndim != 2 or loading.shape[1] != len(amplitude):
        raise ValueError("loading and factor amplitude shapes differ")
    if loading.shape[0] != len(baseline):
        raise ValueError("loading and baseline gene dimensions differ")
    if np.any(loading < 0) or np.any(amplitude <= 0) or np.any(baseline <= 0):
        raise ValueError("mNSF parameters must be positive")
    if (nonspatial_component is None) != (nonspatial_loading is None):
        raise ValueError("nonspatial component and loading must be supplied together")
    if nonspatial_component is not None and nonspatial_loading is not None:
        nuisance = np.asarray(nonspatial_component, dtype=float)
        nuisance_loading = np.asarray(nonspatial_loading, dtype=float)
        if nuisance.ndim != 2 or nuisance_loading.ndim != 2:
            raise ValueError("nonspatial component and loading must be matrices")
        if nuisance_loading.shape[1] != nuisance.shape[1]:
            raise ValueError("nonspatial component dimensions differ")
        if np.any(nuisance <= 0) or np.any(nuisance_loading < 0):
            raise ValueError("nonspatial parameters must be non-negative")
    else:
        nuisance = nuisance_loading = None
    outputs: list[np.ndarray] = []
    offset = 0
    for section, group in zip(sections, field_groups):
        if len(group) != loading.shape[1]:
            raise ValueError("field group factor dimension differs from loading")
        fields = (
            np.column_stack([np.asarray(fit.field_mean, dtype=float) for fit in group])
            if group else np.empty((len(section.coords), 0), dtype=float)
        )
        library = (
            np.asarray(section.library_size, dtype=float)
            if section.library_size is not None
            else np.asarray(
                section.counts.toarray() if hasattr(section.counts, "toarray") else section.counts
            ).sum(axis=1).astype(float)
        )
        components: list[np.ndarray] = []
        if loading.shape[1]:
            components.append(logsumexp(
                np.log(loading)[None, :, :] + np.log(amplitude)[None, None, :] + fields[:, None, :],
                axis=2,
            ))
        if nuisance is not None and nuisance_loading is not None:
            nuisance_slice = nuisance[offset : offset + len(section.coords)]
            components.append(np.log(np.maximum(nuisance_slice @ nuisance_loading.T, 1e-8)))
        baseline_component = np.broadcast_to(
            np.log(baseline)[None, :], (len(section.coords), len(baseline))
        )
        components.append(baseline_component)
        log_mu = np.log(np.maximum(library, 1.0))[:, None] + logsumexp(
            np.stack(components, axis=2), axis=2
        )
        outputs.append(log_mu)
        offset += len(section.coords)
    return np.concatenate(outputs, axis=0)


class MNSFEstimator:
    model_id = "mNSF_NB_VI"

    def __init__(self, config: MNSFConfig = MNSFConfig()) -> None:
        self.config = config
        self.loading_: np.ndarray | None = None
        self.fields_: list[list[FieldFit]] = []
        self.diagnostics_: dict[str, object] = {}
        self.gene_id_: tuple[str, ...] | None = None
        self.factor_amplitude_: np.ndarray | None = None
        self.gene_baseline_: np.ndarray | None = None
        self.dispersion_: np.ndarray | None = None
        self.nonspatial_loading_: np.ndarray | None = None
        self.fit_input_hash_: str = ""
        self.fit_config_hash_: str = ""
        self.fit_environment_hash_: str = ""
        self.inference_nonspatial_component_: np.ndarray | None = None
        self.inference_diagnostics_: dict[str, object] = {}

    def fit(self, sections: Sequence[SectionData]) -> "MNSFEstimator":
        if not sections:
            raise ValueError("at least one section is required")
        tf, tfp = _require_tf()
        execution_mode = _resolve_execution_mode(tf, self.config.execution_mode)
        seed_everything(self.config.seed)
        config_hash = self.config.config_hash or config_fingerprint(
            self.config,
            exclude={
                "checkpoint_dir",
                "resume_checkpoint_dir",
                "checkpoint_steps",
                "config_hash",
                "environment_hash",
                "execution_mode",
            },
        )
        runtime_environment_hash = config_fingerprint({
            "python": tuple(__import__("sys").version_info[:3]),
            "tensorflow": getattr(tf, "__version__", "unknown"),
            "tensorflow_probability": getattr(tfp, "__version__", "unknown"),
        })
        environment_hash = self.config.environment_hash or runtime_environment_hash
        genes = sections[0].gene_id
        if any(section.gene_id != genes for section in sections):
            raise ValueError("all sections must share the frozen gene universe")
        input_hash = sections_content_hash(sections)
        counts = np.concatenate([section.counts.toarray() if hasattr(section.counts, "toarray") else np.asarray(section.counts) for section in sections], axis=0).astype(np.float32)
        if np.any(counts != np.floor(counts)) or np.any(counts < 0):
            raise ValueError("mNSF requires non-negative integer counts")
        if any(section.library_size is not None for section in sections) and not all(
            section.library_size is not None for section in sections
        ):
            raise ValueError("library_size metadata must be present for every section or none")
        library = (
            np.concatenate([np.asarray(section.library_size, dtype=np.float32) for section in sections])
            if all(section.library_size is not None for section in sections)
            else counts.sum(axis=1).astype(np.float32)
        )
        library = np.maximum(library, 1.0)
        log_library = np.log(library)
        basis_list: list[np.ndarray] = []
        k_inv_list: list[np.ndarray] = []
        logdet_prior_list: list[float] = []
        section_ranges: list[tuple[int, int]] = []
        offset = 0
        for section in sections:
            coords, _ = normalize_coordinates(section.coords)
            m = min(self.config.inducing_points, len(coords))
            inducing = _inducing(coords, m)
            basis, prior_inverse, logdet_prior = _inducing_projection(
                coords, inducing, self.config.lengthscale, self.config.ridge
            )
            basis_list.append(basis.astype(np.float32))
            k_inv_list.append(prior_inverse.astype(np.float32))
            logdet_prior_list.append(logdet_prior)
            section_ranges.append((offset, offset + len(coords)))
            offset += len(coords)
        basis_by_section = [tf.constant(value) for value in basis_list]
        prior_inverse_by_section = [tf.constant(value) for value in k_inv_list]
        logdet_prior_by_section = [
            tf.constant(value, dtype=tf.float32) for value in logdet_prior_list
        ]
        y = tf.constant(counts)
        log_l = tf.constant(log_library)
        n, g = counts.shape
        k = self.config.factors
        q = self.config.nonspatial_rank
        if k < 0 or q < 0:
            raise ValueError("factors must be non-negative and nonspatial_rank non-negative")
        if self.config.gp_parameterization != "inducing_values_v2":
            raise ValueError("unsupported GP parameterization")
        if self.config.initialization == "rate_k_invariant":
            loading0, amplitude0, baseline0, dispersion0 = _rate_initial_components(
                counts, library, k, seed=self.config.seed
            )
            raw_w = tf.Variable(np.log(loading0).astype(np.float32))
            raw_a = tf.Variable(_softplus_inverse(amplitude0).astype(np.float32))
            raw_b = tf.Variable(_softplus_inverse(baseline0).astype(np.float32))
            raw_theta = tf.Variable(_softplus_inverse(dispersion0).astype(np.float32))
        elif self.config.initialization == "legacy":
            raw_w = tf.Variable(tf.random.normal((g, k), stddev=0.05, seed=self.config.seed))
            raw_a = tf.Variable(tf.zeros((k,)))
            raw_b = tf.Variable(tf.fill((g,), -3.0))
            raw_theta = tf.Variable(tf.fill((g,), 3.0))
        else:
            raise ValueError("unsupported mNSF initialization")
        raw_v = tf.Variable(tf.random.normal((g, q), stddev=0.05, seed=self.config.seed + 1)) if q else None
        raw_h = tf.Variable(tf.fill((n, q), -4.0)) if q else None
        q_loc: list[tf.Variable] = []
        q_scale: list[tf.Variable] = []
        for index, section in enumerate(sections):
            m = basis_list[index].shape[1]
            q_loc.append(tf.Variable(tf.zeros((m, k))))
            q_scale.append(tf.Variable(tf.fill((m, k), -2.0)))
        if self.config.learning_rate <= 0:
            raise ValueError("learning_rate must be positive")
        if self.config.learning_rate_final is not None and self.config.learning_rate_final <= 0:
            raise ValueError("learning_rate_final must be positive")
        if self.config.learning_rate_final is not None and (
            self.config.learning_rate_decay_steps is None
            or self.config.learning_rate_decay_steps < 1
        ):
            raise ValueError(
                "learning_rate_decay_steps is required when learning_rate_final is set"
            )
        if self.config.diagnostic_interval < 0:
            raise ValueError("diagnostic_interval must be non-negative")
        if self.config.evaluation_interval < 0:
            raise ValueError("evaluation_interval must be non-negative")
        if self.config.evaluation_mc_draws < 1:
            raise ValueError("evaluation_mc_draws must be positive")
        if self.config.optimization_schedule not in {"joint", "staged_shared_first"}:
            raise ValueError(
                "optimization_schedule must be 'joint' or 'staged_shared_first'"
            )
        if self.config.shared_steps < 0:
            raise ValueError("shared_steps must be non-negative")
        optimizer = tf.keras.optimizers.Adam(float(self.config.learning_rate))

        def scheduled_learning_rate(step: int) -> float:
            if self.config.learning_rate_final is None:
                return float(self.config.learning_rate)
            decay_steps = int(self.config.learning_rate_decay_steps)
            fraction = min(max(float(step) / decay_steps, 0.0), 1.0)
            return float(
                self.config.learning_rate
                + fraction * (self.config.learning_rate_final - self.config.learning_rate)
            )

        global_step = tf.Variable(0, dtype=tf.int64, trainable=False, name="global_step")
        variables = [raw_w, raw_a, raw_b, raw_theta, *q_loc, *q_scale]
        if raw_v is not None and raw_h is not None:
            variables.extend([raw_v, raw_h])
        optimizer_variables = [
            variable for variable in variables
            if variable.shape.num_elements() not in (None, 0)
        ]
        # Keras 3 rejects a variable that appears for the first time after
        # the shared-only stage. Build all Adam slots up front so the later
        # transition to joint updates is a scheduling choice, not an
        # optimizer-state mutation.
        optimizer.build(optimizer_variables)

        def deterministic_objective() -> tuple[float, float]:
            """Dense, fixed-seed objective used only for convergence monitoring."""
            eval_losses: list[float] = []
            loading_eval = (
                tf.nn.softmax(raw_w, axis=0)
                if k else tf.zeros((g, 0), dtype=tf.float32)
            )
            amplitude_eval = tf.nn.softplus(raw_a) + 1e-5
            baseline_eval = tf.nn.softplus(raw_b) + 1e-5
            theta_eval = tf.nn.softplus(raw_theta) + 1e-3
            for draw in range(self.config.evaluation_mc_draws):
                kl_eval = tf.constant(0.0, dtype=tf.float32)
                samples_eval = []
                for index, loc in enumerate(q_loc):
                    scale = tf.nn.softplus(q_scale[index]) + 1e-4
                    epsilon = tf.random.stateless_normal(
                        tf.shape(loc),
                        seed=[self.config.seed + 7000, draw * 1000 + index + 1],
                    )
                    samples_eval.append(loc + scale * epsilon)
                    kl_eval += _diagonal_gp_kl(
                        tf,
                        loc,
                        scale,
                        prior_inverse_by_section[index],
                        logdet_prior_by_section[index],
                    )
                field_parts_eval = [
                    basis_by_section[index] @ samples_eval[index]
                    for index in range(len(sections))
                ]
                field_eval = tf.concat([
                    part - tf.reduce_mean(part, axis=0, keepdims=True)
                    for part in field_parts_eval
                ], axis=0)
                components_eval = []
                if k:
                    components_eval.append(
                        tf.reduce_logsumexp(
                            tf.math.log(loading_eval[None, :, :] + 1e-8)
                            + tf.math.log(amplitude_eval[None, None, :])
                            + field_eval[:, None, :],
                            axis=2,
                        )
                    )
                if q and raw_v is not None and raw_h is not None:
                    v_eval = tf.nn.softmax(raw_v, axis=0)
                    h_eval = tf.nn.softplus(raw_h) + 1e-5
                    components_eval.append(
                        tf.math.log(h_eval @ tf.transpose(v_eval) + 1e-8)
                    )
                components_eval.append(
                    tf.broadcast_to(
                        tf.math.log(baseline_eval[None, :]), [n, g]
                    )
                )
                log_mu_eval = log_l[:, None] + tf.reduce_logsumexp(
                    tf.stack(components_eval, axis=2), axis=2
                )
                distribution_eval = tfp.distributions.NegativeBinomial(
                    total_count=theta_eval[None, :],
                    logits=log_mu_eval - tf.math.log(theta_eval[None, :]),
                )
                likelihood_eval = gene_batch_sum_objective(
                    -distribution_eval.log_prob(y), g, g, tf_module=tf
                )
                eval_losses.append(float(
                    (likelihood_eval + 1e-4 * kl_eval / float(n)).numpy()
                ))
            return float(np.mean(eval_losses)), float(
                np.std(eval_losses, ddof=1) if len(eval_losses) > 1 else 0.0
            )
        # The first stage deliberately leaves the spatially varying blocks
        # untouched.  Global factor amplitudes remain active because they are
        # shared across spots; the objective and K=0 parameterization are
        # unchanged.  With K=0 this set is empty, so the schedule is exactly
        # the joint baseline.
        spatial_variable_ids = {
            id(variable)
            for variable in [raw_w, *q_loc, *q_scale]
            if variable.shape.num_elements() not in (None, 0)
        }

        def optimization_stage(step: int) -> str:
            if (
                self.config.optimization_schedule == "staged_shared_first"
                and k > 0
                and step < self.config.shared_steps
            ):
                return "shared_first"
            return "joint"

        best_loss_variable = tf.Variable(
            np.inf, dtype=tf.float32, trainable=False, name="best_loss"
        )
        best_step_variable = tf.Variable(
            -1, dtype=tf.int64, trainable=False, name="best_step"
        )
        best_variables = [
            tf.Variable(
                variable.numpy(),
                dtype=variable.dtype,
                trainable=False,
                name=f"best_variable_{index}",
            )
            for index, variable in enumerate(variables)
        ]
        manager = None
        resume_checkpoint_dir = (
            Path(self.config.resume_checkpoint_dir)
            if self.config.resume_checkpoint_dir is not None
            else None
        )
        if resume_checkpoint_dir is not None and self.config.checkpoint_dir is None:
            raise ValueError("resume_checkpoint_dir requires checkpoint_dir")
        if self.config.checkpoint_dir:
            checkpoint_dir = Path(self.config.checkpoint_dir)
            checkpoint_dir.mkdir(parents=True, exist_ok=True)
            metadata_path = checkpoint_dir / "checkpoint_metadata.json"
            metadata = {
                "checkpoint_schema": "r04.mnsf_checkpoint.v2",
                "input_hash": input_hash,
                "config_hash": config_hash,
                "environment_hash": environment_hash,
            }
            if metadata_path.exists():
                import json
                previous = json.loads(metadata_path.read_text(encoding="utf-8"))
                if previous != metadata:
                    raise RuntimeError("checkpoint hash mismatch; start a new R04 run")
            else:
                atomic_json(metadata_path, metadata)
            checkpoint_values = {
                "optimizer": optimizer,
                "global_step": global_step,
                "best_loss": best_loss_variable,
                "best_step": best_step_variable,
                "raw_w": raw_w,
                "raw_a": raw_a,
                "raw_b": raw_b,
                "raw_theta": raw_theta,
            }
            if raw_v is not None and raw_h is not None:
                checkpoint_values.update({"raw_v": raw_v, "raw_h": raw_h})
            checkpoint_values.update({f"q_loc_{i}": value for i, value in enumerate(q_loc)})
            checkpoint_values.update({f"q_scale_{i}": value for i, value in enumerate(q_scale)})
            checkpoint_values.update({
                f"best_variable_{index}": value
                for index, value in enumerate(best_variables)
            })
            checkpoint = tf.train.Checkpoint(**checkpoint_values)
            manager = tf.train.CheckpointManager(checkpoint, str(checkpoint_dir), max_to_keep=2)
            if resume_checkpoint_dir is not None:
                source_metadata_path = resume_checkpoint_dir / "checkpoint_metadata.json"
                if not source_metadata_path.exists():
                    raise RuntimeError("resume checkpoint metadata is missing")
                import json
                source_metadata = json.loads(source_metadata_path.read_text(encoding="utf-8"))
                if source_metadata.get("checkpoint_schema") != "r04.mnsf_checkpoint.v2":
                    raise RuntimeError("resume checkpoint schema mismatch")
                if source_metadata.get("input_hash") != input_hash:
                    raise RuntimeError("resume checkpoint input hash mismatch")
                if source_metadata.get("environment_hash") != environment_hash:
                    raise RuntimeError("resume checkpoint environment hash mismatch")
                source_checkpoint = tf.train.latest_checkpoint(str(resume_checkpoint_dir))
                if source_checkpoint is None:
                    raise RuntimeError("resume checkpoint is missing")
                checkpoint.restore(source_checkpoint).expect_partial()
            elif manager.latest_checkpoint:
                checkpoint.restore(manager.latest_checkpoint).expect_partial()
        best = float(best_loss_variable.numpy())
        losses: list[float] = []
        gradient_norm_trace: list[dict[str, object]] = []
        parameter_norm_trace: list[dict[str, object]] = []
        learning_rate_trace: list[dict[str, float]] = []
        optimization_stage_trace: list[dict[str, object]] = []
        evaluation_trace: list[dict[str, object]] = []
        start_step = int(global_step.numpy())

        def run_compiled_fit() -> None:
            """Run the same logical steps with TensorFlow graph dispatch.

            The Python loop remains outside the graph so checkpoint and
            deterministic-objective boundaries stay unchanged.  Each graph
            step computes the same pre-update loss/gradients as the eager
            reference, while best-state capture happens before the caller
            applies Adam, matching the historical checkpoint semantics.
            """
            trainable_variables = [
                variable
                for variable in variables
                if variable.shape.num_elements() not in (None, 0)
            ]
            spatial_ids = spatial_variable_ids
            group_specs = {
                "loading": [raw_w],
                "factor_amplitude": [raw_a],
                "baseline": [raw_b],
                "dispersion": [raw_theta],
                "gp_loc": q_loc,
                "gp_scale": q_scale,
                "nonspatial_loading": [raw_v] if raw_v is not None else [],
                "nonspatial_field": [raw_h] if raw_h is not None else [],
            }
            best_assignments = [
                (best_variable, variable)
                for best_variable, variable in zip(best_variables, variables)
                if variable.shape.num_elements() not in (None, 0)
            ]
            step_functions: dict[tuple[str, bool], object] = {}

            def make_step(stage: str, with_diagnostics: bool):
                @tf.function(reduce_retracing=True)
                def compiled_step(gene_index_tf, step_tensor):
                    with tf.GradientTape() as tape:
                        samples = []
                        kl = tf.constant(0.0, dtype=tf.float32)
                        for index, loc in enumerate(q_loc):
                            scale = tf.nn.softplus(q_scale[index]) + 1e-4
                            seed = tf.stack([
                                tf.constant(self.config.seed, dtype=tf.int32),
                                tf.cast(
                                    step_tensor + tf.cast(index + 1, tf.int64),
                                    tf.int32,
                                ),
                            ])
                            eps = tf.random.stateless_normal(tf.shape(loc), seed=seed)
                            samples.append(loc + scale * eps)
                            kl += _diagonal_gp_kl(
                                tf,
                                loc,
                                scale,
                                prior_inverse_by_section[index],
                                logdet_prior_by_section[index],
                            )
                        field_parts = [
                            basis_by_section[index] @ samples[index]
                            for index in range(len(sections))
                        ]
                        field = tf.concat([
                            part - tf.reduce_mean(part, axis=0, keepdims=True)
                            for part in field_parts
                        ], axis=0)
                        loading = (
                            tf.nn.softmax(raw_w, axis=0)
                            if k else tf.zeros((g, 0), dtype=tf.float32)
                        )
                        loading_batch = tf.gather(loading, gene_index_tf, axis=0)
                        raw_b_batch = tf.gather(raw_b, gene_index_tf)
                        raw_theta_batch = tf.gather(raw_theta, gene_index_tf)
                        components = []
                        if k:
                            log_components = (
                                tf.math.log(loading_batch[None, :, :] + 1e-8)
                                + tf.math.log(
                                    tf.nn.softplus(raw_a)[None, None, :] + 1e-5
                                )
                                + field[:, None, :]
                            )
                            components.append(
                                tf.reduce_logsumexp(log_components, axis=2)
                            )
                        if q and raw_v is not None and raw_h is not None:
                            v = tf.nn.softmax(raw_v, axis=0)
                            v_batch = tf.gather(v, gene_index_tf, axis=0)
                            components.append(
                                tf.math.log(
                                    (tf.nn.softplus(raw_h) + 1e-5)
                                    @ tf.transpose(v_batch)
                                    + 1e-8
                                )
                            )
                        components.append(
                            tf.broadcast_to(
                                tf.math.log(
                                    tf.nn.softplus(raw_b_batch)[None, :] + 1e-5
                                ),
                                [n, batch_size],
                            )
                        )
                        log_mu = log_l[:, None] + tf.reduce_logsumexp(
                            tf.stack(components, axis=2), axis=2
                        )
                        theta = tf.nn.softplus(raw_theta_batch)[None, :] + 1e-3
                        dist = tfp.distributions.NegativeBinomial(
                            total_count=theta,
                            logits=log_mu - tf.math.log(theta),
                        )
                        likelihood = gene_batch_sum_objective(
                            -dist.log_prob(tf.gather(y, gene_index_tf, axis=1)),
                            g,
                            batch_size,
                            tf_module=tf,
                        )
                        loss = likelihood + 1e-4 * kl / float(n)
                    gradients = tape.gradient(loss, trainable_variables)
                    if any(gradient is None for gradient in gradients):
                        raise RuntimeError(
                            "compiled mNSF step produced a disconnected gradient"
                        )

                    def persist_best():
                        best_loss_variable.assign(loss)
                        best_step_variable.assign(step_tensor + 1)
                        for best_variable, variable in best_assignments:
                            best_variable.assign(variable)
                        return tf.constant(0, dtype=tf.int32)

                    tf.cond(
                        loss < best_loss_variable,
                        persist_best,
                        lambda: tf.constant(0, dtype=tf.int32),
                    )
                    if with_diagnostics:
                        gradient_by_id = {
                            id(variable): gradient
                            for variable, gradient in zip(
                                trainable_variables, gradients
                            )
                        }

                        def squared_norm(values, source):
                            terms = [
                                tf.reduce_sum(tf.square(source[id(variable)]))
                                for variable in values
                                if id(variable) in source
                            ]
                            return (
                                tf.sqrt(tf.add_n(terms))
                                if terms
                                else tf.constant(0.0, dtype=tf.float32)
                            )

                        parameter_by_id = {
                            id(variable): variable for variable in variables
                        }
                        gradient_norms = tf.stack([
                            squared_norm(values, gradient_by_id)
                            for values in group_specs.values()
                        ])
                        parameter_norms = tf.stack([
                            squared_norm(values, parameter_by_id)
                            for values in group_specs.values()
                        ])
                        return loss, gradient_norms, parameter_norms, gradients
                    return loss, gradients

                return compiled_step

            loss_tensors = []
            for step in range(start_step, self.config.steps):
                optimizer.learning_rate.assign(scheduled_learning_rate(step))
                batch_size = (
                    g
                    if self.config.gene_batch_size is None
                    else int(self.config.gene_batch_size)
                )
                if batch_size < 1:
                    raise ValueError("gene_batch_size must be positive")
                batch_size = min(batch_size, g)
                if batch_size == g:
                    gene_index = np.arange(g, dtype=np.int32)
                else:
                    gene_index = np.random.default_rng(
                        self.config.seed + 9001 + step
                    ).choice(g, size=batch_size, replace=False).astype(np.int32)
                gene_index_tf = tf.constant(gene_index, dtype=tf.int32)
                stage = optimization_stage(step)
                diagnostic_due = (
                    self.config.diagnostic_interval > 0
                    and (
                        (step + 1) % self.config.diagnostic_interval == 0
                        or step == start_step
                        or step + 1 == self.config.steps
                    )
                )
                key = (stage, diagnostic_due)
                if key not in step_functions:
                    step_functions[key] = make_step(stage, diagnostic_due)
                step_result = step_functions[key](
                    gene_index_tf,
                    tf.constant(step, dtype=tf.int64),
                )
                if diagnostic_due:
                    value_tensor, gradient_norms, parameter_norms, gradients = (
                        step_result
                    )
                    gradient_values = gradient_norms.numpy()
                    parameter_values = parameter_norms.numpy()
                    gradient_norm_trace.append({
                        "step": step + 1,
                        "optimization_stage": stage,
                        **dict(zip(group_specs, map(float, gradient_values))),
                    })
                    parameter_norm_trace.append({
                        "step": step + 1,
                        "optimization_stage": stage,
                        **dict(zip(group_specs, map(float, parameter_values))),
                    })
                    learning_rate_trace.append({
                        "step": float(step + 1),
                        "learning_rate": float(optimizer.learning_rate.numpy()),
                    })
                    optimization_stage_trace.append({
                        "step": step + 1,
                        "optimization_stage": stage,
                    })
                else:
                    value_tensor, gradients = step_result
                if evaluation_due := (
                    self.config.evaluation_interval > 0
                    and (
                        (step + 1) % self.config.evaluation_interval == 0
                        or step == start_step
                        or step + 1 == self.config.steps
                    )
                ):
                    evaluation_mean, evaluation_sd = deterministic_objective()
                    evaluation_trace.append({
                        "step": step + 1,
                        "objective_mean": evaluation_mean,
                        "objective_sd": evaluation_sd,
                        "mc_draws": self.config.evaluation_mc_draws,
                        "objective_kind": "dense_full_panel_fixed_stateless_mc",
                    })
                active_pairs = [
                    (gradient, variable)
                    for gradient, variable in zip(gradients, trainable_variables)
                    if gradient is not None
                    and not (
                        stage == "shared_first" and id(variable) in spatial_ids
                    )
                ]
                optimizer.apply_gradients(active_pairs)
                global_step.assign(step + 1)
                loss_tensors.append(value_tensor)
                if (
                    manager is not None
                    and (step + 1) % max(1, self.config.checkpoint_steps) == 0
                ):
                    manager.save(checkpoint_number=step + 1)
            if loss_tensors:
                losses.extend(
                    np.asarray(tf.stack(loss_tensors).numpy(), dtype=float).tolist()
                )

        optimizer_started = time.perf_counter()
        if execution_mode == "compiled":
            run_compiled_fit()
            best = float(best_loss_variable.numpy())
        for step in (
            range(start_step, self.config.steps)
            if execution_mode == "eager"
            else ()
        ):
            optimizer.learning_rate.assign(scheduled_learning_rate(step))
            batch_size = g if self.config.gene_batch_size is None else int(self.config.gene_batch_size)
            if batch_size < 1:
                raise ValueError("gene_batch_size must be positive")
            batch_size = min(batch_size, g)
            if batch_size == g:
                gene_index = np.arange(g, dtype=np.int32)
            else:
                # A stateless per-step sampler is reproducible after a
                # restart and makes the sampled objective independent of
                # TensorFlow's global RNG state.
                gene_index = np.random.default_rng(self.config.seed + 9001 + step).choice(
                    g, size=batch_size, replace=False
                ).astype(np.int32)
            gene_index_tf = tf.constant(gene_index, dtype=tf.int32)
            with tf.GradientTape() as tape:
                samples = []
                kl = tf.constant(0.0, dtype=tf.float32)
                for index, loc in enumerate(q_loc):
                    scale = tf.nn.softplus(q_scale[index]) + 1e-4
                    eps = tf.random.stateless_normal(tf.shape(loc), seed=[self.config.seed, step + index + 1])
                    samples.append(loc + scale * eps)
                    inv = prior_inverse_by_section[index]
                    logdet = logdet_prior_by_section[index]
                    kl += _diagonal_gp_kl(tf, loc, scale, inv, logdet)
                field_parts = [
                    basis_by_section[index] @ samples[index]
                    for index in range(len(sections))
                ]
                field = tf.concat(field_parts, axis=0)
                field = tf.concat([
                    part - tf.reduce_mean(part, axis=0, keepdims=True)
                    for part in field_parts
                ], axis=0)
                loading = (
                    tf.nn.softmax(raw_w, axis=0)
                    if k else tf.zeros((g, 0), dtype=tf.float32)
                )
                loading_batch = tf.gather(loading, gene_index_tf, axis=0)
                raw_b_batch = tf.gather(raw_b, gene_index_tf)
                raw_theta_batch = tf.gather(raw_theta, gene_index_tf)
                components = []
                if k:
                    log_components = tf.math.log(loading_batch[None, :, :] + 1e-8) + tf.math.log(tf.nn.softplus(raw_a)[None, None, :] + 1e-5) + field[:, None, :]
                    components.append(tf.reduce_logsumexp(log_components, axis=2))
                if q and raw_v is not None and raw_h is not None:
                    v = tf.nn.softmax(raw_v, axis=0)
                    v_batch = tf.gather(v, gene_index_tf, axis=0)
                    components.append(tf.math.log((tf.nn.softplus(raw_h) + 1e-5) @ tf.transpose(v_batch) + 1e-8))
                components.append(tf.broadcast_to(tf.math.log(tf.nn.softplus(raw_b_batch)[None, :] + 1e-5), [n, batch_size]))
                log_mu = log_l[:, None] + tf.reduce_logsumexp(tf.stack(components, axis=2), axis=2)
                theta = tf.nn.softplus(raw_theta_batch)[None, :] + 1e-3
                dist = tfp.distributions.NegativeBinomial(total_count=theta, logits=log_mu - tf.math.log(theta))
                # The likelihood is a sum over genes and an average over
                # spots.  Summing over the sampled genes before applying the
                # Horvitz--Thompson factor is essential; scaling a joint
                # (spot, gene) mean would add an extra batch-size factor.
                likelihood = gene_batch_sum_objective(
                    -dist.log_prob(tf.gather(y, gene_index_tf, axis=1)),
                    g,
                    batch_size,
                    tf_module=tf,
                )
                loss = likelihood + 1e-4 * kl / float(n)
            variables = [raw_w, raw_a, raw_b, raw_theta, *q_loc, *q_scale]
            if raw_v is not None and raw_h is not None:
                variables.extend([raw_v, raw_h])
            trainable_variables = [
                variable for variable in variables
                if variable.shape.num_elements() not in (None, 0)
            ]
            gradients = tape.gradient(loss, trainable_variables)
            value = float(loss.numpy())
            stage = optimization_stage(step)
            diagnostic_due = (
                self.config.diagnostic_interval > 0
                and (
                    (step + 1) % self.config.diagnostic_interval == 0
                    or step == start_step
                    or step + 1 == self.config.steps
                )
            )
            if diagnostic_due:
                groups = {
                    "loading": [raw_w],
                    "factor_amplitude": [raw_a],
                    "baseline": [raw_b],
                    "dispersion": [raw_theta],
                    "gp_loc": q_loc,
                    "gp_scale": q_scale,
                    "nonspatial_loading": [raw_v] if raw_v is not None else [],
                    "nonspatial_field": [raw_h] if raw_h is not None else [],
                }
                gradient_by_variable = {
                    id(variable): gradient
                    for variable, gradient in zip(trainable_variables, gradients)
                }

                def squared_norm(values: list[object], *, gradient: bool) -> float:
                    total = 0.0
                    for variable in values:
                        value_tensor = (
                            gradient_by_variable.get(id(variable))
                            if gradient
                            else variable
                        )
                        if value_tensor is not None:
                            total += float(tf.reduce_sum(tf.square(value_tensor)).numpy())
                    return float(np.sqrt(total))

                gradient_norm_trace.append({
                    "step": step + 1,
                    "optimization_stage": stage,
                    **{
                        name: squared_norm(values, gradient=True)
                        for name, values in groups.items()
                    },
                })
                parameter_norm_trace.append({
                    "step": step + 1,
                    "optimization_stage": stage,
                    **{
                        name: squared_norm(values, gradient=False)
                        for name, values in groups.items()
                    },
                })
                learning_rate_trace.append({
                    "step": float(step + 1),
                    "learning_rate": float(optimizer.learning_rate.numpy()),
                })
                optimization_stage_trace.append({
                    "step": step + 1,
                    "optimization_stage": stage,
                })
            evaluation_due = (
                self.config.evaluation_interval > 0
                and (
                    (step + 1) % self.config.evaluation_interval == 0
                    or step == start_step
                    or step + 1 == self.config.steps
                )
            )
            if evaluation_due:
                evaluation_mean, evaluation_sd = deterministic_objective()
                evaluation_trace.append({
                    "step": step + 1,
                    "objective_mean": evaluation_mean,
                    "objective_sd": evaluation_sd,
                    "mc_draws": self.config.evaluation_mc_draws,
                    "objective_kind": "dense_full_panel_fixed_stateless_mc",
                })
            is_best = value < best
            if is_best:
                # ``loss`` was evaluated at the pre-update parameters.  Copy
                # that exact state into TensorFlow shadows before Adam changes
                # it so loss_best and the checkpoint always identify one model.
                best = value
                best_loss_variable.assign(value)
                best_step_variable.assign(step + 1)
                for best_variable, variable in zip(best_variables, variables):
                    best_variable.assign(variable)
            active_pairs = [
                (gradient, variable)
                for gradient, variable in zip(gradients, trainable_variables)
                if gradient is not None
                and not (
                    stage == "shared_first"
                    and id(variable) in spatial_variable_ids
                )
            ]
            optimizer.apply_gradients(active_pairs)
            global_step.assign(step + 1)
            losses.append(value)
            if manager is not None and (step + 1) % max(1, self.config.checkpoint_steps) == 0:
                manager.save(checkpoint_number=step + 1)
        optimizer_wall_seconds = time.perf_counter() - optimizer_started
        if not np.isfinite(best):
            raise RuntimeError("mNSF optimizer produced no persisted best state")
        best_state = [variable.numpy().copy() for variable in best_variables]
        for variable, value in zip(variables, best_state):
            variable.assign(value)
        if manager is not None:
            manager.save(checkpoint_number=int(global_step.numpy()))
        loading = (
            tf.nn.softmax(raw_w, axis=0).numpy()
            if k else np.zeros((g, 0), dtype=float)
        )
        self.loading_ = loading
        self.fields_ = []
        for index, (r0, r1) in enumerate(section_ranges):
            loc = q_loc[index].numpy()
            scale = tf.nn.softplus(q_scale[index]).numpy() + 1e-4
            basis_np = basis_list[index]
            mean = basis_np @ loc
            mean -= mean.mean(axis=0, keepdims=True)
            rng = np.random.default_rng(self.config.seed + index)
            draws = np.stack([basis_np @ (loc + scale * rng.normal(size=loc.shape)) for _ in range(self.config.posterior_draws)], axis=0)
            draws -= draws.mean(axis=1, keepdims=True)
            self.fields_.append([
                FieldFit(self.model_id, f"MNSF_{factor+1:02d}", sections[index].section_id,
                         loading[:, factor], mean[:, factor], draws[:, :, factor].std(axis=0),
                         self.config.lengthscale, float(np.var(mean[:, factor])),
                         input_hash=input_hash, diagnostics={"uncertainty_kind": "variational_inducing_weights", "steps": self.config.steps, "patient_id": sections[index].patient_id})
                for factor in range(k)
            ])
        field_matrix = (
            np.column_stack([
                np.concatenate([group[factor].field_mean for group in self.fields_])
                for factor in range(k)
            ])
            if k else np.empty((len(counts), 0), dtype=float)
        )
        self.diagnostics_ = {
            "loss_final": losses[-1] if losses else None,
            "loss_best": best if np.isfinite(best) else None,
            "best_step": int(best_step_variable.numpy()),
            "loss_trace": [float(value) for value in losses],
            "steps": self.config.steps,
            "start_step": start_step,
            "optimizer_steps_this_call": self.config.steps - start_step,
            "optimizer_wall_seconds": optimizer_wall_seconds,
            "execution_mode": execution_mode,
            "learning_rate_final": self.config.learning_rate_final,
            "learning_rate_decay_steps": self.config.learning_rate_decay_steps,
            "learning_rate_trace": learning_rate_trace,
            "optimization_schedule": self.config.optimization_schedule,
            "shared_steps": self.config.shared_steps,
            "optimization_stage_trace": optimization_stage_trace,
            "evaluation_interval": self.config.evaluation_interval,
            "evaluation_mc_draws": self.config.evaluation_mc_draws,
            "evaluation_trace": evaluation_trace,
            "deterministic_objective_platform": deterministic_objective_platform_summary(
                evaluation_trace
            ),
            "diagnostic_interval": self.config.diagnostic_interval,
            "gradient_norm_trace": gradient_norm_trace,
            "parameter_norm_trace": parameter_norm_trace,
            "backend": "tensorflow_probability", "lengthscale_status": "FIXED_CONFIG_NOT_LEARNED",
            "uncertainty_status": "INDUCING_WEIGHTS_ONLY",
            "config_hash": config_hash,
            "environment_hash": environment_hash,
            "runtime_environment_hash": runtime_environment_hash,
            "environment_hash_override": self.config.environment_hash is not None,
            "gene_batch_size": self.config.gene_batch_size,
            "gene_batch_scaling": "uniform_without_replacement_ht",
            "global_spatial_kl_scaling": "unscaled_once_per_step",
            "gene_sampler": "stateless_per_step_seeded_by_seed_plus_step",
            "resume_checkpoint_dir": (
                str(resume_checkpoint_dir) if resume_checkpoint_dir is not None else None
            ),
            **(
                continuous_factor_diagnostics(loading, field_matrix, signed=False)
                if k else {"n_factors": 0, "effective_rank_entropy": 0.0,
                           "effective_rank_participation": 0.0,
                           "factor_energy": [], "loading_cosine": [],
                           "field_cosine": []}
            ),
        }
        self.gene_id_ = tuple(genes)
        self.factor_amplitude_ = tf.nn.softplus(raw_a).numpy() + 1e-5
        self.gene_baseline_ = tf.nn.softplus(raw_b).numpy() + 1e-5
        self.dispersion_ = tf.nn.softplus(raw_theta).numpy() + 1e-3
        self.nonspatial_loading_ = (
            tf.nn.softmax(raw_v, axis=0).numpy() if raw_v is not None else None
        )
        self.fit_input_hash_ = input_hash
        self.fit_config_hash_ = config_hash
        self.fit_environment_hash_ = environment_hash
        return self

    def _require_fitted(self) -> None:
        if (
            self.loading_ is None
            or self.gene_id_ is None
            or self.factor_amplitude_ is None
            or self.gene_baseline_ is None
            or self.dispersion_ is None
        ):
            raise RuntimeError("mNSF must be fitted before frozen inference")

    def frozen_state(self) -> dict[str, object]:
        """Return only the parameters needed for fixed-model inference."""
        self._require_fitted()
        return {
            "schema": "r04.frozen_model.v1",
            "model_id": self.model_id,
            "config": asdict(self.config),
            "gene_id": list(self.gene_id_ or ()),
            "fit_input_hash": self.fit_input_hash_,
            "fit_config_hash": self.fit_config_hash_,
            "fit_environment_hash": self.fit_environment_hash_,
            "parameters": {
                "loading": self.loading_.tolist(),
                "factor_amplitude": self.factor_amplitude_.tolist(),
                "gene_baseline": self.gene_baseline_.tolist(),
                "dispersion": self.dispersion_.tolist(),
                "nonspatial_loading": (
                    self.nonspatial_loading_.tolist()
                    if self.nonspatial_loading_ is not None else None
                ),
            },
            "diagnostics": self.diagnostics_,
        }

    @classmethod
    def from_frozen_state(cls, state: dict[str, object]) -> "MNSFEstimator":
        if state.get("schema") != "r04.frozen_model.v1" or state.get("model_id") != cls.model_id:
            raise ValueError("unsupported mNSF frozen model artifact")
        config = dict(state["config"])
        if config.get("gp_parameterization") != "inducing_values_v2":
            raise ValueError("legacy mNSF frozen model lacks the D-050 GP parameterization")
        estimator = cls(MNSFConfig(**config))
        parameters = dict(state["parameters"])
        estimator.gene_id_ = tuple(str(value) for value in state["gene_id"])
        estimator.loading_ = np.asarray(parameters["loading"], dtype=float)
        estimator.factor_amplitude_ = np.asarray(parameters["factor_amplitude"], dtype=float)
        estimator.gene_baseline_ = np.asarray(parameters["gene_baseline"], dtype=float)
        estimator.dispersion_ = np.asarray(parameters["dispersion"], dtype=float)
        nonspatial = parameters.get("nonspatial_loading")
        estimator.nonspatial_loading_ = (
            np.asarray(nonspatial, dtype=float) if nonspatial is not None else None
        )
        if (
            estimator.loading_.ndim != 2
            or estimator.loading_.shape[0] != len(estimator.gene_id_)
            or estimator.factor_amplitude_.shape != (estimator.loading_.shape[1],)
            or estimator.gene_baseline_.shape != (len(estimator.gene_id_),)
            or estimator.dispersion_.shape != (len(estimator.gene_id_),)
            or (
                estimator.config.nonspatial_rank
                and (
                    estimator.nonspatial_loading_ is None
                    or estimator.nonspatial_loading_.shape[0] != len(estimator.gene_id_)
                )
            )
        ):
            raise ValueError("mNSF frozen parameter shapes do not match the frozen gene universe")
        estimator.fit_input_hash_ = str(state.get("fit_input_hash", ""))
        estimator.fit_config_hash_ = str(state.get("fit_config_hash", ""))
        estimator.fit_environment_hash_ = str(state.get("fit_environment_hash", ""))
        estimator.diagnostics_ = dict(state.get("diagnostics", {}))
        return estimator

    def subset_to_genes(self, genes: Sequence[str]) -> "MNSFEstimator":
        """Project a frozen model onto genes observed in every validation section."""
        self._require_fitted()
        requested = tuple(genes)
        if not requested or len(set(requested)) != len(requested):
            raise ValueError("conditional gene panel must be non-empty and unique")
        positions = {gene: index for index, gene in enumerate(self.gene_id_ or ())}
        if any(gene not in positions for gene in requested):
            raise ValueError("conditional panel contains genes outside the frozen model")
        indices = np.asarray([positions[gene] for gene in requested], dtype=int)
        projected = MNSFEstimator(self.config)
        projected.gene_id_ = requested
        projected.loading_ = self.loading_[indices, :].copy()
        projected.factor_amplitude_ = self.factor_amplitude_.copy()
        projected.gene_baseline_ = self.gene_baseline_[indices].copy()
        projected.dispersion_ = self.dispersion_[indices].copy()
        projected.nonspatial_loading_ = (
            self.nonspatial_loading_[indices, :].copy()
            if self.nonspatial_loading_ is not None else None
        )
        projected.fit_input_hash_ = self.fit_input_hash_
        projected.fit_config_hash_ = self.fit_config_hash_
        projected.fit_environment_hash_ = self.fit_environment_hash_
        full_energy = np.sum(self.loading_ * self.loading_, axis=0)
        retained_energy = np.sum(projected.loading_ * projected.loading_, axis=0) / np.maximum(full_energy, 1e-12)
        full_rank = (
            int(np.linalg.matrix_rank(self.loading_))
            if self.loading_.shape[1] else 0
        )
        panel_rank = (
            int(np.linalg.matrix_rank(projected.loading_))
            if projected.loading_.shape[1] else 0
        )
        projected.diagnostics_ = {
            **self.diagnostics_,
            "conditional_panel": True,
            "conditional_panel_genes": len(requested),
            "retained_loading_energy": retained_energy.tolist(),
            "full_loading_rank": full_rank,
            "conditional_loading_rank": panel_rank,
            "conditional_rank_loss": panel_rank < full_rank,
        }
        projected.fields_ = []
        return projected

    def infer(
        self,
        sections: Sequence[SectionData],
        *,
        steps: int | None = None,
        posterior_draws: int | None = None,
    ) -> list[list[FieldFit]]:
        """Infer continuous fields for new sections with fitted parameters fixed."""
        if not sections:
            raise ValueError("at least one section is required")
        self._require_fitted()
        tf, tfp = _require_tf()
        execution_mode = _resolve_execution_mode(tf, self.config.execution_mode)
        if any(tuple(section.gene_id) != self.gene_id_ for section in sections):
            raise ValueError("new sections must use the frozen gene universe")
        self.inference_nonspatial_component_ = None
        self.inference_diagnostics_ = {}
        counts_by_section = [
            section.counts.toarray() if hasattr(section.counts, "toarray") else np.asarray(section.counts)
            for section in sections
        ]
        counts = np.concatenate(counts_by_section, axis=0).astype(np.float32)
        if np.any(counts != np.floor(counts)) or np.any(counts < 0):
            raise ValueError("frozen mNSF inference requires non-negative integer counts")
        if any(section.library_size is not None for section in sections) and not all(
            section.library_size is not None for section in sections
        ):
            raise ValueError("library_size metadata must be present for every section or none")
        library = np.maximum(
            np.concatenate([np.asarray(section.library_size, dtype=np.float32) for section in sections])
            if all(section.library_size is not None for section in sections)
            else counts.sum(axis=1),
            1.0,
        )
        basis_list: list[np.ndarray] = []
        k_inv_list: list[np.ndarray] = []
        logdet_prior_list: list[float] = []
        section_ranges: list[tuple[int, int]] = []
        offset = 0
        for section in sections:
            coords, _ = normalize_coordinates(section.coords)
            m = min(self.config.inducing_points, len(coords))
            inducing = _inducing(coords, m)
            basis, prior_inverse, logdet_prior = _inducing_projection(
                coords, inducing, self.config.lengthscale, self.config.ridge
            )
            basis_list.append(basis.astype(np.float32))
            k_inv_list.append(prior_inverse.astype(np.float32))
            logdet_prior_list.append(logdet_prior)
            section_ranges.append((offset, offset + len(coords)))
            offset += len(coords)
        basis_by_section = [tf.constant(value) for value in basis_list]
        prior_inverse_by_section = [tf.constant(value) for value in k_inv_list]
        logdet_prior_by_section = [
            tf.constant(value, dtype=tf.float32) for value in logdet_prior_list
        ]
        k = self.loading_.shape[1]
        g = counts.shape[1]
        q = self.config.nonspatial_rank
        q_loc = [tf.Variable(tf.zeros((basis.shape[1], k))) for basis in basis_list]
        q_scale = [tf.Variable(tf.fill((basis.shape[1], k), -2.0)) for basis in basis_list]
        raw_h = tf.Variable(tf.zeros((len(counts), q))) if q else None
        optimizer = tf.keras.optimizers.Adam(self.config.learning_rate)
        n_steps = self.config.steps if steps is None else int(steps)
        if n_steps < 1:
            raise ValueError("inference steps must be positive")
        seed_everything(self.config.seed + 5000)
        y = tf.constant(counts)
        log_l = tf.constant(np.log(library), dtype=tf.float32)
        loading = tf.constant(self.loading_, dtype=tf.float32)
        factor_amplitude = tf.constant(self.factor_amplitude_, dtype=tf.float32)
        gene_baseline = tf.constant(self.gene_baseline_, dtype=tf.float32)
        dispersion = tf.constant(self.dispersion_, dtype=tf.float32)
        nonspatial_loading = (
            tf.constant(self.nonspatial_loading_, dtype=tf.float32)
            if self.nonspatial_loading_ is not None else None
        )
        variables = [*q_loc, *q_scale]
        if raw_h is not None:
            variables.append(raw_h)
        best = np.inf
        best_state = None
        losses: list[float] = []

        def run_compiled_inference() -> None:
            """Run inference steps in a graph while preserving step order."""
            trainable_variables = [
                variable
                for variable in variables
                if variable.shape.num_elements() not in (None, 0)
            ]
            optimizer.build(trainable_variables)
            best_loss_variable = tf.Variable(
                np.inf, dtype=tf.float32, trainable=False, name="inference_best_loss"
            )
            best_variables = [
                tf.Variable(
                    variable.numpy(),
                    dtype=variable.dtype,
                    trainable=False,
                    name=f"inference_best_variable_{index}",
                )
                for index, variable in enumerate(variables)
            ]
            best_assignments = [
                (best_variable, variable)
                for best_variable, variable in zip(
                    best_variables,
                    variables,
                )
                if variable.shape.num_elements() not in (None, 0)
            ]

            @tf.function(reduce_retracing=True)
            def compiled_step(step_tensor):
                with tf.GradientTape() as tape:
                    samples = []
                    kl = tf.constant(0.0, dtype=tf.float32)
                    for index, loc in enumerate(q_loc):
                        scale = tf.nn.softplus(q_scale[index]) + 1e-4
                        seed = tf.stack([
                            tf.constant(self.config.seed + 5000, dtype=tf.int32),
                            tf.cast(
                                step_tensor + tf.cast(index + 1, tf.int64),
                                tf.int32,
                            ),
                        ])
                        eps = tf.random.stateless_normal(tf.shape(loc), seed=seed)
                        samples.append(loc + scale * eps)
                        kl += _diagonal_gp_kl(
                            tf,
                            loc,
                            scale,
                            prior_inverse_by_section[index],
                            logdet_prior_by_section[index],
                        )
                    field_parts = [
                        basis_by_section[index] @ samples[index]
                        for index in range(len(sections))
                    ]
                    centered = [
                        part - tf.reduce_mean(part, axis=0, keepdims=True)
                        for part in field_parts
                    ]
                    field = tf.concat(centered, axis=0)
                    components = []
                    if k:
                        log_components = (
                            tf.math.log(loading[None, :, :] + 1e-8)
                            + tf.math.log(factor_amplitude[None, None, :])
                            + field[:, None, :]
                        )
                        components.append(tf.reduce_logsumexp(log_components, axis=2))
                    if q and raw_h is not None and nonspatial_loading is not None:
                        h = tf.nn.softplus(raw_h) + 1e-5
                        components.append(
                            tf.math.log(
                                h @ tf.transpose(nonspatial_loading) + 1e-8
                            )
                        )
                    components.append(
                        tf.broadcast_to(
                            tf.math.log(gene_baseline[None, :]), tf.shape(y)
                        )
                    )
                    log_mu = log_l[:, None] + tf.reduce_logsumexp(
                        tf.stack(components, axis=2), axis=2
                    )
                    dist = tfp.distributions.NegativeBinomial(
                        total_count=dispersion[None, :],
                        logits=log_mu - tf.math.log(dispersion[None, :]),
                    )
                    likelihood = gene_batch_sum_objective(
                        -dist.log_prob(y), g, g, tf_module=tf
                    )
                    loss = likelihood + 1e-4 * kl / float(max(len(counts), 1))
                gradients = tape.gradient(loss, trainable_variables)
                if any(gradient is None for gradient in gradients):
                    raise RuntimeError(
                        "compiled mNSF inference produced a disconnected gradient"
                    )

                optimizer.apply_gradients(zip(gradients, trainable_variables))

                def persist_best():
                    best_loss_variable.assign(loss)
                    # The reference inference path records the state after
                    # Adam's update, even though the comparison uses the
                    # pre-update loss.  Preserve that ordering exactly.
                    for best_variable, variable in best_assignments:
                        best_variable.assign(variable)
                    return tf.constant(0, dtype=tf.int32)

                tf.cond(
                    loss < best_loss_variable,
                    persist_best,
                    lambda: tf.constant(0, dtype=tf.int32),
                )
                return loss

            loss_tensors = [
                compiled_step(tf.constant(step, dtype=tf.int64))
                for step in range(n_steps)
            ]
            if loss_tensors:
                losses.extend(
                    np.asarray(tf.stack(loss_tensors).numpy(), dtype=float).tolist()
                )
            nonlocal best, best_state
            best = float(best_loss_variable.numpy())
            best_state = [
                variable.numpy().copy()
                for variable in best_variables
            ]

        if execution_mode == "compiled":
            run_compiled_inference()
        for step in range(n_steps) if execution_mode == "eager" else ():
            with tf.GradientTape() as tape:
                samples = []
                kl = tf.constant(0.0, dtype=tf.float32)
                for index, loc in enumerate(q_loc):
                    scale = tf.nn.softplus(q_scale[index]) + 1e-4
                    eps = tf.random.stateless_normal(tf.shape(loc), seed=[self.config.seed + 5000, step + index + 1])
                    samples.append(loc + scale * eps)
                    inv = prior_inverse_by_section[index]
                    logdet = logdet_prior_by_section[index]
                    kl += _diagonal_gp_kl(tf, loc, scale, inv, logdet)
                field_parts = [
                    basis_by_section[index] @ samples[index]
                    for index in range(len(sections))
                ]
                centered = [part - tf.reduce_mean(part, axis=0, keepdims=True) for part in field_parts]
                field = tf.concat(centered, axis=0)
                components = []
                if k:
                    log_components = (
                        tf.math.log(loading[None, :, :] + 1e-8)
                        + tf.math.log(factor_amplitude[None, None, :])
                        + field[:, None, :]
                    )
                    components.append(tf.reduce_logsumexp(log_components, axis=2))
                if q and raw_h is not None and nonspatial_loading is not None:
                    h = tf.nn.softplus(raw_h) + 1e-5
                    components.append(tf.math.log(h @ tf.transpose(nonspatial_loading) + 1e-8))
                components.append(tf.broadcast_to(tf.math.log(gene_baseline[None, :]), tf.shape(y)))
                log_mu = log_l[:, None] + tf.reduce_logsumexp(tf.stack(components, axis=2), axis=2)
                dist = tfp.distributions.NegativeBinomial(
                    total_count=dispersion[None, :],
                    logits=log_mu - tf.math.log(dispersion[None, :]),
                )
                likelihood = gene_batch_sum_objective(
                    -dist.log_prob(y), g, g, tf_module=tf
                )
                loss = likelihood + 1e-4 * kl / float(max(len(counts), 1))
            gradients = tape.gradient(loss, variables)
            optimizer.apply_gradients(zip(gradients, variables))
            value = float(loss.numpy())
            losses.append(value)
            if value < best:
                best = value
                best_state = [variable.numpy().copy() for variable in variables]
        if best_state is None:
            raise RuntimeError("frozen mNSF inference produced no state")
        for variable, value in zip(variables, best_state):
            variable.assign(value)
        self.inference_nonspatial_component_ = (
            (tf.nn.softplus(raw_h) + 1e-5).numpy()
            if raw_h is not None else None
        )
        self.inference_diagnostics_ = {
            "loss_best": float(best),
            "loss_trace": [float(value) for value in losses],
            "steps": n_steps,
            "nonspatial_rank": q,
            "execution_mode": execution_mode,
        }
        draws_requested = self.config.posterior_draws if posterior_draws is None else int(posterior_draws)
        if draws_requested < 1:
            raise ValueError("posterior_draws must be positive")
        input_hash = sections_content_hash(sections)
        result: list[list[FieldFit]] = []
        for index, section in enumerate(sections):
            loc = q_loc[index].numpy()
            scale = tf.nn.softplus(q_scale[index]).numpy() + 1e-4
            basis_np = basis_list[index]
            mean = basis_np @ loc
            mean -= mean.mean(axis=0, keepdims=True)
            rng = np.random.default_rng(self.config.seed + 5100 + index)
            draws = np.stack([
                basis_np @ (loc + scale * rng.normal(size=loc.shape))
                for _ in range(draws_requested)
            ], axis=0)
            draws -= draws.mean(axis=1, keepdims=True)
            result.append([
                FieldFit(
                    self.model_id,
                    f"MNSF_{factor + 1:02d}",
                    section.section_id,
                    self.loading_[:, factor],
                    mean[:, factor],
                    draws[:, :, factor].std(axis=0),
                    self.config.lengthscale,
                    float(np.var(mean[:, factor])),
                    input_hash=input_hash,
                    diagnostics={
                        "uncertainty_kind": "frozen_variational_inducing_weights",
                        "fit_input_hash": self.fit_input_hash_,
                        "inference_steps": n_steps,
                        "nonspatial_component": "optimized" if q else "absent",
                        "loss_best": best,
                        "loss_trace": [float(value) for value in losses],
                        "patient_id": section.patient_id,
                    },
                )
                for factor in range(k)
            ])
        return result
