"""Signed spatial GP on a negative-binomial nuisance residual."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from ..runtime import atomic_json, config_fingerprint, seed_everything
from ..io_contract import sections_content_hash
from ..spatial import normalize_coordinates
from ..types import FieldFit, SectionData
from ..diagnostics import continuous_factor_diagnostics
from ..objectives import gene_batch_sum_objective
from .mnsf import _diagonal_gp_kl, _inducing_projection


@dataclass(frozen=True)
class SignedResidualGPConfig:
    factors: int = 4
    inducing_points: int = 64
    lengthscale: float = 2.0
    steps: int = 250
    learning_rate: float = 0.01
    posterior_draws: int = 200
    seed: int = 20260807
    ridge: float = 1e-4
    checkpoint_dir: str | None = None
    checkpoint_steps: int = 500
    config_hash: str = ""
    environment_hash: str = ""
    gene_batch_size: int | None = None
    gp_parameterization: str = "inducing_values_v2"


class SignedResidualGPEstimator:
    model_id = "SIGNED_NB_RESIDUAL_GP"

    def __init__(self, config: SignedResidualGPConfig = SignedResidualGPConfig()) -> None:
        self.config = config
        self.loading_: np.ndarray | None = None
        self.fields_: list[list[FieldFit]] = []
        self.diagnostics_: dict[str, object] = {}
        self.reference_profile_: np.ndarray | None = None
        self.fit_input_hash_: str = ""
        self.fit_config_hash_: str = ""
        self.fit_environment_hash_: str = ""

    @staticmethod
    def _pearson_residual(
        counts: np.ndarray,
        section_sizes: Sequence[int] | None = None,
        profile: np.ndarray | None = None,
    ) -> np.ndarray:
        library = np.maximum(
            counts.sum(axis=1) if section_sizes is None else np.asarray(section_sizes, dtype=float),
            1.0,
        )
        normalized = counts / library[:, None]
        reference = normalized.mean(axis=0) if profile is None else np.asarray(profile, dtype=float)
        expected = reference[None, :] * library[:, None]
        # Conservative fixed overdispersion is only the nuisance stage; the
        # signed field likelihood below estimates residual noise.
        variance = expected + expected * expected / 20.0
        return (counts - expected) / np.sqrt(variance + 1e-6)

    def fit(self, sections: Sequence[SectionData]) -> "SignedResidualGPEstimator":
        if not sections:
            raise ValueError("at least one section is required")
        try:
            import tensorflow as tf
        except ImportError as exc:
            raise RuntimeError("TensorFlow is required for the signed residual GP") from exc
        seed_everything(self.config.seed)
        config_hash = self.config.config_hash or config_fingerprint(
            self.config, exclude={"checkpoint_dir", "checkpoint_steps", "config_hash", "environment_hash"}
        )
        environment_hash = self.config.environment_hash or config_fingerprint({
            "python": tuple(__import__("sys").version_info[:3]),
            "tensorflow": getattr(tf, "__version__", "unknown"),
        })
        genes = sections[0].gene_id
        if any(section.gene_id != genes for section in sections):
            raise ValueError("all sections must share the frozen gene universe")
        input_hash = sections_content_hash(sections)
        counts = np.concatenate([section.counts.toarray() if hasattr(section.counts, "toarray") else np.asarray(section.counts) for section in sections], axis=0).astype(np.float32)
        if any(section.library_size is not None for section in sections) and not all(
            section.library_size is not None for section in sections
        ):
            raise ValueError("library_size metadata must be present for every section or none")
        library = (
            np.concatenate([np.asarray(section.library_size, dtype=np.float32) for section in sections])
            if all(section.library_size is not None for section in sections)
            else counts.sum(axis=1)
        )
        self.reference_profile_ = (
            counts / np.maximum(library, 1.0)[:, None]
        ).mean(axis=0)
        residual = self._pearson_residual(counts, section_sizes=library, profile=self.reference_profile_).astype(np.float32)
        basis_list: list[np.ndarray] = []
        k_inv_list: list[np.ndarray] = []
        logdet_prior_list: list[float] = []
        section_ranges: list[tuple[int, int]] = []
        offset = 0
        for section in sections:
            coords, _ = normalize_coordinates(section.coords)
            m = min(self.config.inducing_points, len(coords))
            order = np.lexsort((coords[:, 1], coords[:, 0]))
            inducing = coords[order[np.linspace(0, len(coords) - 1, m).round().astype(int)]]
            basis, prior_inverse, logdet_prior = _inducing_projection(
                coords, inducing, self.config.lengthscale, self.config.ridge
            )
            basis_list.append(basis.astype(np.float32))
            k_inv_list.append(prior_inverse.astype(np.float32))
            logdet_prior_list.append(logdet_prior)
            section_ranges.append((offset, offset + len(coords)))
            offset += len(coords)
        x = tf.constant(residual)
        g = residual.shape[1]
        k = self.config.factors
        if self.config.gp_parameterization != "inducing_values_v2":
            raise ValueError("unsupported GP parameterization")
        raw_lambda = tf.Variable(tf.random.normal((g, k), stddev=0.05, seed=self.config.seed))
        q_loc = [tf.Variable(tf.zeros((basis.shape[1], k))) for basis in basis_list]
        q_scale = [tf.Variable(tf.fill((basis.shape[1], k), -2.0)) for basis in basis_list]
        optimizer = tf.keras.optimizers.Adam(self.config.learning_rate)
        global_step = tf.Variable(0, dtype=tf.int64, trainable=False, name="global_step")
        manager = None
        if self.config.checkpoint_dir:
            checkpoint_dir = Path(self.config.checkpoint_dir)
            checkpoint_dir.mkdir(parents=True, exist_ok=True)
            metadata_path = checkpoint_dir / "checkpoint_metadata.json"
            metadata = {"input_hash": input_hash, "config_hash": config_hash, "environment_hash": environment_hash}
            if metadata_path.exists():
                import json
                if json.loads(metadata_path.read_text(encoding="utf-8")) != metadata:
                    raise RuntimeError("checkpoint hash mismatch; start a new R04 run")
            else:
                atomic_json(metadata_path, metadata)
            checkpoint_values = {"optimizer": optimizer, "global_step": global_step, "raw_lambda": raw_lambda}
            checkpoint_values.update({f"q_loc_{i}": value for i, value in enumerate(q_loc)})
            checkpoint_values.update({f"q_scale_{i}": value for i, value in enumerate(q_scale)})
            checkpoint = tf.train.Checkpoint(**checkpoint_values)
            manager = tf.train.CheckpointManager(checkpoint, str(checkpoint_dir), max_to_keep=2)
            if manager.latest_checkpoint:
                checkpoint.restore(manager.latest_checkpoint).expect_partial()
        losses: list[float] = []
        variables = [raw_lambda, *q_loc, *q_scale]
        start_step = int(global_step.numpy())
        for step in range(start_step, self.config.steps):
            batch_size = g if self.config.gene_batch_size is None else int(self.config.gene_batch_size)
            if batch_size < 1:
                raise ValueError("gene_batch_size must be positive")
            batch_size = min(batch_size, g)
            if batch_size == g:
                gene_index = np.arange(g, dtype=np.int32)
            else:
                gene_index = np.random.default_rng(self.config.seed + 9001 + step).choice(
                    g, size=batch_size, replace=False
                ).astype(np.int32)
            gene_index_tf = tf.constant(gene_index, dtype=tf.int32)
            with tf.GradientTape() as tape:
                fields = []
                kl = tf.constant(0.0, dtype=tf.float32)
                for index, basis in enumerate(basis_list):
                    b = tf.constant(basis)
                    scale = tf.nn.softplus(q_scale[index]) + 1e-4
                    eps = tf.random.stateless_normal(tf.shape(q_loc[index]), seed=[self.config.seed + 7, step + index + 1])
                    u = q_loc[index] + scale * eps
                    value = b @ u
                    value -= tf.reduce_mean(value, axis=0, keepdims=True)
                    fields.append(value)
                    inv = tf.constant(k_inv_list[index])
                    logdet = tf.constant(logdet_prior_list[index], dtype=tf.float32)
                    kl += _diagonal_gp_kl(tf, q_loc[index], scale, inv, logdet)
                field = tf.concat(fields, axis=0)
                loading = tf.linalg.l2_normalize(raw_lambda, axis=0)
                loading_batch = tf.gather(loading, gene_index_tf, axis=0)
                predicted = field @ tf.transpose(loading_batch)
                residual_batch = tf.gather(x, gene_index_tf, axis=1)
                likelihood = gene_batch_sum_objective(
                    (residual_batch - predicted) ** 2,
                    g,
                    batch_size,
                    tf_module=tf,
                )
                loss = likelihood + 1e-4 * kl / float(max(len(residual), 1))
            gradients = tape.gradient(loss, variables)
            optimizer.apply_gradients(zip(gradients, variables))
            global_step.assign(step + 1)
            losses.append(float(loss.numpy()))
            if manager is not None and (step + 1) % max(1, self.config.checkpoint_steps) == 0:
                manager.save(checkpoint_number=step + 1)
        if manager is not None:
            manager.save(checkpoint_number=int(global_step.numpy()))
        loading = tf.linalg.l2_normalize(raw_lambda, axis=0).numpy()
        self.loading_ = loading
        self.fields_ = []
        for index, section in enumerate(sections):
            loc = q_loc[index].numpy()
            scale = tf.nn.softplus(q_scale[index]).numpy() + 1e-4
            mean = basis_list[index] @ loc
            mean -= mean.mean(axis=0, keepdims=True)
            rng = np.random.default_rng(self.config.seed + 100 + index)
            draws = np.stack([basis_list[index] @ (loc + scale * rng.normal(size=loc.shape)) for _ in range(self.config.posterior_draws)], axis=0)
            draws -= draws.mean(axis=1, keepdims=True)
            self.fields_.append([
                FieldFit(self.model_id, f"SIGNED_{factor+1:02d}", section.section_id,
                         loading[:, factor], mean[:, factor], draws[:, :, factor].std(axis=0),
                         self.config.lengthscale, float(np.var(mean[:, factor])),
                         input_hash=input_hash, diagnostics={"uncertainty_kind": "variational_inducing_weights", "nuisance": "NB_PEARSON", "steps": self.config.steps, "patient_id": section.patient_id})
                for factor in range(k)
            ])
        field_matrix = np.column_stack([
            np.concatenate([group[factor].field_mean for group in self.fields_])
            for factor in range(k)
        ])
        self.diagnostics_ = {
            "backend": "tensorflow", "steps": self.config.steps, "nuisance": "NB_PEARSON",
            "likelihood": "MSE_ON_NB_PEARSON_RESIDUAL", "nuisance_overdispersion": 20.0,
            "lengthscale_status": "FIXED_CONFIG_NOT_LEARNED",
            "config_hash": config_hash, "environment_hash": environment_hash,
            "loss_final": losses[-1] if losses else None,
            "loss_best": min(losses) if losses else None,
            "loss_trace": [float(value) for value in losses],
            "gene_batch_size": self.config.gene_batch_size,
            "gene_batch_scaling": "uniform_without_replacement_ht",
            "global_spatial_kl_scaling": "unscaled_once_per_step",
            "gp_kl": "diagonal_variational_inducing_weights",
            "gene_sampler": "stateless_per_step_seeded_by_seed_plus_step",
            **continuous_factor_diagnostics(loading, field_matrix, signed=True),
        }
        self.fit_input_hash_ = input_hash
        self.fit_config_hash_ = config_hash
        self.fit_environment_hash_ = environment_hash
        self.gene_id_ = tuple(genes)
        return self

    def _require_fitted(self) -> None:
        if self.loading_ is None or self.reference_profile_ is None:
            raise RuntimeError("signed residual GP must be fitted before frozen inference")

    def frozen_state(self) -> dict[str, object]:
        self._require_fitted()
        return {
            "schema": "r04.frozen_model.v1",
            "model_id": self.model_id,
            "config": asdict(self.config),
            "gene_id": list(self.gene_id_) if hasattr(self, "gene_id_") else [],
            "fit_input_hash": self.fit_input_hash_,
            "fit_config_hash": self.fit_config_hash_,
            "fit_environment_hash": self.fit_environment_hash_,
            "parameters": {
                "loading": self.loading_.tolist(),
                "reference_profile": self.reference_profile_.tolist(),
            },
            "diagnostics": self.diagnostics_,
        }

    @classmethod
    def from_frozen_state(cls, state: dict[str, object]) -> "SignedResidualGPEstimator":
        if state.get("schema") != "r04.frozen_model.v1" or state.get("model_id") != cls.model_id:
            raise ValueError("unsupported signed residual GP frozen model artifact")
        config = dict(state["config"])
        if config.get("gp_parameterization") != "inducing_values_v2":
            raise ValueError("legacy signed frozen model lacks the D-050 GP parameterization")
        estimator = cls(SignedResidualGPConfig(**config))
        estimator.gene_id_ = tuple(str(value) for value in state.get("gene_id", ()))
        parameters = dict(state["parameters"])
        estimator.loading_ = np.asarray(parameters["loading"], dtype=float)
        estimator.reference_profile_ = np.asarray(parameters["reference_profile"], dtype=float)
        if (
            estimator.loading_.ndim != 2
            or estimator.loading_.shape[0] != len(estimator.gene_id_)
            or estimator.reference_profile_.shape != (len(estimator.gene_id_),)
        ):
            raise ValueError("signed frozen parameter shapes do not match the frozen gene universe")
        estimator.fit_input_hash_ = str(state.get("fit_input_hash", ""))
        estimator.fit_config_hash_ = str(state.get("fit_config_hash", ""))
        estimator.fit_environment_hash_ = str(state.get("fit_environment_hash", ""))
        estimator.diagnostics_ = dict(state.get("diagnostics", {}))
        return estimator

    def subset_to_genes(self, genes: Sequence[str]) -> "SignedResidualGPEstimator":
        """Project a frozen signed model onto a common observed validation panel."""
        self._require_fitted()
        requested = tuple(genes)
        if not requested or len(set(requested)) != len(requested):
            raise ValueError("conditional gene panel must be non-empty and unique")
        positions = {gene: index for index, gene in enumerate(getattr(self, "gene_id_", ())) }
        if any(gene not in positions for gene in requested):
            raise ValueError("conditional panel contains genes outside the frozen model")
        indices = np.asarray([positions[gene] for gene in requested], dtype=int)
        projected = SignedResidualGPEstimator(self.config)
        projected.gene_id_ = requested
        projected.loading_ = self.loading_[indices, :].copy()
        projected.reference_profile_ = self.reference_profile_[indices].copy()
        projected.fit_input_hash_ = self.fit_input_hash_
        projected.fit_config_hash_ = self.fit_config_hash_
        projected.fit_environment_hash_ = self.fit_environment_hash_
        full_energy = np.sum(self.loading_ * self.loading_, axis=0)
        retained_energy = np.sum(projected.loading_ * projected.loading_, axis=0) / np.maximum(full_energy, 1e-12)
        full_rank = int(np.linalg.matrix_rank(self.loading_))
        panel_rank = int(np.linalg.matrix_rank(projected.loading_))
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
        """Infer signed continuous residual fields with loadings fixed."""
        if not sections:
            raise ValueError("at least one section is required")
        self._require_fitted()
        try:
            import tensorflow as tf
        except ImportError as exc:
            raise RuntimeError("TensorFlow is required for frozen signed residual inference") from exc
        if getattr(self, "gene_id_", ()) and any(tuple(section.gene_id) != self.gene_id_ for section in sections):
            raise ValueError("new sections must use the frozen gene universe")
        counts = np.concatenate([
            section.counts.toarray() if hasattr(section.counts, "toarray") else np.asarray(section.counts)
            for section in sections
        ], axis=0).astype(np.float32)
        if np.any(counts != np.floor(counts)) or np.any(counts < 0):
            raise ValueError("frozen signed inference requires non-negative integer counts")
        if any(section.library_size is not None for section in sections) and not all(
            section.library_size is not None for section in sections
        ):
            raise ValueError("library_size metadata must be present for every section or none")
        library = (
            np.concatenate([np.asarray(section.library_size, dtype=np.float32) for section in sections])
            if all(section.library_size is not None for section in sections)
            else counts.sum(axis=1)
        )
        residual = self._pearson_residual(
            counts, section_sizes=library, profile=self.reference_profile_
        ).astype(np.float32)
        basis_list: list[np.ndarray] = []
        k_inv_list: list[np.ndarray] = []
        logdet_prior_list: list[float] = []
        for section in sections:
            coords, _ = normalize_coordinates(section.coords)
            m = min(self.config.inducing_points, len(coords))
            order = np.lexsort((coords[:, 1], coords[:, 0]))
            inducing = coords[order[np.linspace(0, len(coords) - 1, m).round().astype(int)]]
            basis, prior_inverse, logdet_prior = _inducing_projection(
                coords, inducing, self.config.lengthscale, self.config.ridge
            )
            basis_list.append(basis.astype(np.float32))
            k_inv_list.append(prior_inverse.astype(np.float32))
            logdet_prior_list.append(logdet_prior)
        k = self.loading_.shape[1]
        q_loc = [tf.Variable(tf.zeros((basis.shape[1], k))) for basis in basis_list]
        q_scale = [tf.Variable(tf.fill((basis.shape[1], k), -2.0)) for basis in basis_list]
        optimizer = tf.keras.optimizers.Adam(self.config.learning_rate)
        n_steps = self.config.steps if steps is None else int(steps)
        if n_steps < 1:
            raise ValueError("inference steps must be positive")
        seed_everything(self.config.seed + 6000)
        x = tf.constant(residual)
        loading = tf.constant(self.loading_, dtype=tf.float32)
        variables = [*q_loc, *q_scale]
        best = np.inf
        best_state = None
        for step in range(n_steps):
            with tf.GradientTape() as tape:
                fields = []
                prior = tf.constant(0.0, dtype=tf.float32)
                for index, basis in enumerate(basis_list):
                    loc = q_loc[index]
                    scale = tf.nn.softplus(q_scale[index]) + 1e-4
                    eps = tf.random.stateless_normal(tf.shape(loc), seed=[self.config.seed + 6000, step + index + 1])
                    u = loc + scale * eps
                    value = tf.constant(basis) @ u
                    value -= tf.reduce_mean(value, axis=0, keepdims=True)
                    fields.append(value)
                    inv = tf.constant(k_inv_list[index])
                    logdet = tf.constant(logdet_prior_list[index], dtype=tf.float32)
                    prior += _diagonal_gp_kl(tf, loc, scale, inv, logdet)
                field = tf.concat(fields, axis=0)
                predicted = field @ tf.transpose(loading)
                likelihood = gene_batch_sum_objective(
                    (x - predicted) ** 2,
                    counts.shape[1],
                    counts.shape[1],
                    tf_module=tf,
                )
                loss = likelihood + 1e-4 * prior / float(max(len(counts), 1))
            gradients = tape.gradient(loss, variables)
            optimizer.apply_gradients(zip(gradients, variables))
            value = float(loss.numpy())
            if value < best:
                best = value
                best_state = [variable.numpy().copy() for variable in variables]
        if best_state is None:
            raise RuntimeError("frozen signed inference produced no state")
        for variable, value in zip(variables, best_state):
            variable.assign(value)
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
            rng = np.random.default_rng(self.config.seed + 6100 + index)
            draws = np.stack([
                basis_np @ (loc + scale * rng.normal(size=loc.shape))
                for _ in range(draws_requested)
            ], axis=0)
            draws -= draws.mean(axis=1, keepdims=True)
            result.append([
                FieldFit(
                    self.model_id,
                    f"SIGNED_{factor + 1:02d}",
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
                        "nuisance": "NB_PEARSON_FROZEN_PROFILE",
                        "loss_best": best,
                        "patient_id": section.patient_id,
                    },
                )
                for factor in range(k)
            ])
        return result
