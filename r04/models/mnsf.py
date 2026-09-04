"""A compact variational multi-section NSF-style count model (PyTorch).

This is a PyTorch port of the historical TensorFlow/TFP implementation.
Science is preserved: continuous overlapping irregular latent fields,
NB count objective with offset, nonspatial nuisance, gene-minibatch
Horvitz-Thompson scaling, inducing-point GP, variational KL, K=0 fair
baseline, fold/gene split, seed semantics and checkpoint/best-state.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import time
import json
from typing import Sequence

import numpy as np
from scipy.special import logsumexp
import torch
import torch.nn.functional as F

from ..runtime import atomic_json, config_fingerprint, seed_everything
from ..io_contract import sections_content_hash, sections_objective_input_hash
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
    execution_mode: str = "auto"


class MissingTorchBackend(RuntimeError):
    """Raised when PyTorch is not available."""


def _get_device() -> torch.device:
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def _resolve_torch_execution_mode(requested: str) -> str:
    """Resolve the supported Torch execution mode without changing the math.

    The refactor currently has one implemented execution path: native eager
    Torch.  ``auto`` is retained as a protocol-friendly alias and may select
    CUDA through :func:`_get_device`, but it does not claim to use
    ``torch.compile``.  A future compiled path must be added deliberately with
    its own equivalence tests before it is accepted here.
    """
    if requested not in {"auto", "eager"}:
        raise ValueError(
            "execution_mode must be 'auto' or 'eager'; 'compiled' is not implemented"
        )
    return "eager"


def _diagonal_gp_kl_torch(
    loc: torch.Tensor,
    scale: torch.Tensor,
    prior_inverse: torch.Tensor,
    logdet_prior: float,
) -> torch.Tensor:
    """KL[q(u)||N(0,K)] for diagonal q covariance, per factor jointly (torch)."""
    # loc, scale: (m, k), prior_inverse: (m, m)
    # quadratic = sum( loc^T K^{-1} loc ) per factor summed
    # loc: (m,k), prior_inverse @ loc : (m,k)
    quadratic = torch.sum((prior_inverse @ loc) * loc)
    # trace term: sum_i diag(K^{-1})_{ii} * scale_{i,f}^2
    diag = torch.diag(prior_inverse)  # (m,)
    trace = torch.sum((scale * scale) * diag[:, None])
    logdet_q = torch.sum(torch.log(scale * scale))
    factors = float(loc.shape[1]) if loc.shape[1] > 0 else 0.0
    dimension = float(loc.numel())
    return 0.5 * (quadratic + trace - dimension + factors * float(logdet_prior) - logdet_q)


_TORCH_INTERRUPT_COUNTER = 0
_TORCH_INTERRUPT_TRIGGER: int | None = None


def _maybe_raise_torch_interrupt() -> None:
    """Hook for interrupt tests that monkeypatch ``torch.optim.Adam.step``."""
    global _TORCH_INTERRUPT_COUNTER
    if _TORCH_INTERRUPT_TRIGGER is None:
        return
    _TORCH_INTERRUPT_COUNTER += 1
    if _TORCH_INTERRUPT_COUNTER == _TORCH_INTERRUPT_TRIGGER:
        raise RuntimeError("simulated interruption")


def _maybe_raise_tf_interrupt() -> None:
    """Backward-compat alias for ``_maybe_raise_torch_interrupt``."""
    _maybe_raise_torch_interrupt()


def _diagonal_gp_kl(tf_module, loc, scale, prior_inverse, logdet):  # type: ignore[no-untyped-def]
    """Backward-compat shim for tests that import the TF-named helper.

    The torch port is :func:`_diagonal_gp_kl_torch`; this wrapper preserves the
    historical ``_diagonal_gp_kl(tf, ...)`` call signature used in
    ``tests/test_r04_models.py`` and returns an object with ``.numpy()``.
    """

    def _to_numpy(x):  # type: ignore[no-untyped-def]
        if hasattr(x, "numpy"):
            try:
                return np.asarray(x.numpy())
            except Exception:
                pass
        return np.asarray(x)

    loc_np = _to_numpy(loc).astype(float)
    scale_np = _to_numpy(scale).astype(float)
    inv_np = _to_numpy(prior_inverse).astype(float)
    logdet_f = float(_to_numpy(logdet).squeeze() if hasattr(_to_numpy(logdet), "squeeze") else _to_numpy(logdet)) if np.asarray(_to_numpy(logdet)).size == 1 else float(np.asarray(_to_numpy(logdet)).flat[0])
    # numpy KL: same formula as torch version
    quadratic = float(np.sum((inv_np @ loc_np) * loc_np)) if loc_np.size else 0.0
    diag = np.diag(inv_np) if inv_np.size else np.array([], dtype=float)
    trace = float(np.sum((scale_np * scale_np) * diag[:, None])) if scale_np.size and diag.size else 0.0
    logdet_q = float(np.sum(np.log(scale_np * scale_np))) if scale_np.size else 0.0
    factors = float(loc_np.shape[1]) if loc_np.ndim == 2 and loc_np.shape[1] > 0 else 0.0
    dimension = float(loc_np.size)
    val = 0.5 * (quadratic + trace - dimension + factors * logdet_f - logdet_q)

    class _Wrap:
        def __init__(self, v: float) -> None:
            self._v = float(v)

        def numpy(self) -> float:  # type: ignore[no-untyped-def]
            return self._v

        def __float__(self) -> float:
            return float(self._v)

    return _Wrap(val)


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
    counts = np.asarray(counts, dtype=float)
    library = np.asarray(library, dtype=float)
    if counts.ndim != 2 or library.shape != (len(counts),) or factors < 0:
        raise ValueError("counts, library and factors have incompatible shapes")
    if np.any(counts < 0) or np.any(library <= 0):
        raise ValueError("counts and library must be non-negative and positive")
    rate = counts / library[:, None]
    gene_mean = np.maximum(rate.mean(axis=0), 1e-4)
    if factors == 0:
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


def _nb_log_prob(
    y: torch.Tensor,
    mu: torch.Tensor,
    theta: torch.Tensor,
) -> torch.Tensor:
    """NB log prob with mean mu and dispersion theta (total_count)."""
    # y, mu, theta broadcastable; theta shape (1, g) or (g,)
    # mu >0, theta>0, y >=0 integer
    # log_prob = lgamma(y+theta)-lgamma(theta)-lgamma(y+1) + theta*log(theta/(theta+mu)) + y*log(mu/(theta+mu))
    # Use lgamma for stability
    y_f = y.float()
    theta_f = theta.float()
    mu_f = mu.float().clamp(min=1e-8)
    theta_f = theta_f.clamp(min=1e-8)
    term = torch.lgamma(y_f + theta_f) - torch.lgamma(theta_f) - torch.lgamma(y_f + 1)
    term = term + theta_f * (torch.log(theta_f) - torch.log(theta_f + mu_f))
    term = term + y_f * (torch.log(mu_f) - torch.log(theta_f + mu_f))
    return term


def _stateless_normal(shape: torch.Size, seed: int, device: torch.device, dtype: torch.dtype) -> torch.Tensor:
    gen = torch.Generator(device=device)
    # torch.Generator seed is 64-bit
    gen.manual_seed(int(seed) & 0xFFFFFFFFFFFFFFFF)
    return torch.randn(shape, generator=gen, device=device, dtype=dtype)


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
        self.execution_mode_requested_legacy_: str | None = None
        self.execution_mode_migration_: str | None = None

    def fit(self, sections: Sequence[SectionData]) -> "MNSFEstimator":
        if not sections:
            raise ValueError("at least one section is required")
        device = _get_device()
        execution_mode = _resolve_torch_execution_mode(self.config.execution_mode)
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
            "torch": getattr(torch, "__version__", "unknown"),
        })
        environment_hash = self.config.environment_hash or runtime_environment_hash
        genes = sections[0].gene_id
        if any(section.gene_id != genes for section in sections):
            raise ValueError("all sections must share the frozen gene universe")
        input_hash = sections_content_hash(sections)
        objective_input_hash = sections_objective_input_hash(sections)
        counts_np = np.concatenate([section.counts.toarray() if hasattr(section.counts, "toarray") else np.asarray(section.counts) for section in sections], axis=0).astype(np.float32)
        if np.any(counts_np != np.floor(counts_np)) or np.any(counts_np < 0):
            raise ValueError("mNSF requires non-negative integer counts")
        if any(section.library_size is not None for section in sections) and not all(
            section.library_size is not None for section in sections
        ):
            raise ValueError("library_size metadata must be present for every section or none")
        library_np = (
            np.concatenate([np.asarray(section.library_size, dtype=np.float32) for section in sections])
            if all(section.library_size is not None for section in sections)
            else counts_np.sum(axis=1).astype(np.float32)
        )
        library_np = np.maximum(library_np, 1.0)
        log_library_np = np.log(library_np)
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
        # Torch tensors
        basis_torch = [torch.tensor(b, dtype=torch.float32, device=device) for b in basis_list]
        prior_inv_torch = [torch.tensor(v, dtype=torch.float32, device=device) for v in k_inv_list]
        y_torch = torch.tensor(counts_np, dtype=torch.float32, device=device)
        log_l_torch = torch.tensor(log_library_np, dtype=torch.float32, device=device)
        n, g = counts_np.shape
        k = self.config.factors
        q = self.config.nonspatial_rank
        if k < 0 or q < 0:
            raise ValueError("factors must be non-negative and nonspatial_rank non-negative")
        if self.config.gp_parameterization != "inducing_values_v2":
            raise ValueError("unsupported GP parameterization")
        if self.config.initialization == "rate_k_invariant":
            loading0, amplitude0, baseline0, dispersion0 = _rate_initial_components(
                counts_np, library_np, k, seed=self.config.seed
            )
            raw_w_t = torch.tensor(np.log(np.maximum(loading0, 1e-8)).astype(np.float32) if k > 0 else np.zeros((g, 0), dtype=np.float32), device=device)
            raw_a_t = torch.tensor(_softplus_inverse(amplitude0).astype(np.float32) if k > 0 else np.zeros((0,), dtype=np.float32), device=device)
            raw_b_t = torch.tensor(_softplus_inverse(baseline0).astype(np.float32), device=device)
            raw_theta_t = torch.tensor(_softplus_inverse(dispersion0).astype(np.float32), device=device)
        elif self.config.initialization == "legacy":
            gen = torch.Generator(device=device)
            gen.manual_seed(self.config.seed)
            raw_w_t = torch.randn((g, k), generator=gen, device=device, dtype=torch.float32) * 0.05 if k > 0 else torch.zeros((g, 0), device=device, dtype=torch.float32)
            raw_a_t = torch.zeros((k,), device=device, dtype=torch.float32) if k > 0 else torch.zeros((0,), device=device, dtype=torch.float32)
            raw_b_t = torch.full((g,), -3.0, device=device, dtype=torch.float32)
            raw_theta_t = torch.full((g,), 3.0, device=device, dtype=torch.float32)
        else:
            raise ValueError("unsupported mNSF initialization")
        # Ensure requires_grad
        params_to_optimize: list[torch.nn.Parameter] = []
        def _make_param(t: torch.Tensor) -> torch.nn.Parameter:
            p = torch.nn.Parameter(t)
            return p
        raw_w = _make_param(raw_w_t) if raw_w_t.numel() > 0 else torch.nn.Parameter(torch.zeros((g, 0), device=device, dtype=torch.float32))
        # For k=0, raw_w has shape (g,0) but still a Parameter; optimizer will handle empty?
        # We'll keep it but filter empty numel later.
        raw_a = _make_param(raw_a_t) if raw_a_t.numel() > 0 else torch.nn.Parameter(torch.zeros((0,), device=device, dtype=torch.float32))
        raw_b = _make_param(raw_b_t)
        raw_theta = _make_param(raw_theta_t)
        if q:
            gen2 = torch.Generator(device=device)
            gen2.manual_seed(self.config.seed + 1)
            raw_v_t = torch.randn((g, q), generator=gen2, device=device, dtype=torch.float32) * 0.05
            raw_h_t = torch.full((n, q), -4.0, device=device, dtype=torch.float32)
            raw_v = _make_param(raw_v_t)
            raw_h = _make_param(raw_h_t)
        else:
            raw_v = None
            raw_h = None
        q_loc: list[torch.nn.Parameter] = []
        q_scale_raw: list[torch.nn.Parameter] = []
        for index, section in enumerate(sections):
            m = basis_list[index].shape[1]
            q_loc.append(torch.nn.Parameter(torch.zeros((m, k), device=device, dtype=torch.float32)))
            q_scale_raw.append(torch.nn.Parameter(torch.full((m, k), -2.0, device=device, dtype=torch.float32)))
        if self.config.learning_rate <= 0:
            raise ValueError("learning_rate must be positive")
        if self.config.learning_rate_final is not None and self.config.learning_rate_final <= 0:
            raise ValueError("learning_rate_final must be positive")
        if self.config.learning_rate_final is not None and (
            self.config.learning_rate_decay_steps is None
            or self.config.learning_rate_decay_steps < 1
        ):
            raise ValueError("learning_rate_decay_steps is required when learning_rate_final is set")
        if self.config.diagnostic_interval < 0:
            raise ValueError("diagnostic_interval must be non-negative")
        if self.config.evaluation_interval < 0:
            raise ValueError("evaluation_interval must be non-negative")
        if self.config.evaluation_mc_draws < 1:
            raise ValueError("evaluation_mc_draws must be positive")
        if self.config.optimization_schedule not in {"joint", "staged_shared_first"}:
            raise ValueError("optimization_schedule must be 'joint' or 'staged_shared_first'")
        if self.config.shared_steps < 0:
            raise ValueError("shared_steps must be non-negative")

        def scheduled_lr(step: int) -> float:
            if self.config.learning_rate_final is None:
                return float(self.config.learning_rate)
            decay_steps = int(self.config.learning_rate_decay_steps)
            fraction = min(max(float(step) / decay_steps, 0.0), 1.0)
            return float(self.config.learning_rate + fraction * (self.config.learning_rate_final - self.config.learning_rate))

        # Build optimizer with non-empty params
        all_params: list[torch.nn.Parameter] = [raw_w, raw_a, raw_b, raw_theta, *q_loc, *q_scale_raw]
        if raw_v is not None and raw_h is not None:
            all_params.extend([raw_v, raw_h])
        optim_params = [p for p in all_params if p.numel() > 0]
        optimizer = torch.optim.Adam(optim_params, lr=float(self.config.learning_rate))

        # Keep the compact p_* representation used by resume, while making
        # the positional contract explicit for independent checkpoint audits.
        q_loc_keys = [f"p_{4 + index}" for index in range(len(q_loc))]
        q_scale_start = 4 + len(q_loc)
        q_scale_keys = [f"p_{q_scale_start + index}" for index in range(len(q_scale_raw))]
        nuisance_start = q_scale_start + len(q_scale_raw)
        parameter_layout: dict[str, object] = {
            "schema": "r04.mnsf_parameter_layout.v1",
            "raw_w": "p_0",
            "raw_a": "p_1",
            "raw_b": "p_2",
            "raw_theta": "p_3",
            "q_loc": q_loc_keys,
            "q_scale": q_scale_keys,
            "raw_v": f"p_{nuisance_start}" if raw_v is not None else None,
            "raw_h": f"p_{nuisance_start + 1}" if raw_h is not None else None,
        }

        def checkpoint_payload() -> dict[str, object]:
            return {
                "checkpoint_schema": "r04.mnsf_checkpoint.v3_torch",
                "backend": "torch",
                "input_hash": input_hash,
                "objective_input_hash": objective_input_hash,
                "config_hash": config_hash,
                "environment_hash": environment_hash,
                "seed": self.config.seed,
                "factors": k,
                "global_step": global_step,
                "best_loss": best_loss_val,
                "best_step": best_step,
                "parameter_layout": parameter_layout,
                "best_state": best_state,
                "params": {
                    f"p_{idx}": p.detach().cpu().clone()
                    for idx, p in enumerate(all_params)
                },
                "optimizer_state": optimizer.state_dict(),
            }

        def _validate_param_map(param_map: object, *, label: str) -> None:
            if not isinstance(param_map, dict):
                raise RuntimeError(f"{label} parameter map is missing or invalid")
            expected = {f"p_{idx}": parameter for idx, parameter in enumerate(all_params)}
            missing = sorted(set(expected) - set(param_map))
            if missing:
                raise RuntimeError(
                    f"{label} is missing parameter tensors: {', '.join(missing)}"
                )
            for key, parameter in expected.items():
                value = param_map[key]
                if not isinstance(value, torch.Tensor):
                    raise RuntimeError(f"{label} tensor {key} is not a Torch tensor")
                if tuple(value.shape) != tuple(parameter.shape):
                    raise RuntimeError(
                        f"{label} tensor {key} shape mismatch: "
                        f"{tuple(value.shape)} != {tuple(parameter.shape)}"
                    )
                if not bool(torch.isfinite(value).all().item()):
                    raise RuntimeError(f"{label} tensor {key} contains non-finite values")

        def _validate_resume_payload(
            payload: object,
            *,
            source: Path,
            require_config_hash: bool,
        ) -> None:
            if not isinstance(payload, dict):
                raise RuntimeError(f"resume checkpoint is not a mapping: {source}")
            if (
                payload.get("checkpoint_schema") != "r04.mnsf_checkpoint.v3_torch"
                or payload.get("backend") != "torch"
            ):
                raise RuntimeError("resume checkpoint schema/backend mismatch")
            if payload.get("input_hash") != input_hash:
                raise RuntimeError("resume checkpoint input hash mismatch")
            if payload.get("objective_input_hash") != objective_input_hash:
                raise RuntimeError("resume checkpoint objective input hash mismatch")
            if require_config_hash and payload.get("config_hash") != config_hash:
                raise RuntimeError("resume checkpoint config hash mismatch")
            if payload.get("environment_hash") != environment_hash:
                raise RuntimeError("resume checkpoint environment hash mismatch")
            _validate_param_map(payload.get("params"), label="resume checkpoint params")
            best_state = payload.get("best_state")
            if best_state is None:
                raise RuntimeError("resume checkpoint best_state is missing")
            _validate_param_map(best_state, label="resume checkpoint best_state")
            if not isinstance(payload.get("optimizer_state"), dict):
                raise RuntimeError("resume checkpoint optimizer_state is missing or invalid")

        # Checkpoint handling - torch.save dict
        global_step = 0
        best_loss_val = float("inf")
        best_step = -1
        best_state: dict[str, torch.Tensor] | None = None

        # Helper to snapshot best state
        def _snapshot_best():
            snap: dict[str, torch.Tensor] = {}
            for idx, p in enumerate(all_params):
                snap[f"p_{idx}"] = p.detach().cpu().clone()
            return snap

        def _restore_best(snap: dict[str, torch.Tensor]):
            for idx, p in enumerate(all_params):
                if f"p_{idx}" in snap:
                    p.data.copy_(snap[f"p_{idx}"].to(p.device))

        manager_dir: Path | None = None
        resume_dir: Path | None = None
        if self.config.resume_checkpoint_dir is not None and self.config.checkpoint_dir is None:
            raise ValueError("resume_checkpoint_dir requires checkpoint_dir")
        if self.config.checkpoint_dir:
            manager_dir = Path(self.config.checkpoint_dir)
            manager_dir.mkdir(parents=True, exist_ok=True)
            metadata_path = manager_dir / "checkpoint_metadata.json"
            metadata = {
                "checkpoint_schema": "r04.mnsf_checkpoint.v3_torch",
                "backend": "torch",
                "input_hash": input_hash,
                "objective_input_hash": objective_input_hash,
                "config_hash": config_hash,
                "environment_hash": environment_hash,
                "seed": self.config.seed,
                "factors": k,
            }
            if metadata_path.exists():
                previous = json.loads(metadata_path.read_text(encoding="utf-8"))
                if previous != metadata:
                    raise RuntimeError("checkpoint hash mismatch; start a new R04 run")
            else:
                atomic_json(metadata_path, metadata)
            # Try to load existing checkpoint for continuation
            ckpt_path = manager_dir / "checkpoint.pt"
            # Resume from explicit dir
            if self.config.resume_checkpoint_dir is not None:
                resume_dir = Path(self.config.resume_checkpoint_dir)
                resume_meta = resume_dir / "checkpoint_metadata.json"
                if not resume_meta.exists():
                    raise RuntimeError("resume checkpoint metadata is missing")
                src_meta = json.loads(resume_meta.read_text(encoding="utf-8"))
                if src_meta.get("checkpoint_schema") != "r04.mnsf_checkpoint.v3_torch":
                    raise RuntimeError("resume checkpoint schema mismatch")
                if src_meta.get("backend") != "torch":
                    raise RuntimeError("resume checkpoint backend mismatch")
                if src_meta.get("input_hash") != input_hash:
                    raise RuntimeError("resume checkpoint input hash mismatch")
                if src_meta.get("objective_input_hash") != objective_input_hash:
                    raise RuntimeError("resume checkpoint objective input hash mismatch")
                if src_meta.get("environment_hash") != environment_hash:
                    raise RuntimeError("resume checkpoint environment hash mismatch")
                src_ckpt = resume_dir / "checkpoint.pt"
                if not src_ckpt.exists():
                    raise RuntimeError("resume checkpoint is missing")
                payload = torch.load(str(src_ckpt), map_location=device, weights_only=True)
                _validate_resume_payload(
                    payload, source=src_ckpt, require_config_hash=False
                )
                # Restore params and optimizer.  The parameter order is part of
                # the checkpoint contract; missing tensors fail closed.
                for idx, p in enumerate(all_params):
                    p.data.copy_(payload["params"][f"p_{idx}"].to(device))
                optimizer.load_state_dict(payload["optimizer_state"])
                # Move optimizer state to device
                for state in optimizer.state.values():
                    for kk, vv in state.items():
                        if isinstance(vv, torch.Tensor):
                            state[kk] = vv.to(device)
                global_step = int(payload.get("global_step", 0))
                best_loss_val = float(payload.get("best_loss", float("inf")))
                best_step = int(payload.get("best_step", -1))
                best_state = payload.get("best_state")
                if best_state is not None:
                    # best_state tensors are on cpu, keep as is
                    pass
            elif ckpt_path.exists():
                payload = torch.load(str(ckpt_path), map_location=device, weights_only=True)
                # Fail-closed on schema/backend
                if payload.get("checkpoint_schema") != "r04.mnsf_checkpoint.v3_torch" or payload.get("backend") != "torch":
                    raise RuntimeError("checkpoint schema/backend mismatch; start a new R04 run")
                if payload.get("input_hash") != input_hash or payload.get("config_hash") != config_hash or payload.get("environment_hash") != environment_hash:
                    raise RuntimeError("checkpoint hash mismatch; start a new R04 run")
                _validate_resume_payload(
                    payload, source=ckpt_path, require_config_hash=True
                )
                for idx, p in enumerate(all_params):
                    p.data.copy_(payload["params"][f"p_{idx}"].to(device))
                optimizer.load_state_dict(payload["optimizer_state"])
                for state in optimizer.state.values():
                    for kk, vv in state.items():
                        if isinstance(vv, torch.Tensor):
                            state[kk] = vv.to(device)
                global_step = int(payload.get("global_step", 0))
                best_loss_val = float(payload.get("best_loss", float("inf")))
                best_step = int(payload.get("best_step", -1))
                best_state = payload.get("best_state")

        # Diagnostics traces
        losses: list[float] = []
        gradient_norm_trace: list[dict[str, object]] = []
        parameter_norm_trace: list[dict[str, object]] = []
        learning_rate_trace: list[dict[str, float]] = []
        optimization_stage_trace: list[dict[str, object]] = []
        evaluation_trace: list[dict[str, object]] = []
        start_step = int(global_step)

        def optimization_stage(step: int) -> str:
            if (
                self.config.optimization_schedule == "staged_shared_first"
                and k > 0
                and step < self.config.shared_steps
            ):
                return "shared_first"
            return "joint"

        # Deterministic objective helper (dense, fixed seed)
        def deterministic_objective() -> tuple[float, float]:
            eval_losses: list[float] = []
            # no grad
            with torch.no_grad():
                for draw in range(self.config.evaluation_mc_draws):
                    kl_eval = torch.tensor(0.0, device=device)
                    samples_eval: list[torch.Tensor] = []
                    for idx, loc in enumerate(q_loc):
                        scale = F.softplus(q_scale_raw[idx]) + 1e-4
                        eps = _stateless_normal(loc.shape, self.config.seed + 7000 + draw * 1000 + idx + 1, device, loc.dtype)
                        samples_eval.append(loc + scale * eps)
                        kl_eval = kl_eval + _diagonal_gp_kl_torch(loc, scale, prior_inv_torch[idx], logdet_prior_list[idx])
                    field_parts = [basis_torch[idx] @ samples_eval[idx] for idx in range(len(sections))]
                    field_centered = [p - p.mean(dim=0, keepdim=True) for p in field_parts]
                    field_eval = torch.cat(field_centered, dim=0) if field_centered and field_centered[0].shape[1] > 0 else torch.zeros((n, 0), device=device)
                    # Build components
                    components: list[torch.Tensor] = []
                    if k:
                        loading_eval = F.softmax(raw_w, dim=0)
                        amp_eval = F.softplus(raw_a) + 1e-5
                        # log loading + log amp + field
                        # field_eval: (n,k), loading: (g,k)
                        log_loading = torch.log(loading_eval.clamp(min=1e-8))  # (g,k)
                        log_amp = torch.log(amp_eval.clamp(min=1e-8))  # (k,)
                        # broadcast: (n,g,k)
                        log_comp = log_loading[None, :, :] + log_amp[None, None, :] + field_eval[:, None, :]
                        components.append(torch.logsumexp(log_comp, dim=2))
                    if q and raw_v is not None and raw_h is not None:
                        v_eval = F.softmax(raw_v, dim=0)
                        h_eval = F.softplus(raw_h) + 1e-5
                        components.append(torch.log((h_eval @ v_eval.T).clamp(min=1e-8)))
                    baseline_eval = F.softplus(raw_b) + 1e-5
                    components.append(torch.log(baseline_eval.clamp(min=1e-8))[None, :].expand(n, g))
                    log_mu_eval = log_l_torch[:, None] + torch.logsumexp(torch.stack(components, dim=2), dim=2)
                    theta_eval = (F.softplus(raw_theta) + 1e-3)[None, :]
                    mu_eval = torch.exp(log_mu_eval)
                    # NB log prob dense
                    lp = _nb_log_prob(y_torch, mu_eval, theta_eval)
                    likelihood_eval = gene_batch_sum_objective(-lp, g, g)
                    loss_eval = likelihood_eval + 1e-4 * kl_eval / float(n)
                    eval_losses.append(float(loss_eval.item()))
            return float(np.mean(eval_losses)), float(np.std(eval_losses, ddof=1) if len(eval_losses) > 1 else 0.0)

        # Identify spatial params for staged schedule
        spatial_ids = {id(raw_w), *[id(p) for p in q_loc], *[id(p) for p in q_scale_raw]}
        # Remove empty params from spatial set (raw_w/raw_a when k=0)
        spatial_ids = {i for i in spatial_ids if i is not None}

        optimizer_started = time.perf_counter()
        for step in range(start_step, self.config.steps):
            # scheduled lr
            lr = scheduled_lr(step)
            for pg in optimizer.param_groups:
                pg["lr"] = lr
            batch_size = g if self.config.gene_batch_size is None else int(self.config.gene_batch_size)
            if batch_size < 1:
                raise ValueError("gene_batch_size must be positive")
            batch_size = min(batch_size, g)
            if batch_size == g:
                gene_index = np.arange(g, dtype=np.int64)
            else:
                gene_index = np.random.default_rng(self.config.seed + 9001 + step).choice(g, size=batch_size, replace=False).astype(np.int64)
            gene_index_t = torch.tensor(gene_index, dtype=torch.long, device=device)
            stage = optimization_stage(step)
            # Forward + loss
            optimizer.zero_grad()
            # Sample GP
            samples: list[torch.Tensor] = []
            kl = torch.tensor(0.0, device=device)
            for idx, loc in enumerate(q_loc):
                scale = F.softplus(q_scale_raw[idx]) + 1e-4
                eps = _stateless_normal(loc.shape, self.config.seed + (step * 1000 + idx + 1) * 917, device, loc.dtype)
                # Use same seed scheme as before but deterministic; use step+index+1 pattern
                # Keep compatible with simple deterministic seeding
                samples.append(loc + scale * eps)
                kl = kl + _diagonal_gp_kl_torch(loc, scale, prior_inv_torch[idx], logdet_prior_list[idx])
            field_parts = [basis_torch[idx] @ samples[idx] for idx in range(len(sections))] if k > 0 else [torch.zeros((basis_torch[idx].shape[0], 0), device=device) for idx in range(len(sections))]
            field_centered = [p - p.mean(dim=0, keepdim=True) if p.shape[1] > 0 else p for p in field_parts]
            field = torch.cat(field_centered, dim=0) if field_centered else torch.zeros((n, 0), device=device)
            # Loading etc
            if k:
                loading = F.softmax(raw_w, dim=0)
                loading_batch = loading[gene_index_t, :]
                amp = F.softplus(raw_a) + 1e-5
                log_comp = torch.log(loading_batch.clamp(min=1e-8))[None, :, :] + torch.log(amp.clamp(min=1e-8))[None, None, :] + field[:, None, :]
            # Components
            components: list[torch.Tensor] = []
            if k:
                components.append(torch.logsumexp(log_comp, dim=2))
            if q and raw_v is not None and raw_h is not None:
                v = F.softmax(raw_v, dim=0)
                v_batch = v[gene_index_t, :]
                h = F.softplus(raw_h) + 1e-5
                components.append(torch.log((h @ v_batch.T).clamp(min=1e-8)))
            baseline_batch = (F.softplus(raw_b) + 1e-5)[gene_index_t]
            components.append(torch.log(baseline_batch.clamp(min=1e-8))[None, :].expand(n, batch_size))
            log_mu = log_l_torch[:, None] + torch.logsumexp(torch.stack(components, dim=2), dim=2)
            theta_batch = (F.softplus(raw_theta) + 1e-3)[gene_index_t][None, :]
            y_batch = y_torch[:, gene_index_t]
            mu = torch.exp(log_mu)
            lp = _nb_log_prob(y_batch, mu, theta_batch)
            likelihood = gene_batch_sum_objective(-lp, g, batch_size)
            loss = likelihood + 1e-4 * kl / float(n)
            value = float(loss.item())
            # Best state capture before optimizer step (pre-update loss)
            if value < best_loss_val:
                best_loss_val = value
                best_step = step + 1
                best_state = _snapshot_best()
            # Diagnostic due
            diagnostic_due = (
                self.config.diagnostic_interval > 0
                and (
                    (step + 1) % self.config.diagnostic_interval == 0
                    or step == start_step
                    or step + 1 == self.config.steps
                )
            )
            if diagnostic_due:
                loss.backward()
                # Compute gradient norms per group
                grad_by_id = {id(p): p.grad for p in all_params if p.grad is not None}
                groups = {
                    "loading": [raw_w],
                    "factor_amplitude": [raw_a],
                    "baseline": [raw_b],
                    "dispersion": [raw_theta],
                    "gp_loc": q_loc,
                    "gp_scale": q_scale_raw,
                    "nonspatial_loading": [raw_v] if raw_v is not None else [],
                    "nonspatial_field": [raw_h] if raw_h is not None else [],
                }
                def _norm(vals, use_grad: bool) -> float:
                    total = 0.0
                    for v in vals:
                        t = grad_by_id.get(id(v)) if use_grad else v
                        if t is not None:
                            total += float((t.detach() ** 2).sum().item())
                    return float(np.sqrt(total))
                gradient_norm_trace.append({
                    "step": step + 1,
                    "optimization_stage": stage,
                    **{name: _norm(vals, True) for name, vals in groups.items()},
                })
                parameter_norm_trace.append({
                    "step": step + 1,
                    "optimization_stage": stage,
                    **{name: _norm(vals, False) for name, vals in groups.items()},
                })
                learning_rate_trace.append({"step": float(step + 1), "learning_rate": float(lr)})
                optimization_stage_trace.append({"step": step + 1, "optimization_stage": stage})
                # Handle staged masking before step
                if stage == "shared_first":
                    for p in all_params:
                        if id(p) in spatial_ids and p.grad is not None:
                            p.grad = None
                _maybe_raise_torch_interrupt()
                optimizer.step()
                optimizer.zero_grad()
            else:
                loss.backward()
                if stage == "shared_first":
                    for p in all_params:
                        if id(p) in spatial_ids and p.grad is not None:
                            p.grad = None
                _maybe_raise_torch_interrupt()
                optimizer.step()
                optimizer.zero_grad()
            # Evaluation trace
            evaluation_due = (
                self.config.evaluation_interval > 0
                and (
                    (step + 1) % self.config.evaluation_interval == 0
                    or step == start_step
                    or step + 1 == self.config.steps
                )
            )
            if evaluation_due:
                em, es = deterministic_objective()
                evaluation_trace.append({
                    "step": step + 1,
                    "objective_mean": em,
                    "objective_sd": es,
                    "mc_draws": self.config.evaluation_mc_draws,
                    "objective_kind": "dense_full_panel_fixed_stateless_mc",
                })
            global_step = step + 1
            losses.append(value)
            if manager_dir is not None and (step + 1) % max(1, self.config.checkpoint_steps) == 0:
                torch.save(checkpoint_payload(), str(manager_dir / "checkpoint.pt"))
        optimizer_wall_seconds = time.perf_counter() - optimizer_started
        if not np.isfinite(best_loss_val) or best_state is None:
            raise RuntimeError("mNSF optimizer produced no persisted best state")
        # Restore best state
        _restore_best(best_state)
        if manager_dir is not None:
            torch.save(checkpoint_payload(), str(manager_dir / "checkpoint.pt"))
        # Extract final parameters as numpy
        with torch.no_grad():
            if k:
                loading_np = F.softmax(raw_w, dim=0).cpu().numpy()
            else:
                loading_np = np.zeros((g, 0), dtype=float)
            amp_np = (F.softplus(raw_a) + 1e-5).cpu().numpy() if raw_a.numel() > 0 else np.zeros((0,), dtype=float)
            baseline_np = (F.softplus(raw_b) + 1e-5).cpu().numpy()
            disp_np = (F.softplus(raw_theta) + 1e-3).cpu().numpy()
            nonspatial_loading_np = F.softmax(raw_v, dim=0).cpu().numpy() if raw_v is not None else None
        self.loading_ = loading_np
        self.factor_amplitude_ = amp_np
        self.gene_baseline_ = baseline_np
        self.dispersion_ = disp_np
        self.nonspatial_loading_ = nonspatial_loading_np
        # Build fields with posterior draws (numpy)
        self.fields_ = []
        for index, (r0, r1) in enumerate(section_ranges):
            loc_np = q_loc[index].detach().cpu().numpy()
            scale_np = F.softplus(q_scale_raw[index]).detach().cpu().numpy() + 1e-4
            basis_np = basis_list[index]
            mean = basis_np @ loc_np
            mean -= mean.mean(axis=0, keepdims=True)
            rng = np.random.default_rng(self.config.seed + index)
            draws = np.stack([basis_np @ (loc_np + scale_np * rng.normal(size=loc_np.shape)) for _ in range(self.config.posterior_draws)], axis=0)
            draws -= draws.mean(axis=1, keepdims=True)
            self.fields_.append([
                FieldFit(self.model_id, f"MNSF_{factor+1:02d}", sections[index].section_id,
                         loading_np[:, factor], mean[:, factor], draws[:, :, factor].std(axis=0),
                         self.config.lengthscale, float(np.var(mean[:, factor])),
                         input_hash=input_hash, diagnostics={"uncertainty_kind": "variational_inducing_weights", "steps": self.config.steps, "patient_id": sections[index].patient_id})
                for factor in range(k)
            ])
        field_matrix = (
            np.column_stack([
                np.concatenate([group[factor].field_mean for group in self.fields_])
                for factor in range(k)
            ])
            if k else np.empty((len(counts_np), 0), dtype=float)
        )
        self.diagnostics_ = {
            "loss_final": losses[-1] if losses else None,
            "loss_best": best_loss_val if np.isfinite(best_loss_val) else None,
            "best_step": int(best_step),
            "loss_trace": [float(v) for v in losses],
            "steps": self.config.steps,
            "start_step": start_step,
            "optimizer_steps_this_call": self.config.steps - start_step,
            "optimizer_wall_seconds": optimizer_wall_seconds,
            "execution_mode": execution_mode,
            "execution_mode_effective": execution_mode,
            "execution_mode_requested": self.config.execution_mode,
            "learning_rate_final": self.config.learning_rate_final,
            "learning_rate_decay_steps": self.config.learning_rate_decay_steps,
            "learning_rate_trace": learning_rate_trace,
            "optimization_schedule": self.config.optimization_schedule,
            "shared_steps": self.config.shared_steps,
            "optimization_stage_trace": optimization_stage_trace,
            "evaluation_interval": self.config.evaluation_interval,
            "evaluation_mc_draws": self.config.evaluation_mc_draws,
            "evaluation_trace": evaluation_trace,
            "deterministic_objective_platform": deterministic_objective_platform_summary(evaluation_trace),
            "diagnostic_interval": self.config.diagnostic_interval,
            "gradient_norm_trace": gradient_norm_trace,
            "parameter_norm_trace": parameter_norm_trace,
            "backend": "torch",
            "lengthscale_status": "FIXED_CONFIG_NOT_LEARNED",
            "uncertainty_status": "INDUCING_WEIGHTS_ONLY",
            "config_hash": config_hash,
            "environment_hash": environment_hash,
            "input_hash": input_hash,
            "objective_input_hash": objective_input_hash,
            "runtime_environment_hash": runtime_environment_hash,
            "environment_hash_override": self.config.environment_hash is not None,
            "gene_batch_size": self.config.gene_batch_size,
            "gene_batch_scaling": "uniform_without_replacement_ht",
            "global_spatial_kl_scaling": "unscaled_once_per_step",
            "gp_kl": "diagonal_variational_inducing_weights",
            "gene_sampler": "stateless_per_step_seeded_by_seed_plus_step",
            "resume_checkpoint_dir": str(resume_dir) if resume_dir is not None else None,
            **(
                continuous_factor_diagnostics(loading_np, field_matrix, signed=False)
                if k else {"n_factors": 0, "effective_rank_entropy": 0.0,
                           "effective_rank_participation": 0.0,
                           "factor_energy": [], "loading_cosine": [],
                           "field_cosine": []}
            ),
        }
        self.gene_id_ = tuple(genes)
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
        self._require_fitted()
        return {
            "schema": "r04.frozen_model.v2_torch",
            "backend": "torch",
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
        if state.get("schema") != "r04.frozen_model.v2_torch" or state.get("model_id") != cls.model_id:
            raise ValueError("unsupported mNSF frozen model artifact")
        if state.get("backend") != "torch":
            raise ValueError("unsupported frozen model backend")
        config = dict(state["config"])  # type: ignore[arg-type]
        if config.get("gp_parameterization") != "inducing_values_v2":
            raise ValueError("legacy mNSF frozen model lacks the D-050 GP parameterization")
        legacy_execution_mode = config.get("execution_mode")
        if legacy_execution_mode == "compiled":
            # Historical Torch artifacts used this label before a compiled
            # execution path existed.  Preserve the provenance while mapping
            # the false label to the only implemented, mathematically exact
            # path.  New callers requesting ``compiled`` still fail closed.
            config["execution_mode"] = "eager"
        estimator = cls(MNSFConfig(**config))
        if legacy_execution_mode == "compiled":
            estimator.execution_mode_requested_legacy_ = "compiled"
            estimator.execution_mode_migration_ = "LEGACY_FALSE_LABEL_TO_EAGER"
        parameters = dict(state["parameters"])  # type: ignore[arg-type]
        estimator.gene_id_ = tuple(str(v) for v in state["gene_id"])  # type: ignore[arg-type]
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
        estimator.diagnostics_ = dict(state.get("diagnostics", {}))  # type: ignore[arg-type]
        if estimator.execution_mode_migration_ is not None:
            estimator.diagnostics_.update({
                "execution_mode_requested_legacy": estimator.execution_mode_requested_legacy_,
                "execution_mode_migration": estimator.execution_mode_migration_,
                "execution_mode": "eager",
                "execution_mode_effective": "eager",
            })
        return estimator

    def subset_to_genes(self, genes: Sequence[str]) -> "MNSFEstimator":
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
        full_rank = int(np.linalg.matrix_rank(self.loading_)) if self.loading_.shape[1] else 0
        panel_rank = int(np.linalg.matrix_rank(projected.loading_)) if projected.loading_.shape[1] else 0
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
        if not sections:
            raise ValueError("at least one section is required")
        self._require_fitted()
        device = _get_device()
        execution_mode = _resolve_torch_execution_mode(self.config.execution_mode)
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
        basis_torch = [torch.tensor(b, dtype=torch.float32, device=device) for b in basis_list]
        prior_inv_torch = [torch.tensor(v, dtype=torch.float32, device=device) for v in k_inv_list]
        k = self.loading_.shape[1]
        g = counts.shape[1]
        q = self.config.nonspatial_rank
        q_loc = [torch.nn.Parameter(torch.zeros((b.shape[1], k), device=device, dtype=torch.float32)) for b in basis_list]
        q_scale_raw = [torch.nn.Parameter(torch.full((b.shape[1], k), -2.0, device=device, dtype=torch.float32)) for b in basis_list]
        raw_h = torch.nn.Parameter(torch.zeros((len(counts), q), device=device, dtype=torch.float32)) if q else None
        # Optimizer for inference
        infer_params: list[torch.nn.Parameter] = [*q_loc, *q_scale_raw]
        if raw_h is not None:
            infer_params.append(raw_h)
        infer_params = [p for p in infer_params if p.numel() > 0]
        optimizer = torch.optim.Adam(infer_params, lr=float(self.config.learning_rate)) if infer_params else None
        n_steps = self.config.steps if steps is None else int(steps)
        if n_steps < 1:
            raise ValueError("inference steps must be positive")
        seed_everything(self.config.seed + 5000)
        y_torch = torch.tensor(counts, dtype=torch.float32, device=device)
        log_l_torch = torch.tensor(np.log(library), dtype=torch.float32, device=device)
        loading_t = torch.tensor(self.loading_, dtype=torch.float32, device=device)
        amp_t = torch.tensor(self.factor_amplitude_, dtype=torch.float32, device=device) if k else torch.zeros((0,), device=device, dtype=torch.float32)
        baseline_t = torch.tensor(self.gene_baseline_, dtype=torch.float32, device=device)
        dispersion_t = torch.tensor(self.dispersion_, dtype=torch.float32, device=device)
        nonspatial_loading_t = (
            torch.tensor(self.nonspatial_loading_, dtype=torch.float32, device=device)
            if self.nonspatial_loading_ is not None else None
        )
        best = float("inf")
        best_state: list[torch.Tensor] | None = None
        losses: list[float] = []
        # Use all params list for snapshot
        all_infer_params = [*q_loc, *q_scale_raw]
        if raw_h is not None:
            all_infer_params.append(raw_h)
        for step in range(n_steps):
            if optimizer is not None:
                optimizer.zero_grad()
            # forward
            samples: list[torch.Tensor] = []
            kl = torch.tensor(0.0, device=device)
            for idx, loc in enumerate(q_loc):
                scale = F.softplus(q_scale_raw[idx]) + 1e-4
                eps = _stateless_normal(loc.shape, self.config.seed + 5000 + step * 1000 + idx + 1, device, loc.dtype)
                samples.append(loc + scale * eps)
                kl = kl + _diagonal_gp_kl_torch(loc, scale, prior_inv_torch[idx], logdet_prior_list[idx])
            field_parts = [basis_torch[idx] @ samples[idx] for idx in range(len(sections))] if k > 0 else [torch.zeros((basis_torch[idx].shape[0], 0), device=device) for idx in range(len(sections))]
            centered = [p - p.mean(dim=0, keepdim=True) if p.shape[1] > 0 else p for p in field_parts]
            field = torch.cat(centered, dim=0) if centered else torch.zeros((len(counts), 0), device=device)
            components: list[torch.Tensor] = []
            if k:
                log_comp = torch.log(loading_t.clamp(min=1e-8))[None, :, :] + torch.log(amp_t.clamp(min=1e-8))[None, None, :] + field[:, None, :]
                components.append(torch.logsumexp(log_comp, dim=2))
            if q and raw_h is not None and nonspatial_loading_t is not None:
                h = F.softplus(raw_h) + 1e-5
                components.append(torch.log((h @ nonspatial_loading_t.T).clamp(min=1e-8)))
            components.append(torch.log(baseline_t.clamp(min=1e-8))[None, :].expand(len(counts), g))
            log_mu = log_l_torch[:, None] + torch.logsumexp(torch.stack(components, dim=2), dim=2)
            mu = torch.exp(log_mu)
            theta_exp = dispersion_t[None, :]
            lp = _nb_log_prob(y_torch, mu, theta_exp)
            likelihood = gene_batch_sum_objective(-lp, g, g)
            loss = likelihood + 1e-4 * kl / float(max(len(counts), 1))
            value = float(loss.item())
            losses.append(value)
            if value < best:
                best = value
                best_state = [p.detach().cpu().clone() for p in all_infer_params]
            if optimizer is not None:
                loss.backward()
                _maybe_raise_torch_interrupt()
                optimizer.step()
        if best_state is None:
            raise RuntimeError("frozen mNSF inference produced no state")
        for p, b in zip(all_infer_params, best_state):
            p.data.copy_(b.to(p.device))
        self.inference_nonspatial_component_ = (
            (F.softplus(raw_h) + 1e-5).detach().cpu().numpy()
            if raw_h is not None else None
        )
        self.inference_diagnostics_ = {
            "loss_best": float(best),
            "loss_trace": [float(v) for v in losses],
            "steps": n_steps,
            "nonspatial_rank": q,
            "execution_mode": execution_mode,
            "execution_mode_effective": execution_mode,
            "execution_mode_requested": self.config.execution_mode,
        }
        draws_requested = self.config.posterior_draws if posterior_draws is None else int(posterior_draws)
        if draws_requested < 1:
            raise ValueError("posterior_draws must be positive")
        input_hash = sections_content_hash(sections)
        result: list[list[FieldFit]] = []
        for index, section in enumerate(sections):
            loc_np = q_loc[index].detach().cpu().numpy()
            scale_np = F.softplus(q_scale_raw[index]).detach().cpu().numpy() + 1e-4
            basis_np = basis_list[index]
            mean = basis_np @ loc_np
            mean -= mean.mean(axis=0, keepdims=True)
            rng = np.random.default_rng(self.config.seed + 5100 + index)
            draws = np.stack([
                basis_np @ (loc_np + scale_np * rng.normal(size=loc_np.shape))
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
                        "loss_trace": [float(v) for v in losses],
                        "patient_id": section.patient_id,
                    },
                )
                for factor in range(k)
            ])
        return result
