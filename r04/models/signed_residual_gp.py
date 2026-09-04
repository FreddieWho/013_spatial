"""Signed spatial GP on a negative-binomial nuisance residual (PyTorch)."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import json
from typing import Sequence

import numpy as np
import torch
import torch.nn.functional as F

from ..runtime import atomic_json, config_fingerprint, seed_everything
from ..io_contract import sections_content_hash
from ..spatial import normalize_coordinates
from ..types import FieldFit, SectionData
from ..diagnostics import continuous_factor_diagnostics
from ..objectives import gene_batch_sum_objective
from .mnsf import _diagonal_gp_kl_torch, _inducing_projection, _get_device, _maybe_raise_torch_interrupt, _stateless_normal


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
    resume_checkpoint_dir: str | None = None
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
        device = _get_device()
        seed_everything(self.config.seed)
        config_hash = self.config.config_hash or config_fingerprint(
            self.config, exclude={"checkpoint_dir", "resume_checkpoint_dir", "checkpoint_steps", "config_hash", "environment_hash"}
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
        library = np.maximum(library, 1.0)
        self.reference_profile_ = (
            counts / library[:, None]
        ).mean(axis=0).astype(float)
        residual_np = self._pearson_residual(counts, section_sizes=library, profile=self.reference_profile_).astype(np.float32)
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
        # torch tensors
        basis_torch = [torch.tensor(b, dtype=torch.float32, device=device) for b in basis_list]
        prior_inv_torch = [torch.tensor(v, dtype=torch.float32, device=device) for v in k_inv_list]
        residual_torch = torch.tensor(residual_np, dtype=torch.float32, device=device)
        n, g = residual_np.shape
        k = self.config.factors
        if self.config.gp_parameterization != "inducing_values_v2":
            raise ValueError("unsupported GP parameterization")
        # Parameters: raw_lambda (g,k) init N(0,0.05) seeded, q_loc/q_scale_raw per section (m,k) zeros/-2.0
        gen = torch.Generator(device=device)
        gen.manual_seed(int(self.config.seed) & 0xFFFFFFFFFFFFFFFF)
        raw_lambda_init = torch.randn((g, k), generator=gen, device=device, dtype=torch.float32) * 0.05 if k > 0 and g > 0 and k * g > 0 else torch.zeros((g, k), device=device, dtype=torch.float32)
        # For k==0 edge, produce (g,0) zeros
        if k == 0:
            raw_lambda_init = torch.zeros((g, 0), device=device, dtype=torch.float32)
        raw_lambda = torch.nn.Parameter(raw_lambda_init)
        q_loc: list[torch.nn.Parameter] = []
        q_scale_raw: list[torch.nn.Parameter] = []
        for idx in range(len(basis_list)):
            m = basis_list[idx].shape[1]
            q_loc.append(torch.nn.Parameter(torch.zeros((m, k), device=device, dtype=torch.float32)))
            q_scale_raw.append(torch.nn.Parameter(torch.full((m, k), -2.0, device=device, dtype=torch.float32)))
        all_params: list[torch.nn.Parameter] = [raw_lambda, *q_loc, *q_scale_raw]
        optim_params = [p for p in all_params if p.numel() > 0]
        optimizer = torch.optim.Adam(optim_params, lr=float(self.config.learning_rate))

        global_step = 0
        best_loss_val = float("inf")
        best_step = -1
        best_state: dict[str, torch.Tensor] | None = None

        def _snapshot_best() -> dict[str, torch.Tensor]:
            snap: dict[str, torch.Tensor] = {}
            for idx_p, p in enumerate(all_params):
                snap[f"p_{idx_p}"] = p.detach().cpu().clone()
            return snap

        def _restore_best(snap: dict[str, torch.Tensor]) -> None:
            for idx_p, p in enumerate(all_params):
                key = f"p_{idx_p}"
                if key in snap:
                    p.data.copy_(snap[key].to(p.device))

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
            ckpt_path = manager_dir / "checkpoint.pt"
            if self.config.resume_checkpoint_dir is not None:
                resume_dir = Path(self.config.resume_checkpoint_dir)
                resume_meta = resume_dir / "checkpoint_metadata.json"
                if not resume_meta.exists():
                    raise RuntimeError("resume checkpoint metadata is missing")
                src_meta = json.loads(resume_meta.read_text(encoding="utf-8"))
                if src_meta.get("checkpoint_schema") != "r04.mnsf_checkpoint.v3_torch":
                    raise RuntimeError("resume checkpoint schema mismatch")
                if src_meta.get("input_hash") != input_hash:
                    raise RuntimeError("resume checkpoint input hash mismatch")
                if src_meta.get("environment_hash") != environment_hash:
                    raise RuntimeError("resume checkpoint environment hash mismatch")
                src_ckpt = resume_dir / "checkpoint.pt"
                if not src_ckpt.exists():
                    raise RuntimeError("resume checkpoint is missing")
                payload = torch.load(str(src_ckpt), map_location=device)
                for idx_p, p in enumerate(all_params):
                    key = f"p_{idx_p}"
                    if key in payload.get("params", {}):
                        p.data.copy_(payload["params"][key].to(device))
                optimizer.load_state_dict(payload["optimizer_state"])
                for state in optimizer.state.values():
                    for kk, vv in state.items():
                        if isinstance(vv, torch.Tensor):
                            state[kk] = vv.to(device)
                global_step = int(payload.get("global_step", 0))
                best_loss_val = float(payload.get("best_loss", float("inf")))
                best_step = int(payload.get("best_step", -1))
                best_state = payload.get("best_state")
            elif ckpt_path.exists():
                payload = torch.load(str(ckpt_path), map_location=device)
                if payload.get("checkpoint_schema") != "r04.mnsf_checkpoint.v3_torch" or payload.get("backend") != "torch":
                    raise RuntimeError("checkpoint schema/backend mismatch; start a new R04 run")
                if payload.get("input_hash") != input_hash or payload.get("config_hash") != config_hash or payload.get("environment_hash") != environment_hash:
                    raise RuntimeError("checkpoint hash mismatch; start a new R04 run")
                for idx_p, p in enumerate(all_params):
                    key = f"p_{idx_p}"
                    if key in payload.get("params", {}):
                        p.data.copy_(payload["params"][key].to(device))
                optimizer.load_state_dict(payload["optimizer_state"])
                for state in optimizer.state.values():
                    for kk, vv in state.items():
                        if isinstance(vv, torch.Tensor):
                            state[kk] = vv.to(device)
                global_step = int(payload.get("global_step", 0))
                best_loss_val = float(payload.get("best_loss", float("inf")))
                best_step = int(payload.get("best_step", -1))
                best_state = payload.get("best_state")

        losses: list[float] = []
        start_step = int(global_step)
        for step in range(start_step, self.config.steps):
            batch_size = g if self.config.gene_batch_size is None else int(self.config.gene_batch_size)
            if batch_size < 1:
                raise ValueError("gene_batch_size must be positive")
            batch_size = min(batch_size, g)
            if batch_size == g:
                gene_index = np.arange(g, dtype=np.int64)
            else:
                gene_index = np.random.default_rng(self.config.seed + 9001 + step).choice(
                    g, size=batch_size, replace=False
                ).astype(np.int64)
            gene_index_t = torch.tensor(gene_index, dtype=torch.long, device=device)
            optimizer.zero_grad()
            # sample GP inducing values
            field_parts: list[torch.Tensor] = []
            kl = torch.tensor(0.0, device=device)
            for index, basis in enumerate(basis_torch):
                loc = q_loc[index]
                scale = F.softplus(q_scale_raw[index]) + 1e-4
                eps = _stateless_normal(loc.shape, self.config.seed + 7 + step * 1000 + index + 1, device, loc.dtype)
                u = loc + scale * eps
                value = basis @ u
                value = value - value.mean(dim=0, keepdim=True)
                field_parts.append(value)
                kl = kl + _diagonal_gp_kl_torch(loc, scale, prior_inv_torch[index], logdet_prior_list[index])
            field = torch.cat(field_parts, dim=0) if field_parts else torch.zeros((n, 0), device=device)
            loading = F.normalize(raw_lambda, p=2, dim=0) if raw_lambda.numel() > 0 and k > 0 else raw_lambda
            # loading shape (g,k)
            loading_batch = loading[gene_index_t, :] if k > 0 else loading
            predicted = field @ loading_batch.T if k > 0 else torch.zeros((n, batch_size), device=device, dtype=torch.float32)
            residual_batch = residual_torch[:, gene_index_t] if batch_size != g else residual_torch
            # residual_batch shape (n, batch_size), predicted same
            likelihood = gene_batch_sum_objective(
                (residual_batch - predicted) ** 2,
                g,
                batch_size,
            )
            loss = likelihood + 1e-4 * kl / float(max(n, 1))
            value = float(loss.item())
            if value < best_loss_val:
                best_loss_val = value
                best_step = step + 1
                best_state = _snapshot_best()
            losses.append(value)
            loss.backward()
            _maybe_raise_torch_interrupt()
            optimizer.step()
            global_step = step + 1
            if manager_dir is not None and (step + 1) % max(1, self.config.checkpoint_steps) == 0:
                payload = {
                    "checkpoint_schema": "r04.mnsf_checkpoint.v3_torch",
                    "backend": "torch",
                    "input_hash": input_hash,
                    "config_hash": config_hash,
                    "environment_hash": environment_hash,
                    "seed": self.config.seed,
                    "factors": k,
                    "global_step": global_step,
                    "best_loss": best_loss_val,
                    "best_step": best_step,
                    "best_state": best_state,
                    "params": {f"p_{idx}": p.detach().cpu().clone() for idx, p in enumerate(all_params)},
                    "optimizer_state": optimizer.state_dict(),
                }
                torch.save(payload, str(manager_dir / "checkpoint.pt"))
        if manager_dir is not None:
            payload = {
                "checkpoint_schema": "r04.mnsf_checkpoint.v3_torch",
                "backend": "torch",
                "input_hash": input_hash,
                "config_hash": config_hash,
                "environment_hash": environment_hash,
                "seed": self.config.seed,
                "factors": k,
                "global_step": global_step,
                "best_loss": best_loss_val,
                "best_step": best_step,
                "best_state": best_state,
                "params": {f"p_{idx}": p.detach().cpu().clone() for idx, p in enumerate(all_params)},
                "optimizer_state": optimizer.state_dict(),
            }
            torch.save(payload, str(manager_dir / "checkpoint.pt"))
        if not np.isfinite(best_loss_val) or best_state is None:
            raise RuntimeError("signed residual GP optimizer produced no persisted best state")
        _restore_best(best_state)
        with torch.no_grad():
            loading_np = F.normalize(raw_lambda, p=2, dim=0).detach().cpu().numpy() if k > 0 else np.zeros((g, 0), dtype=float)
        self.loading_ = loading_np
        self.fields_ = []
        for index, section in enumerate(sections):
            loc = q_loc[index].detach().cpu().numpy()
            scale = F.softplus(q_scale_raw[index]).detach().cpu().numpy() + 1e-4
            mean = basis_list[index] @ loc
            mean -= mean.mean(axis=0, keepdims=True)
            rng = np.random.default_rng(self.config.seed + 100 + index)
            draws = np.stack([basis_list[index] @ (loc + scale * rng.normal(size=loc.shape)) for _ in range(self.config.posterior_draws)], axis=0)
            draws -= draws.mean(axis=1, keepdims=True)
            self.fields_.append([
                FieldFit(self.model_id, f"SIGNED_{factor+1:02d}", section.section_id,
                         loading_np[:, factor], mean[:, factor], draws[:, :, factor].std(axis=0),
                         self.config.lengthscale, float(np.var(mean[:, factor])),
                         input_hash=input_hash, diagnostics={"uncertainty_kind": "variational_inducing_weights", "nuisance": "NB_PEARSON", "steps": self.config.steps, "patient_id": section.patient_id})
                for factor in range(k)
            ])
        field_matrix = np.column_stack([
            np.concatenate([group[factor].field_mean for group in self.fields_])
            for factor in range(k)
        ]) if k > 0 else np.empty((n, 0), dtype=float)
        # Determine best loss trace semantics: best is min(losses) but now best_state restores best
        loss_best = best_loss_val if np.isfinite(best_loss_val) else (min(losses) if losses else None)
        self.diagnostics_ = {
            "backend": "torch", "steps": self.config.steps, "nuisance": "NB_PEARSON",
            "likelihood": "MSE_ON_NB_PEARSON_RESIDUAL", "nuisance_overdispersion": 20.0,
            "lengthscale_status": "FIXED_CONFIG_NOT_LEARNED",
            "config_hash": config_hash, "environment_hash": environment_hash,
            "loss_final": losses[-1] if losses else None,
            "loss_best": loss_best,
            "best_step": int(best_step) if best_step >= 0 else None,
            "loss_trace": [float(value) for value in losses],
            "gene_batch_size": self.config.gene_batch_size,
            "gene_batch_scaling": "uniform_without_replacement_ht",
            "global_spatial_kl_scaling": "unscaled_once_per_step",
            "gp_kl": "diagonal_variational_inducing_weights",
            "gene_sampler": "stateless_per_step_seeded_by_seed_plus_step",
            "start_step": start_step,
            "resume_checkpoint_dir": str(resume_dir) if resume_dir is not None else None,
            **continuous_factor_diagnostics(loading_np, field_matrix, signed=True),
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
            "schema": "r04.frozen_model.v2_torch",
            "backend": "torch",
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
        schema = state.get("schema")
        if schema == "r04.frozen_model.v1":
            raise ValueError(
                "legacy TensorFlow frozen model (r04.frozen_model.v1) is not readable by the PyTorch runtime; "
                "use historical evidence only"
            )
        if schema != "r04.frozen_model.v2_torch" or state.get("model_id") != cls.model_id:
            raise ValueError("unsupported signed residual GP frozen model artifact")
        if state.get("backend") is not None and state.get("backend") != "torch":
            raise ValueError("unsupported frozen model backend")
        config = dict(state["config"])  # type: ignore[arg-type]
        if config.get("gp_parameterization") != "inducing_values_v2":
            raise ValueError("legacy signed frozen model lacks the D-050 GP parameterization")
        estimator = cls(SignedResidualGPConfig(**config))
        estimator.gene_id_ = tuple(str(value) for value in state.get("gene_id", ()))
        parameters = dict(state["parameters"])  # type: ignore[arg-type]
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
        estimator.diagnostics_ = dict(state.get("diagnostics", {}))  # type: ignore[arg-type]
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
        device = _get_device()
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
        library = np.maximum(library, 1.0)
        residual = self._pearson_residual(
            counts, section_sizes=library, profile=self.reference_profile_
        ).astype(np.float32)
        residual_torch = torch.tensor(residual, dtype=torch.float32, device=device)
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
        basis_torch = [torch.tensor(b, dtype=torch.float32, device=device) for b in basis_list]
        prior_inv_torch = [torch.tensor(v, dtype=torch.float32, device=device) for v in k_inv_list]
        k = self.loading_.shape[1]
        g = counts.shape[1]
        n = counts.shape[0]
        loading_torch = torch.tensor(self.loading_, dtype=torch.float32, device=device)
        q_loc = [torch.nn.Parameter(torch.zeros((basis.shape[1], k), device=device, dtype=torch.float32)) for basis in basis_list]
        q_scale_raw = [torch.nn.Parameter(torch.full((basis.shape[1], k), -2.0, device=device, dtype=torch.float32)) for basis in basis_list]
        all_infer_params: list[torch.nn.Parameter] = [*q_loc, *q_scale_raw]
        optim_params = [p for p in all_infer_params if p.numel() > 0]
        optimizer = torch.optim.Adam(optim_params, lr=float(self.config.learning_rate)) if optim_params else None
        n_steps = self.config.steps if steps is None else int(steps)
        if n_steps < 1:
            raise ValueError("inference steps must be positive")
        seed_everything(self.config.seed + 6000)
        best = float("inf")
        best_state: list[torch.Tensor] | None = None
        losses: list[float] = []
        for step in range(n_steps):
            if optimizer is not None:
                optimizer.zero_grad()
            fields: list[torch.Tensor] = []
            kl = torch.tensor(0.0, device=device)
            for index, basis in enumerate(basis_torch):
                loc = q_loc[index]
                scale = F.softplus(q_scale_raw[index]) + 1e-4
                eps = _stateless_normal(loc.shape, self.config.seed + 6000 + step * 1000 + index + 1, device, loc.dtype)
                u = loc + scale * eps
                value = basis @ u
                value = value - value.mean(dim=0, keepdim=True)
                fields.append(value)
                kl = kl + _diagonal_gp_kl_torch(loc, scale, prior_inv_torch[index], logdet_prior_list[index])
            field = torch.cat(fields, dim=0) if fields else torch.zeros((n, 0), device=device)
            predicted = field @ loading_torch.T if k > 0 else torch.zeros((n, g), device=device, dtype=torch.float32)
            likelihood = gene_batch_sum_objective(
                (residual_torch - predicted) ** 2,
                g,
                g,
            )
            loss = likelihood + 1e-4 * kl / float(max(n, 1))
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
            raise RuntimeError("frozen signed inference produced no state")
        for p, b in zip(all_infer_params, best_state):
            p.data.copy_(b.to(p.device))
        draws_requested = self.config.posterior_draws if posterior_draws is None else int(posterior_draws)
        if draws_requested < 1:
            raise ValueError("posterior_draws must be positive")
        input_hash = sections_content_hash(sections)
        result: list[list[FieldFit]] = []
        for index, section in enumerate(sections):
            loc = q_loc[index].detach().cpu().numpy()
            scale = F.softplus(q_scale_raw[index]).detach().cpu().numpy() + 1e-4
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
                        "loss_trace": [float(v) for v in losses],
                        "patient_id": section.patient_id,
                    },
                )
                for factor in range(k)
            ])
        return result
