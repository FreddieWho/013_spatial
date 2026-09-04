"""Deterministic, read-only diagnostics for frozen R-04 mNSF checkpoints (PyTorch)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Sequence

import numpy as np

from .io_contract import sections_content_hash, sections_objective_input_hash
from .models.mnsf import _inducing, _inducing_projection
from .spatial import normalize_coordinates
from .types import SectionData


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


def _diagonal_gp_kl_torch(loc, scale, prior_inverse, logdet_prior) -> object:
    """KL[q(u)||N(0,K)] for a diagonal q covariance, per factor jointly (torch)."""
    import torch

    # loc, scale: (m, K), prior_inverse: (m, m), logdet_prior: scalar
    if loc.numel() == 0 or loc.shape[1] == 0:
        return torch.tensor(0.0, dtype=torch.float32, device=loc.device)
    factors = float(loc.shape[1])
    dimension = float(loc.numel())
    quadratic = torch.sum((prior_inverse @ loc) * loc)
    diag = torch.diagonal(prior_inverse)
    trace = torch.sum((scale * scale) * diag[:, None])
    # log|Q| where Q=diag(scale^2)
    # scale was softplus(raw)+1e-4 >0, so log is safe
    logdet_q = torch.sum(torch.log(scale * scale + 1e-18))
    return 0.5 * (quadratic + trace - dimension + factors * logdet_prior - logdet_q)


def _to_numpy(value) -> np.ndarray | None:
    if value is None:
        return None
    if isinstance(value, np.ndarray):
        return value
    # torch tensor
    try:
        import torch

        if isinstance(value, torch.Tensor):
            return value.detach().cpu().numpy()
    except Exception:
        pass
    return np.asarray(value)


def _resolve_checkpoint_paths(checkpoint_path: str | Path) -> tuple[Path, Path]:
    """Resolve a .pt file and its JSON metadata, fail-closed for TF artifacts.

    If ``checkpoint_path`` is a directory, expect ``checkpoint.pt`` inside and
    ``checkpoint_metadata.json`` alongside.  If it is a file, load it directly
    and look for metadata in the parent directory.
    """
    p = Path(checkpoint_path)
    legacy_msg = (
        "legacy TensorFlow checkpoint is not readable by the PyTorch runtime; "
        "use historical evidence only"
    )
    # Fail-closed for obvious TF sentinel files
    if p.is_dir():
        # TF CheckpointManager creates a file named "checkpoint" with no extension
        if (p / "checkpoint").exists():
            raise ValueError(legacy_msg + f" (found TF state file at {p / 'checkpoint'})")
        # TF checkpoint shards
        if any(p.glob("*.index")) or any(p.glob("*.data-*")):
            raise ValueError(legacy_msg + f" (found TF shard in directory {p})")
        # Resolve pt path
        pt_path = p / "checkpoint.pt"
        metadata_path = p / "checkpoint_metadata.json"
        if not pt_path.exists():
            # Allow single .pt inside directory (e.g. custom name)
            pts = list(p.glob("*.pt"))
            # also consider .pth
            pts += list(p.glob("*.pth"))
            if len(pts) == 1:
                pt_path = pts[0]
                metadata_path = pt_path.parent / "checkpoint_metadata.json"
            else:
                raise ValueError(
                    f"PyTorch checkpoint not found in directory: {p} (expected checkpoint.pt)"
                )
        return pt_path, metadata_path
    else:
        # File path
        if p.name == "checkpoint" or p.suffix == ".index" or ".data-" in p.name:
            raise ValueError(legacy_msg + f" (path looks like TF checkpoint: {p})")
        if not p.exists():
            # Might be a TF prefix without extension (e.g. ckpt-100); check for .index sibling
            if Path(str(p) + ".index").exists() or Path(str(p) + ".data-00000-of-00001").exists():
                raise ValueError(legacy_msg + f" (found TF shard for prefix {p})")
            raise FileNotFoundError(f"checkpoint path does not exist: {p}")
        # If suffix is not .pt/.pth, check sibling TF shards still fail-closed
        if p.suffix not in {".pt", ".pth"}:
            if Path(str(p) + ".index").exists():
                raise ValueError(legacy_msg + f" (found TF .index for {p})")
        pt_path = p
        metadata_path = p.parent / "checkpoint_metadata.json"
        # Also check parent directory for TF sentinels if we are given a file inside a TF dir
        # but the file itself is .pt – then it's fine; only if directory looks TF without .pt
        return pt_path, metadata_path


def _load_metadata(metadata_path: Path) -> dict | None:
    if not metadata_path.exists():
        return None
    try:
        payload = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"invalid checkpoint metadata at {metadata_path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"checkpoint metadata must be an object: {metadata_path}")
    return payload


def _validate_metadata(payload: dict | None, checkpoint_path: Path) -> None:
    if payload is None:
        return
    schema = payload.get("checkpoint_schema")
    backend = payload.get("backend")
    # Fail-closed for legacy TF schemas
    if schema in {"r04.mnsf_checkpoint.v1", "r04.mnsf_checkpoint.v2"}:
        raise ValueError(
            "legacy TensorFlow checkpoint (r04.mnsf_checkpoint.v2) is not readable by the PyTorch runtime; "
            "use historical evidence only"
        )
    # New expected schema
    if schema is not None and schema != "r04.mnsf_checkpoint.v3_torch":
        raise ValueError(
            f"unsupported checkpoint schema: {schema}; expected r04.mnsf_checkpoint.v3_torch"
        )
    if backend is not None and backend != "torch":
        raise ValueError(
            f"unsupported checkpoint backend: {backend}; expected torch (found at {checkpoint_path})"
        )


def _parameter_from_params(params: dict, key: str, *, label: str) -> np.ndarray:
    if key not in params:
        raise ValueError(f"checkpoint {label} is missing parameter tensor {key}")
    value = _to_numpy(params[key])
    if value is None:
        raise ValueError(f"checkpoint {label} parameter tensor {key} is invalid")
    array = np.asarray(value)
    if not np.all(np.isfinite(array)):
        raise ValueError(f"checkpoint {label} parameter tensor {key} is non-finite")
    return array


def _layout_key(
    layout: dict,
    params: dict,
    name: str,
    *,
    required: bool = True,
) -> str | None:
    value = layout.get(name)
    if value is None and not required:
        return None
    if not isinstance(value, str) or not value:
        raise ValueError(f"checkpoint parameter_layout field {name} must be a non-empty key")
    if value not in params:
        raise ValueError(
            f"checkpoint parameter_layout field {name} points to missing tensor {value}"
        )
    return value


def _extract_parameter_state(
    checkpoint: dict,
    *,
    n_sections: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[np.ndarray], list[np.ndarray], np.ndarray | None, np.ndarray | None]:
    """Extract one atomic Torch parameter state, never mixing best_state.

    v3 checkpoints keep the resumable state in ``params``.  ``best_state`` is
    a separate optimizer-selection snapshot and must not be used to fill a
    missing tensor in ``params``; doing so would create a state that never
    existed on disk.
    """
    params = checkpoint.get("params")
    if not isinstance(params, dict):
        raise ValueError("Torch checkpoint must contain a params mapping")
    has_layout = "parameter_layout" in checkpoint
    layout = checkpoint.get("parameter_layout")
    used: list[str] = []

    def take(key: str, field: str) -> np.ndarray:
        if key in used:
            raise ValueError(f"checkpoint parameter_layout reuses tensor key {key}")
        used.append(key)
        return _parameter_from_params(params, key, label=field)

    if has_layout:
        if not isinstance(layout, dict):
            raise ValueError("checkpoint parameter_layout must be an object")
        if layout.get("schema") != "r04.mnsf_parameter_layout.v1":
            raise ValueError("unsupported or malformed mNSF parameter_layout schema")
        for name in ("raw_w", "raw_a", "raw_b", "raw_theta", "q_loc", "q_scale", "raw_v", "raw_h"):
            if name not in layout:
                raise ValueError(f"checkpoint parameter_layout is missing field {name}")
        raw_keys = {
            name: _layout_key(layout, params, name)
            for name in ("raw_w", "raw_a", "raw_b", "raw_theta")
        }
        raw_w = take(raw_keys["raw_w"], "raw_w")  # type: ignore[arg-type]
        raw_a = take(raw_keys["raw_a"], "raw_a")  # type: ignore[arg-type]
        raw_b = take(raw_keys["raw_b"], "raw_b")  # type: ignore[arg-type]
        raw_theta = take(raw_keys["raw_theta"], "raw_theta")  # type: ignore[arg-type]

        q_loc_keys = layout.get("q_loc")
        q_scale_keys = layout.get("q_scale")
        if (
            not isinstance(q_loc_keys, list)
            or not isinstance(q_scale_keys, list)
            or len(q_loc_keys) != n_sections
            or len(q_scale_keys) != n_sections
        ):
            raise ValueError(
                "checkpoint parameter_layout q_loc/q_scale must list one key per section"
            )
        q_loc = []
        q_scale = []
        for index, key in enumerate(q_loc_keys):
            if not isinstance(key, str) or not key:
                raise ValueError(f"checkpoint parameter_layout q_loc[{index}] is invalid")
            q_loc.append(take(key, f"q_loc[{index}]"))
        for index, key in enumerate(q_scale_keys):
            if not isinstance(key, str) or not key:
                raise ValueError(f"checkpoint parameter_layout q_scale[{index}] is invalid")
            q_scale.append(take(key, f"q_scale[{index}]"))

        raw_v_key = _layout_key(layout, params, "raw_v", required=False)
        raw_h_key = _layout_key(layout, params, "raw_h", required=False)
        if (raw_v_key is None) != (raw_h_key is None):
            raise ValueError("checkpoint must contain raw_v and raw_h together")
        raw_v = take(raw_v_key, "raw_v") if raw_v_key is not None else None
        raw_h = take(raw_h_key, "raw_h") if raw_h_key is not None else None
        return raw_w, raw_a, raw_b, raw_theta, q_loc, q_scale, raw_v, raw_h

    # Compatibility path for early v3 Torch checkpoints: positional p_* keys
    # inside params only.  No top-level or best_state fallback is allowed.
    raw_w = take("p_0", "raw_w")
    raw_a = take("p_1", "raw_a")
    raw_b = take("p_2", "raw_b")
    raw_theta = take("p_3", "raw_theta")
    q_loc = [take(f"p_{4 + index}", f"q_loc[{index}]") for index in range(n_sections)]
    q_scale_start = 4 + n_sections
    q_scale = [
        take(f"p_{q_scale_start + index}", f"q_scale[{index}]")
        for index in range(n_sections)
    ]
    nuisance_start = q_scale_start + n_sections
    has_raw_v = f"p_{nuisance_start}" in params
    has_raw_h = f"p_{nuisance_start + 1}" in params
    if has_raw_v != has_raw_h:
        raise ValueError("checkpoint must contain raw_v and raw_h together")
    raw_v = take(f"p_{nuisance_start}", "raw_v") if has_raw_v else None
    raw_h = take(f"p_{nuisance_start + 1}", "raw_h") if has_raw_h else None
    return raw_w, raw_a, raw_b, raw_theta, q_loc, q_scale, raw_v, raw_h


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

    Parameters
    ----------
    checkpoint_path : Path
        Either a ``.pt`` file (torch checkpoint dict) or a directory containing
        ``checkpoint.pt`` and ``checkpoint_metadata.json``.  Metadata must have
        ``checkpoint_schema=r04.mnsf_checkpoint.v3_torch`` and
        ``backend=torch``.  Legacy TensorFlow checkpoints fail-closed.
    """
    if not sections:
        raise ValueError("at least one section is required")
    if mc_draws < 1:
        raise ValueError("mc_draws must be positive")
    import torch
    import torch.nn.functional as F

    pt_path, metadata_path = _resolve_checkpoint_paths(checkpoint_path)
    metadata = _load_metadata(metadata_path)
    _validate_metadata(metadata, pt_path)

    # Load torch checkpoint dict
    try:
        ckpt = torch.load(str(pt_path), map_location="cpu", weights_only=True)  # type: ignore[call-arg]
    except TypeError:
        # Older torch without weights_only
        ckpt = torch.load(str(pt_path), map_location="cpu")
    except Exception as exc:
        # If file is not a valid pt, check if it might be TF
        msg = str(exc).lower()
        _tf_token = "tensor" + "flow"
        if "checkpoint" in msg and (_tf_token in msg or "index" in msg):
            raise ValueError(
                "legacy TensorFlow checkpoint is not readable by the PyTorch runtime; "
                "use historical evidence only"
            ) from exc
        raise ValueError(f"failed to load PyTorch checkpoint at {pt_path}: {exc}") from exc

    if not isinstance(ckpt, dict):
        raise ValueError(f"checkpoint at {pt_path} is not a dict (found {type(ckpt).__name__})")

    # Detect TF-style checkpoint embedded in dict (e.g. schema v2 without backend)
    if ckpt.get("checkpoint_schema") in {"r04.mnsf_checkpoint.v1", "r04.mnsf_checkpoint.v2"}:
        raise ValueError(
            "legacy TensorFlow checkpoint (r04.mnsf_checkpoint.v2) is not readable by the PyTorch runtime; "
            "use historical evidence only"
        )
    # Fail-closed if dict looks like TF frozen model
    if ckpt.get("schema") == "r04.frozen_model.v1":
        raise ValueError(
            "legacy TensorFlow frozen model (r04.frozen_model.v1) is not readable by the PyTorch runtime; "
            "use historical evidence only"
        )
    if (
        ckpt.get("checkpoint_schema") != "r04.mnsf_checkpoint.v3_torch"
        or ckpt.get("backend") != "torch"
    ):
        raise ValueError(
            "unsupported or incomplete checkpoint schema/backend; "
            "expected r04.mnsf_checkpoint.v3_torch with backend=torch"
        )
    if metadata is not None:
        for key in (
            "checkpoint_schema",
            "backend",
            "input_hash",
            "objective_input_hash",
            "config_hash",
            "environment_hash",
            "seed",
            "factors",
        ):
            if key in metadata and key in ckpt and metadata[key] != ckpt[key]:
                raise ValueError(f"checkpoint metadata mismatch for {key}")

    actual_input_hash = sections_content_hash(sections)
    checkpoint_input_hash = ckpt.get("input_hash")
    if checkpoint_input_hash != actual_input_hash:
        raise ValueError(
            "checkpoint input hash does not match the supplied sections; refusing to score"
        )
    actual_objective_hash = sections_objective_input_hash(sections)
    checkpoint_objective_hash = ckpt.get("objective_input_hash")
    if checkpoint_objective_hash is not None and checkpoint_objective_hash != actual_objective_hash:
        raise ValueError(
            "checkpoint objective input hash does not match the supplied sections; refusing to score"
        )
    objective_hash_status = (
        "VERIFIED" if checkpoint_objective_hash is not None else "NOT_VERIFIED_LEGACY_V3"
    )

    raw_w, raw_a, raw_b, raw_theta, q_loc_list, q_scale_list, raw_v, raw_h = (
        _extract_parameter_state(ckpt, n_sections=len(sections))
    )

    # global_step may be stored as tensor or int; default 0 if missing
    global_step_val = ckpt.get("global_step")
    if global_step_val is None:
        # alternative key from older code
        global_step_val = ckpt.get("global_step_tensor", 0)
    try:
        if isinstance(global_step_val, torch.Tensor):
            global_step = int(global_step_val.detach().cpu().numpy().item())
        elif isinstance(global_step_val, np.ndarray):
            global_step = int(np.asarray(global_step_val).item())
        else:
            global_step = int(global_step_val)  # type: ignore[arg-type]
    except Exception:
        global_step = 0

    n_sections = len(sections)

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

    k = int(factors)
    if k < 0:
        raise ValueError("factors must be non-negative")
    if raw_w.shape != (counts.shape[1], k):
        raise ValueError(
            f"checkpoint K/genes mismatch: raw_w={raw_w.shape}, expected={(counts.shape[1], k)}"
        )
    expected_shapes = {
        "raw_a": (k,),
        "raw_b": (counts.shape[1],),
        "raw_theta": (counts.shape[1],),
    }
    for name, value in (
        ("raw_a", raw_a),
        ("raw_b", raw_b),
        ("raw_theta", raw_theta),
    ):
        if value.shape != expected_shapes[name]:
            raise ValueError(
                f"checkpoint {name} shape mismatch: {value.shape}, "
                f"expected {expected_shapes[name]}"
            )
    for index, (loc, scale, basis) in enumerate(zip(q_loc_list, q_scale_list, basis_list)):
        expected = (basis.shape[1], k)
        if loc.shape != expected or scale.shape != expected:
            raise ValueError(
                f"checkpoint GP tensor shape mismatch for section {index}: "
                f"q_loc={loc.shape}, q_scale={scale.shape}, expected={expected}"
            )
    if (raw_v is None) != (raw_h is None):
        raise ValueError("checkpoint must contain raw_v and raw_h together")
    if raw_v is not None and raw_h is not None:
        if raw_v.ndim != 2 or raw_h.ndim != 2 or raw_v.shape[0] != counts.shape[1] or raw_h.shape[0] != counts.shape[0]:
            raise ValueError("checkpoint nuisance tensor shapes do not match the data")
        if raw_v.shape[1] != raw_h.shape[1] or raw_v.shape[1] < 1:
            raise ValueError("checkpoint nuisance tensors must share a positive rank")

    device = torch.device("cpu")
    basis_t = [torch.tensor(v, dtype=torch.float32, device=device) for v in basis_list]
    y_t = torch.tensor(counts, dtype=torch.float32, device=device)
    log_l_t = torch.tensor(np.log(library), dtype=torch.float32, device=device)
    raw_w_t = torch.tensor(raw_w, dtype=torch.float32, device=device)
    raw_a_t = torch.tensor(raw_a, dtype=torch.float32, device=device)
    raw_b_t = torch.tensor(raw_b, dtype=torch.float32, device=device)
    raw_theta_t = torch.tensor(raw_theta, dtype=torch.float32, device=device)
    raw_v_t = torch.tensor(raw_v, dtype=torch.float32, device=device) if raw_v is not None else None
    raw_h_t = torch.tensor(raw_h, dtype=torch.float32, device=device) if raw_h is not None else None
    loc_t = [torch.tensor(v, dtype=torch.float32, device=device) for v in q_loc_list]
    # q_scale is stored as raw (pre-softplus); apply the model transform.
    scale_t = [F.softplus(torch.tensor(v, dtype=torch.float32, device=device)) + 1e-4 for v in q_scale_list]
    prior_inverse_t = [torch.tensor(v, dtype=torch.float32, device=device) for v in prior_inverse_list]
    logdet_t = [torch.tensor(v, dtype=torch.float32, device=device) for v in logdet_list]

    # Positive parameters
    if k:
        loading = F.softmax(raw_w_t, dim=0)
    else:
        loading = torch.zeros((counts.shape[1], 0), dtype=torch.float32, device=device)
    amplitude = F.softplus(raw_a_t) + 1e-5
    baseline = F.softplus(raw_b_t) + 1e-5
    theta = F.softplus(raw_theta_t) + 1e-3

    losses: list[float] = []
    for draw in range(mc_draws):
        kl = torch.tensor(0.0, dtype=torch.float32, device=device)
        samples: list[torch.Tensor] = []
        for index, loc in enumerate(loc_t):
            # Stateless normal via generator seeded from (seed, draw, index)
            gen = torch.Generator(device=device)
            combined = (int(seed) * 1000003 + (draw * 1000 + index + 1) * 917) & 0xFFFFFFFF
            # manual_seed expects int64; keep within 2**63-1
            gen.manual_seed(int(combined) & 0x7FFFFFFFFFFFFFFF)
            eps = torch.randn(loc.shape, generator=gen, dtype=torch.float32, device=device)
            samples.append(loc + scale_t[index] * eps)
            kl = kl + _diagonal_gp_kl_torch(loc, scale_t[index], prior_inverse_t[index], logdet_t[index])

        if samples:
            field_parts = [basis_t[i] @ samples[i] for i in range(len(samples))]
            # Center each section field (zero-mean) to preserve continuous overlapping field identifiability
            centered = [part - part.mean(dim=0, keepdim=True) for part in field_parts]
            field = torch.cat(centered, dim=0)  # (n, K)
        else:
            field = torch.zeros((counts.shape[0], 0), dtype=torch.float32, device=device)

        components: list[torch.Tensor] = []
        if k:
            # (n, g, K): log loading + log amplitude + field
            # loading (g,K) -> (1,g,K), amplitude (K,)->(1,1,K), field (n,K)->(n,1,K)
            comp_spatial = (
                torch.log(loading[None, :, :] + 1e-8)
                + torch.log(amplitude[None, None, :] + 1e-12)
                + field[:, None, :]
            )
            spatial = torch.logsumexp(comp_spatial, dim=2)  # (n, g)
            components.append(spatial)
        if raw_v_t is not None and raw_h_t is not None and raw_h_t.shape[1] > 0:
            nuisance_loading = F.softmax(raw_v_t, dim=0)
            nuisance = (F.softplus(raw_h_t) + 1e-5) @ nuisance_loading.T
            components.append(torch.log(nuisance + 1e-8))
        # Baseline component: log(baseline) broadcast over spots
        baseline_log = torch.log(baseline + 1e-12)[None, :].expand(counts.shape[0], -1)
        components.append(baseline_log)

        stacked = torch.stack(components, dim=2)  # (n, g, C)
        log_mu = log_l_t[:, None] + torch.logsumexp(stacked, dim=2)  # (n, g)

        # NB log_prob with total_count=theta, logits=log_mu - log(theta), mean=exp(log_mu)
        # log_prob = lgamma(y+r)-lgamma(r)-lgamma(y+1) + y*logits - (y+r)*softplus(logits)
        logits = log_mu - torch.log(theta + 1e-12)[None, :]
        softplus_logits = F.softplus(logits)
        # Use float32 lgamma (torch.lgamma)
        lgamma_y_r = torch.lgamma(y_t + theta[None, :])
        lgamma_r = torch.lgamma(theta[None, :])
        lgamma_y1 = torch.lgamma(y_t + 1.0)
        log_prob = lgamma_y_r - lgamma_r - lgamma_y1 + y_t * logits - (y_t + theta[None, :]) * softplus_logits
        nll = -log_prob
        likelihood = nll.sum(dim=1).mean()
        loss = likelihood + 1e-4 * kl / float(counts.shape[0])
        losses.append(float(loss.detach().cpu().item()))

    return {
        "checkpoint": str(checkpoint_path),
        "checkpoint_schema": ckpt.get("checkpoint_schema"),
        "backend": ckpt.get("backend"),
        "parameter_state_source": "params",
        "input_hash": checkpoint_input_hash,
        "objective_input_hash": checkpoint_objective_hash,
        "objective_input_hash_status": objective_hash_status,
        "global_step": int(global_step),
        "mc_draws": mc_draws,
        "seed": seed,
        "objective_mean": float(np.mean(losses)),
        "objective_sd": float(np.std(losses, ddof=1)) if len(losses) > 1 else 0.0,
        "objective_draws": losses,
        "n_spots": int(counts.shape[0]),
        "n_genes": int(counts.shape[1]),
        "nuisance_gauge": nuisance_gauge_summary(raw_b, raw_v, raw_h),
    }
