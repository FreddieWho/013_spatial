"""Stable JSON/TSV serialisation for R-04 control-plane artifacts."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Iterable

from .candidates import candidate_registry_hash
from .runtime import atomic_json
from .types import CandidateField, FieldFit


def write_candidate_registry(path: Path, candidates: Iterable[CandidateField]) -> str:
    values = list(candidates)
    registry_hash = candidate_registry_hash(values)
    atomic_json(path, {
        "schema": "r04.candidate_registry.v1",
        "candidate_registry_hash": registry_hash,
        "candidates": [asdict(value) for value in values],
    })
    return registry_hash


def write_field_fits(path: Path, fits: Iterable[FieldFit]) -> None:
    values = []
    for fit in fits:
        value = asdict(fit)
        value["loading"] = fit.loading.tolist()
        value["field_mean"] = fit.field_mean.tolist()
        value["field_sd"] = fit.field_sd.tolist()
        values.append(value)
    atomic_json(path, {"schema": "r04.field_fit.v1", "fields": values})


def write_frozen_model(path: Path, estimator: object) -> None:
    """Persist a model's fixed parameters for later section-level inference."""
    frozen_state = getattr(estimator, "frozen_state", None)
    if frozen_state is None or not callable(frozen_state):
        raise TypeError("estimator does not expose a frozen_state method")
    atomic_json(path, frozen_state())


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_frozen_model(path: Path) -> object:
    """Load a frozen R-04 model without rebuilding or refitting its parameters."""
    payload = read_json(path)
    if payload.get("schema") != "r04.frozen_model.v1":
        raise ValueError("frozen model schema is missing or unsupported")
    from .models import MNSFEstimator, SignedResidualGPEstimator

    model_id = payload.get("model_id")
    if model_id == MNSFEstimator.model_id:
        return MNSFEstimator.from_frozen_state(payload)
    if model_id == SignedResidualGPEstimator.model_id:
        return SignedResidualGPEstimator.from_frozen_state(payload)
    raise ValueError(f"unsupported frozen model: {model_id}")


def read_candidate_registry(path: Path) -> tuple[str, list[CandidateField]]:
    payload = read_json(path)
    if payload.get("schema") != "r04.candidate_registry.v1":
        raise ValueError("candidate registry schema is missing or unsupported")
    raw_candidates = payload.get("candidates")
    if not isinstance(raw_candidates, list):
        raise ValueError("candidate registry candidates must be a list")
    candidates = []
    for value in raw_candidates:
        if not isinstance(value, dict):
            raise ValueError("candidate registry contains a non-object candidate")
        value = dict(value)
        value["member_factor_ids"] = tuple(value.get("member_factor_ids", ()))
        candidates.append(CandidateField(**value))
    stored = str(payload.get("candidate_registry_hash", ""))
    actual = candidate_registry_hash(candidates)
    if not stored or stored != actual:
        raise ValueError("candidate registry hash is missing or invalid")
    return stored, candidates
