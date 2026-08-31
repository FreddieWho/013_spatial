"""Deterministic execution, resource and checkpoint helpers."""

from __future__ import annotations

import hashlib
import json
import os
import random
import tempfile
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np


DEFAULT_STORAGE_RESERVE_BYTES = 1_200_000_000_000
STORAGE_RESERVE_ENV = "R04_STORAGE_RESERVE_BYTES"


def derived_seed(master_seed: int, *parts: object) -> int:
    payload = "|".join([str(master_seed), *(str(part) for part in parts)]).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "big") % (2**31 - 1)


def seed_everything(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    try:
        import tensorflow as tf
        tf.keras.utils.set_random_seed(seed)
        tf.config.experimental.enable_op_determinism()
    except ImportError:
        pass


def config_fingerprint(value: object, *, exclude: Iterable[str] = ()) -> str:
    """Hash serialisable model settings while excluding runtime-only fields."""
    payload = asdict(value) if is_dataclass(value) else dict(value) if isinstance(value, dict) else value
    if isinstance(payload, dict):
        excluded = set(exclude)
        payload = {key: item for key, item in payload.items() if key not in excluded}
    encoded = json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def atomic_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        handle.write(payload)
        temporary = Path(handle.name)
    temporary.replace(path)


def available_bytes(path: Path = Path(".")) -> int:
    return int(os.statvfs(path).f_bavail * os.statvfs(path).f_frsize)


def configured_storage_reserve_bytes() -> int:
    """Return the default reserve, with an explicit execution-only override."""
    raw = os.environ.get(STORAGE_RESERVE_ENV)
    if raw is None:
        return DEFAULT_STORAGE_RESERVE_BYTES
    try:
        value = int(raw)
    except ValueError as exc:
        raise ValueError(f"{STORAGE_RESERVE_ENV} must be an integer") from exc
    if value < 0:
        raise ValueError(f"{STORAGE_RESERVE_ENV} must be non-negative")
    return value


def resource_status(
    path: Path = Path("."), reserve_bytes: int | None = None
) -> str:
    effective_reserve = (
        configured_storage_reserve_bytes()
        if reserve_bytes is None
        else reserve_bytes
    )
    if effective_reserve < 0:
        raise ValueError("reserve_bytes must be non-negative")
    return "BLOCKED_STORAGE" if available_bytes(path) < effective_reserve else "OK"


def checkpoint_is_compatible(path: Path, *, input_hash: str, config_hash: str, environment_hash: str) -> bool:
    if not path.exists():
        return False
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return value.get("input_hash") == input_hash and value.get("config_hash") == config_hash and value.get("environment_hash") == environment_hash
