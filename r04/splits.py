"""Patient/block grouped and buffered spatial split helpers."""

from __future__ import annotations

import hashlib
from typing import Iterable

import numpy as np


def _stable_group_order(groups: Iterable[str], seed: int) -> list[str]:
    unique = sorted({str(group) for group in groups})
    return sorted(unique, key=lambda group: hashlib.sha256(f"{seed}|{group}".encode()).hexdigest())


def grouped_fold_ids(groups: Iterable[str], n_folds: int = 5, seed: int = 20260807) -> np.ndarray:
    values = np.asarray([str(group) for group in groups])
    unique = _stable_group_order(values, seed)
    if len(unique) < 2:
        raise ValueError("at least two biological groups are required")
    fold = {group: index % min(n_folds, len(unique)) for index, group in enumerate(unique)}
    return np.asarray([fold[group] for group in values], dtype=int)


def assert_group_disjoint(train_groups: Iterable[str], test_groups: Iterable[str]) -> None:
    overlap = set(map(str, train_groups)) & set(map(str, test_groups))
    if overlap:
        raise ValueError(f"biological groups cross split boundary: {sorted(overlap)}")


def buffered_spatial_train_mask(coords: np.ndarray, holdout: np.ndarray, buffer: float) -> np.ndarray:
    """Exclude all training spots within ``buffer`` of held-out spots."""
    coordinates = np.asarray(coords, dtype=float)
    holdout = np.asarray(holdout, dtype=bool)
    if len(coordinates) != len(holdout) or buffer < 0:
        raise ValueError("invalid buffered split inputs")
    if not holdout.any() or holdout.all():
        raise ValueError("holdout must contain both train and test spots")
    distance = np.sqrt(((coordinates[:, None, :] - coordinates[holdout][None, :, :]) ** 2).sum(axis=2))
    return (distance.min(axis=1) > buffer) & ~holdout
