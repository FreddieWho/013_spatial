"""Coordinate geometry used by R-04 models and nulls."""

from __future__ import annotations

import numpy as np


def median_nearest_neighbor_distance(coords: np.ndarray) -> float:
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or coords.shape[1] != 2 or len(coords) < 2:
        raise ValueError("at least two 2-D coordinates are required")
    delta = coords[:, None, :] - coords[None, :, :]
    distance = np.sqrt(np.sum(delta * delta, axis=2))
    np.fill_diagonal(distance, np.inf)
    nearest = distance.min(axis=1)
    value = float(np.median(nearest))
    if not np.isfinite(value) or value <= 0:
        raise ValueError("coordinates have no positive nearest-neighbor scale")
    return value


def normalize_coordinates(coords: np.ndarray) -> tuple[np.ndarray, float]:
    coords = np.asarray(coords, dtype=float)
    scale = median_nearest_neighbor_distance(coords)
    center = np.median(coords, axis=0)
    return (coords - center) / scale, scale


def matern32_kernel(coords_a: np.ndarray, coords_b: np.ndarray, lengthscale: float) -> np.ndarray:
    """Matérn-3/2 correlation kernel over 2-D coordinates."""
    if lengthscale <= 0:
        raise ValueError("lengthscale must be positive")
    a = np.asarray(coords_a, dtype=float)
    b = np.asarray(coords_b, dtype=float)
    distance = np.sqrt(np.maximum(0.0, ((a[:, None, :] - b[None, :, :]) ** 2).sum(axis=2)))
    root3 = np.sqrt(3.0) * distance / lengthscale
    return (1.0 + root3) * np.exp(-root3)


def spatial_lengthscale_is_interior(lengthscale: float, lower: float = 0.75, upper: float = 8.0) -> bool:
    return lower < float(lengthscale) < upper
