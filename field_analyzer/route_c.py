"""Route C (minimal): cost-based routing between visible endpoints.

Scope of this file is deliberately narrow: it implements the **C0-style reference**
(a path algorithm on one support field) and measures what such a path actually means
on S3:
  * does it follow the true bridge when the bridge is present?
  * does it still return a path when the bridge is broken (i.e. can it refuse)?
  * does it connect the wrong endpoint?

A full C1 (two-field common support with a refusal rule) is not implemented here.
The point of this node is to show what "a path exists" is worth before building it.
"""

from __future__ import annotations

import heapq
from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass
class PathResult:
    path: np.ndarray  # (m, 2) coordinates in grid units
    cost: float
    nodes: int

    def to_dict(self) -> dict[str, Any]:
        return {"cost": float(self.cost), "nodes": int(self.nodes)}


def dijkstra(
    n: int,
    node_cost: np.ndarray,
    start: int,
    goal: int,
) -> PathResult:
    """4-neighbour least-cost path on an n×n lattice with per-node cost."""
    node_cost = np.asarray(node_cost, dtype=np.float64)
    if node_cost.shape != (n * n,):
        raise ValueError("node_cost must have length n*n")
    dist = np.full(n * n, np.inf)
    prev = np.full(n * n, -1, dtype=np.int64)
    dist[start] = 0.0
    heap = [(0.0, int(start))]
    while heap:
        d, u = heapq.heappop(heap)
        if d > dist[u]:
            continue
        if u == goal:
            break
        i, j = divmod(u, n)
        for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ni, nj = i + di, j + dj
            if not (0 <= ni < n and 0 <= nj < n):
                continue
            v = ni * n + nj
            step = 0.5 * (node_cost[u] + node_cost[v])
            nd = d + step
            if nd < dist[v]:
                dist[v] = nd
                prev[v] = u
                heapq.heappush(heap, (nd, v))
    if not np.isfinite(dist[goal]):
        return PathResult(path=np.empty((0, 2)), cost=float("inf"), nodes=0)
    nodes = []
    cur = int(goal)
    while cur != -1:
        nodes.append(cur)
        if cur == start:
            break
        cur = int(prev[cur])
    nodes.reverse()
    coords = np.array([[idx % n, idx // n] for idx in nodes], dtype=float)
    return PathResult(path=coords, cost=float(dist[goal]), nodes=len(nodes))


def support_cost(field: np.ndarray, *, eps_quantile: float = 0.5) -> np.ndarray:
    """Cost = 1 / (eps + normalised support); high support is cheap."""
    f = np.asarray(field, dtype=np.float64)
    f = np.where(np.isfinite(f), f, 0.0)
    positive = f[f > 0]
    eps = float(np.quantile(positive, 0.1)) if len(positive) else 1e-6
    scale = float(np.quantile(positive, 0.9)) if len(positive) else 1.0
    u = np.clip(f / max(scale, 1e-9), 0.0, 1.0)
    return 1.0 / (max(eps, 1e-6) + u)


def uniform_cost(n: int) -> np.ndarray:
    return np.ones(n * n, dtype=np.float64)


def path_metrics(
    path: np.ndarray,
    n: int,
    truth_bridge: np.ndarray,
    *,
    support_field: np.ndarray | None = None,
    support_threshold: float = 0.3,
    gap_center: float = 32.0,
    gap_half_width: float = 4.0,
) -> dict[str, Any]:
    """What does the returned path actually traverse?"""
    if len(path) == 0:
        return {"path_length": 0.0, "refused": True}
    idx = (path[:, 1].astype(int) * n + path[:, 0].astype(int))
    values = truth_bridge[idx]
    on_support = values >= support_threshold
    own = np.asarray(support_field, dtype=np.float64)[idx] if support_field is not None else None
    xs = path[:, 0]
    in_gap = np.abs(xs - gap_center) <= gap_half_width
    length = float(np.sum(np.linalg.norm(np.diff(path, axis=0), axis=1))) if len(path) > 1 else 0.0
    out = {
        "path_length": length,
        "refused": False,
        "mean_truth_support_along_path": float(np.mean(values)),
        "fraction_on_true_bridge": float(np.mean(on_support)),
        "fraction_in_gap_band": float(np.mean(in_gap)),
        "max_truth_along_path": float(np.max(values)),
        "min_truth_along_path": float(np.min(values)),
    }
    if own is not None:
        out["variant_support_mean_along_path"] = float(np.mean(own))
        out["variant_support_min_along_path"] = float(np.min(own))
        out["variant_support_below_0.05_fraction"] = float(np.mean(own < 0.05))
    return out
