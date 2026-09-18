"""T1 PC 匹配修复：abs-cosine 分组 + medoid 参照符号对齐 + 定向平均。

只修 PC 语义（A02）。不对通用生物学签名取绝对值（03 §5 禁止）。
"""

from __future__ import annotations

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform


def abs_cosine_matrix(L: np.ndarray) -> np.ndarray:
    """Pairwise |cosine| — PC 符号不唯一，匹配必须符号不变。"""
    nrm = L / np.linalg.norm(L, axis=1, keepdims=True)
    return np.abs(nrm @ nrm.T)


def match_groups(L: np.ndarray, cutoff: float = 0.75) -> np.ndarray:
    """Average-linkage clustering on 1-|cos| distance."""
    sim = abs_cosine_matrix(L)
    np.fill_diagonal(sim, 1.0)
    Z = linkage(squareform(1 - np.clip(sim, 0, 1), checks=False), method="average")
    return fcluster(Z, t=1 - cutoff, criterion="distance")


def oriented_mean(L: np.ndarray, mem: list[int]) -> tuple[np.ndarray, list[int], float]:
    """Medoid-referenced sign alignment, then mean.

    Returns (mean_loading, flips, conflict) where flips[i] = ±1 applied to
    member i, conflict = fraction of members disagreeing with medoid sign
    pattern (opposing directions are NOT silently averaged away).
    """
    if len(mem) == 1:
        return L[mem[0]].copy(), [1], 0.0
    sub = np.abs((L[mem] / np.linalg.norm(L[mem], axis=1, keepdims=True))
                 @ (L[mem] / np.linalg.norm(L[mem], axis=1, keepdims=True)).T)
    np.fill_diagonal(sub, 0.0)
    medoid_pos = int(np.argmax(sub.sum(axis=1)))
    ref = L[mem[medoid_pos]]
    flips, aligned = [], []
    for i in mem:
        s = 1 if float(L[i] @ ref) >= 0 else -1
        flips.append(s)
        aligned.append(s * L[i])
    A = np.stack(aligned)
    mean = A.mean(axis=0)
    # conflict: fraction of members whose RAW direction opposed the medoid
    # reference (i.e. needed a -1 flip). v/-v pair -> 0.5 (expected,
    # resolved); all-aligned group -> 0.0. Orthogonal strangers grouped by a
    # loose cutoff also show high conflict: do not trust their mean.
    n_flip = sum(1 for s in flips if s < 0)
    conflict = float(n_flip / len(flips))
    return mean, flips, conflict
