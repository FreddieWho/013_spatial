"""T1 回归测试：PC 符号翻转必须不改变分组（A02，能失败的测试）。

旧实现用 signed cosine (nrm @ nrm.T) 直接聚类：同一根轴 v -> -v
会被分成两组；组内 raw 平均还会互相抵消（toy 审计：mean_norm=0.0）。
修复后必须用 |cosine| + 定向平均。
"""

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

from r16.recovery import pc_matching as M


def _legacy_groups(L: np.ndarray) -> np.ndarray:
    nrm = L / np.linalg.norm(L, axis=1, keepdims=True)
    sim = nrm @ nrm.T
    np.fill_diagonal(sim, 1.0)
    Z = linkage(squareform(1 - np.clip(sim, -1, 1), checks=False), method="average")
    return fcluster(Z, t=1 - 0.75, criterion="distance")


def test_legacy_signed_matching_splits_flipped_axis():
    """旧代码的失败演示：v 与 -v 被分成两组（A02 实锤）。"""
    rng = np.random.default_rng(7)
    v = rng.normal(size=500)
    v /= np.linalg.norm(v)
    L = np.stack([v, -v,
                  rng.normal(size=500) / np.sqrt(500),
                  rng.normal(size=500) / np.sqrt(500)])
    lab = _legacy_groups(L)
    assert lab[0] != lab[1], "legacy unexpectedly sign-invariant?! re-check A02"


def test_fixed_abs_matching_merges_flipped_axis():
    """修复后：v 与 -v 必须同组，且定向平均不抵消。"""
    rng = np.random.default_rng(7)
    v = rng.normal(size=500)
    v /= np.linalg.norm(v)
    L = np.stack([v, -v,
                  rng.normal(size=500) / np.sqrt(500),
                  rng.normal(size=500) / np.sqrt(500)])
    lab = M.match_groups(L, cutoff=0.75)
    assert lab[0] == lab[1]
    mean, flips, conflict = M.oriented_mean(L, [0, 1])
    assert np.linalg.norm(mean) > 0.9, f"oriented mean collapsed: {np.linalg.norm(mean)}"
    # flips record the applied sign correction; conflict counts members whose
    # raw direction opposed the medoid (here the -v member: expected, resolved)
    assert flips == [1, -1] and conflict == 0.5


def test_oriented_mean_conflict_counts_flips():
    """conflict = 需要翻转的成员比例：全同向组 0.0，含反向成员组 >0。"""
    rng = np.random.default_rng(11)
    a = rng.normal(size=300); a /= np.linalg.norm(a)
    b = rng.normal(size=300); b /= np.linalg.norm(b)
    c = -(a + 0.01 * rng.normal(size=300))
    c /= np.linalg.norm(c)
    _, _, conflict_clean = M.oriented_mean(np.stack([a, b]), [0, 1])
    assert conflict_clean == 0.0  # 同向（点积>0），无需翻转
    mean, flips, conflict = M.oriented_mean(np.stack([a, c]), [0, 1])
    assert conflict == 0.5 and flips.count(-1) == 1
    assert np.linalg.norm(mean) > 0.9  # 定向后不抵消
