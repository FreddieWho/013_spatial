"""T5 遮蔽任务公共模块（03 §7 防泄漏合同）。

窗生成与 GT 无关；隐藏派生不可见；泄漏探针可测。
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree


def generate_windows(coords: np.ndarray, n_windows: int, window_size: float,
                     seed: int) -> list[dict]:
    """与答案无关的规则窗口：网格中心 + 固定半窗。GT 不得作为输入。"""
    rng = np.random.default_rng(seed)
    mins = coords.min(axis=0)
    span = coords.max(axis=0) - mins
    n_side = int(np.ceil(np.sqrt(n_windows)))
    xs = np.linspace(mins[0], mins[0] + span[0], n_side + 2)[1:-1]
    ys = np.linspace(mins[1], mins[1] + span[1], n_side + 2)[1:-1]
    grid = [(x, y) for x in xs for y in ys]
    pick = rng.choice(len(grid), size=min(n_windows, len(grid)), replace=False)
    return [{"center": grid[i], "half": window_size / 2} for i in pick]


def window_mask(coords: np.ndarray, win: dict) -> np.ndarray:
    cx, cy = win["center"]
    h = win["half"]
    return ((np.abs(coords[:, 0] - cx) <= h) & (np.abs(coords[:, 1] - cy) <= h))


def label_windows(wins: list[dict], coords: np.ndarray, gt: np.ndarray,
                  pos_frac: float = 0.1) -> tuple[list[int], list[int]]:
    """事后套 GT：窗内 GT 覆盖>=pos_frac 为阳性；零覆盖为阴性；部分覆盖单列."""
    pos, neg, partial = [], [], []
    for i, w in enumerate(wins):
        m = window_mask(coords, w)
        if m.sum() == 0:
            continue
        f = gt[m].mean()
        if f >= pos_frac:
            pos.append(i)
        elif f == 0:
            neg.append(i)
        else:
            partial.append(i)
    return pos, neg


def visible_features(coords: np.ndarray, expr: np.ndarray, win: dict,
                     k: int = 8) -> tuple[np.ndarray, np.ndarray]:
    """仅可见点的 kNN 均值特征；隐藏点不参与任何消息。返回 (feat, visible)."""
    vis = ~window_mask(coords, win)
    tree = cKDTree(coords[vis])
    _, idx = tree.query(coords, k=min(k, vis.sum()))
    # 全片特征：隐藏点的邻域只能看可见点（03 §7：可见点图不得用隐藏分子值；
    # 反向——隐藏点预测时用可见邻居——是合法的场预测输入）。
    feat = expr[vis][idx].mean(axis=1)
    return feat, vis


def score_windows(coords: np.ndarray, expr: np.ndarray, wins: list[dict],
                  gt: np.ndarray) -> dict:
    pos, neg = label_windows(wins, coords, gt)
    return {"n_pos": len(pos), "n_neg": len(neg), "n_windows": len(wins)}
