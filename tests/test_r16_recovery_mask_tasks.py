"""T5 遮蔽任务防泄漏回归测试（03 §7）。

硬约束（fail-closed）：
1. 遮蔽窗生成与 GT 无关（随机/规则窗口，事后套 GT 定阴阳）。
2. 隐藏窗内表达及一切派生（score/mask/图特征/归一化统计）不可见；
   可见点空间图不得用隐藏分子值。
3. 隐藏值改成任意数不改变输入特征与预测（泄漏探针）。
4. 阴性窗与阳性窗同样生成；部分覆盖单列。
"""

import numpy as np

from r16.recovery import mask_tasks as T


def _toy(n=400, seed=0):
    rng = np.random.default_rng(seed)
    coords = rng.uniform(0, 10, size=(n, 2))
    expr = rng.normal(size=n)
    return coords, expr


def test_windows_independent_of_gt():
    coords, _ = _toy()
    gt = np.zeros(len(coords), dtype=int)
    gt[:20] = 1  # 假 GT 簇
    wins = T.generate_windows(coords, n_windows=10, window_size=2.0, seed=1)
    # 窗生成不拿 gt 作输入（签名检查：函数不接受 gt 参数）
    import inspect
    assert "gt" not in inspect.signature(T.generate_windows).parameters
    pos, neg = T.label_windows(wins, coords, gt)
    # 部分覆盖窗单列（既非阳性亦非阴性），故 pos+neg<=windows
    assert len(pos) + len(neg) <= 10 and len(neg) > 0


def test_hidden_values_do_not_leak_into_features():
    coords, expr = _toy()
    wins = T.generate_windows(coords, n_windows=4, window_size=2.0, seed=2)
    feat1, _ = T.visible_features(coords, expr, wins[0])
    # 隐藏值改成任意数，输入特征与预测必须不变
    expr2 = expr.copy()
    hid = T.window_mask(coords, wins[0])
    expr2[hid] = 999.0
    feat2, _ = T.visible_features(coords, expr2, wins[0])
    assert np.array_equal(feat1, feat2, equal_nan=True)


def test_negative_windows_scored():
    coords, expr = _toy()
    gt = np.zeros(len(coords), dtype=int)
    gt[:10] = 1
    wins = T.generate_windows(coords, n_windows=8, window_size=2.0, seed=3)
    pos, neg = T.label_windows(wins, coords, gt)
    assert len(neg) > 0, "必须有阴性窗参与评分"
    res = T.score_windows(coords, expr, wins, gt)
    assert "n_pos" in res and "n_neg" in res
