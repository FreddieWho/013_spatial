"""T2 Lane A 统计核心：唯一 patient-first 汇总 + 同算子统计量（A05/A06）。

观察与所有 surrogate 必须调用同一函数。患者等权；有效箱记录点数。
"""

from __future__ import annotations

import numpy as np


def summarize_patient_profiles(by_patient: dict[str, list[np.ndarray]]) -> dict[str, np.ndarray]:
    """患者内先汇总：patient_profile[p,b] = median_s profile[p,s,b]。

    输入：patient -> 该患者的切片 profile 列表（每片 n_bins 向量）。
    输出：patient -> median profile。NaN 感知（空箱不算零）。
    """
    out = {}
    for p, secs in by_patient.items():
        out[p] = np.nanmedian(np.stack(secs), axis=0)
    return out


def cohort_max_statistic(patient_profiles: dict[str, np.ndarray]) -> float:
    """T_obs = max_b median_p patient_profile[p,b]（03 §3 定义）。

    诊断用峰高统计量；正式主指标优先固定对比/留出改善（T2 任务要求）。
    """
    M = np.nanmedian(np.stack(list(patient_profiles.values())), axis=0)
    return float(np.nanmax(M))


def cohort_contrast_statistic(patient_profiles: dict[str, np.ndarray],
                              inside=(-4, -3, -2), outside=(2, 3, 4)) -> float:
    """固定对比统计量：区内均值 - 区外均值（cohort median profile 上）。

    不追 argmax，避免峰位投票的自由度；正负凹陷对称可检。
    """
    M = np.nanmedian(np.stack(list(patient_profiles.values())), axis=0)
    nb = len(M)
    idx_in = [b + 4 for b in inside if -4 <= b <= 8 and b + 4 < nb]
    idx_out = [b + 4 for b in outside if -4 <= b <= 8 and b + 4 < nb]
    return float(np.nanmean(M[idx_in]) - np.nanmean(M[idx_out]))


def matched_null_p(t_obs: float, null_stats: list[float]) -> float:
    """p = (1 + #{T_null >= T_obs}) / (B+1)；B=联合重复次数。"""
    null_stats = np.asarray(null_stats, dtype=float)
    return float((np.sum(null_stats >= t_obs) + 1) / (len(null_stats) + 1))
