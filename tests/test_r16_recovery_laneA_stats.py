"""T2 回归测试：患者先汇总 + 观察/null 同算子 + 患者权重（A05/A06）。

旧代码的失败演示必须保留为测试：cohort-median-peak vs
section-draw pooled peaks 是两个不同的统计量；section 投票给多片患者多票。
"""

import numpy as np

from r16.recovery import laneA_stats as S


def _toy_profiles():
    rng = np.random.default_rng(3)
    # 2 patients x 2 sections x 5 bins; patient A has a real bump at bin 2
    profs = {
        ("pA", "s1"): np.array([0, 0, 2.0, 0, 0]) + rng.normal(0, 0.1, 5),
        ("pA", "s2"): np.array([0, 0, 1.8, 0, 0]) + rng.normal(0, 0.1, 5),
        ("pB", "s1"): rng.normal(0, 0.1, 5),
        ("pB", "s2"): rng.normal(0, 0.1, 5),
    }
    by_pat = {"pA": [profs[("pA", "s1")], profs[("pA", "s2")]],
              "pB": [profs[("pB", "s1")], profs[("pB", "s2")]]}
    return by_pat


def test_patient_first_then_cohort_statistic():
    by_pat = _toy_profiles()
    pp = S.summarize_patient_profiles(by_pat)
    assert set(pp) == {"pA", "pB"}
    # patient A keeps its bump after within-patient median
    assert int(np.argmax(pp["pA"])) == 2
    t = S.cohort_max_statistic(pp)
    assert t > 0.8, f"cohort stat should see the bump, got {t}"


def test_section_votes_overweight_multisection_patient():
    """多片患者在 section 投票下票多；patient-first 下每人一票。"""
    by_pat = {"pA": [np.array([0, 0, 5.0, 0, 0])] * 4,
              "pB": [np.array([5.0, 0, 0, 0, 0])],
              "pC": [np.array([5.0, 0, 0, 0, 0])]}
    # legacy section-vote agree: bin2 gets 4/6 votes
    votes = []
    for p, ss in by_pat.items():
        for v in ss:
            votes.append(int(np.argmax(v)))
    legacy_agree = max(np.unique(votes, return_counts=True)[1]) / len(votes)
    assert legacy_agree > 0.5  # 旧口径下 pA 一人说了算
    # patient-first: 3 patients, modal bin0 has 2/3
    pp = S.summarize_patient_profiles(by_pat)
    peaks = {p: int(np.argmax(v)) for p, v in pp.items()}
    assert peaks == {"pA": 2, "pB": 0, "pC": 0}


def test_null_uses_identical_operator():
    """null 的每次 draw 必须走同一 patient-first + max 算子（A05）。"""
    by_pat = _toy_profiles()
    t_obs = S.cohort_max_statistic(S.summarize_patient_profiles(by_pat))
    rng = np.random.default_rng(0)
    # under null (pure noise), matched-operator p should usually be large
    draws = [S.cohort_max_statistic(
        S.summarize_patient_profiles(
            {p: [rng.normal(0, 0.1, 5) for _ in ss] for p, ss in by_pat.items()}))
        for _ in range(50)]
    p = (sum(d >= t_obs for d in draws) + 1) / (len(draws) + 1)
    assert p <= 0.1, f"bump should beat noise draws, p={p}"
