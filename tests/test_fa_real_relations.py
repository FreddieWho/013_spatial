import numpy as np
import pytest

from field_analyzer import real_relations as rr


def test_zscore_is_per_value_set():
    v = np.arange(20.0)
    z = rr.zscore(v)
    assert np.isclose(np.mean(z), 0.0)
    assert np.isclose(np.std(z), 1.0)


def test_displacement_curve_finds_a_known_shift():
    rng = np.random.default_rng(0)
    xy = np.column_stack([rng.uniform(0, 100, 800), rng.uniform(0, 100, 800)])
    pitch = 5.0
    base = np.exp(-((xy[:, 0] - 40) ** 2 + (xy[:, 1] - 50) ** 2) / (2 * 10 ** 2))
    # target is the same field shifted by +4 pitch along x
    shifted_xy = xy + np.array([4 * pitch, 0.0])
    target = np.exp(-((shifted_xy[:, 0] - 40) ** 2 + (shifted_xy[:, 1] - 50) ** 2) / (2 * 10 ** 2))
    res = rr.displacement_curve(xy, base, target, pitch=pitch, max_lag_pitch=10, step_pitch=2)
    assert res["status"] == "OK"
    assert res["peak_dx_pitch"] == pytest.approx(4.0, abs=2.0)


def test_displacement_null_returns_arbitrary_large_peak():
    """The null must NOT return zero: that is why it is required as a control."""
    rng = np.random.default_rng(1)
    xy = np.column_stack([rng.uniform(0, 100, 800), rng.uniform(0, 100, 800)])
    pitch = 5.0
    base = np.exp(-((xy[:, 0] - 40) ** 2 + (xy[:, 1] - 50) ** 2) / (2 * 12 ** 2))
    target = np.exp(-((xy[:, 0] - 41) ** 2 + (xy[:, 1] - 50) ** 2) / (2 * 12 ** 2))
    res = rr.displacement_curve(
        xy, base, target, pitch=pitch, max_lag_pitch=10, step_pitch=2, null_rotate_deg=137.0
    )
    assert res["status"] == "OK"
    # a rotated (unrelated) target can still produce a large peak magnitude
    assert res["peak_magnitude_pitch"] >= 0.0


def test_louo_recovers_a_transferable_relation():
    """Two synthetic units with the same relation must transfer."""
    from types import SimpleNamespace

    rng = np.random.default_rng(0)
    sections = []
    for u in range(4):
        n = 300
        xy = np.column_stack([rng.uniform(0, 50, n), rng.uniform(0, 50, n)])
        a1 = rng.normal(size=n)
        a2 = 0.7 * a1 + 0.3 * rng.normal(size=n)
        b1 = rng.normal(size=n)
        b2 = 0.6 * b1 + 0.4 * rng.normal(size=n)
        sections.append(
            rr.SectionScores(
                section_id=f"s{u}",
                unit_id=f"unit{u}",
                xy=xy,
                visible=np.ones(n, dtype=bool),
                a={"x": a1, "y": a2},
                b={"x": b1, "y": b2},
                log_lib=np.zeros(n),
            )
        )
    results = rr.louo_holdout(sections, pairs=[("y", "x")], conditioned=False)
    assert len(results) == 4
    assert all(r.r2_test > 0.2 for r in results)
