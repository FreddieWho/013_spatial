import numpy as np
from scipy import sparse

from field_analyzer import geometry
from field_analyzer.real_relations import SectionScores
from field_analyzer.spatial_lag import compare_pairs, holdout_mask, misalign_values, summarise


def _section(section_id, unit_id, xy, predictor, target):
    n = len(xy)
    return SectionScores(
        section_id=section_id,
        unit_id=unit_id,
        xy=xy,
        visible=np.ones(n, dtype=bool),
        a={"p": predictor, "t": np.zeros(n)},
        b={"p": np.zeros(n), "t": target},
        log_lib=np.zeros(n),
    )


def test_misalign_keeps_missing_and_changes_alignment():
    xy = np.column_stack([np.linspace(0, 10, 30), np.zeros(30)])
    values = np.linspace(0, 1, 30)
    values[0] = np.nan
    moved = misalign_values(xy, values, 137.0)
    assert np.isnan(moved[0])
    assert np.isfinite(moved[1:]).all()
    assert not np.allclose(moved[1:], values[1:])


def test_holdout_mask_is_stable_and_empty_at_zero_fraction():
    a = holdout_mask(100, 0.1, 0, "ZEN36")
    b = holdout_mask(100, 0.1, 0, "ZEN36")
    c = holdout_mask(100, 0.1, 0, "ZEN40")
    assert np.array_equal(a, b)
    assert not np.array_equal(a, c)
    assert holdout_mask(20, 0.0, 0, "ZEN36").sum() == 0


def test_masked_neighbour_ignores_hidden_values():
    xy = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0]])
    graph = geometry.build_adjacency(xy, np.ones(3, dtype=bool), k=1, max_edge_pitch=None).graph
    values = np.array([0.0, 999.0, 0.0])
    observed = np.array([True, False, True])
    feat = geometry.neighbor_sum(graph, values, observed)
    assert feat[1] == 0.0
    assert 999.0 not in feat


def test_cross_neighbour_beats_same_neighbour_when_target_is_the_lag():
    rng = np.random.default_rng(0)
    sections = []
    for u in range(4):
        n = 180
        xy = rng.uniform(0, 40, size=(n, 2))
        predictor = rng.normal(size=n)
        graph = geometry.build_adjacency(xy, np.ones(n, dtype=bool), k=6).graph
        target = geometry.neighbor_sum(graph, predictor, np.ones(n, dtype=bool))
        target = target + 0.05 * rng.normal(size=n)
        sections.append(_section(f"s{u}", f"u{u}", xy, predictor, target))
    rows = compare_pairs(
        sections, pairs=[("t", "p")], conditioned=False, holdout_fraction=0.0,
        protocol="full", angles=(137.0, 90.0),
    )
    assert len(rows) == 4
    assert all(r["r2_cross_neighbour"] > r["r2_same_neighbour"] for r in rows)
    assert all(r["r2_cross_neighbour"] > r["r2_rot_137"] for r in rows)
    assert all(r["passes"] for r in rows)


def test_independent_predictor_does_not_upgrade():
    rng = np.random.default_rng(1)
    sections = []
    for u in range(4):
        n = 160
        xy = rng.uniform(0, 30, size=(n, 2))
        predictor = rng.normal(size=n)
        target = np.sin(xy[:, 0] / 3.0)
        sections.append(_section(f"s{u}", f"u{u}", xy, predictor, target))
    rows = compare_pairs(
        sections, pairs=[("t", "p")], conditioned=False, holdout_fraction=0.0, protocol="full",
    )
    summary = summarise(rows + [{**r, "protocol": "masked"} for r in rows])
    assert summary["upgrade_af01"] is False
    assert all(r["r2_cross_neighbour"] < r["r2_same_neighbour"] for r in rows)
