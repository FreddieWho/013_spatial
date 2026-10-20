import numpy as np
from scipy import sparse
from pathlib import Path

from field_analyzer import poisoning
from field_analyzer.pipeline import InputSpec


def synthetic_case(n=80, g=30, hidden=12, seed=0):
    rng = np.random.default_rng(seed)
    counts = sparse.csr_matrix(rng.poisson(3.0, size=(n, g)).astype(float))
    genes = np.array([f"g{i}" for i in range(g)])
    xy = np.column_stack([np.repeat(np.arange(10), 8), np.tile(np.arange(8), 10)]).astype(float)
    support = np.ones(n, dtype=bool)
    visible = np.ones(n, dtype=bool)
    visible[:hidden] = False
    witness = frozenset({f"g{i}" for i in range(20, 25)})
    spec = InputSpec(
        fields={"f1": [f"g{i}" for i in range(0, 9)], "f2": [f"g{i}" for i in range(9, 18)]},
        witness_genes=witness,
        k_neighbors=4,
        max_edge_pitch=2.0,
        min_library_genes=5,
    )
    return counts, genes, xy, support, visible, spec


def test_run_poisoning_passes_on_clean_pipeline():
    counts, genes, xy, support, visible, spec = synthetic_case()
    report = poisoning.run_poisoning(
        counts,
        genes,
        xy,
        support,
        visible,
        spec,
        package_dir=Path(poisoning.__file__).parent,
        section_id="synthetic",
    )
    assert report["passed"], report
    assert all(s["identical"] for s in report["scenarios"])
    assert report["evaluation"]["metric_changed"] is True
    assert report["evaluation"]["prediction_max_abs_delta"] == 0.0


def test_compare_inputs_detects_real_difference():
    a = {"x": np.array([1.0, 2.0]), "flag": np.array([True, False])}
    b = {"x": np.array([1.0, 2.0001]), "flag": np.array([True, False])}
    res = poisoning.compare_inputs(a, b, "s", "d")
    assert res.identical is False
    assert res.per_array_delta["x"] > 0

    c = {"x": np.array([1.0, np.nan]), "flag": np.array([True, False])}
    d = {"x": np.array([1.0, 5.0]), "flag": np.array([True, False])}
    assert poisoning.compare_inputs(c, d, "s", "d").identical is False


def test_leaky_slide_statistic_would_be_caught():
    """Negative control: pre-mask whole-slide statistics differ across poisonings."""
    counts, genes, xy, support, visible, spec = synthetic_case()
    hidden = np.flatnonzero(~visible)
    base_stats = np.asarray(counts.todense()).mean(axis=0)  # the bug: before masking
    poisoned = poisoning.poison_hidden_extreme(counts, hidden)
    poisoned_stats = np.asarray(poisoned.todense()).mean(axis=0)
    res = poisoning.compare_inputs(
        {"slide_mean": base_stats}, {"slide_mean": poisoned_stats}, "leak", "leaky control"
    )
    assert res.identical is False, "a pre-mask statistic must be detected as a leak"

    honest_base = np.asarray(counts[visible].todense()).mean(axis=0)
    honest_poisoned = np.asarray(poisoned[visible].todense()).mean(axis=0)
    res2 = poisoning.compare_inputs(
        {"slide_mean": honest_base}, {"slide_mean": honest_poisoned}, "ok", "masked control"
    )
    assert res2.identical is True


def test_static_source_check_is_clean():
    report = poisoning.static_source_check(Path(poisoning.__file__).parent)
    assert report["clean"], report


def test_all_missing_hidden_truth_cannot_satisfy_the_reverse_assertion():
    """N1/F1: exercise the *production* decision functions, not a copy of the rule."""
    all_nan = np.full(5, np.nan)
    summary = poisoning.truth_delta(all_nan, all_nan)
    assert summary["truth_max_abs_delta"] is None
    assert summary["n_finite_hidden_truth"] == 0
    assert poisoning.metric_change({"f": summary}) is False

    unchanged = poisoning.truth_delta(np.array([1.0, 2.0]), np.array([1.0, 2.0]))
    assert unchanged["truth_max_abs_delta"] == 0.0
    assert poisoning.metric_change({"f": unchanged}) is False

    moved = poisoning.truth_delta(np.array([1.0, 2.0]), np.array([1.0, 2.5]))
    assert poisoning.metric_change({"f": moved}) is True

    # an inf sentinel from the old code path must not count as a change
    assert poisoning.metric_change({"f": {"truth_max_abs_delta": float("inf")}}) is False
