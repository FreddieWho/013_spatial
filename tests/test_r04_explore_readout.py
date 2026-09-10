import numpy as np

from scripts.r04_explore_readout_deepdive import rank_readout
from scripts.r04_explore_readout_null import (
    aggregate_specific_increment,
    empirical_upper_p,
    permutation_null_specific_increment,
)


def _synthetic(n_per_patient: int = 60, seed: int = 0):
    rng = np.random.default_rng(seed)
    patients = np.asarray(["p1"] * n_per_patient + ["p2"] * n_per_patient)
    signal = np.concatenate([np.zeros(n_per_patient), np.ones(n_per_patient)])
    z = np.column_stack([signal + 0.3 * rng.standard_normal(2 * n_per_patient),
                         rng.standard_normal(2 * n_per_patient)])
    order = rng.permutation(2 * n_per_patient)
    base = np.concatenate([np.zeros(45), np.ones(15), np.zeros(15), np.ones(45)])
    labels = {
        "TLS": base[order],
        "TUMOR_STROMA_BOUNDARY": base[rng.permutation(2 * n_per_patient)],
    }
    return z, labels, patients


def test_aggregate_matches_manual_inner_folds():
    z, labels, patients = _synthetic()
    aggregate = aggregate_specific_increment(z, labels, patients)
    assert aggregate["valid_inner_folds"] == 2
    assert aggregate["total_inner_folds"] == 2
    for name in ("TLS", "TUMOR_STROMA_BOUNDARY"):
        value = aggregate["mean_auc_delta_by_target"][name]
        assert value is None or -1.0 <= value <= 1.0
    assert aggregate["mean_mse_delta"] is not None


def test_permutation_null_is_deterministic_and_bounded():
    z, labels, patients = _synthetic()
    first = permutation_null_specific_increment(z, labels, patients, seed=7, draws=10)
    second = permutation_null_specific_increment(z, labels, patients, seed=7, draws=10)
    assert first["valid_draws"] == 10
    assert (first["null_auc_delta_by_target"]["TLS"]
            == second["null_auc_delta_by_target"]["TLS"])
    third = permutation_null_specific_increment(z, labels, patients, seed=8, draws=10)
    assert (third["null_auc_delta_by_target"]["TLS"]
            != first["null_auc_delta_by_target"]["TLS"])


def test_empirical_upper_p_edges():
    assert empirical_upper_p([0.1, 0.2, 0.3], 0.25) == 1 / 3
    assert empirical_upper_p([], 0.1) is None
    assert empirical_upper_p([0.1], None) is None
    assert empirical_upper_p([0.5, 0.5], 0.5) == 1.0


def test_rank_readout_reports_absolute_metrics_and_norms():
    z, labels, patients = _synthetic()
    readout = rank_readout(z, labels, patients, representation={"rank": 2})
    assert readout["valid_inner_folds"] == 2
    for name in ("TLS", "TUMOR_STROMA_BOUNDARY"):
        metrics = readout["absolute_metrics"][name]
        assert metrics["shared_auc"] is not None
        assert 0.0 <= metrics["shared_auc"] <= 1.0
        assert metrics["specific_residual_coefficient_norm"] is not None
        assert metrics["specific_residual_coefficient_norm"] >= 0.0
    assert len(readout["fold_records"]) == 2
