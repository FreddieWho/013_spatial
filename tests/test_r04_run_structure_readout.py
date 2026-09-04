import numpy as np

from scripts.r04_run_structure_readout import _aggregate_inner_readout, _inner_patient_folds


def test_inner_patient_folds_are_disjoint_and_cover_all_patients():
    patients = np.asarray(["p3", "p1", "p2", "p1", "p3", "p2", "p4"])
    folds = _inner_patient_folds(patients)
    assert len(folds) == 3
    seen = set()
    for train, test, test_patients in folds:
        assert not np.any(train & test)
        assert np.all(train | test)
        assert set(patients[test].tolist()) == set(test_patients)
        seen.update(test_patients)
    assert seen == {"p1", "p2", "p3", "p4"}


def test_inner_readout_aggregation_keeps_invalid_fold_audit():
    computed = {
        "inner_fold": 0,
        "status": "COMPUTED",
        "readout": {
            "shared_only": {"mean_mse": 2.0},
            "shared_plus_specific_residual": {"mean_mse": 1.0},
            "specific_increment": {
                "mean_mse_delta_specific_minus_shared": -1.0,
                "auc_delta_by_target": {"TLS": 0.1, "TUMOR_STROMA_BOUNDARY": None},
            },
        },
    }
    invalid = {"inner_fold": 1, "status": "NOT_TESTABLE_NO_TARGET_VARIATION"}
    result = _aggregate_inner_readout([computed, invalid])
    assert result["valid_inner_folds"] == 1
    assert result["total_inner_folds"] == 2
    assert result["specific_minus_shared_mean_mse_delta"] == -1.0
    assert len(result["fold_results"]) == 2
    assert result["fold_results"][1]["status"] == "NOT_TESTABLE_NO_TARGET_VARIATION"
