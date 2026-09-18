"""Only tests this package's toy demonstrations, not the repository pipeline."""
import numpy as np
import pytest
from audit_reproductions import run, groups, cohort_stat, matched_null_stats, monte_carlo_p

@pytest.fixture(scope='module')
def result():
    return run()

def test_scope_is_explicit(result):
    assert not result['patient_data_used'] and not result['repository_training_run']

def test_sign_flip_splits_legacy_axes(result):
    assert result['pc_sign']['legacy_groups'] == 2
    assert result['pc_sign']['sign_invariant_groups'] == 1

def test_sign_alignment_required_before_averaging(result):
    assert result['pc_sign']['legacy_mean_norm'] == pytest.approx(0)
    assert result['pc_sign']['oriented_mean_norm'] == pytest.approx(1)

def test_rotation_is_not_individual_axis_matching(result):
    assert result['pc_rotation']['best_individual_axis_cosine'] < .75
    np.testing.assert_allclose(result['pc_rotation']['canonical_correlations'], [1,1])

def test_individually_novel_can_be_wholly_known(result):
    assert result['known_axis_mixture']['max_individual_axis_correlation'] < .7
    assert result['known_axis_mixture']['joint_axis_r_squared'] == pytest.approx(1)

def test_nested_contour_stability_is_rank_geometry(result):
    c = result['contour_stability']
    np.testing.assert_allclose(c['ordered_jaccards'], [3/4, 2/3])
    np.testing.assert_allclose(c['ordered_jaccards'], c['shuffled_jaccards'])
    assert c['mean_jaccard'] == pytest.approx(17/24)

def test_null_recomputes_the_identical_statistic():
    z = np.arange(4*3*5, dtype=float).reshape(4,3,5)
    np.testing.assert_allclose(matched_null_stats(z), [cohort_stat(x) for x in z])

def test_pooled_and_cohort_nulls_are_not_interchangeable(result):
    r = result['laneA_statistic_mismatch']
    assert r['legacy_null_count'] == 30*r['correct_null_count']
    assert r['legacy_p'] > r['matched_statistic_p']

def test_patient_votes_not_section_votes(result):
    r = result['patient_weighting']
    assert r['original_section_median'] != r['repeated_section_median']
    assert r['patient_first_median'] == r['original_section_median']

def test_nearest_neighbor_map_is_not_guaranteed_permutation(result):
    assert not result['nearest_neighbor_remap']['is_bijection']

def test_empirical_p_has_plus_one_correction():
    assert monte_carlo_p(2., np.array([0.,1.])) == pytest.approx(1/3)

def test_reject_zero_axes():
    with pytest.raises(ValueError):
        groups(np.zeros((2,4)), sign_invariant=True)
