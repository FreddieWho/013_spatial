import numpy as np
import pytest

from field_analyzer import synthetic as syn


def test_seed_roles_are_disjoint_and_complete():
    seen = []
    for role, seeds in syn.SEED_ROLES.items():
        assert len(set(seeds)) == len(seeds), role
        seen.extend(seeds)
    assert len(seen) == len(set(seen))
    assert syn.role_of_seed(0) == "train"
    assert syn.role_of_seed(9) == "development"
    assert syn.role_of_seed(20) == "evaluation"
    with pytest.raises(ValueError):
        syn.role_of_seed(99)


def test_s1_pair_has_matched_radial_profile():
    xy = syn.grid_xy()
    s = syn.build_s1(0)
    a = s.fields["circular"]
    b = s.fields["sided"]
    pa = syn._radial_profile(a, xy)
    pb = syn._radial_profile(b, xy)
    assert np.max(np.abs(pa - pb)) < 0.01, "the S1 control should share the radial mean"
    assert np.max(np.abs(a - b)) > 0.1, "but the angular shape must differ"


def test_s1_bimodal_keeps_centroid_and_mass():
    xy = syn.grid_xy()
    s = syn.build_s1(1)
    u, v = s.fields["unimodal"], s.fields["bimodal"]
    cu = np.sum(xy * u[:, None], axis=0) / u.sum()
    cv = np.sum(xy * v[:, None], axis=0) / v.sum()
    assert np.linalg.norm(cu - cv) < 0.05
    assert abs(u.sum() - v.sum()) / u.sum() < 0.02


def test_s0_fields_are_structureless_relative_to_each_other():
    s = syn.build_s0(0)
    a, b = s.fields["f0"], s.fields["f1"]
    assert abs(np.corrcoef(a, b)[0, 1]) < 0.6


def test_s3_variants_differ_as_designed():
    xy = syn.grid_xy()
    s = syn.build_s3(0)
    cont = s.fields["bridge_continuous"]
    broken = s.fields["bridge_broken"]
    mid = np.abs(xy[:, 0] - 32.0) < 3.5
    assert cont[mid].max() > 0.5
    assert broken[mid].max() == pytest.approx(0.0, abs=1e-9)
    assert "endpoint_wrong" in s.fields
    wrong = s.fields["endpoint_wrong"]
    left = s.fields["endpoint_left"]
    assert wrong.argmax() != left.argmax()


def test_s4_keeps_true_bridge_while_adding_fields():
    base = syn.build_s4(0, n_distractors=0)
    many = syn.build_s4(0, n_distractors=16)
    assert len(many.fields) == len(base.fields) + 16
    assert np.array_equal(base.fields["bridge_continuous"], many.fields["bridge_continuous"])


def test_s5_visible_values_are_bit_identical():
    a, b = syn.make_s5_pair(0)
    va, vb = a.counts.copy(), b.counts.copy()
    mask = a.visible_mask
    assert np.array_equal(va[mask], vb[mask])
    assert not np.array_equal(a.fields["bridge"][~mask], b.fields["bridge"][~mask])


def test_counts_are_nonnegative_integers():
    s = syn.build_s1(2)
    counts = syn.counts_from_fields(s, seed=2)
    assert counts.min() >= 0
    assert np.allclose(counts, np.round(counts))
    assert counts.shape == (syn.GRID * syn.GRID, len(s.gene_names()))


def test_materialise_separates_truth_from_inputs(tmp_path):
    a, b = syn.make_s5_pair(0)  # S5 has a real hidden disc
    base = syn.materialise(a, tmp_path)
    assert (base / "observed_inputs" / "counts_visible.npz").exists()
    assert (base / "latent_truth" / "fields.npz").exists()
    assert (base / "evaluation_only" / "answers.npz").exists()
    obs = np.load(base / "observed_inputs" / "counts_visible.npz")
    truth = np.load(base / "latent_truth" / "fields.npz")
    mask = obs["visible_mask"]
    assert np.all(np.isnan(obs["counts"][~mask]))
    assert truth["bridge"][~mask].max() > 0.5


def test_g2_mechanism_is_not_the_geometric_generator():
    f = syn.diffuse_field(conductance_x=1.0, conductance_y=0.2)
    grid = f.reshape(syn.GRID, syn.GRID)
    axis = np.arange(syn.GRID) - syn.GRID // 2
    mass_x = grid.sum(axis=0)
    mass_y = grid.sum(axis=1)
    spread_x = float(np.sum(mass_x * axis ** 2) / mass_x.sum())
    spread_y = float(np.sum(mass_y * axis ** 2) / mass_y.sum())
    # anisotropy comes from conductance, not from a gaussian template
    assert spread_x > 1.15 * spread_y
    iso = syn.diffuse_field(conductance_x=1.0, conductance_y=1.0).reshape(syn.GRID, syn.GRID)
    ix = float(np.sum(iso.sum(axis=0) * axis ** 2) / iso.sum())
    iy = float(np.sum(iso.sum(axis=1) * axis ** 2) / iso.sum())
    assert abs(ix - iy) / max(ix, iy) < 0.05
