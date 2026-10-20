import numpy as np
import pytest
from scipy import sparse

from field_analyzer import observation as obs
from field_analyzer.pipeline import InputSpec, build_inputs


def make_roles(genes, witness):
    roles = {}
    for g in genes:
        roles[g] = obs.ROLE_B_WITNESS if g in witness else obs.ROLE_BACKGROUND
    return obs.GeneRoles.from_mapping(roles)


def test_hidden_must_be_nan():
    genes = np.array(["g1", "g2"])
    roles = make_roles(genes, set())
    values = np.array([[1.0, 2.0], [3.0, 4.0]])
    mask = np.array([True, False])
    visible = values.copy()
    visible[1] = np.nan
    bundle = obs.ObservationBundle(
        sample_id="s",
        unit_id="u",
        source_id="src",
        patient_status="known_within_study",
        barcodes=np.array(["b1", "b2"]),
        source_rows=np.array([0, 1]),
        xy=np.array([[0.0, 0.0], [1.0, 0.0]]),
        coordinate_unit="synthetic_pitch",
        scale_provenance={"source": "test"},
        tissue_support=np.array([True, True]),
        visible_mask=mask,
        genes=genes,
        gene_roles=roles,
        values_visible=visible,
    )
    assert bundle.n_visible == 1

    leaked = values.copy()
    with pytest.raises(obs.MaskViolation):
        obs.ObservationBundle(**{**bundle.__dict__, "values_visible": leaked, "visible_mask": mask})


def test_witness_gene_on_axis_is_rejected():
    genes = np.array(["g1", "witness"])
    roles = make_roles(genes, {"witness"})
    with pytest.raises(obs.MaskViolation):
        obs.ObservationBundle(
            sample_id="s",
            unit_id="u",
            source_id="src",
            patient_status="unknown",
            barcodes=np.array(["b1"]),
            source_rows=np.array([0]),
            xy=np.array([[0.0, 0.0]]),
            coordinate_unit="synthetic_pitch",
            scale_provenance={},
            tissue_support=np.array([True]),
            visible_mask=np.array([True]),
            genes=np.array(["g1"]),
            gene_roles=roles,
            values_visible=np.array([[1.0]]),
        ).gene_roles.assert_no_witness(genes, "test")


def test_visible_stats_ignore_hidden_rows():
    values = np.array([[1.0], [2.0], [1e6], [1e6]])
    mask = np.array([True, True, False, False])
    visible = values.copy()
    visible[~mask] = np.nan
    mean, std = obs.visible_mean_std(visible, mask)
    assert mean[0] == pytest.approx(1.5)
    assert std[0] == pytest.approx(0.5)


def test_library_excludes_witness_and_hidden():
    counts = sparse.csr_matrix(
        np.array(
            [
                [10.0, 1.0, 5.0],
                [2.0, 7.0, 3.0],
                [100.0, 100.0, 100.0],
            ]
        )
    )
    mask = np.array([True, True, False])
    gene_mask = np.array([True, True, False])  # third gene is a witness
    lib = obs.visible_library_size(counts, mask, gene_mask)
    assert lib[0] == pytest.approx(11.0)
    assert lib[1] == pytest.approx(9.0)
    assert np.isnan(lib[2])


def test_save_load_roundtrip(tmp_path):
    genes = np.array(["g1", "g2"])
    roles = make_roles(genes, set())
    bundle = obs.ObservationBundle(
        sample_id="s",
        unit_id="u",
        source_id="src",
        patient_status="unknown",
        barcodes=np.array(["b1", "b2"]),
        source_rows=np.array([0, 1]),
        xy=np.array([[0.0, 0.0], [1.0, 1.0]]),
        coordinate_unit="synthetic_pitch",
        scale_provenance={"source": "test"},
        tissue_support=np.array([True, True]),
        visible_mask=np.array([True, False]),
        genes=genes,
        gene_roles=roles,
        values_visible=np.array([[1.0, 2.0], [np.nan, np.nan]]),
        field_ids=["f1"],
        gene_sets={"f1": ["g1"]},
    )
    obs.save_observation(bundle, tmp_path)
    loaded, sidecar = obs.load_observation(tmp_path)
    assert loaded.values_visible.shape == bundle.values_visible.shape
    assert sidecar["contains_hidden_values"] is False
    assert obs.fingerprint_np(loaded.values_visible) == obs.fingerprint_np(bundle.values_visible)


def test_load_refuses_hidden_payload(tmp_path):
    np.savez(tmp_path / "observation.npz", hidden_values=np.zeros(3))
    (tmp_path / "observation.json").write_text("{}")
    with pytest.raises(obs.MaskViolation):
        obs.load_observation(tmp_path)


def test_split_for_evaluation_returns_separate_object():
    values = np.arange(12.0).reshape(4, 3)
    mask = np.array([True, False, True, False])
    visible, truth = obs.split_for_evaluation(values, mask)
    assert np.isnan(visible[~mask]).all()
    assert np.array_equal(truth.rows, np.array([1, 3]))
    assert truth.values.shape == (2, 3)


def test_pipeline_inputs_are_hidden_invariant():
    rng = np.random.default_rng(0)
    n, g = 120, 40
    counts = sparse.csr_matrix(rng.poisson(2.0, size=(n, g)).astype(np.float64))
    genes = np.array([f"g{i}" for i in range(g)])
    xy = np.column_stack([np.repeat(np.arange(10), 12), np.tile(np.arange(12), 10)]).astype(float)
    support = np.ones(n, dtype=bool)
    mask = np.ones(n, dtype=bool)
    mask[:20] = False  # a hidden block

    witness = {f"g{i}" for i in range(20, 25)}
    spec = InputSpec(fields={"f1": [f"g{i}" for i in range(0, 5)], "f2": [f"g{i}" for i in range(5, 10)]},
                     witness_genes=frozenset(witness), k_neighbors=4, min_library_genes=5)

    base = build_inputs(counts, genes, xy, support, mask, spec)

    poisoned = counts.copy().tolil()
    poisoned[:20] = np.asarray(poisoned[:20].todense()) + 1e6
    poisoned = poisoned.tocsr()
    other = build_inputs(poisoned, genes, xy, support, mask, spec)

    for key in base.arrays:
        a, b = base.arrays[key], other.arrays[key]
        if a.dtype.kind in "UO":
            assert np.array_equal(a, b), key
        else:
            assert np.array_equal(a, b, equal_nan=True), key


def test_leaky_builder_would_be_caught():
    """Negative control: a builder that uses pre-mask statistics must differ."""
    rng = np.random.default_rng(1)
    n, g = 60, 20
    counts = rng.poisson(3.0, size=(n, g)).astype(float)
    mask = np.ones(n, dtype=bool)
    mask[:10] = False
    leaky_before = counts.mean(axis=0)  # computed before masking (the bug)
    counts2 = counts.copy()
    counts2[:10] += 1e5
    leaky_after = counts2.mean(axis=0)
    assert not np.allclose(leaky_before, leaky_after)

    honest_before = counts[mask].mean(axis=0)
    honest_after = counts2[mask].mean(axis=0)
    assert np.allclose(honest_before, honest_after)


def test_external_library_excludes_all_field_genes():
    rng = np.random.default_rng(3)
    n, g = 80, 400
    counts = sparse.csr_matrix(rng.poisson(5.0, size=(n, g)).astype(np.float64))
    genes = np.array([f"g{i}" for i in range(g)])
    xy = np.column_stack([np.repeat(np.arange(10), 8), np.tile(np.arange(8), 10)]).astype(float)
    support = np.ones(n, dtype=bool)
    mask = np.ones(n, dtype=bool)
    fields = {"f1": [f"g{i}" for i in range(10)], "f2": [f"g{i}" for i in range(10, 20)]}
    witness = frozenset({f"g{i}" for i in range(20, 30)})
    spec = InputSpec(fields=fields, witness_genes=witness, k_neighbors=4)
    base = build_inputs(counts, genes, xy, support, mask, spec)
    assert base.meta["library_scope"] == "all_genes_excluding_all_analyzed_field_genes"
    assert base.meta["n_library_genes"] == g - 30

    # perturbing a field (A) gene must not move the denominator
    perturbed = counts.copy().tolil()
    perturbed[:, 0] = 0.0
    changed = build_inputs(perturbed.tocsr(), genes, xy, support, mask, spec)
    assert np.array_equal(base.arrays["library"], changed.arrays["library"])
