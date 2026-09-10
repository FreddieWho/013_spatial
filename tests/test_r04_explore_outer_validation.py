import numpy as np
import pytest

from scripts.r04_explore_outer_validation import (
    apply_frozen_shared,
    assemble_full_effect,
    fit_shared_on_training,
    load_heldout_gt_catalog,
    permutation_p,
    reconstruct_split_effect,
)


def _fake_split(n_spots=12, n_genes=6, k=2, seed=0, gene_offset=0):
    rng = np.random.default_rng(seed)
    return {
        "field_mean": rng.standard_normal((n_spots, k)),
        "evaluation_loading": np.abs(rng.standard_normal((n_genes, k))) + 0.1,
        "factor_amplitude": np.abs(rng.standard_normal(k)) + 0.5,
        "section_ids": np.asarray(["s1"] * (n_spots // 2) + ["s2"] * (n_spots - n_spots // 2)),
        "evaluation_indices": np.arange(gene_offset, gene_offset + n_genes),
        "evaluation_gene_ids": [f"g{i}" for i in range(gene_offset, gene_offset + n_genes)],
        "frozen_gene_ids": [f"g{i}" for i in range(12)],
        "barcodes": np.asarray([f"bc{i}" for i in range(n_spots)]),
        "patient_ids": np.asarray(["p1"] * (n_spots // 2) + ["p2"] * (n_spots - n_spots // 2)),
    }


def test_reconstruct_matches_manual_rate_math():
    split = _fake_split()
    effect = reconstruct_split_effect(
        split["field_mean"], split["evaluation_loading"],
        split["factor_amplitude"], split["section_ids"])
    assert effect.shape == (12, 6)
    assert np.isfinite(effect).all()
    # manual spot check on section s1 rows
    rows = np.flatnonzero(split["section_ids"] == "s1")
    field = np.exp(np.clip(split["field_mean"][rows], -30.0, 30.0))
    centered = field - field.mean(axis=0, keepdims=True)
    expected = centered @ (split["evaluation_loading"]
                           * split["factor_amplitude"][None, :]).T
    assert np.allclose(effect[rows], expected)


def test_assemble_places_genes_and_rejects_mismatch():
    first = _fake_split(seed=1, gene_offset=0)
    second = _fake_split(seed=2, gene_offset=6)
    full, genes, sids, bcs = assemble_full_effect([first, second])
    assert full.shape == (12, 12)
    assert genes == [f"g{i}" for i in range(12)]
    assert np.isfinite(full).all()
    bad = dict(second)
    bad["section_ids"] = np.asarray(["sX"] * 12)
    with pytest.raises(ValueError):
        assemble_full_effect([first, bad])


def test_fit_apply_roundtrip_and_null_shape():
    rng = np.random.default_rng(3)
    n = 120
    patients = np.asarray(["p1"] * 60 + ["p2"] * 60)
    signal = np.concatenate([np.zeros(60), np.ones(60)])
    z = np.column_stack([signal + 0.3 * rng.standard_normal(n),
                         rng.standard_normal(n)])
    order = rng.permutation(n)
    base = np.concatenate([np.zeros(45), np.ones(15), np.zeros(15), np.ones(45)])
    labels = {"TLS": base[order],
              "TUMOR_STROMA_BOUNDARY": base[rng.permutation(n)]}
    w = np.full(n, 1.0 / n)
    frozen = fit_shared_on_training(z, labels, w)
    assert frozen["beta"].shape == (2,)
    assert set(frozen["mappings"]) == {"TLS", "TUMOR_STROMA_BOUNDARY"}
    out = apply_frozen_shared(z, frozen, labels, w)
    for name in labels:
        assert 0.0 <= out["auc_by_structure"][name] <= 1.0
    null = permutation_p(z, frozen, labels, w, seed=5, draws=10,
                         group_ids=np.asarray(["s1"] * 60 + ["s2"] * 60))
    assert len(null["TLS"]) == 10


def test_heldout_catalog_keeps_only_confirmatory(tmp_path):
    annot = tmp_path / "annot.csv"
    annot.write_text("barcode,pathology_annotation\nbc1,TLS positive\nbc2,stroma\n", encoding="utf-8")
    registry = tmp_path / "registry.tsv"
    registry.write_text(
        "structure_id\tphysical_unit_id\tconfirmation_status\tallowed_use\tgt_geometry_locator\n"
        f"TLS\tsecA\tCONFIRMATORY\texploratory readout\t{annot}#pathology_annotation=TLS positive\n"
        f"TLS\tsecB\tNONCONFIRMATORY_CONTEXT_NORMAL_MUCOSA\tdiscovery only\t{annot}#pathology_annotation=TLS positive\n",
        encoding="utf-8")
    catalog = load_heldout_gt_catalog(registry, ["secA", "secB"])
    assert set(catalog["structures"]["TLS"]["sections"]) == {"secA"}
    assert catalog["skipped_nonconfirmatory"] == 1


def test_apply_frozen_shared_ignores_unannotated_spots():
    from scripts.r04_explore_outer_validation import apply_frozen_shared
    rng = np.random.default_rng(11)
    n = 80
    z = rng.standard_normal((n, 2))
    frozen = {
        "beta": np.array([1.0, 0.0]),
        "feature_mean": np.zeros(2),
        "feature_scale": np.ones(2),
        "intercept": 0.0,
        "mappings": {
            "TLS": {"intercept": 0.0, "slope": 1.0},
            "TUMOR_STROMA_BOUNDARY": {"intercept": 0.0, "slope": 1.0},
        },
    }
    signal = (z[:, 0] > 0).astype(float)
    labels = {"TLS": signal.copy(), "TUMOR_STROMA_BOUNDARY": signal.copy()}
    labels["TLS"][::2] = np.nan  # half unannotated must not poison the AUC
    out = apply_frozen_shared(z, frozen, labels, np.ones(n))
    assert out["auc_by_structure"]["TLS"] is not None
    assert out["auc_by_structure"]["TLS"] > 0.9


def test_permutation_p_keeps_annotated_set_fixed():
    from scripts.r04_explore_outer_validation import permutation_p
    rng = np.random.default_rng(12)
    n = 60
    z = rng.standard_normal((n, 2))
    frozen = {
        "beta": np.array([1.0, 0.0]),
        "feature_mean": np.zeros(2),
        "feature_scale": np.ones(2),
        "intercept": 0.0,
        "mappings": {
            "TLS": {"intercept": 0.0, "slope": 1.0},
            "TUMOR_STROMA_BOUNDARY": {"intercept": 0.0, "slope": 1.0},
        },
    }
    labels = {"TLS": np.full(n, np.nan),
              "TUMOR_STROMA_BOUNDARY": (rng.random(n) > 0.5).astype(float)}
    labels["TLS"][:30] = (rng.random(30) > 0.5).astype(float)
    null = permutation_p(z, frozen, labels, np.ones(n), seed=3, draws=8,
                         group_ids=np.asarray(["s1"] * 30 + ["s2"] * 30))
    assert len(null["TLS"]) == 8
    assert len(null["TUMOR_STROMA_BOUNDARY"]) == 8
