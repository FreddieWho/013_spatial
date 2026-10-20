import numpy as np
import pytest

from field_analyzer import roles as role_mod
from field_analyzer.observation import ContractError


def test_split_is_order_independent_and_disjoint():
    members = [f"g{i}" for i in range(30)]
    a1, b1 = role_mod.split_field_genes(members, seed=0)
    a2, b2 = role_mod.split_field_genes(list(reversed(members)), seed=0)
    assert a1 == a2 and b1 == b2
    assert not set(a1) & set(b1)
    assert set(a1) | set(b1) == set(members)
    assert len(a1) == 20 and len(b1) == 10


def test_split_changes_with_seed():
    members = [f"g{i}" for i in range(40)]
    a0, _ = role_mod.split_field_genes(members, seed=0)
    a1, _ = role_mod.split_field_genes(members, seed=1)
    assert a0 != a1


def test_witness_union_wins_over_feature():
    fields = {"f1": [f"g{i}" for i in range(9)], "f2": [f"g{i}" for i in range(9, 18)]}
    splits = role_mod.build_field_splits(fields, [f"g{i}" for i in range(18)], seed=0)
    genes = [f"g{i}" for i in range(18)] + ["bg1"]
    roles = role_mod.build_gene_roles(splits, genes)
    witness = set().union(*(set(s.b_genes) for s in splits.values()))
    for g in witness:
        assert roles.role_of(g) == "B_WITNESS"
    assert roles.role_of("bg1") == "BACKGROUND"


def test_feature_axis_excludes_witness_union():
    fields = {"f1": [f"g{i}" for i in range(9)]}
    splits = role_mod.build_field_splits(fields, [f"g{i}" for i in range(9)], seed=0)
    scored = role_mod.feature_genes(fields, splits)["f1"]
    assert not set(scored) & set(splits["f1"].b_genes)
    with pytest.raises(ContractError):
        role_mod.assert_witness_union_excluded(scored + splits["f1"].b_genes, splits)


def test_field_with_too_few_present_genes_is_reported():
    fields = {"f1": ["g1", "g2", "g3"]}
    splits = role_mod.build_field_splits(fields, ["g1"], seed=0)
    assert splits["f1"].n_present == 1
    assert len(splits["f1"].a_genes) <= 1
