"""Deterministic A/B gene splitting and the gene-role registry.

Rules (05_COMMON_PROTOCOL §4, addendum §2):
  * split by gene *name* with a fixed seed, before any expression is read;
  * the witness (B) set is the union over all analysed fields and is excluded from
    every input layer, including the normalization denominator;
  * a field whose A half has fewer than ``min_measurable`` genes is not testable;
  * aliasing two names for one signature does not create a second field.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .observation import (
    ROLE_A_FEATURE,
    ROLE_BACKGROUND,
    ROLE_B_WITNESS,
    ContractError,
    GeneRoles,
)


def _stable_key(gene: str, seed: int) -> str:
    return hashlib.sha256(f"{seed}::{gene}".encode()).hexdigest()


def split_field_genes(
    members: Sequence[str],
    *,
    seed: int = 0,
    a_fraction: float = 2.0 / 3.0,
    min_measurable: int = 3,
) -> tuple[list[str], list[str]]:
    """Split one field's genes into (A features, B witness) by name only."""
    unique = sorted({str(g) for g in members})
    if not unique:
        raise ContractError("empty field definition")
    ordered = sorted(unique, key=lambda g: (_stable_key(g, seed), g))
    n_a = int(round(len(ordered) * a_fraction))
    n_a = max(1, min(len(ordered) - 1, n_a)) if len(ordered) > 1 else len(ordered)
    a_genes = sorted(ordered[:n_a])
    b_genes = sorted(ordered[n_a:])
    if len(ordered) >= 2 * min_measurable and len(a_genes) < min_measurable:
        raise ContractError(
            f"field A half has {len(a_genes)} genes (< min_measurable={min_measurable})"
        )
    return a_genes, b_genes


@dataclass
class FieldSplit:
    field_id: str
    a_genes: list[str]
    b_genes: list[str]
    n_defined: int
    n_present: int

    def to_dict(self) -> dict:
        return {
            "field_id": self.field_id,
            "n_defined": self.n_defined,
            "n_present": self.n_present,
            "n_a": len(self.a_genes),
            "n_b": len(self.b_genes),
        }


def build_field_splits(
    fields: Mapping[str, Sequence[str]],
    genes_present: Sequence[str],
    *,
    seed: int = 0,
    a_fraction: float = 2.0 / 3.0,
    min_measurable: int = 3,
) -> dict[str, FieldSplit]:
    present = set(map(str, genes_present))
    splits: dict[str, FieldSplit] = {}
    for field_id in sorted(fields):
        members = [g for g in fields[field_id]]
        a_genes, b_genes = split_field_genes(
            members, seed=seed, a_fraction=a_fraction, min_measurable=min_measurable
        )
        splits[field_id] = FieldSplit(
            field_id=field_id,
            a_genes=[g for g in a_genes if g in present],
            b_genes=[g for g in b_genes if g in present],
            n_defined=len(set(members)),
            n_present=len({g for g in members if g in present}),
        )
    return splits


def build_gene_roles(splits: Mapping[str, FieldSplit], genes: Sequence[str]) -> GeneRoles:
    """Union of all witness genes wins over feature status."""
    a_union: set[str] = set()
    b_union: set[str] = set()
    for split in splits.values():
        a_union.update(split.a_genes)
        b_union.update(split.b_genes)
    roles: dict[str, str] = {}
    for gene in map(str, genes):
        if gene in b_union:
            roles[gene] = ROLE_B_WITNESS
        elif gene in a_union:
            roles[gene] = ROLE_A_FEATURE
        else:
            roles[gene] = ROLE_BACKGROUND
    return GeneRoles.from_mapping(roles)


def feature_genes(fields: Mapping[str, Sequence[str]], splits: Mapping[str, FieldSplit]) -> dict[str, list[str]]:
    """Field definitions restricted to their A half (what the pipeline scores)."""
    return {fid: list(splits[fid].a_genes) for fid in fields}


def assert_witness_union_excluded(
    genes: Sequence[str], splits: Mapping[str, FieldSplit], context: str = ""
) -> None:
    present = set(map(str, genes))
    witness = set().union(*(set(s.b_genes) for s in splits.values())) if splits else set()
    leaked = sorted(witness.intersection(present))
    if leaked:
        where = f" ({context})" if context else ""
        raise ContractError(f"witness genes present in an input axis{where}: {leaked[:8]}")
