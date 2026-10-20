"""Poisoning harness: prove that hidden values cannot reach model inputs.

Runtime checks (added to a real or synthetic section):
  P1 hidden expression -> extreme values
  P2 hidden expression -> random permutation
  P3 witness (B) genes -> permuted across spots
Plus a static check that this package contains no label/TLS reader.

Every run also asserts the *reverse* direction: the evaluation quantity must
change, otherwise the test would be vacuous (evaluation never read the hidden
region).
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy import sparse

from .observation import ContractError, fingerprint_arrays, max_abs_delta
from .pipeline import InputSpec, build_inputs, build_truth, ridge_predict

# reader/loader patterns only: a metadata string naming a legacy asset is not a
# reader, so the scan looks for imports and openers, not for the bare word.
FORBIDDEN_SOURCE_PATTERNS = (
    r"import\s+tls",
    r"from\s+tls",
    r"load_tls",
    r"read_tls",
    r"tls_label",
    r"human_label",
    r"annotation_tsv",
    r"hooknet",
)
SELF_EXCLUDED = {"poisoning.py"}


def truth_delta(t0: np.ndarray, t1: np.ndarray) -> dict[str, Any]:
    """Finite statistics of the evaluation change; None when nothing is finite."""
    diff = np.asarray(t0, dtype=np.float64) - np.asarray(t1, dtype=np.float64)
    finite = np.isfinite(diff)
    return {
        "n_finite_hidden_truth": int(finite.sum()),
        "truth_max_abs_delta": float(np.max(np.abs(diff[finite]))) if finite.any() else None,
        "base_mean_hidden": float(np.nanmean(t0)) if np.isfinite(t0).any() else None,
        "poisoned_mean_hidden": float(np.nanmean(t1)) if np.isfinite(t1).any() else None,
    }


def metric_change(summaries: Mapping[str, Mapping[str, Any]]) -> bool:
    """The reverse assertion: at least one FINITE, NON-ZERO evaluation change.

    An all-missing hidden truth yields ``None`` and must fail: otherwise a vacuous
    test would pass.
    """
    deltas = [s.get("truth_max_abs_delta") for s in summaries.values()]
    return any(d is not None and np.isfinite(d) and d > 0 for d in deltas)


@dataclass
class PoisonResult:
    scenario: str
    description: str
    arrays_compared: int
    max_abs_delta: float
    per_array_delta: dict[str, float]
    identical: bool

    def to_dict(self) -> dict[str, Any]:
        return {
            "scenario": self.scenario,
            "description": self.description,
            "arrays_compared": self.arrays_compared,
            "max_abs_delta": self.max_abs_delta,
            "identical": self.identical,
            "per_array_delta": self.per_array_delta,
            "changed_arrays": [k for k, v in self.per_array_delta.items() if v != 0.0],
        }


def poison_hidden_extreme(counts: sparse.csr_matrix, hidden_rows: np.ndarray) -> sparse.csr_matrix:
    out = counts.copy().tolil()
    if len(hidden_rows):
        block = np.asarray(out[hidden_rows].todense())
        scale = max(float(block.max()) if block.size else 0.0, 1.0)
        out[hidden_rows] = np.asarray(block) + 1000.0 * scale
    return out.tocsr()


def poison_hidden_permute(counts: sparse.csr_matrix, hidden_rows: np.ndarray, seed: int = 0) -> sparse.csr_matrix:
    out = counts.copy().tolil()
    if len(hidden_rows) > 1:
        rng = np.random.default_rng(seed)
        block = np.asarray(out[hidden_rows].todense())
        perm = rng.permutation(len(hidden_rows))
        out[hidden_rows] = block[perm]
    return out.tocsr()


def poison_witness(counts: sparse.csr_matrix, genes: np.ndarray, witness: Sequence[str], seed: int = 0) -> sparse.csr_matrix:
    """Permute witness-gene values across spots, preserving each column's values."""
    cols = np.flatnonzero(np.isin(np.asarray(genes).astype(str), list(witness)))
    if not len(cols):
        return counts.copy()
    coo = counts.tocoo()
    rng = np.random.default_rng(seed)
    rows = coo.row.copy()
    for col in cols:
        selected = coo.col == col
        n = int(selected.sum())
        if n > 1:
            rows[selected] = rng.permutation(coo.row[selected])
    return sparse.coo_matrix((coo.data, (rows, coo.col)), shape=counts.shape).tocsr()


def compare_inputs(base: Mapping[str, np.ndarray], other: Mapping[str, np.ndarray], scenario: str, description: str) -> PoisonResult:
    if set(base) != set(other):
        missing = sorted(set(base) ^ set(other))
        raise ContractError(f"input key mismatch after poisoning: {missing[:5]}")
    deltas: dict[str, float] = {}
    for key in sorted(base):
        a, b = base[key], other[key]
        if a.dtype.kind in "US" or b.dtype.kind in "US" or a.dtype == object:
            deltas[key] = 0.0 if np.array_equal(a, b) else float("inf")
        else:
            deltas[key] = max_abs_delta(a, b)
    worst = max(deltas.values()) if deltas else 0.0
    return PoisonResult(
        scenario=scenario,
        description=description,
        arrays_compared=len(deltas),
        max_abs_delta=float(worst),
        per_array_delta=deltas,
        identical=bool(worst == 0.0),
    )


def static_source_check(package_dir: Path) -> dict[str, Any]:
    """Scan our own package for label/TLS readers (we should read none)."""
    hits: dict[str, list[str]] = {}
    for path in sorted(Path(package_dir).rglob("*.py")):
        if path.name in SELF_EXCLUDED:
            continue
        text = path.read_text()
        found = [p for p in FORBIDDEN_SOURCE_PATTERNS if re.search(p, text, flags=re.IGNORECASE)]
        if found:
            hits[str(path)] = found
    return {"forbidden_patterns": list(FORBIDDEN_SOURCE_PATTERNS), "hits": hits, "clean": not hits}


def run_poisoning(
    counts: sparse.csr_matrix,
    genes: np.ndarray,
    xy: np.ndarray,
    tissue_support: np.ndarray,
    visible_mask: np.ndarray,
    spec: InputSpec,
    *,
    package_dir: Path,
    section_id: str,
    seed: int = 0,
) -> dict[str, Any]:
    hidden_rows = np.flatnonzero(~np.asarray(visible_mask, dtype=bool))
    base = build_inputs(counts, genes, xy, tissue_support, visible_mask, spec)
    base_truth = build_truth(counts, genes, hidden_rows, spec)
    base_pred = {fid: ridge_predict(base, fid) for fid in spec.fields}

    results: list[PoisonResult] = []

    poisoned_extreme = poison_hidden_extreme(counts, hidden_rows)
    extreme_inputs = build_inputs(poisoned_extreme, genes, xy, tissue_support, visible_mask, spec)
    results.append(
        compare_inputs(
            base.arrays,
            extreme_inputs.arrays,
            "P1_hidden_extreme",
            "hidden expression set to +1000x max; inputs must be bit-identical",
        )
    )

    poisoned_perm = poison_hidden_permute(counts, hidden_rows, seed=seed)
    perm_inputs = build_inputs(poisoned_perm, genes, xy, tissue_support, visible_mask, spec)
    results.append(
        compare_inputs(
            base.arrays,
            perm_inputs.arrays,
            "P2_hidden_permute",
            "hidden rows permuted among themselves; inputs must be bit-identical",
        )
    )

    poisoned_witness = poison_witness(counts, genes, sorted(spec.witness_genes), seed=seed)
    witness_inputs = build_inputs(poisoned_witness, genes, xy, tissue_support, visible_mask, spec)
    results.append(
        compare_inputs(
            base.arrays,
            witness_inputs.arrays,
            "P3_witness_permute",
            "witness (B) genes permuted across spots; inputs must be bit-identical",
        )
    )

    # reverse direction: evaluation must move
    extreme_truth = build_truth(poisoned_extreme, genes, hidden_rows, spec)
    prediction_delta = 0.0
    for fid in spec.fields:
        pred_extreme = ridge_predict(extreme_inputs, fid)
        prediction_delta = max(prediction_delta, max_abs_delta(base_pred[fid], pred_extreme))
    metric = {}
    for fid in sorted(spec.fields):
        t0, t1 = base_truth[fid], extreme_truth[fid]
        if len(t0) == 0:
            continue
        metric[fid] = truth_delta(t0, t1)
    metric_changed = metric_change(metric)

    field_genes_present = set(map(str, genes))
    analysed = set(map(str, spec.witness_genes))
    for members in spec.fields.values():
        analysed.update(map(str, members))
    library_info = {
        "n_field_gene_union_present": int(len(analysed & field_genes_present)),
        "n_library_genes": int(base.meta.get("n_library_genes", -1)),
        "library_scope": base.meta.get("library_scope"),
        "detection_scope": base.meta.get("detection_scope"),
    }
    entry_points = {
        "library_size": "library" in base.arrays,
        "depth": "depth" in base.arrays,
        "detection": "detection" in base.arrays,
        "whole_slide_stats": "slide_mean" in base.arrays and "slide_std" in base.arrays,
        "expression_weighted_adjacency": any(k.startswith("neighbor::") for k in base.arrays),
        "structural_adjacency": any(k.startswith("neighbor_struct::") for k in base.arrays),
        "query_eligibility": "query_eligible" in base.arrays,
    }
    covered = all(entry_points.values())

    passed = (
        all(r.identical for r in results)
        and metric_changed
        and prediction_delta == 0.0
        and covered
        and static_source_check(package_dir)["clean"]
    )

    return {
        "section_id": section_id,
        "n_obs": int(counts.shape[0]),
        "n_visible": int(np.asarray(visible_mask, dtype=bool).sum()),
        "n_hidden": int(len(hidden_rows)),
        "scenarios": [r.to_dict() for r in results],
        "input_fingerprints": fingerprint_arrays(base.arrays),
        "entry_points": entry_points,
        "library": library_info,
        "evaluation": {
            "prediction_max_abs_delta": float(prediction_delta),
            "truth_under_poisoning": metric,
            "metric_changed": bool(metric_changed),
        },
        "static_source_check": static_source_check(package_dir),
        "passed": bool(passed),
        "notes": [
            "input arrays are compared bit-for-bit (dtype+shape+bytes)",
            "the reverse assertion requires the evaluation quantity to change, "
            "otherwise the hidden region was never read by the metric",
        ],
    }
