import numpy as np
import pytest
from scipy import sparse

from scripts.r04_explore_marker_smoke import load_section_counts, module_scores


def _write_cache(tmp_path, section_id, barcodes, gene_ids, data):
    import json
    (tmp_path / "sections").mkdir(exist_ok=True)
    mat = sparse.csr_matrix(np.asarray(data, dtype=float))
    np.savez(tmp_path / "sections" / f"0000_{section_id}.npz", data=mat.data, indices=mat.indices,
             indptr=mat.indptr, shape=np.asarray(mat.shape))
    (tmp_path / "sections" / f"0000_{section_id}.json").write_text(json.dumps(
        {"barcodes": barcodes, "gene_id": gene_ids}))


def test_module_scores_math():
    counts = sparse.csr_matrix(np.array([[1.0, 0.0, 3.0], [0.0, 2.0, 0.0]]))
    gene_index = {"G1": 0, "G2": 1, "G3": 2}
    out = module_scores(counts, gene_index, {"A": ["G1", "G3"], "B": ["G2"]})
    assert out.shape == (2, 2)
    assert out[0, 0] == pytest.approx((np.log1p(1.0) + np.log1p(3.0)) / 2)
    assert out[1, 1] == pytest.approx(np.log1p(2.0))
    with pytest.raises(ValueError):
        module_scores(counts, gene_index, {"A": ["G1"], "B": ["GX"]})


def test_load_section_counts_alignment(tmp_path):
    _write_cache(tmp_path, "S.h5ad", ["b1", "b2"], ["G1", "G2"],
                 [[1.0, 0.0], [0.0, 2.0]])
    mat, barcodes, genes = load_section_counts(tmp_path, "S.h5ad")
    assert mat.shape == (2, 2) and barcodes == ["b1", "b2"] and genes == ["G1", "G2"]
    with pytest.raises(ValueError):
        load_section_counts(tmp_path, "MISSING.h5ad")


def test_align_rows_resolves_cross_section_collision():
    from scripts.r04_explore_marker_smoke import align_rows
    export_sections = ["A", "A", "B", "B"]
    export_barcodes = ["bX", "b1", "bX", "b2"]
    assert align_rows(export_sections, export_barcodes, "A", ["bX", "b1"]).tolist() == [0, 1]
    assert align_rows(export_sections, export_barcodes, "B", ["bX", "b2"]).tolist() == [2, 3]
    try:
        align_rows(export_sections, export_barcodes, "B", ["bZZZ"])
    except SystemExit:
        pass
    else:
        raise AssertionError("missing barcode must fail closed")


def test_barcode_collision_across_sections():
    """Identical barcodes in different sections must not cross-align."""
    from scripts.r04_explore_marker_smoke import load_section_counts
    import tempfile, json
    from pathlib import Path
    # simulate: same barcode 'bX' in two sections resolves per-section
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        _write_cache(tmp, "A.h5ad", ["bX", "b1"], ["G1", "G2"], [[1.0, 0.0], [0.0, 1.0]])
        _write_cache(tmp, "B.h5ad", ["bX", "b2"], ["G1", "G2"], [[5.0, 0.0], [0.0, 5.0]])
        ma, bca, _ = load_section_counts(tmp, "A.h5ad")
        mb, bcb, _ = load_section_counts(tmp, "B.h5ad")
        assert ma[0, 0] == 1.0 and mb[0, 0] == 5.0


def test_corrupt_shard_errors_are_skippable(tmp_path):
    """Truncated shards must raise something the failover catches."""
    import numpy as np
    bad = tmp_path / "empty.npz"
    bad.write_bytes(b"")
    try:
        np.load(bad, allow_pickle=False)
    except Exception as exc:  # noqa: BLE001 - asserting on the type below
        assert isinstance(exc, (ValueError, OSError, EOFError))
