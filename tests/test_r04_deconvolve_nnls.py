import json
import numpy as np
import pytest
from scipy import sparse


def _write_ref(tmp_path):
    # 4 genes x 18 cells covering all six Majors; gene order: G1 G2 G3 G4
    majors = ["T", "B", "Mye", "ILC", "Epi", "Stromal"]
    rng = np.random.default_rng(0)
    counts = np.abs(rng.standard_normal((4, 18))) + 0.5
    counts[0, :3] += 10.0  # G1 high in T cells
    mat = sparse.csr_matrix(counts)
    np.savez(tmp_path / "reference_counts.npz",
             counts_data=mat.data, counts_indices=mat.indices,
             counts_indptr=mat.indptr, counts_shape=np.asarray(mat.shape),
             panel_symbols=np.asarray(["G1", "G2", "G3", "G4"]),
             panel_ensg=np.asarray(["E1", "E2", "E3", "E4"]),
             barcodes=np.asarray([f"c{i}" for i in range(18)]),
             major=np.asarray([m for m in majors for _ in range(3)]),
             sub=np.asarray(["a"] * 18),
             patient=np.asarray(["p"] * 18),
             tissue=np.asarray(["Tumor"] * 18))
    (tmp_path / "provenance.json").write_text(json.dumps(
        {"schema": "r04.composition_reference.v1"}))


def _write_counts(tmp_path, genes):
    # 3 spots x len(genes); spot0 looks like class A
    rng = np.random.default_rng(1)
    n = len(genes)
    X = np.abs(rng.standard_normal((3, n))) + 0.5
    X[0, 0] += 20.0
    mat = sparse.csr_matrix(X)
    np.savez(tmp_path / "counts.npz", data=mat.data, indices=mat.indices,
             indptr=mat.indptr, shape=np.asarray(mat.shape))
    (tmp_path / "genes.txt").write_text("\n".join(genes) + "\n")


def test_positional_mapping_with_gaps_and_duplicates(tmp_path, monkeypatch):
    import sys
    from scripts import r04_deconvolve_nnls as mod
    _write_ref(tmp_path)
    # panel order has an empty slot and a duplicated symbol
    genes = ["G1", "", "G2", "G3", "G1", "G4"]
    _write_counts(tmp_path, ["G1", "Gx", "G2", "G3", "G1", "G4"])
    argv = ["prog", "--ref-dir", str(tmp_path),
            "--counts", str(tmp_path / "counts.npz"),
            "--genes", str(tmp_path / "genes.txt"),
            "--output", str(tmp_path / "out.npz"), "--min-shared-genes", "2"]
    monkeypatch.setattr(sys, "argv", argv)
    assert mod.main() == 0
    d = np.load(tmp_path / "out.npz")
    assert d["proportions"].shape == (3, 6)
    assert np.allclose(d["proportions"].sum(axis=1), 1.0)
    # spot0 enriched in G1 -> class A should dominate
    assert d["proportions"][0, 0] > 0.5


def test_differential_markers_recovers_spiked_gene():
    from scripts.r04_deconvolve_nnls import differential_markers
    import numpy as np
    from scipy import sparse
    rng = np.random.default_rng(0)
    # 20 genes x 60 cells; gene 0 high only in class T
    counts = np.abs(rng.standard_normal((20, 60))) + 0.5
    counts[0, :30] += 8.0
    ref = {
        "matrix": sparse.csr_matrix(counts),
        "genes": [f"G{i}" for i in range(20)],
        "major": np.asarray(["T"] * 30 + ["B"] * 30),
    }
    picked = differential_markers(ref, top_n=3)
    assert "G0" in picked
    assert len(picked) <= 18  # 6 classes x 3, fixture has 2 classes present
