import h5py
import numpy as np
import pytest
from scipy import sparse

from field_analyzer.adapters import AxisSpec, read_h5ad_counts, scan_h5ad_header
from field_analyzer.observation import ContractError, MaskViolation


def write_h5ad(path, counts, genes, barcodes, xy, *, encoding="csr_matrix", obs_extra=None):
    counts = sparse.csr_matrix(counts)
    with h5py.File(path, "w") as f:
        if encoding is None:
            grp = f.create_group("X")
            grp.attrs["shape"] = np.array(counts.shape)
            grp["data"] = counts.data
            grp["indices"] = counts.indices
            grp["indptr"] = counts.indptr
        else:
            grp = f.create_group("X")
            grp.attrs["encoding-type"] = encoding
            grp.attrs["encoding-version"] = "0.1.0"
            grp.attrs["shape"] = np.array(counts.shape)
            grp["data"] = counts.data
            grp["indices"] = counts.indices
            grp["indptr"] = counts.indptr
        obs = f.create_group("obs")
        obs["_index"] = np.array(barcodes, dtype="S")
        for key, values in (obs_extra or {}).items():
            obs[key] = values
        var = f.create_group("var")
        var["_index"] = np.array(genes, dtype="S")
        obsm = f.create_group("obsm")
        obsm["spatial"] = np.asarray(xy, dtype=np.float64)
        f.create_group("uns")


@pytest.fixture()
def simple(tmp_path):
    path = tmp_path / "s.h5ad"
    counts = np.array([[1, 0, 2], [0, 3, 0], [4, 0, 0], [0, 0, 5]])
    write_h5ad(
        path,
        counts,
        ["g1", "g2", "g3"],
        ["b1", "b2", "b3", "b4"],
        [[0, 0], [1, 0], [0, 1], [1, 1]],
        obs_extra={
            "in_tissue": np.array([1, 1, 0, 1]),
            "n_genes_by_counts": np.array([2, 1, 1, 1]),
            "total_counts": np.array([3, 3, 4, 5]),
        },
    )
    return path


def test_reads_csr_with_explicit_axis(simple):
    table = read_h5ad_counts(simple)
    assert table.counts.shape == (4, 3)
    assert table.report.axis_order == "obs_by_gene"
    assert table.report.aux_kept == ["in_tissue"]
    assert "n_genes_by_counts" in table.report.aux_dropped
    assert table.xy.shape == (4, 2)


def test_transposed_file_is_rejected_not_quietly_fixed(tmp_path):
    path = tmp_path / "t.h5ad"
    counts = np.array([[1, 0], [0, 2], [3, 0]])  # 3 genes x 2 obs
    write_h5ad(path, counts, ["g1", "g2", "g3"], ["b1", "b2"], [[0, 0], [1, 1]])
    with pytest.raises(ContractError) as err:
        read_h5ad_counts(path)
    assert "transposed" in str(err.value)

    table = read_h5ad_counts(path, axis=AxisSpec(axis_order="gene_by_obs"))
    assert table.counts.shape == (2, 3)
    assert table.report.axis_order == "gene_by_obs"


def test_declared_counts_must_match(simple):
    with pytest.raises(ContractError):
        read_h5ad_counts(simple, axis=AxisSpec(n_obs=99))


def test_duplicate_symbols(tmp_path):
    path = tmp_path / "d.h5ad"
    counts = np.array([[1, 2, 3], [4, 5, 6]])
    write_h5ad(path, counts, ["g1", "g1", "g2"], ["b1", "b2"], [[0, 0], [1, 1]])
    with pytest.raises(ContractError):
        read_h5ad_counts(path)
    table = read_h5ad_counts(path, on_duplicate="sum")
    assert table.n_genes == 2
    assert table.counts[0, 0] == 3


def test_zero_depth_is_reported_not_dropped(tmp_path):
    path = tmp_path / "z.h5ad"
    counts = np.array([[0, 0], [1, 1], [2, 0]])
    write_h5ad(path, counts, ["g1", "g2"], ["b1", "b2", "b3"], [[0, 0], [1, 0], [2, 0]])
    table = read_h5ad_counts(path)
    assert table.report.zero_depth_spots == 1
    assert table.n_obs == 3


def test_missing_spatial_is_refused(tmp_path):
    path = tmp_path / "n.h5ad"
    with h5py.File(path, "w") as f:
        grp = f.create_group("X")
        grp.attrs["encoding-type"] = "csr_matrix"
        grp.attrs["shape"] = np.array([1, 1])
        grp["data"] = np.array([1.0])
        grp["indices"] = np.array([0])
        grp["indptr"] = np.array([0, 1])
        obs = f.create_group("obs")
        obs["_index"] = np.array([b"b1"])
        var = f.create_group("var")
        var["_index"] = np.array([b"g1"])
        f.create_group("obsm")
    with pytest.raises(MaskViolation):
        read_h5ad_counts(path)


def test_header_scan_matches_full_read(simple):
    header = scan_h5ad_header(simple)
    table = read_h5ad_counts(simple)
    assert header["shape"] == table.counts.shape
    assert header["n_barcodes"] == table.n_obs
    assert header["declared_axis_consistent"] is True


def test_square_matrix_is_rejected_without_explicit_declaration(tmp_path):
    """F2 regression: a square matrix has two equally valid shapes."""
    path = tmp_path / "sq.h5ad"
    counts = np.array([[1.0, 0.0], [0.0, 2.0]])
    write_h5ad(path, counts, ["g1", "g2"], ["b1", "b2"], [[0, 0], [1, 1]])
    with pytest.raises(ContractError) as err:
        read_h5ad_counts(path)
    assert "square" in str(err.value)

    table = read_h5ad_counts(
        path, axis=AxisSpec(allow_square=True, square_reason="synthetic square test")
    )
    assert table.counts.shape == (2, 2)
    assert any("square matrix accepted" in n for n in table.report.notes)


def test_detection_denominator_excludes_field_genes():
    import numpy as np
    from scipy import sparse

    from field_analyzer.pipeline import InputSpec, build_inputs

    rng = np.random.default_rng(0)
    n, g = 40, 300
    counts = sparse.csr_matrix(rng.poisson(3.0, size=(n, g)).astype(float))
    genes = np.array([f"g{i}" for i in range(g)])
    xy = np.column_stack([np.repeat(np.arange(8), 5), np.tile(np.arange(5), 8)]).astype(float)
    mask = np.ones(n, dtype=bool)
    fields = {"f1": [f"g{i}" for i in range(10)]}
    spec = InputSpec(fields=fields, witness_genes=frozenset({f"g{i}" for i in range(10, 20)}), k_neighbors=4)
    inputs = build_inputs(counts, genes, xy, mask, mask, spec)
    assert inputs.meta["detection_scope"] == inputs.meta["library_scope"]
    assert inputs.meta["n_library_genes"] == g - 20
