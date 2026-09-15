from pathlib import Path

import numpy as np

from r16.section_io import usz_tls_labels


def test_usz_tls_labels_align_and_exclude_unassigned():
    root = Path("/home/huyudi/013_spatial")
    y = usz_tls_labels("TLS_VISIUM_USZ::LC3", ("NOPE-1",), root)
    if y is None:
        return
    assert y.shape == (1,) and np.isnan(y[0])
    # real barcodes from the h5ad should yield TLS or 0, not all-nan if file present
    import h5py
    from r04.io_contract import _decode_strings
    path = root / "data/other_sources/zenodo_usz_tls_visium/TLS_VISIUM_USZ/h5ad_preprocessed/LC3.h5ad"
    if not path.exists():
        return
    with h5py.File(path, "r") as handle:
        bcs = _decode_strings(handle["obs"]["_index"][:20])
    y2 = usz_tls_labels("TLS_VISIUM_USZ::LC3", bcs, root)
    assert y2 is not None and y2.shape == (20,)
    assert np.isfinite(y2).any()
    assert set(np.unique(y2[np.isfinite(y2)])).issubset({0.0, 1.0})
