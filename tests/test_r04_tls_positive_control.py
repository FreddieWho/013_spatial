import numpy as np
import pytest

from scripts.r04_tls_positive_control import _uauc, PLASMA_GENES, SINGLE_GENES


def test_uauc_perfect_and_null():
    y = np.array([0, 0, 1, 1])
    assert _uauc(y, np.array([0.1, 0.2, 0.8, 0.9])) == pytest.approx(1.0)
    assert _uauc(y, np.array([0.5, 0.5, 0.5, 0.5])) == pytest.approx(0.5)
    assert _uauc(y, np.array([0.9, 0.8, 0.2, 0.1])) == pytest.approx(0.0)


def test_gene_sets_are_panel_compatible():
    # plasma/single genes were verified present in the frozen panel; guard set
    assert len(set(PLASMA_GENES)) == len(PLASMA_GENES)
    assert len(set(SINGLE_GENES)) == len(SINGLE_GENES)
    assert not set(PLASMA_GENES) & set(SINGLE_GENES)
