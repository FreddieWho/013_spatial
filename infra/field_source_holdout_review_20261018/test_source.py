import sys
from pathlib import Path
import numpy as np
from scipy import sparse
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import field_source_holdout_20261018 as run


def test_source_unknowns_excluded_and_units_preserved():
    rows=[dict(section_id=f's{i}',unit_id=f'u{i//2}',fold=i//2,source_key=('https://study/a' if i<2 else 'https://study/b' if i<4 else 'title:unknown'),organ='Bowel') for i in range(10)]
    ss=run.splits(rows)
    assert len(ss)==7
    for s in ss[5:]:
        assert set(s['train']+s['test'])=={'s0','s1','s2','s3'}
        assert not set(s['train'])&set(s['test'])


def test_input_normalization_ignores_hidden_counts():
    x=sparse.csr_matrix([[2.,8.,10.],[0.,0.,100.]])
    z,qc=run.follow.gene_input(x,['a','b','hidden'],['a','b'])
    x2=sparse.csr_matrix([[2.,8.,999999.],[0.,0.,1.]])
    z2,qc2=run.follow.gene_input(x2,['a','b','hidden'],['a','b'])
    np.testing.assert_array_equal(z,z2);np.testing.assert_array_equal(qc,qc2)
    np.testing.assert_allclose(np.expm1(z[0]),[2000.,8000.],rtol=1e-6)
    assert not z[1].any()
