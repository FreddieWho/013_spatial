import importlib.util
from pathlib import Path
import numpy as np
from scipy import sparse

spec=importlib.util.spec_from_file_location('peripheral',Path(__file__).resolve().parents[1]/'scripts/tls_peripheral.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_hidden_counts_cannot_change_visible_scores():
    rng=np.random.default_rng(12);X=rng.poisson(2,size=(12,100)).astype(float)
    genes=np.array([f'g{i}' for i in range(100)]);panel={'a':list(genes[:20]),'b':list(genes[20:40])}
    visible=np.arange(12)>=3
    a=m.visible_scores(sparse.csr_matrix(X),genes,visible,panel)
    X[~visible]=rng.uniform(1e6,1e9,size=(3,100))
    b=m.visible_scores(sparse.csr_matrix(X),genes,visible,panel)
    np.testing.assert_array_equal(a,b)


def test_mask_edges_and_distance_bands():
    xy=np.column_stack([np.array([0,100,130,200,230,300,400,700,999,1000]),np.zeros(10)])
    labels=np.array(['TLS']+['']*9)
    d,v,n,f=m.geometry(xy,labels,130)
    assert np.flatnonzero(n).tolist()==[3,4,5]
    assert np.flatnonzero(f).tolist()==[7,8]
    _,_,n2,_=m.geometry(xy,labels,230)
    assert np.flatnonzero(n2).tolist()==[5]
    assert not v[2]


def test_gaussian_generator_preserves_expected_variance():
    # Ensemble, not spatial demeaning: the latter would shrink long-range variance.
    coords=np.array([[0,0],[0,2],[1,1]])
    Y=m.gaussian_lattice(coords,150,1000,np.random.default_rng(1))
    assert np.all((Y.var(axis=1)>.85)&(Y.var(axis=1)<1.15))
    assert .65<np.corrcoef(Y)[0,1]<.92


def test_null_contrast_variance_includes_covariance_and_constant_cancels():
    xy=np.array([[0.,0],[100.,0],[200.,0],[400.,0]])
    w=np.array([.5,.5,-.5,-.5])
    *_,var=m.fit_setup(xy,w,1)
    assert var[0]==1
    assert var[-1]<var[1]
    # Constant offset must not change the empirical variogram or statistic.
    Y=np.random.default_rng(5).normal(size=(4,20))
    setup=m.fit_setup(xy,w,1)
    p,e=m.fit_null(Y,setup,w);p2,e2=m.fit_null(Y+50,setup,w)
    np.testing.assert_allclose(p,p2,atol=1e-12)
