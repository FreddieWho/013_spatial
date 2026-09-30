import numpy as np
from scipy.spatial.distance import pdist

from r16.recovery.corrected import (
    NaturalSpline, classify_curve, hex_coordinates, classify_window,
    masked_predict, feature_block, fit_from_moments,
)


def test_natural_spline_represents_line_and_boundary_conditions():
    x = np.linspace(.1, 1, 91)
    spline = NaturalSpline(x)
    X = spline(x)
    beta = np.linalg.lstsq(X, 2 + 3*x, rcond=None)[0]
    q = np.linspace(-.2, 1.2, 120)
    assert np.max(abs(spline(q)@beta-(2+3*q))) < 1e-10
    assert np.max(abs(spline.second_derivative_at_boundaries())) < 1e-10


def test_shape_is_translation_invariant():
    v = np.linspace(-.02, .02, 101)
    assert classify_curve(v,.02)==classify_curve(v+.5,.02)=='low-to-high'


def test_unknown_is_not_negative():
    assert classify_window(np.array([0., np.nan]))=='UNKNOWN'
    assert classify_window(np.zeros(20))=='NEGATIVE'
    assert classify_window(np.r_[np.ones(4),np.zeros(14),[np.nan]*2])=='POSITIVE'


def test_hex_neighbors_are_equidistant():
    xy=hex_coordinates(np.array([[0,0],[0,2],[1,1]]))
    assert np.allclose(pdist(xy),1)


def setup_toy():
    rng=np.random.default_rng(5)
    counts=rng.poisson(5,(100,15)).astype(float)
    xy=np.array([(i,j) for i in range(10) for j in range(10)],float)
    genes=[f'g{i}' for i in range(15)]
    axes={f'a{i}':[f'g{i+4}'] for i in range(6)}
    return counts,xy,genes,axes


def test_readout_is_not_in_any_feature_or_depth():
    counts,xy,genes,axes=setup_toy()
    axes['a0'].append('g2')
    before=feature_block(counts,genes,['g0','g1'],['g2','g3'],axes)
    counts[:,2:4]=1e6
    after=feature_block(counts,genes,['g0','g1'],['g2','g3'],axes)
    assert np.array_equal(before['X'],after['X'])
    assert not np.array_equal(before['y'],after['y'])


def test_complete_mask_prediction_ignores_hidden_expression():
    counts,xy,genes,axes=setup_toy()
    hidden=(abs(xy[:,0]-4.5)<2)&(abs(xy[:,1]-4.5)<2)
    before=masked_predict(counts,xy,genes,['g0','g1'],['g2','g3'],axes,hidden)
    counts[hidden]=1e9
    after=masked_predict(counts,xy,genes,['g0','g1'],['g2','g3'],axes,hidden)
    for k in before:
        assert np.array_equal(before[k],after[k]),k


def test_moments_solver_matches_independent_lstsq():
    rng=np.random.default_rng(10)
    X=np.c_[np.ones(100),rng.normal(size=(100,5))]
    y=X@np.arange(6)+rng.normal(size=100)
    w=rng.uniform(.1,2,100)
    beta=fit_from_moments(X.T@(w[:,None]*X),X.T@(w*y))
    oracle=np.linalg.lstsq(X*np.sqrt(w[:,None]),y*np.sqrt(w),rcond=None)[0]
    assert np.allclose(beta,oracle,atol=1e-9)


def test_batched_real_runner_matches_direct_feature_contract(monkeypatch):
    from r16.recovery import repair_pipeline as pipeline
    from r16.recovery.corrected import neighbor_values
    counts,xy,genes,axes=setup_toy()
    axes['a0'] += ['g0','g2']
    monkeypatch.setattr(pipeline,'axes',lambda:axes)
    defs={'toy':{'input_genes':['g0','g1'],'readout_genes':['g2','g3']}}
    ids,X,y,bad=next(pipeline.raw_blocks(counts,xy,genes,defs))
    direct=feature_block(counts,genes,['g0','g1'],['g2','g3'],axes)
    assert ids==['toy'] and not bad
    assert np.allclose(X[:,0,:10],direct['X'])
    assert np.allclose(y[:,0],direct['y'])
    assert np.allclose(X[:,0,10],neighbor_values(xy,direct['X'][:,-1],exclude_self=True))


def test_hex_foci_do_not_bridge_across_missing_spots():
    from r16.recovery.corrected import tls_components,signed_tls_distance
    xy=hex_coordinates(np.array([[0,0],[0,2],[0,6],[1,1]]))
    labels=np.array([1.,1.,1.,0.])
    assert sorted(map(len,tls_components(xy,labels)))==[1,2]
    d=signed_tls_distance(xy,labels)
    assert np.all(d[:3]<0) and d[3]>0


def test_padded_unmeasured_gene_is_not_a_measured_zero(monkeypatch):
    from r16.recovery import repair_pipeline as p
    counts,xy,genes,axes=setup_toy();counts[:,2]=0
    monkeypatch.setattr(p,'axes',lambda:axes)
    defs={'toy':{'input_genes':['g0','g1'],'readout_genes':['g2','g3']}}
    ids,_,_,bad=next(p.raw_blocks(counts,xy,genes,defs,unmeasured=['g2']))
    assert ids==[] and bad==[('toy','MISSING_PROGRAM_GENES')]


def test_structure_score_ignores_hidden_expression():
    from r16.recovery.repair_pipeline import reconstructed_scores
    counts,xy,genes,_=setup_toy();hidden=(abs(xy[:,0]-4.5)<2)&(abs(xy[:,1]-4.5)<2)
    first=reconstructed_scores(counts,xy,genes,[['g0','g1'],['g2']],hidden)
    counts[hidden]=1e6
    second=reconstructed_scores(counts,xy,genes,[['g0','g1'],['g2']],hidden)
    assert np.array_equal(first,second)
