import importlib.util
from pathlib import Path
import numpy as np
from scipy import sparse

ROOT=Path(__file__).resolve().parents[1]
def module(name):
 spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def test_tumor_boundary_never_uses_unknown_or_tls_as_non_tumor():
 m=module('tls_multiscale_gobp')
 xy=np.array([[0,0],[100,0],[200,0],[300,0]],float)
 _,_,_,edge=m.build_endpoints(xy,np.array(['TUM','UNASSIGNED','TLS','TUM']))
 assert not edge.any()
 _,_,_,edge=m.build_endpoints(xy,np.array(['TUM','NOR','TLS','TUM']))
 assert edge.tolist()==[True,False,False,False]


def test_disconnected_graph_cannot_create_a_bridge():
 m=module('tls_multiscale_forward')
 g=sparse.csr_matrix([[0,100,0],[100,0,0],[0,0,0]])
 from scipy.sparse.csgraph import dijkstra
 _,pred=dijkstra(g,directed=False,indices=0,return_predecessors=True)
 assert m.trace_path(pred,0,2) is None


def test_corridor_contrast_constant_zero_and_correct_injected_bridge():
 m=module('tls_multiscale_forward')
 xx,yy=np.meshgrid(np.arange(0,2001,50),np.arange(-500,501,50));xy=np.c_[xx.ravel(),yy.ravel()]
 along,side,length=m.path_coordinates(xy,np.array([0,0]),np.array([2000,0]))
 visible=(np.linalg.norm(xy,axis=1)>130)&(np.linalg.norm(xy-[2000,0],axis=1)>130)
 result,status,_=m.corridor_weights(along,side,length,visible,100)
 assert status=='ELIGIBLE';W,V,center=result
 np.testing.assert_allclose(W.sum(axis=1),0,atol=1e-14)
 np.testing.assert_allclose(W@center.astype(float),1,atol=1e-14)
 assert np.all(W[:,~visible]==0)
 assert np.all(V[:,~visible]==0)


def test_endpoint_adjustment_matches_explicit_fit_and_ignores_hidden():
 m=module('tls_multiscale_forward');rng=np.random.default_rng(11)
 Z=np.c_[np.ones(100),rng.normal(size=(100,4))];Y=rng.normal(size=(100,8));train=np.arange(100)<70
 w=np.zeros(100);w[70:80]=.1;w[80:90]=-.1
 a=m.endpoint_adjust(w,Z,train)
 beta=np.linalg.lstsq(Z[train].T@Z[train]+np.eye(5)*1e-6,Z[train].T@Y[train],rcond=None)[0]
 np.testing.assert_allclose(a@Y,w@(Y-Z@beta),atol=1e-12)
 assert np.all(a[90:]==0)
 Y[90:]=1e9
 np.testing.assert_allclose(a@Y,w@(Y-Z@beta),atol=1e-12)


def test_short_scale_not_coerced_into_millimeter_bucket():
 m=module('tls_multiscale_forward')
 assert m.scale_id(250)==0
 assert m.scale_id(500)==1
 assert m.scale_id(1000)==2


def test_reverse_features_ignore_all_hidden_expression_and_depth():
 m=module('tls_multiscale_reverse')
 xx,yy=np.meshgrid(np.arange(-1500,1501,100),np.arange(-1500,1501,100));xy=np.c_[xx.ravel(),yy.ravel()]
 centre=int(np.argmin(np.sum(xy*xy,axis=1)));q=[{'spot_index':centre}]
 rng=np.random.default_rng(18);S=rng.normal(size=(len(xy),4)).astype(np.float32);lib=rng.integers(100,10000,size=len(xy)).astype(float)
 a,b=m.pool_query_features(xy,S,lib,q,500)
 hidden=np.linalg.norm(xy-xy[centre],axis=1)<=500
 S[hidden]=rng.uniform(1e4,1e9,size=(hidden.sum(),4));lib[hidden]=1e12
 aa,bb=m.pool_query_features(xy,S,lib,q,500)
 np.testing.assert_array_equal(a,aa);np.testing.assert_array_equal(b,bb)


def test_reverse_moment_fit_matches_explicit_weighted_ridge():
 m=module('tls_multiscale_reverse_models');rng=np.random.default_rng(29)
 Fs=[rng.normal(size=(20,3,5)),rng.normal(size=(25,3,5))]
 ys=[np.r_[np.ones(6),np.zeros(14)],np.r_[np.ones(9),np.zeros(16)]]
 st=[m.moments(F,y,np.ones(3,bool)) for F,y in zip(Fs,ys)]
 fit=m.fit_moments(sum(t[0] for t in st),sum(t[1] for t in st),sum(t[2].astype(int) for t in st))
 weights=[np.where(y==1,.5/2/sum(y==1),.5/2/sum(y==0)) for y in ys]
 F=np.concatenate(Fs);y=np.concatenate(ys);w=np.concatenate(weights)
 for j in range(3):
  X=F[:,j];mu=w@X;sd=np.sqrt(w@((X-mu)**2));Z=(X-mu)/sd
  beta=np.linalg.solve(Z.T@(w[:,None]*Z)+np.eye(5),Z.T@(w*(y-.5)))
  np.testing.assert_allclose(m.predict(F,fit)[:,j],.5+Z@beta,atol=1e-12)


def test_auc_and_ap_are_tie_aware():
 m=module('tls_multiscale_reverse_models')
 from sklearn.metrics import roc_auc_score,average_precision_score
 y=np.array([0,1,1,0,0,1]);p=np.array([[0,1],[0,1],[.2,1],[.2,1],[.5,1],[.5,1]],float)
 auc,ap=m.metrics(y,p)
 for j in range(2):
  np.testing.assert_allclose(auc[j],roc_auc_score(y,p[:,j]))
  np.testing.assert_allclose(ap[j],average_precision_score(y,p[:,j]))

def test_reference_low_rank_residual_matches_full_operator():
 m=module('tls_multiscale_controls');rng=np.random.default_rng(37)
 Z=np.c_[np.ones(200),rng.normal(size=(200,4))];Y=rng.normal(size=(200,9))
 visible=np.arange(200)>=20;centre=(np.arange(200)>=100)&(np.arange(200)<130)
 w=np.zeros(200);w[centre]=1/30;w[40:70]=-1/30
 coef,local=m.residual_parts(w,Z,visible,centre)
 exact=m.F.endpoint_adjust(w,Z,visible&~centre)
 low=(w+local)@Y-coef@(Z[visible].T@Y[visible])
 np.testing.assert_allclose(low,exact@Y,atol=1e-12)
 Y[~visible]=1e12
 np.testing.assert_allclose((w+local)@Y-coef@(Z[visible].T@Y[visible]),low,atol=1e-12)

def test_fixed_span_zero_rotation_reproduces_true_corridor():
 m=module('tls_multiscale_controls')
 xx,yy=np.meshgrid(np.arange(-500,2501,50),np.arange(-700,701,50));xy=np.c_[xx.ravel(),yy.ravel()]
 points=np.array([[0,0],[2000,0]],float);visible=np.ones(len(xy),bool)
 along,side,length=m.F.path_coordinates(xy,*points)
 true,status,counts=m.F.corridor_weights(along,side,length,visible,100)
 assert status=='ELIGIBLE'
 W,_,centre=true;span=(along[centre].min(),along[centre].max());labels=np.full(len(xy),'NOR')
 ref,code,_=m.reference_geometry(xy,labels,points,length,visible,100,counts,m.composition(W,labels),span,0)
 assert code==0
 np.testing.assert_array_equal(ref[0],W)
 np.testing.assert_array_equal(ref[1],centre)

def test_complete_search_rank_is_conservative_for_ties_and_rotation_equivariant():
 m=module('tls_multiscale_bridge_calibration');rng=np.random.default_rng(51)
 E=rng.normal(size=(7,40,9));positive=np.where(E>0,5,0);negative=np.where(E<0,5,0)
 T,_=m.selected_statistics(E,positive,negative)
 TT,_=m.selected_statistics(np.roll(E,7,axis=1),np.roll(positive,7,axis=1),np.roll(negative,7,axis=1))
 np.testing.assert_allclose(TT,np.roll(T,7,axis=0),atol=1e-12)
 np.testing.assert_array_equal(m.diagnostic_rank(np.zeros((40,9))),np.ones(9))
 planted=np.ones((1,40,1));planted[0,0]=8
 score,_=m.selected_statistics(planted,np.full_like(planted,5),np.zeros_like(planted))
 np.testing.assert_allclose(m.diagnostic_rank(score),1/40)
 assert m.bound(0,400,False)==0 and m.bound(400,400,True)==1

def test_reference_quadrant_support_does_not_accept_concentrated_angles():
 m=module('tls_multiscale_controls_summary')
 assert m.reference_support(np.arange(4.5,180,4.5))[0]
 assert not m.reference_support(np.linspace(5,40,25))[0]

def test_halo_edge_references_do_not_read_hidden_means_or_variances():
 m=module('tls_multiscale_halo_edge');rng=np.random.default_rng(67)
 W=np.zeros((4,200));W[0,20:80]=1/60;W[1,80:140]=1/60;W[2,40:100]=1/60;W[3,100:160]=1/60
 V=np.zeros((4,200));V[:2,20:]=1;V[2:,40:]=1
 op=dict(status=np.zeros((2,2),np.uint8),W=sparse.csr_matrix(W),V=sparse.csr_matrix(V),cases=[dict(a=0,b=1,mask_um=130)])
 S=rng.normal(size=(200,3));a=m.evaluate(S,op);S[:20]=1e10;b=m.evaluate(S,op)
 for x,y in zip(a,b):np.testing.assert_array_equal(x,y)
 S[20:40]=1e9;c=m.evaluate(S,op)
 np.testing.assert_array_equal(b[0][1],c[0][1]);np.testing.assert_array_equal(b[1][:,1],c[1][:,1])

def test_distance_contrasts_do_not_jump_across_unsupported_bins():
 m=module('tls_multiscale_halo_edge')
 rows=[dict(mask_um='130',lo_um=str(lo),hi_um=str(hi),status=status) for lo,hi,status in [(130,230,'EVALUABLE'),(230,400,'NOT_TESTABLE_LT30'),(400,800,'EVALUABLE')]]
 cases=m.contrasts(rows,'halo');assert [(c['a'],c['b']) for c in cases]==[(0,1),(1,2)]
 assert all(c['status']=='NOT_TESTABLE_BIN_SUPPORT' for c in cases)
 shifts=m.translations();assert shifts.shape==(85,2) and len(np.unique(shifts,axis=0))==85
