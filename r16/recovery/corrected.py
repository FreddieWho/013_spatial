"""Audited inference primitives (D-140). No I/O, implicit fitting, or GT inputs.

Program features exclude input/readout genes from both composition and QC.
Mask prediction constructs every molecular feature from visible rows only.
"""
from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.interpolate import BSpline
from scipy.linalg import null_space
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import connected_components


class NotTestable(ValueError):
    pass


class NaturalSpline:
    def __init__(self, x):
        x=np.asarray(x,float)
        if not np.isfinite(x).all() or np.unique(x).size < 4:
            raise NotTestable('INSUFFICIENT_DISTINCT_DISTANCES')
        q=np.quantile(x,[0,1/3,2/3,1])
        if np.min(np.diff(q)) <= 1e-10:
            q=np.linspace(x.min(),x.max(),4)
        self.bounds=q[[0,-1]]
        t=np.r_[[q[0]]*4,q[1:3],[q[-1]]*4]
        self.bs=BSpline(t,np.eye(6),3)
        self.transform=null_space(self.bs.derivative(2)(self.bounds))

    def __call__(self,x):
        x=np.asarray(x,float)
        clipped=np.clip(x,*self.bounds)
        # Natural splines continue linearly outside the training support.
        return (self.bs(clipped)+(x-clipped)[...,None]*self.bs.derivative(1)(clipped))@self.transform

    def second_derivative_at_boundaries(self):
        return self.bs.derivative(2)(self.bounds)@self.transform


def ns_basis(x,df=3):
    if df!=3:
        raise ValueError('only the frozen df=3 basis is supported')
    return NaturalSpline(x)(x)


def classify_curve(curve, observed_sd, floor=.2):
    v=np.asarray(curve,float)
    if len(v)<4 or not np.isfinite(v).all() or observed_sd<=1e-12:
        return 'NOT_TESTABLE'
    if np.ptp(v)/observed_sd < floor:
        return 'flat'
    diff=np.diff(v); tol=max(np.ptp(v),observed_sd)*1e-7
    if np.mean(diff>=-tol)>=.95 and (v[-1]-v[0])/observed_sd>=floor:
        return 'low-to-high'
    if np.mean(diff<=tol)>=.95 and (v[0]-v[-1])/observed_sd>=floor:
        return 'high-to-low'
    return 'nonmonotonic'


def hex_coordinates(array_coords):
    xy=np.asarray(array_coords,float)
    if xy.ndim!=2 or xy.shape[1]!=2 or not np.isfinite(xy).all():
        raise ValueError('array coordinates must be finite n by 2')
    if not np.allclose(xy,np.round(xy)):
        raise ValueError('array coordinates must be explicit integer indices')
    return xy*np.array([np.sqrt(3)/2,.5])


def tls_components(coords, labels):
    ix=np.flatnonzero(np.asarray(labels)==1)
    if not len(ix):
        return []
    pairs=cKDTree(np.asarray(coords)[ix]).query_pairs(1.001,output_type='ndarray')
    w=sparse.csr_matrix((np.ones(len(pairs)),(pairs[:,0],pairs[:,1])),shape=(len(ix),len(ix)))
    nc,lab=connected_components(w.maximum(w.T),directed=False)
    return [ix[lab==k] for k in range(nc)]


def signed_tls_distance(coords,labels):
    inside=np.asarray(labels)==1
    if not inside.any() or inside.all():
        raise NotTestable('NO_TLS_OR_NO_OUTSIDE')
    dist=cKDTree(coords[inside]).query(coords)[0]
    dist[inside]=-cKDTree(coords[~inside]).query(coords[inside])[0]
    return dist


def classify_window(labels, min_known=.8, pos_fraction=.1):
    y=np.asarray(labels,float)
    known=np.isfinite(y)
    if not len(y) or known.mean()<min_known:
        return 'UNKNOWN'
    if np.mean(y[known]==1)>=pos_fraction:
        return 'POSITIVE'
    if known.all() and np.all(y==0):
        return 'NEGATIVE'
    return 'PARTIAL' if np.any(y[known]==1) else 'UNKNOWN'


def feature_block(counts,genes,input_genes,readout_genes,axes,mode='raw'):
    counts=np.asarray(counts,float)
    gx={g:i for i,g in enumerate(genes)}
    if set(input_genes)&set(readout_genes):
        raise ValueError('input/readout overlap')
    if not set(input_genes+readout_genes)<=set(gx):
        raise NotTestable('MISSING_PROGRAM_GENES')
    exclude=set(input_genes+readout_genes)
    bg=[i for i,g in enumerate(genes) if g not in exclude]
    if not bg:
        raise NotTestable('EMPTY_QC_BACKGROUND')
    lib=counts[:,bg].sum(1)
    detected=(counts[:,bg]>0).sum(1)
    scale=(1e4/np.maximum(lib,1))[:,None] if mode=='depth_normalized' else 1
    if mode not in {'raw','depth_normalized'}:
        raise ValueError(mode)
    def score(gs):
        cols=[gx[g] for g in gs if g in gx]
        if not cols:
            raise NotTestable('EMPTY_COMPOSITION_AXIS')
        return np.log1p(counts[:,cols]*scale).mean(1)
    clean={a:[g for g in gs if g in gx and g not in exclude] for a,gs in axes.items()}
    cols=[np.ones(len(counts)),np.log1p(lib),np.log1p(detected)]
    cols.extend(score(gs) for gs in clean.values())
    cols.append(score(input_genes))
    return {'X':np.column_stack(cols),'y':score(readout_genes),
            'clean_axes':clean,'background_library':lib}


def fit_from_moments(gram,rhs):
    """Batched weighted least squares with training-only feature scaling."""
    g=np.asarray(gram,float); b=np.asarray(rhs,float)
    mass=g[...,0,0]
    mu=g[...,0,1:]/mass[...,None]
    cov=g[...,1:,1:]-mass[...,None,None]*mu[...,None,:]*mu[..., :,None]
    sd=np.sqrt(np.maximum(np.diagonal(cov,axis1=-2,axis2=-1)/mass[...,None],0))
    sd=np.where(sd>1e-8,sd,1)
    stdgram=cov/(sd[...,None,:]*sd[..., :,None])
    cross=(b[...,1:]-mu*b[...,0,None])/sd
    coef=np.einsum('...ij,...j->...i',np.linalg.pinv(stdgram,rcond=1e-10,hermitian=True),cross)/sd
    intercept=b[...,0]/mass-np.sum(mu*coef,axis=-1)
    return np.concatenate([intercept[...,None],coef],axis=-1)


def fit_linear(X,y,weights=None):
    X=np.asarray(X,float); y=np.asarray(y,float)
    w=np.ones(len(y)) if weights is None else np.asarray(weights,float)
    return fit_from_moments(X.T@(w[:,None]*X),X.T@(w*y))


def neighbor_values(coords,values,query=None,k=8,exclude_self=False):
    coords=np.asarray(coords,float)
    query=coords if query is None else np.asarray(query,float)
    if exclude_self and len(coords)<2:
        raise NotTestable('INSUFFICIENT_VISIBLE_POINTS')
    use=min(k+(1 if exclude_self else 0),len(coords))
    d,ix=cKDTree(coords).query(query,k=use)
    if use==1:
        d,ix=d[:,None],ix[:,None]
    if exclude_self:
        d,ix=d[:,1:],ix[:,1:]
    w=1/np.maximum(d,1e-6)
    w/=w.sum(1,keepdims=True)
    vals=np.asarray(values)[ix]
    return np.sum(vals*w[(...,)+(None,)*(vals.ndim-2)],axis=1)


def masked_predict(counts,coords,genes,input_genes,readout_genes,axes,hidden,mode='raw'):
    """No hidden molecular values are accessed, including for normalization."""
    hidden=np.asarray(hidden,bool); vis=~hidden
    if hidden.sum()==0 or vis.sum()<20:
        raise NotTestable('INSUFFICIENT_VISIBLE_OR_HIDDEN')
    visible=np.asarray(counts)[vis]
    f=feature_block(visible,genes,input_genes,readout_genes,axes,mode)
    x=f['X']; y=f['y']; cv=np.asarray(coords)[vis]; ch=np.asarray(coords)[hidden]
    nearest=cKDTree(cv).query(ch)[1]
    xte=x[nearest]
    ntr=neighbor_values(cv,x[:,-1],exclude_self=True)
    nte=neighbor_values(cv,x[:,-1],ch)
    designs=[(x[:,:-1],xte[:,:-1]),(x,xte),(np.c_[x,ntr],np.c_[xte,nte])]
    result={f'M{i}':b@fit_linear(a,y) for i,(a,b) in enumerate(designs)}
    result['NN']=y[nearest]
    result['KNN8']=neighbor_values(cv,y,ch)
    return result


def patient_interval(values,draws=5000,seed=20260921):
    v=np.asarray(values,float); v=v[np.isfinite(v)]
    if len(v)<2:
        return [None,None]
    rng=np.random.default_rng(seed)
    boot=np.median(v[rng.integers(len(v),size=(draws,len(v)))],axis=1)
    return np.quantile(boot,[.025,.975]).tolist()
