#!/usr/bin/env python3
"""Full-program section-held-out linear prediction, train-only program selection."""
from pathlib import Path
import argparse,importlib.util,json
import numpy as np
from scipy.stats import rankdata

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'infra/tls_multiscale_gobp_20261002/reverse'
def mod(name):
 spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/f'{name}.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
M=mod('tls_multiscale_gobp');RADII=[230,500,1000]


def moments(F,y,valid):
    P,D=F.shape[1:];mean=np.zeros((2,P,D));second=np.zeros((2,P,D,D));present=np.zeros((2,P),bool)
    for cls in [0,1]:
        use=y==cls
        if not use.any():continue
        for a in range(0,P,64):
            b=min(a+64,P);X=np.ascontiguousarray(F[use,a:b].transpose(1,0,2),dtype=float)
            mean[cls,a:b]=X.mean(axis=1);second[cls,a:b]=np.matmul(X.transpose(0,2,1),X)/X.shape[1]
        present[cls]=valid
    mean*=present[:,:,None];second*=present[:,:,None,None]
    return mean,second,present


def fit_moments(mean_sum,second_sum,count,indices=None):
    if indices is None:indices=np.arange(mean_sum.shape[-1])
    means=mean_sum[:,:,indices]/np.maximum(count[:,:,None],1)
    seconds=second_sum[:,:,indices][:,:,:,indices]/np.maximum(count[:,:,None,None],1)
    mu=(means[0]+means[1])/2
    cov=(seconds[0]+seconds[1])/2-mu[:,:,None]*mu[:,None,:]
    sd=np.sqrt(np.maximum(np.diagonal(cov,axis1=1,axis2=2),1e-12))
    target=.25*(means[1]-means[0])
    Z=cov/sd[:,:,None]/sd[:,None,:];rhs=target/sd
    coef=np.linalg.solve(Z+np.eye(len(indices))[None],rhs)/sd
    intercept=.5-np.sum(mu*coef,axis=1);valid=np.all(count>=2,axis=0)
    coef[~valid]=0;intercept[~valid]=.5
    return coef,intercept,valid


def predict(F,fit,indices=None):
    coef,bias,_=fit
    if indices is not None:F=F[:,:,indices]
    return np.einsum('npd,pd->np',F,coef,optimize=True)+bias


def metrics(y,pred):
    if pred.ndim==1:pred=pred[:,None]
    n1=int(sum(y==1));n0=int(sum(y==0))
    if not n1 or not n0:return np.full(pred.shape[1],np.nan),np.full(pred.shape[1],np.nan)
    rank=rankdata(pred,axis=0,method='average')
    auc=(rank[y==1].sum(axis=0)-n1*(n1+1)/2)/(n1*n0)
    # Exact tie-aware AP, matching threshold-based average_precision_score.
    order=np.argsort(-pred,axis=0,kind='stable');s=np.take_along_axis(pred,order,axis=0);yy=y[order]
    tp=np.cumsum(yy,axis=0);ends=np.r_[s[:-1]!=s[1:],np.ones((1,s.shape[1]),bool)]
    rec=tp/n1;precision=tp/np.arange(1,len(y)+1)[:,None]
    ap=np.zeros(pred.shape[1])
    for j in range(pred.shape[1]):
        ix=np.flatnonzero(ends[:,j]);ap[j]=np.sum(np.diff(np.r_[0,rec[ix,j]])*precision[ix,j])
    return auc,ap


def stats_cache():
    directory=OUT/'moments';directory.mkdir(exist_ok=True)
    for source in sorted(OUT.glob('*_features.npz')):
        target=directory/source.name.replace('_features','_moments')
        if target.exists():continue
        z=np.load(source);F=z['features'];y=z['y'];v=z['valid_program']
        mean,second,present=moments(F,y,v)
        bm,bs,bp=moments(z['baseline'][:,None,:],y,np.array([True]))
        with target.with_suffix('.partial').open('wb') as f:np.savez_compressed(f,mean=mean,second=second,present=present,base_mean=bm,base_second=bs,base_present=bp)
        target.with_suffix('.partial').replace(target);print('moments',source.stem,flush=True)


def load_sources(radius):
    main=json.loads((OUT.parent/'contract.json').read_text())
    data={}
    for info in main['sources']:
        sid=info['section_id'];data[sid]=dict(info=info,z=np.load(OUT/f'{sid}_{radius}_features.npz'),stats=np.load(OUT/'moments'/f'{sid}_{radius}_moments.npz'))
    return data


def combine_stats(data,train,prefix=''):
    arrays=[data[s]['stats'] for s in train]
    return (sum(a[prefix+'mean'] for a in arrays),sum(a[prefix+'second'] for a in arrays),sum(a[prefix+'present'].astype(int) for a in arrays))


def joint_data(data,ids,selected):
    Xs=[];ys=[];sids=[]
    for sid in ids:
        z=data[sid]['z'];F=z['features'];B=z['baseline'];y=z['y']
        X=np.concatenate([F[:,selected,:].reshape(len(y),-1),B],axis=1) if len(y) else np.zeros((0,len(selected)*15+17))
        Xs.append(X);ys.append(y);sids.extend([sid]*len(y))
    return Xs,ys,sids


def fit_joint(Xs,ys):
    stats=[moments(X[:,None,:],y,np.array([True])) for X,y in zip(Xs,ys)]
    return fit_moments(sum(s[0] for s in stats),sum(s[1] for s in stats),sum(s[2].astype(int) for s in stats))


def train(jobs_override=None,dest_override=None,completion_path=None):
    assert (OUT/'features_complete.json').exists()
    dest=Path(dest_override) if dest_override is not None else OUT/'models';dest.mkdir(parents=True,exist_ok=True);summary=[]
    for radius in RADII:
        data=load_sources(radius);names=next(iter(data.values()))['z']['set_ids'];P=len(names)
        jobs=[('leave_section',s,[s]) for s in data]+[('leave_cohort',c,[s for s in data if data[s]['info']['cohort']==c]) for c in ['USZ','GSE175540']]
        if jobs_override is not None:jobs=jobs_override
        for split,job,test in jobs:
            target=dest/f'{radius}_{split}_{job}.npz';receipt=dest/f'{radius}_{split}_{job}.json'
            if receipt.exists():summary.append(json.loads(receipt.read_text()));continue
            train_ids=[s for s in data if s not in test];mean,second,count=combine_stats(data,train_ids)
            full=fit_moments(mean,second,count);radial=fit_moments(mean,second,count,np.arange(5));angular=fit_moments(mean,second,count,np.arange(5,15));base=fit_moments(*combine_stats(data,train_ids,'base_'))
            # Training data alone determines the combined model's GO panel.
            train_auc_sum=np.zeros(P);train_auc_n=np.zeros(P);common=np.ones(P,bool)
            for sid in train_ids:
                z=data[sid]['z'];y=z['y'];valid=z['valid_program'];common&=valid
                if len(np.unique(y))<2:continue
                pred=predict(z['features'],full);auc,_=metrics(y,pred)
                good=np.isfinite(auc)&valid&full[2];train_auc_sum[good]+=auc[good];train_auc_n[good]+=1
            train_auc=np.divide(train_auc_sum,train_auc_n,out=np.full(P,-np.inf),where=train_auc_n>0)
            candidate=np.flatnonzero(common[:1691]&full[2][:1691]&np.isfinite(train_auc[:1691]))
            selected=candidate[np.argsort(-train_auc[candidate],kind='stable')[:10]]
            Xs,ys,_=joint_data(data,train_ids,selected);joint=fit_joint(Xs,ys)
            rows=[];prediction_parts=[];label_parts=[];sid_parts=[];base_parts=[];joint_parts=[]
            for sid in test:
                z=data[sid]['z'];y=z['y'];F=z['features'];n1=int(sum(y==1));n0=int(sum(y==0));yok=n1>0 and n0>0
                preds={}
                for arm,fit,ix in [('full',full,None),('radial',radial,np.arange(5)),('angular',angular,np.arange(5,15))]:
                    pr=predict(F,fit,ix);auc,ap=metrics(y,pr);ok=z['valid_program']&fit[2]
                    for j,name in enumerate(names):
                        status='INTERNAL_SECTION_EVALUATION' if yok and ok[j] else ('NOT_TESTABLE_SINGLE_CLASS' if not yok else 'NOT_TESTABLE_COVERAGE_OR_TRAINING')
                        rows.append(dict(mask_um=radius,split=split,heldout_job=job,section_id=sid,cohort=data[sid]['info']['cohort'],arm=arm,program=str(name),n_positive=n1,n_negative=n0,auc=float(auc[j]) if status=='INTERNAL_SECTION_EVALUATION' else '',ap=float(ap[j]) if status=='INTERNAL_SECTION_EVALUATION' else '',status=status))
                    if arm=='full':prediction_parts.append(pr.astype(np.float32))
                prbase=predict(z['baseline'][:,None,:],base)[:,0]
                Xtest,_,_=joint_data(data,[sid],selected);prjoint=predict(Xtest[0][:,None,:],joint)[:,0]
                for arm,pr,fitok in [('geometry_depth',prbase,bool(base[2][0])),('train_selected_GO',prjoint,bool(joint[2][0] and z['valid_program'][selected].all() and len(selected)))]:
                    auc,ap=metrics(y,pr);status='INTERNAL_SECTION_EVALUATION' if yok and fitok else ('NOT_TESTABLE_SINGLE_CLASS' if not yok else 'NOT_TESTABLE_COVERAGE_OR_TRAINING')
                    rows.append(dict(mask_um=radius,split=split,heldout_job=job,section_id=sid,cohort=data[sid]['info']['cohort'],arm=arm,program=arm,n_positive=n1,n_negative=n0,auc=float(auc[0]) if status=='INTERNAL_SECTION_EVALUATION' else '',ap=float(ap[0]) if status=='INTERNAL_SECTION_EVALUATION' else '',status=status))
                base_parts.append(prbase);joint_parts.append(prjoint);label_parts.append(y);sid_parts.extend([sid]*len(y))
            M.tsv(dest/f'{radius}_{split}_{job}_metrics.tsv',rows)
            np.savez_compressed(target,prediction=np.concatenate(prediction_parts),y=np.concatenate(label_parts),section_ids=np.array(sid_parts),baseline=np.concatenate(base_parts),combined=np.concatenate(joint_parts),set_ids=names,selected_GO=names[selected],training_auc=train_auc)
            rr=dict(status='COMPLETE_INTERNAL_VALIDATION',mask_um=radius,split=split,job=job,train_sections=train_ids,test_sections=test,selected_GO=names[selected].tolist(),n_metric_rows=len(rows),n_predictions=sum(len(v) for v in label_parts),sha256=M.sha(target))
            M.js(receipt,rr);summary.append(rr);print('fit',radius,split,job,'selected',len(selected),flush=True)
    M.js(Path(completion_path) if completion_path is not None else OUT/'model_complete.json',dict(status='COMPLETE_INTERNAL_NOT_PATIENT_CONFIRMATION',folds=len(summary),jobs=summary))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['stats','train']);a=ap.parse_args();{'stats':stats_cache,'train':train}[a.stage]()
