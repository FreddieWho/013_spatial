import importlib.util
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
def module():
    spec=importlib.util.spec_from_file_location('patient_validation',ROOT/'scripts/tls_patient_validation.py')
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def test_group_split_excludes_both_preparations():
    m=module();mapping={'ffpe':'p1','frozen':'p1','other':'p2'}
    train,test=m.split_group(mapping,'p1')
    assert train==['other'] and set(test)=={'ffpe','frozen'}

def test_override_fit_does_not_read_heldout_labels_for_training(tmp_path):
    m=module();r=m.R;r.OUT=tmp_path;(tmp_path/'features_complete.json').write_text('{}');r.RADII=[500]
    rng=np.random.default_rng(42);data={}
    for sid in ['a','b','c','d','ffpe','frozen']:
        F=rng.normal(size=(12,3,15));B=rng.normal(size=(12,17));y=np.arange(12)%2;valid=np.ones(3,bool)
        mean,second,present=r.moments(F,y,valid);bm,bs,bp=r.moments(B[:,None,:],y,np.ones(1,bool))
        data[sid]={'info':{'cohort':'GSE175540'},'z':{'features':F,'baseline':B,'y':y,'valid_program':valid,'set_ids':np.array(['GO1','GO2','B_AXIS'])},'stats':{'mean':mean,'second':second,'present':present,'base_mean':bm,'base_second':bs,'base_present':bp}}
    r.load_sources=lambda radius:data
    for iteration in range(2):
        dest=tmp_path/f'run{iteration}'
        r.train(jobs_override=[('leave_patient','p1',['ffpe','frozen'])],dest_override=dest,completion_path=tmp_path/f'done{iteration}.json')
        if iteration==0:
            for sid in ['ffpe','frozen']:
                data[sid]['z']['y']=1-data[sid]['z']['y']
                # A cached test moment contamination would change training too.
                data[sid]['stats']['mean'][:]=1e6
    a=np.load(tmp_path/'run0/500_leave_patient_p1.npz');b=np.load(tmp_path/'run1/500_leave_patient_p1.npz')
    for key in ['prediction','baseline','combined','training_auc','selected_GO']:np.testing.assert_array_equal(a[key],b[key])
    assert not (tmp_path/'model_complete.json').exists()
