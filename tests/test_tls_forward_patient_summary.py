import importlib.util
from pathlib import Path
import numpy as np
p=Path(__file__).resolve().parents[1]/'scripts/tls_forward_patient_summary.py'
s=importlib.util.spec_from_file_location('forward_group',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
def test_patient_weighting_preserves_missing_and_never_merges_scales():
    meta=[{'section_id':sid,'cohort':'GSE175540','scale':scale} for sid,scale in [('a','short'),('b','short'),('c','short'),('a','long')]]
    values={'effect':np.array([[2.,np.nan],[4.,np.nan],[10.,7.],[20.,8.]])}
    u,v,n,s,a=m.aggregate(meta,values,['scale'],{'a':'p1','b':'p1','c':'p2'})
    assert len(u)==3 and len(s)==2
    np.testing.assert_equal(a['effect_median'][0],[6.5,7.])
    np.testing.assert_equal(a['effect_old_section_median'][0],[4.,7.])
    np.testing.assert_equal(a['effect_n_units'][0],[2,1])
    np.testing.assert_equal(a['effect_n_complete_units'][1],[0,0])
    np.testing.assert_equal(a['effect_median'][1],[20.,8.])
