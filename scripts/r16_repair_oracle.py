#!/usr/bin/env python3
"""Independent explicit-design weighted-lstsq check of repaired patient LOPO."""
import json
from collections import Counter
import numpy as np
from r16.recovery.repair_pipeline import (OUT,sections,load,definitions,axes,
                                         readtable,dump,table)
from r16.recovery.corrected import feature_block,neighbor_values


def main():
    genes=json.load(open(OUT/'genes.json'));d=definitions()['P-mhc2']
    saved={(r['cohort'],r['patient']):r for r in readtable(OUT/'six_raw_prediction.tsv') if r['program']=='P-mhc2'}
    checks=[]
    for cohort in ['ST-CRC','USZ']:
        arrays=[]
        for r in sections([cohort]):
            counts,xy,_=load(r)
            f=feature_block(counts,genes,d['input_genes'],d['readout_genes'],axes())
            X=np.c_[f['X'],neighbor_values(xy,f['X'][:,-1],exclude_self=True)]
            arrays.append((r['patient'],X,f['y']))
        ns=Counter(p for p,_,_ in arrays)
        for patient in ns:
            train=[(p,x,y) for p,x,y in arrays if p!=patient]
            weights=np.concatenate([np.full(len(y),1/(len(y)*ns[p])) for p,x,y in train])
            X=np.concatenate([x for p,x,y in train]);y=np.concatenate([y for p,x,y in train]);rootw=np.sqrt(weights)
            for model,ncols in [('M0',9),('M1',10),('M2',11)]:
                beta=np.linalg.lstsq(X[:,:ncols]*rootw[:,None],y*rootw,rcond=None)[0]
                mse=float(np.mean([np.mean((yy-xx[:,:ncols]@beta)**2) for p,xx,yy in arrays if p==patient]))
                actual=float(saved[(cohort,patient)]['mse_'+model]);error=abs(mse-actual)
                assert np.isclose(mse,actual,rtol=1e-6,atol=1e-9),(cohort,patient,model,mse,actual)
                checks.append(dict(cohort=cohort,patient=patient,model=model,explicit_weighted_lstsq_mse=mse,moment_mse=actual,absolute_error=error))
    table('explicit_lstsq_oracle.tsv',checks)
    dump('explicit_lstsq_oracle.json',dict(status='PASS',checks=len(checks),max_absolute_error=max(r['absolute_error'] for r in checks),independence='Explicit stacked spot designs and weighted np.linalg.lstsq; no moment solver used'))


if __name__=='__main__':main()
