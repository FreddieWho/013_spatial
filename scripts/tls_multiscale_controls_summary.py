#!/usr/bin/env python3
"""Every-program native-reference contrasts; empirical quantiles are not p-values."""
from pathlib import Path
import argparse,importlib.util,json,warnings
import numpy as np

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'infra/tls_multiscale_gobp_20261002'
spec=importlib.util.spec_from_file_location('m',ROOT/'scripts/tls_multiscale_gobp.py')
M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M)


def reference_support(angles):
    counts=np.bincount(np.minimum((np.asarray(angles)/45).astype(int),3),minlength=4)
    return len(angles)>=19 and np.all(counts>=2),counts


def summarize_source(sid):
    source=OUT/'controls'/sid
    if (source/'summary_complete.json').exists():return
    assert (source/'scoring_complete.json').exists()
    info=json.loads((source/'geometry_complete.json').read_text());n=info['rows'];P=1742
    if not n:
        M.js(source/'summary_complete.json',dict(status='NO_EVALUABLE_GEOMETRY',rows=0));return
    rows=M.tab(OUT/'forward'/sid/'effect_rows.tsv')
    references={}
    for path in sorted(source.glob('chunk_*_rows.tsv')):
        for r in M.tab(path):references.setdefault(int(r['effect_row']),[]).append(r)
    actual=np.load(OUT/'forward'/sid/'endpoint_residual.npy',mmap_mode='r')
    samples=np.load(source/'adjusted.npy',mmap_mode='r') if info['matched_references'] else None
    keys=['native_z','native_quantile','reference_mean','reference_sd','tissue_z','tissue_quantile']
    result={key:np.lib.format.open_memmap(source/f'{key}.npy',mode='w+',dtype='float32',shape=(n,P)) for key in keys}
    for value in result.values():value[:]=np.nan
    qualification=[]
    for row in rows:
        i=int(row['effect_row']);refs=references.get(i,[])
        ids=[int(r['reference_id']) for r in refs];angles=[float(r['angle_degrees']) for r in refs]
        tissue=[r for r in refs if r['tissue_matched']=='True']
        usable,quarters=reference_support(angles);tusable,tquarters=reference_support([float(r['angle_degrees']) for r in tissue])
        shape_ok=float(row['length_um'])>=2*float(row['width_um'])
        qualification.append(dict(effect_row=i,n_references=len(refs),n_tissue_references=len(tissue),
            quadrant_counts=json.dumps(quarters.tolist()),tissue_quadrant_counts=json.dumps(tquarters.tolist()),
            reference_qualified=bool(usable),tissue_reference_qualified=bool(tusable),
            bridge_aspect_qualified=shape_ok,status='DIRECTION_REFERENCE_AVAILABLE' if usable else 'INSUFFICIENT_DIRECTION_REFERENCES'))
        if len(refs)<2:continue
        for prefix,ix in [('',ids),('tissue_',[int(r['reference_id']) for r in tissue])]:
            if len(ix)<2:continue
            R=np.asarray(samples[ix],dtype=float);v=np.asarray(actual[i],dtype=float)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning)
                mu=np.nanmean(R,axis=0);sd=np.nanstd(R,axis=0,ddof=1)
            valid=np.isfinite(v)&np.isfinite(mu)&(sd>1e-12)&np.all(np.isfinite(R),axis=0)
            z=np.divide(v-mu,sd,out=np.full(P,np.nan),where=valid)
            q=np.where(valid,(np.sum(R<v,axis=0)+.5*np.sum(R==v,axis=0))/len(ix),np.nan)
            result['tissue_z' if prefix else 'native_z'][i]=z
            result['tissue_quantile' if prefix else 'native_quantile'][i]=q
            if not prefix:
                result['reference_mean'][i]=mu;result['reference_sd'][i]=sd
    for value in result.values():value.flush()
    M.tsv(source/'reference_qualification.tsv',qualification)
    M.js(source/'summary_complete.json',dict(status='FULL_PANEL_NATIVE_REFERENCES_NOT_CALIBRATED',rows=n,
        reference_qualified=sum(r['reference_qualified'] for r in qualification),
        bridge_reference_qualified=sum(r['reference_qualified'] and r['bridge_aspect_qualified'] for r in qualification),
        tissue_reference_qualified=sum(r['tissue_reference_qualified'] for r in qualification),n_gobp=1691,n_controls=51))
    print('native summary',sid,n,flush=True)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--section',default='');args=ap.parse_args()
    sources=json.loads((OUT/'contract.json').read_text())['sources']
    for r in sources:
        if not args.section or r['section_id']==args.section:summarize_source(r['section_id'])
    if not args.section:
        reports=[dict(section_id=r['section_id'],**json.loads((OUT/'controls'/r['section_id']/'summary_complete.json').read_text())) for r in sources]
        M.tsv(OUT/'controls/source_summary.tsv',reports)
        M.js(OUT/'controls/summary_complete.json',dict(status='ALL26_FULL_PANEL_NATIVE_REFERENCES_NOT_CALIBRATED',sources=26,n_gobp=1691,n_controls=51))


if __name__=='__main__':main()
