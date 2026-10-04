#!/usr/bin/env python3
"""Complete array summaries; no selected subset or uncalibrated discovery declaration."""
from pathlib import Path
import importlib.util,json,warnings
import numpy as np

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'infra/tls_multiscale_gobp_20261002'
spec=importlib.util.spec_from_file_location('m',ROOT/'scripts/tls_multiscale_gobp.py');M=importlib.util.module_from_spec(spec);spec.loader.exec_module(M)

def main():
    dest=OUT/'summary';dest.mkdir(exist_ok=True)
    contract=json.loads((OUT/'contract.json').read_text());allmeta=[];allraw=[];alladj=[];allcont=[];allcounts=[];support=[]
    for info in contract['sources']:
        sid=info['section_id'];d=OUT/'forward'/sid;receipt=json.loads((d/'scoring_complete.json').read_text())
        support.append(dict(section_id=sid,cohort=info['cohort'],**receipt))
        if not receipt['evaluable_bridge_geometries']:continue
        assert (d/'continuity_complete.json').exists()
        rows=M.tab(d/'effect_rows.tsv');pairs=M.tab(d/'pairs.tsv')
        pairmap={(r['pair_id'],r['mode']):r for r in pairs}
        raw=np.load(d/'raw_effect.npy',mmap_mode='r');adj=np.load(d/'endpoint_residual.npy',mmap_mode='r');positive=np.load(d/'positive_segments.npy',mmap_mode='r')
        groups={}
        for r in rows:
            p=pairmap[(r['pair_id'],r['mode'])]
            key=tuple(r[k] for k in ['family','mode','mask_um','width_um','length_stratum'])+(p['target_family'],)
            groups.setdefault(key,[]).append((int(r['effect_row']),p['source_endpoint']))
        for key,items in sorted(groups.items()):
            sources=sorted({s for _,s in items});R=[];A=[];C=[];N=[]
            with warnings.catch_warnings():
                warnings.simplefilter('ignore',RuntimeWarning)
                for s in sources:
                    ids=[i for i,ss in items if ss==s];rv=np.asarray(raw[ids]);av=np.asarray(adj[ids]);pc=np.asarray(positive[ids])
                    R.append(np.nanmedian(rv,axis=0));A.append(np.nanmedian(av,axis=0))
                    valid=np.isfinite(rv);den=valid.sum(axis=0);N.append(den)
                    C.append(np.divide(np.sum(valid&(pc>=4),axis=0),den,out=np.full(raw.shape[1],np.nan),where=den>0))
                allraw.append(np.nanmedian(R,axis=0));alladj.append(np.nanmedian(A,axis=0));allcont.append(np.nanmedian(C,axis=0));allcounts.append(np.sum(N,axis=0))
            allmeta.append(dict(section_id=sid,cohort=info['cohort'],family=key[0],mode=key[1],mask_um=key[2],width_um=key[3],length_stratum=key[4],target_family=key[5],n_source_components=len(sources),n_geometry_rows=len(items),status='DESCRIPTIVE_NOT_CALIBRATED'))
        print('summarized',sid,len(groups),'complete geometry groups',flush=True)
    names=np.load(OUT/'scores'/f"{contract['sources'][0]['section_id']}.npz")['set_ids']
    np.savez_compressed(dest/'bridge_all_programs.npz',raw_median=np.array(allraw,dtype=np.float32),endpoint_residual_median=np.array(alladj,dtype=np.float32),positive_continuity_fraction=np.array(allcont,dtype=np.float32),n_valid_geometry=np.array(allcounts,dtype=np.int32),set_ids=names)
    M.tsv(dest/'bridge_groups.tsv',allmeta);M.tsv(dest/'source_completion.tsv',support)
    M.js(dest/'receipt.json',dict(status='FULL_FORWARD_DESCRIPTIVE_SUMMARY',sources=len(support),group_rows=len(allmeta),gobp=1691,controls=51,formal_p='NOT_CALIBRATED',inverse='SEPARATE_PENDING'))
if __name__=='__main__':main()
