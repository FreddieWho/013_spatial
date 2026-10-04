#!/usr/bin/env python3
"""Local metadata and definition audit for D-167 (no expression fitting)."""
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import h5py
import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'infra/conclusion_audit_20261002'
def module(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod

def geometry(row, col, pixels):
    old = np.column_stack([col + .5*(row % 2), row*np.sqrt(3)/2])*100
    new = np.column_stack([col/2, row*np.sqrt(3)/2])*100
    pair = cKDTree(new).query_pairs(100.01, output_type='ndarray')
    valid = np.linalg.norm(new[pair[:, 0]]-new[pair[:, 1]], axis=1)>99.99
    pair = pair[valid]
    d0 = np.linalg.norm(old[pair[:, 0]]-old[pair[:, 1]], axis=1)
    dp = np.linalg.norm(pixels[pair[:, 0]]-pixels[pair[:, 1]], axis=1)
    # Metadata pixels establish that these six directions have equal spacing.
    return {'n_neighbor_pairs':len(pair), 'old_um_unique':np.unique(np.round(d0, 3)).tolist(),
            'pixel_distance_q05_q50_q95': np.quantile(dp,[.05,.5,.95]).tolist(),
            'parity_mismatch':int(np.count_nonzero((row-col)%2))}

records=[]
for p in sorted((ROOT/'data/GEO/GSE175540/raw').glob('*tissue_positions_list.csv.gz')):
    with gzip.open(p,'rt') as f: a=list(csv.reader(f))
    a=np.array([r[1:] for r in a],float)
    rec=geometry(a[:,1],a[:,2],a[:,3:5]); rec['source']=str(p.relative_to(ROOT)); records.append(rec)
for p in sorted((ROOT/'data/other_sources/zenodo_usz_tls_visium/TLS_VISIUM_USZ/h5ad_preprocessed').glob('*.h5ad')):
    with h5py.File(p) as f:
        rec=geometry(f['obs/x_array'][:],f['obs/y_array'][:],f['obsm/spatial'][:])
    rec['source']=str(p.relative_to(ROOT)); records.append(rec)
(OUT/'geometry_metadata.json').write_text(json.dumps(records,indent=2)+'\n')
mp=json.loads((ROOT/'infra/r04/marker_proxy_combined.json').read_text())
axes=module('axes','r16/axes.py')
pool=module('pool','scripts/tls_pool_expand.py')
lab=module('lab','scripts/tls_label_field.py')
covs={c:set(mp['classes'][c]['voted_genes']) for c in ['T','Mye','Epi','Stromal']}
covs['Plasma']=set(axes.PLASMA_GENES)
covs['detection']=set(pool.signature_symbols())
with (ROOT/'infra/gobp_halo_20261001/target_sets_v7.tsv').open() as f:
    readouts=[(r['set_id'],r['genes'].split(',')) for r in csv.DictReader(f,delimiter='\t')]
hm=lab.read_gmt(lab.HALLMARK)
readouts += [(k,hm[k]) for k in ['HALLMARK_INTERFERON_GAMMA_RESPONSE','HALLMARK_INTERFERON_ALPHA_RESPONSE','HALLMARK_G2M_CHECKPOINT','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']]
readouts += [('GENE_'+g.replace('-','_'),[g]) for g in lab.SINGLE_GENES]
rows=[]
for name,genes in readouts:
    for cov,markers in covs.items():
        shared=sorted(set(genes)&markers)
        rows.append(dict(set_id=name,covariate=cov,n_readout=len(set(genes)),n_overlap=len(shared),genes=','.join(shared)))
with (OUT/'covariate_overlap.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
print('geometry sources',len(records),'readouts',len(readouts),'readouts overlapping composition',len({r['set_id'] for r in rows if r['n_overlap'] and r['covariate']!='detection'}))
print('old distance values',records[0]['old_um_unique'],'pixel spacing',records[0]['pixel_distance_q05_q50_q95'])
