"""Bounded CPU diagnostic: eight USZ sections, five existing program readouts.
Replay T4 numbers; compare nuisance definitions without changing production.
Run with the same environment as reproduce.py. No model training or new data.
"""
import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from r16 import axes as A
from r16 import census as C
from r16.section_io import load_section_symbols, usz_tls_labels

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
REC = ROOT/'infra/r16/recovery_20260918'
defs = json.load(open(REC/'programs/program_definitions.json'))
proxy = json.load(open(ROOT/'infra/r04/marker_proxy_combined.json'))
axes = {k:v['voted_genes'] for k,v in proxy['classes'].items() if k not in A.RETIRED_AXES}
axes['Plasma'] = A.PLASMA_GENES
records, sections, diagnostics = [], [], []
go = json.load(open(REC/'programs_go/definitions.json'))
go_rep = json.load(open(REC/'programs_go/id_map.json'))
needed_go = set().union(*(set(go[g]['input_genes']+go[g]['readout_genes']) for g in go_rep))
for row in json.load(open(ROOT/'infra/r04/role_manifests/external_validation_manifest.json'))['rows']:
    sec = load_section_symbols(row)
    mat = sec.counts.tocsr()
    genes = list(sec.gene_id)
    gx = {}
    for i,g in enumerate(genes):
        gx.setdefault(g,i)
    lib = np.asarray(mat.sum(axis=1)).ravel()
    ng = np.asarray((mat>0).sum(axis=1)).ravel()
    Q = np.c_[np.ones(len(lib)),np.log10(lib+1),(ng-ng.mean())/(ng.std()+1e-9)]
    needed = sorted(set().union(*map(set,axes.values())))
    aa = np.log1p(mat[:,[gx[g] for g in needed if g in gx]].toarray().astype(np.float32))
    ag = {g:i for i,g in enumerate([g for g in needed if g in gx])}
    scale = aa.std(0); scale[scale<1e-8]=1
    zz = (aa-aa.mean(0))/scale
    C6 = np.column_stack([A.axis_scores(zz,ag,gs)[0] for gs in axes.values()])
    def score(gs):
        return np.log1p(mat[:,[gx[g] for g in gs if g in gx]].toarray().astype(float)).mean(1)
    vals = {p:(score(d['input_genes']),score(d['readout_genes'])) for p,d in defs.items() if p!='P-plasma'}
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=9).fit(sec.coords)
    _, ix = nn.kneighbors(sec.coords)
    items = {}
    for p,(inp,ro) in vals.items():
        excluded = set(defs[p]['input_genes']+defs[p]['readout_genes'])
        clean_axes = {k:[g for g in gs if g not in excluded] for k,gs in axes.items()}
        cleaned = np.column_stack([A.axis_scores(zz,ag,gs)[0] for gs in clean_axes.values()])
        raw = np.column_stack([score(gs) for gs in axes.values()])
        raw_clean = np.column_stack([score(gs) if gs else np.zeros(len(lib)) for gs in clean_axes.values()])
        items[p] = {'y':ro,'inp':inp,'neigh':inp[ix[:,1:]].mean(1),
                    'legacy':np.c_[Q,C6], 'exclude_program_input_readout':np.c_[Q,cleaned],
                    'mean_log_axes_sensitivity':np.c_[Q,raw],
                    'mean_log_axes_exclude_program':np.c_[Q,raw_clean]}
    sections.append({'patient':sec.patient_id,'section':sec.section_id,'items':items})
    lab = usz_tls_labels(sec.section_id,tuple(sec.barcode),ROOT)
    tls = np.flatnonzero(lab==1)
    old_w = C.spatial_weights(sec.coords,C.KNN_SPATIAL)
    old_foci = connected_components(old_w[tls][:,tls],directed=False)[0]
    # Correct physical embedding for a Visium staggered array index grid.
    xy = np.asarray(sec.coords,dtype=float)*np.array([np.sqrt(3),1.0])
    pairs = cKDTree(xy).query_pairs(2.01,output_type='ndarray')
    ew = sparse.csr_matrix((np.ones(len(pairs)),(pairs[:,0],pairs[:,1])),shape=old_w.shape)
    ew = ew.maximum(ew.T)
    hex_foci = connected_components(ew[tls][:,tls],directed=False)[0]
    diagnostics.append({'section':sec.section_id,'patient':sec.patient_id,'spots':len(lab),
                        'tls':len(tls),'unknown':int(np.isnan(lab).sum()),
                        'legacy_knn_foci':int(old_foci),'hex_neighbor_foci_diagnostic':int(hex_foci),
                        'missing_go_genes':sorted(needed_go-set(gx)),
                        'source_gene_duplicates':len(genes)-len(gx)})
    print(sec.section_id, 'loaded',len(lab),'spots',flush=True)

for p in sections[0]['items']:
    for held in sorted({s['patient'] for s in sections}):
        train = [s['items'][p] for s in sections if s['patient']!=held]
        test = [s['items'][p] for s in sections if s['patient']==held]
        ytr = np.concatenate([s['y'] for s in train]); yte = np.concatenate([s['y'] for s in test])
        for variant in ['legacy','exclude_program_input_readout','mean_log_axes_sensitivity','mean_log_axes_exclude_program']:
            xtr = np.vstack([s[variant] for s in train]); xte = np.vstack([s[variant] for s in test])
            ms = []
            for step in range(3):
                if step:
                    key = 'inp' if step==1 else 'neigh'
                    xtr = np.c_[xtr,np.concatenate([s[key] for s in train])]
                    xte = np.c_[xte,np.concatenate([s[key] for s in test])]
                beta = np.linalg.lstsq(xtr,ytr,rcond=None)[0]
                ms.append(float(np.mean((xte@beta-yte)**2)))
            records.append({'program':p,'patient':held,'variant':variant,
                            'mse_M0':ms[0],'mse_M1':ms[1],'mse_M2':ms[2],
                            'dM1_pct':100*(ms[0]-ms[1])/ms[0],
                            'dM2_pct':100*(ms[1]-ms[2])/ms[1]})
df = pd.DataFrame(records)
df.to_csv(OUT/'usz_t4_diagnostic.tsv',sep='\t',index=False)
old = pd.read_csv(REC/'incremental_prediction.tsv',sep='\t')
old = old[(old.cohort=='USZ') & (old.program!='P-plasma')]
paired = df[df.variant=='legacy'].merge(old,on=['program','patient'],suffixes=('_audit','_stored'))
errors = {key:float(np.max(np.abs(paired[key+'_audit']-paired[key+'_stored']))) for key in ['mse_M0','mse_M1','mse_M2','dM1_pct','dM2_pct']}
summary = df.groupby(['program','variant'])[['dM1_pct','dM2_pct']].median().reset_index()
summary.to_csv(OUT/'usz_t4_diagnostic_summary.tsv',sep='\t',index=False)
(OUT/'usz_data_checks.json').write_text(json.dumps({'sections':diagnostics,'legacy_replay_max_errors':errors,
    'note':'Alternative baselines are sensitivity diagnostics, not accepted replacement results.'},indent=2))
print('replay max errors',errors)
print(summary.to_string(index=False))
