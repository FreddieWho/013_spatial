"""Bounded real-data oracles and hidden-value interventions for the full GO stage."""
import csv
import json
from pathlib import Path

import numpy as np

from r16.recovery.corrected import feature_block, masked_predict
from r16.recovery.gobp_stage import FeatureEngine, geometry, mask_predictions
from r16.recovery.mask_tasks import generate_windows, window_mask
from r16.recovery.repair_pipeline import axes
from scripts.r16_gobp_stage import SOURCE, OUT, digest, atomic_json


def main():
    defs = json.loads((SOURCE / 'go_representatives.json').read_text())
    ordered = sorted(defs, key=lambda g:(len(defs[g]['input_genes']) + len(defs[g]['readout_genes']), g))
    chosen = [ordered[0], ordered[len(ordered)//2], ordered[-1]]
    definitions = {g:defs[g] for g in chosen}
    genes = json.loads((SOURCE / 'genes.json').read_text())
    rows = []
    for s in json.loads((SOURCE / 'sections.json').read_text()):
        if s['cohort'] == 'DISCOVERY':
            continue
        counts = np.load(SOURCE / 'cache' / f"{s['key']}_counts.npy", mmap_mode='r')
        xy = np.load(SOURCE / 'cache' / f"{s['key']}_coords.npy")
        windows = generate_windows(xy, 12, 8, 20260921)
        wi, hidden = next((i, window_mask(xy, w)) for i,w in enumerate(windows) if window_mask(xy,w).any() and (~window_mask(xy,w)).sum()>=20)
        geo = geometry(xy, hidden)
        engine = FeatureEngine(counts, genes, axes(), definitions)
        changed = np.array(counts)
        changed[hidden] = 1e6
        altered_engine = FeatureEngine(changed, genes, axes(), definitions)
        for mode in ['raw','depth_normalized']:
            ids,x,y,bad = next(engine.blocks(mode))
            ai,ax,ay,ab = next(altered_engine.blocks(mode))
            assert not bad and not ab and ids == ai == chosen
            predicted = mask_predictions(x,y,geo)
            altered = mask_predictions(ax,ay,geo)
            assert np.array_equal(predicted, altered), (s['key'],mode,'hidden intervention')
            path = OUT / 'sections' / f"{s['key']}_{mode}.npz"
            saved = np.load(path)
            saved_index = {g:i for i,g in enumerate(saved['ids'])}
            for j,g in enumerate(chosen):
                d = definitions[g]
                f = feature_block(counts, genes, d['input_genes'], d['readout_genes'], axes(), mode)
                feature_delta = max(float(np.max(abs(f['X']-x[:,j]))),float(np.max(abs(f['y']-y[:,j]))))
                assert feature_delta<1e-10
                direct = masked_predict(counts,xy,genes,d['input_genes'],d['readout_genes'],axes(),hidden,mode)
                direct_array = np.column_stack([direct[m] for m in ['M0','M1','M2','NN','KNN8']])
                prediction_delta = float(np.max(abs(direct_array - predicted[:,j])))
                np.testing.assert_allclose(direct_array,predicted[:,j],rtol=1e-7,atol=1e-7)
                mse = ((f['y'][hidden,None] - direct_array)**2).mean(0)
                saved_delta = float(np.max(abs(saved['masked_mse'][saved_index[g],wi]-mse)))
                np.testing.assert_allclose(saved['masked_mse'][saved_index[g],wi],mse,rtol=1e-7,atol=1e-7)
                rows.append(dict(section=s['section'],program=g,mode=mode,window=wi,feature_max_delta=feature_delta,
                                 prediction_max_delta=prediction_delta,saved_mse_max_delta=saved_delta,hidden_change_max_delta=0.))
            saved.close()
        print('oracle',s['key'],flush=True)
    with (OUT / 'real_data_oracles.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    atomic_json(OUT / 'oracle_receipt.json',dict(status='PASS',program_section_mode_interventions=len(rows),
                  predicted_models_per_intervention=5,program_selection='deterministic smallest/median/largest gene sets',
                  max_feature_difference=max(r['feature_max_delta'] for r in rows),
                  max_prediction_difference=max(r['prediction_max_delta'] for r in rows),
                  max_saved_mse_difference=max(r['saved_mse_max_delta'] for r in rows),
                  max_hidden_intervention_difference=0,output_sha256=digest(OUT / 'real_data_oracles.tsv'),code_sha256=digest(Path(__file__))))


if __name__ == '__main__':
    main()
