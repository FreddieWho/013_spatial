"""Known-model real-geometry calibration; never substitutes for a fitted spatial null."""
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.stats import norm

from r16.recovery.corrected import fit_linear
from r16.recovery.gobp_stage import FeatureEngine
from r16.recovery.repair_pipeline import axes
from scripts.r16_gobp_stage import SOURCE, OUT, atomic_json, direction_weights, digest


def smoothing(xy, width):
    if width == 0:
        return sparse.eye(len(xy), format='csr')
    coo = cKDTree(xy).sparse_distance_matrix(cKDTree(xy), 3 * width, output_type='coo_matrix')
    values = np.exp(-.5 * (coo.data / width) ** 2)
    matrix = sparse.csr_matrix((values, (coo.row, coo.col)), shape=(len(xy), len(xy)))
    # Each spot's simulated noise has exactly unit marginal variance.
    return sparse.diags(1 / np.sqrt(np.asarray(matrix.multiply(matrix).sum(1)).ravel())) @ matrix


def main():
    defs = json.loads((SOURCE / 'go_representatives.json').read_text())
    ordered = sorted(defs, key=lambda g: (len(defs[g]['input_genes']) + len(defs[g]['readout_genes']), g))
    chosen = [ordered[0], ordered[len(ordered) // 2], ordered[-1]]
    selected = {g: defs[g] for g in chosen}
    genes = json.loads((SOURCE / 'genes.json').read_text())
    plan = dict(program_selection='Smallest, median-sized, largest frozen covered gene sets, deterministic tie by name',
                programs=chosen, widths_pitch=[0, 1, 3, 6], amplitudes_sd=[0, .1, .25, .5, 1],
                draws=1000, seed=20260922, alpha=.05,
                null='Independent Gaussian innovations smoothed on observed coordinates, row L2 normalization; known covariance',
                injection='Smoothed frozen directional contrast template, RMS 1 before nuisance projection',
                patient_test='Exact two-sided sign flip on three patient-level contrasts',
                comparison='Known-covariance Gaussian z vs misspecified independent-spot Gaussian z',
                boundary='Power/FPR only for stated simulation. No fitted real-data variogram, no real-data spatial p.',
                code_sha256=digest(Path(__file__)))
    atomic_json(OUT / 'calibration_contract.json', plan)
    groups = []
    for section in json.loads((SOURCE / 'sections.json').read_text()):
        if section['cohort'] != 'USZ':
            continue
        xy = np.load(SOURCE / 'cache' / f"{section['key']}_coords.npy")
        labels = np.load(SOURCE / 'cache' / f"{section['key']}_labels.npy")
        _, dw = direction_weights(xy, labels)
        if not dw.shape[1]:
            continue
        counts = np.load(SOURCE / 'cache' / f"{section['key']}_counts.npy", mmap_mode='r')
        engine = FeatureEngine(counts, genes, axes(), selected)
        designs = {}
        for mode in ['raw', 'depth_normalized']:
            ids, x, _, bad = next(engine.blocks(mode))
            assert not bad and ids == chosen
            designs[mode] = x[:, :, :9]
        groups.append((xy, dw.mean(1), designs))
    assert len(groups) == 3
    from itertools import product
    signs = np.array(list(product([-1., 1.], repeat=3)))
    rows = []
    for width in plan['widths_pitch']:
        moments = {mode: [] for mode in ['raw', 'depth_normalized']}
        for xy, w, designs in groups:
            S = smoothing(xy, width)
            signal = S @ w
            signal /= np.sqrt(np.mean(signal * signal))
            for mode, x in designs.items():
                patient = []
                for j in range(len(chosen)):
                    q = x[:, j]
                    v = w - q @ fit_linear(q, w)
                    a = S.T @ v
                    patient.append((float(a @ a), float(v @ v), float(v @ signal)))
                moments[mode].append(patient)
        for mode, m in moments.items():
            m = np.asarray(m)  # patient, program, true variance/iid variance/signal mean
            for j, program in enumerate(chosen):
                rng = np.random.default_rng(plan['seed'] + j + int(width * 11))
                # Exact projected Gaussian distribution, avoiding a needless dense field bank.
                noise = rng.normal(size=(plan['draws'], 3)) * np.sqrt(m[:, j, 0])[None, :]
                for amplitude in plan['amplitudes_sd']:
                    values = noise + amplitude * m[:, j, 2][None, :]
                    stat = abs(values.sum(1))
                    exact_p = np.mean(abs(values @ signs.T) >= stat[:, None] - 1e-12, axis=1)
                    known_p = 2 * norm.sf(stat / np.sqrt(m[:, j, 0].sum()))
                    iid_p = 2 * norm.sf(stat / np.sqrt(m[:, j, 1].sum()))
                    rows.append(dict(program=program, mode=mode, width_pitch=width, amplitude_sd=amplitude,
                                     n_patients=3, draws=plan['draws'], min_exact_p=float(exact_p.min()),
                                     exact_patient_rejection_rate=float(np.mean(exact_p <= .05)),
                                     known_covariance_rejection_rate=float(np.mean(known_p <= .05)),
                                     naive_iid_rejection_rate=float(np.mean(iid_p <= .05))))
        print('calibrated known-model width', width, flush=True)
    import csv
    with (OUT / 'calibration_results.tsv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    atomic_json(OUT / 'calibration_receipt.json', dict(status='COMPLETED', scenarios=len(rows),
                simulated_draws=len(rows) * plan['draws'], output_sha256=digest(OUT / 'calibration_results.tsv'),
                exact_direction_p_floor=.25, real_data_spatial_p='NOT_CALIBRATED', T2P2='NOT_RUN',
                interpretation='Three-patient two-sided sign-flip test has zero rejection power at alpha .05 even without multiple testing.'))


if __name__ == '__main__':
    main()
