# R-04 control plane

R-04 discovers continuous latent spatial fields from molecule-only inputs.
Candidate identity is a cross-run field equivalence class; it is never a spot
cluster. A field output must contain a coordinate-aligned posterior mean, SD,
length scale, loading and input hash.

`K_model` is the number of latent spatial molecular-effect directions allowed by
the representation, not the number of biological structures. `K_eff` is only
assigned after held-out prediction, restart/subspace stability, repeatability
and spatial-null evidence. Structures and latent directions can map many to
many, so factor numbering, sign, rotation and a single factor label are not
valid structure evidence; structure claims use rotation-robust effect
representations and structure-specific readouts.

Torch v3 checkpoints keep the resumable parameter state in `params` and may
also keep a separate `best_state`; audits never mix the two. A
`parameter_layout` is validated when present, and the objective input hash
additionally covers an explicit `library_size` offset. Missing legacy v3
objective hashes are reported as unverified rather than invented. The active
execution modes are `auto` and `eager`; the historical false `compiled` frozen
label is migrated to eager with provenance, while a new compiled request fails
closed.

The discovery process may read raw/filtered 10x counts, gene IDs and spatial
coordinates only. GT, images, target-derived metadata, structure distances and
outcomes belong exclusively to the post-freeze `anchor_eval` process.

## Current implementation status

- `r04/` contains contracts, geometry, composition residualization, grouped
  splits, local reference adapters, single PyTorch model adapters, candidate
  matching and machine-readable gates. Role-pure manifests, conditional
  common-panel inference and panel-cache loaders are included. The active
  compute backend is PyTorch (`environment.lock.json` `backend=torch`,
  CPU/CUDA share one `torch.autograd`/`torch.optim` implementation); legacy
  TensorFlow artifacts are sealed as historical evidence only.
- `scripts/r04_prepare.py` builds the molecule-only input manifest.
- `scripts/r04_smoke.py` runs the deterministic dual-model synthetic smoke.
- `infra/r04/smoke/` is a reproducible pilot artifact, not a scientific result.
- `infra/r04/input_manifest.json` and `infra/r04/role_manifests/` are the
  current 71-section READY input surfaces. The training panel cache is split
  into three auditable parts and combined by
  `infra/r04/panel_cache_training/training_panel_manifest.json`.
- `infra/r04/pilot_htan_run/` is a one-section HTAN loader/fit pilot and is
  explicitly `FIT_COMPLETE_NOT_VALIDATED`; it is not a patient-level result.
- The signed adapter still models a signed GP on a fixed-overdispersion NB
  Pearson nuisance residual, and mNSF uncertainty still covers inducing weights
  only; these remain later scientific limitations rather than substitutes for
  the current count-model calibration.
- The 2026-08-13 v2 calibration uses inducing-value projection, a K-invariant
  rate-scale initialisation, one fit/infer likelihood convention, direction-
  matched section-preserving nulls and held-out-gene K scoring in which every K
  sees the same data, sections, evaluation genes and K=0 baseline.
  `scientific_gate_20260813_v10.json` is `PASS_SCIENTIFIC_CALIBRATION`: 18
  no-field, one-field, independent two-field and collinear controls pass their
  corresponding gates on disconnected, crescent and branch supports.
- `k_calibration_aggregate_v4_20260813.json` selects K=2 on a planted K_eff=2
  control across three independent data seeds; K=4 has zero scientific support
  and is not at the search boundary, so K>4 was not searched.
  `k0_gate_aggregate_v4_20260813.json` independently selects K=0 on no-field
  controls. Standard errors use data seeds, not folds, as independent units.
  The gate replays generator truth, checks all numerical values are finite, and
  re-aggregates every nested K source in addition to checking its hash. Earlier
  aggregates and gates are retained only as withdrawn history because their
  gene splits changed with K or their best-state semantics predated D-055.
  Synthetic selection does not set the final K for patient data.
- mNSF checkpoint v2 persists the historical best loss, step and parameter
  state, and the fit CLI exposes model-isolated durable checkpoint directories.
  The signed residual GP now uses the same per-spot gene-sum objective scale in
  fit and frozen inference; it remains a residual pilot rather than a complete
  signed NB posterior.
- Formal real-data R-04 remains pending, but the K=0/K=3 extreme-fold
  diagnostic and its five-seed frozen inference panel are complete. The
  provenance-checked outputs are in `restart_stability_panel_20260901_frozen.json`,
  `k_semantics_and_downstream_robustness_20260901.json` and the current linked
  handoff `k_semantics_and_downstream_robustness_20260901_final.json`: fold 0 consistently
  favors K=0, fold 4 consistently favors K=3, and the result is diagnostic only.
  Ten K=3 training-role effect exports are now available with zero fit updates,
  complete array/hash checks and validation GT sealed. The descriptive
  shared-versus-specific readout is limited to two paired patients per extreme
  fold, so downstream K robustness remains `not_tested` and factor naming stays
  forbidden. A six-cell K=2 Stage-A manifest is frozen, but its first CPU fit
  stopped before the first checkpoint; the stage is `BLOCKED_CPU_RUNTIME` with
  no K=2 scientific result. Strict training audit remains
  `STOPPED_OFF_PLATFORM`; GPU remains prohibited.
  Legacy TensorFlow checkpoints, JSON/NPZ and audits are sealed and must not be
  used as inputs to new Torch training; new science restarts from Torch.
  The current runtime storage guard is a 1.2 TB hard floor and a 1.4 TB soft
  warning line under D-057; historical 2 TB values remain dated records only.

## Environment

The default CPU environment is the user-level environment recorded in
`environment.lock.json`. GPU is not assumed. A GPU can only change runtime
after the CPU pilot estimates the need and a separate rental approval exists.
