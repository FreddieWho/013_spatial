# 2026-09-21 implementation repair replay

Protocol: `docs/decisions.md` D-140 and D-141. Historical outputs under `infra/r16/recovery_20260918/` are input-only; their former entrypoint source is preserved in `legacy_source/`.

All runs use existing local data and CPU. No GPU, new dataset, corridor scan, K search, or model retraining. The GO phase distributes independent sections across 8 CPU workers; aggregation remains deterministic and patient-based. Other phases run sequentially.

```bash
export LD_LIBRARY_PATH=/opt/anaconda3/lib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/spatial-repair-mpl
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=.
python -m r16.recovery.repair_pipeline prepare
python -m r16.recovery.repair_pipeline six
python -m r16.recovery.repair_pipeline go
python -m r16.recovery.repair_pipeline mask
python -m r16.recovery.repair_pipeline spline
python -m r16.recovery.repair_pipeline pc
python scripts/r04_finalize_phase.py
python scripts/r16_repair_report.py
```

Run from the project root. Outputs replace this repair directory's current artifacts, not historical analyses. `*_receipt.json` with `COMPLETED` is required; logs and `RUNNING` progress files alone do not establish completion. Early incomplete GO attempts were interrupted during source-coverage corrections and CPU execution improvements; they are not scientific results.

The 9012-gene analysis panel is shared by the external sources and the discovery-cache column universe. It is **not** a claim that every discovery source measured all 9012 genes. `discovery_source_coverage.tsv` and per-section `unmeasured_genes` control eligibility. The six original frozen programs are not silently truncated: `six_source_coverage.tsv` distinguishes source availability from matched-panel eligibility. In particular, ST-CRC source data fully cover epi/smmhc, although the cross-cohort matched panel does not.

- `go_aliases.json`: every alias uses its executed representative's exact split/hash.
- `*_prediction.tsv`: leave-one-patient-out errors; sections and patients have equal weight at their respective aggregation levels.
- `masked_prediction_intervention.json`: complete prediction invariance after all hidden counts are replaced; not merely a neighborhood-feature probe.
- `masked_structure_*` and `l009_*`: unknown labels excluded; biological N is patients. L009 centers are selected from known TLS and are targeted controls, not blind localization.
- `spline_results.tsv`: true natural splines; train-only knots for buffered spatial holdout. All association p values remain NOT_CALIBRATED. A noiseless linear injection checks this specific implementation capability only.
- `pc_actual_spot_regression.tsv`: actual molecular scores regressed on composition; loading-space projection is not a biological novelty test.
- `r04_final_gate.json`: specificity and the untested common residual field are distinguished. Operational stop is preserved.

No result here establishes a calibrated hidden-structure posterior, cross-plane prediction, a causal effect, or independent confirmatory pathway discovery. Bootstrap intervals condition on the fitted cross-validation predictions and are exploratory.
