#!/usr/bin/env bash
set -euo pipefail
cd /home/huyudi/013_spatial
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 CUDA_VISIBLE_DEVICES=-1
PY=/opt/anaconda3/bin/python
case "${1:-check}" in
  prepare)   "$PY" scripts/field_raw_gene_spatial_20261019.py prepare ;;
  fit)       "$PY" scripts/field_raw_gene_spatial_20261019.py fit ;;
  evaluate)  "$PY" scripts/field_raw_gene_spatial_20261019.py evaluate ;;
  summarize) "$PY" scripts/field_raw_gene_spatial_20261019.py summarize ;;
  smoke)     "$PY" scripts/field_raw_gene_spatial_20261019.py smoke "${@:2}" ;;
  check)     "$PY" scripts/field_raw_gene_spatial_20261019.py check "${@:2}" ;;
  all)       "$PY" scripts/field_raw_gene_spatial_20261019.py prepare \
             && "$PY" scripts/field_raw_gene_spatial_20261019.py fit \
             && "$PY" scripts/field_raw_gene_spatial_20261019.py evaluate \
             && "$PY" scripts/field_raw_gene_spatial_20261019.py summarize ;;
  *) echo "usage: $0 {prepare|fit|evaluate|summarize|smoke|check|all}" >&2; exit 2 ;;
esac
