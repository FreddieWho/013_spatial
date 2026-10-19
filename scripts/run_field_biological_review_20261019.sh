#!/usr/bin/env bash
set -euo pipefail
cd /home/huyudi/013_spatial
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=-1
PY=/opt/anaconda3/bin/python
case "${1:-check}" in
  check)     "$PY" scripts/field_biological_review_20261019.py check "${@:2}" ;;
  prepare)   "$PY" scripts/field_biological_review_20261019.py prepare ;;
  run)       "$PY" scripts/field_biological_review_20261019.py run ;;
  summarize) "$PY" scripts/field_biological_review_20261019.py summarize ;;
  *) echo "usage: $0 {check|prepare|run|summarize}" >&2; exit 2 ;;
esac
