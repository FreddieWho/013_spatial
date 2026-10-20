#!/usr/bin/env bash
set -euo pipefail
export PYTHONNOUSERSITE=1 CUDA_VISIBLE_DEVICES=-1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /home/huyudi/013_spatial
exec /opt/anaconda3/bin/python scripts/field_followup_20261018.py "$@"
