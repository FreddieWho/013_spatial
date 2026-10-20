#!/usr/bin/env bash
set -euo pipefail
export PYTHONNOUSERSITE=1
export CUDA_VISIBLE_DEVICES=-1
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd /home/huyudi/013_spatial
exec /opt/anaconda3/bin/python scripts/field_discovery_20261018.py "$@"
