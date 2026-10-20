#!/usr/bin/env bash
set -euo pipefail
cd /home/huyudi/013_spatial
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=-1
/opt/anaconda3/bin/python scripts/field_source_holdout_20261018.py prepare > infra/field_source_holdout_20261018/prepare.log 2>&1
/opt/anaconda3/bin/python scripts/field_source_holdout_20261018.py run > infra/field_source_holdout_20261018/run.log 2>&1
