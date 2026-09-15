#!/usr/bin/env bash
# GraphST rental setup (D-121). Ubuntu 22.04 + CUDA 12.x driver image.
# V100 (sm_70): torch 2.1.x cu118 is the safe documented combo
# (forward-compatible driver; sm_70 still supported; NOT newest torch).
# Deliberately NO: DGL / torch-geometric (GraphST core is pure torch),
# R+rpy2 (we cluster with Leiden, not mclust), conda (venv is faster).
# Usage: bash scripts/r16_graphst_setup.sh <workdir>
set -euxo pipefail

WORKDIR="${1:?usage: r16_graphst_setup.sh <workdir>}"
cd "$WORKDIR"

nvidia-smi --query-gpu=name,compute_cap,memory.total --format=csv
python3 --version

python3 -m venv .venv-graphst
# shellcheck disable=SC1091
source .venv-graphst/bin/activate
pip install --upgrade pip wheel setuptools
# torch first (CUDA 11.8 wheels; driver CUDA 12.x runs them fine)
pip install "torch==2.1.2+cu118" --index-url https://download.pytorch.org/whl/cu118
pip install "scanpy==1.9.8" "anndata==0.10.7" "scikit-learn==1.3.2" \
  "leidenalg==0.9.1" "igraph==0.11.4" "pandas==2.1.4" "tqdm" "h5py" "pyarrow"
# GraphST source (no release tags; pin the commit we run)
if [ ! -d GraphST-src ]; then
  git clone --depth 1 https://github.com/JinmiaoChenLab/GraphST.git GraphST-src
fi
git -C GraphST-src rev-parse HEAD | tee graphst_version.txt
export PYTHONPATH="$WORKDIR/GraphST-src:$PYTHONPATH"

python3 - <<'EOF'
import torch
assert torch.cuda.is_available(), "CUDA not visible"
cap = torch.cuda.get_device_capability(0)
print("gpu:", torch.cuda.get_device_name(0), "sm_%d%d" % cap, "torch:", torch.__version__)
assert cap >= (7, 0), "compute capability too old?"
import scanpy, anndata, sklearn, igraph, leidenalg
print("scanpy/anndata/sklearn ok")
import sys
sys.path.insert(0, "GraphST-src")
from GraphST import GraphST
print("GraphST import ok")
EOF
echo "SETUP_OK venv=$WORKDIR/.venv-graphst"
