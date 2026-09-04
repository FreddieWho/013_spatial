#!/usr/bin/env bash
# R-04 fold-4 GPU 运行 — 远端环境安装（D-096/D-098）
# 在实例上以 root 执行；幂等；已在 aiGalaxy ubuntu22_cuda12.4 镜像验证（2026-09-04）。
set -euo pipefail

MINICONDA=/root/miniconda3
PROJ=/home/huyudi/013_spatial

# 1) Miniconda (python 3.12.7)
if [ ! -x "$MINICONDA/bin/python" ]; then
  INSTALLER=/root/Miniconda3-py312_24.9.2-0-Linux-x86_64.sh
  if [ ! -f "$INSTALLER" ]; then
    curl -fL --retry 5 -o "$INSTALLER" \
      https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-py312_24.9.2-0-Linux-x86_64.sh
  fi
  bash "$INSTALLER" -b -p "$MINICONDA"
fi
"$MINICONDA/bin/python" --version

# 2) Python 依赖（--isolated 忽略平台预置 pip 配置；tuna 镜像；PyPI 的 torch 2.5.1 即 cu124 构建）
if ! "$MINICONDA/bin/python" -c "import torch" 2>/dev/null; then
  "$MINICONDA/bin/pip" install --isolated --timeout 120 --retries 10 \
    torch==2.5.1 numpy==2.0.2 scipy==1.14.1 h5py==3.11.0 \
    -i https://pypi.tuna.tsinghua.edu.cn/simple
fi

# 3) CUDA + 确定性自检（CUBLAS_WORKSPACE_CONFIG 必须随运行设置）
export CUBLAS_WORKSPACE_CONFIG=:4096:8
"$MINICONDA/bin/python" - <<'PY'
import sys, numpy, scipy, torch
assert torch.cuda.is_available(), "CUDA not available"
torch.use_deterministic_algorithms(True)
a = torch.randn(512, 512, device="cuda"); _ = a @ a; torch.cuda.synchronize()
print("torch", torch.__version__, "| cuda_build", torch.version.cuda, "| python", sys.version.split()[0])
print("numpy", numpy.__version__, "| scipy", scipy.__version__)
print("device", torch.cuda.get_device_name(0),
      "| vram_GB", round(torch.cuda.get_device_properties(0).total_memory / 1e9, 1))
print("deterministic_matmul_ok")
PY

# 4) 环境记录（本次运行的 env lock；runtime_environment_hash 由代码自动进入 cell.json）
GPU_NAME=$("$MINICONDA/bin/python" -c "import torch; print(torch.cuda.get_device_name(0))")
cat > "$PROJ/infra/r04/environment.lock.gpu_$(date +%Y%m%d).json" <<EOF
{
  "schema": "r04.environment.v1",
  "python": "3.12.7",
  "environment_path": "$MINICONDA",
  "backend": "torch",
  "packages": {"numpy": "2.0.2", "scipy": "1.14.1", "h5py": "3.11.0", "torch": "2.5.1+cu124"},
  "gpu_devices_at_capture": ["$GPU_NAME"],
  "note": "Rented GPU instance per D-096/D-098. Runs set R04_STORAGE_RESERVE_BYTES=100GB (D-086 precedent) and CUBLAS_WORKSPACE_CONFIG=:4096:8.",
  "captured_at": "$(date +%F)"
}
EOF

echo SETUP_DONE
