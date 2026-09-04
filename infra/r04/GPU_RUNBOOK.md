# R-04 GPU 资源档案与租用选型（实测，2026-09-04/05 两次会话）

本文件固化实测峰值与开销，供以后租用直接选型。数据来自两台同规格实例（RTX 4090 24GB / 16 核 / 47GB 内存 / 149GB 系统盘，¥1.65/h，aiGalaxy）。

## 实测峰值（代表性 fold 配对跑：fit 8400×2 + 推断 2400×2×2）

| 资源 | 峰值 | 出现位置 |
| --- | --- | --- |
| 显存 | **13.6 GB**（预检 100 步+推断 50 步；正式配对跑 fold-0 10.3GB、fold-4 9.1GB） | 峰值由批形状决定，不随步数增长；预检值即为上界 |
| 内存 | ~16 GB RSS（47 段面板 fit 的历史实测） | 数据加载 + fit |
| 磁盘 | ~10 GB（环境 ~8GB + 数据 0.65GB + 产物 <10MB） | / |
| CPU | 16 核够用（瓶颈在 GPU 与 I/O） | |

## 时间与费用实测

| 项目 | 时间 | 说明 |
| --- | --- | --- |
| 环境安装（miniconda + torch 2.5.1 等） | ~5 min | tuna 镜像；幂等脚本 `gpu_fold4_20260905/setup_remote.sh` |
| 数据上传 651MB | ~6 min @1.8MB/s | rsync，32Mbps 带宽 |
| 预检（2-step + hash 锚点校验） | ~2 min | fail-closed |
| 单 fold 配对跑（K=0+K=3 fit + 推断） | **~50–60 min** | fold-4 optimizer wall：K=0 336s、K=3 562s（90,812 spot） |
| 端点审计 ×2 | ~3 min | |
| **单 fold 全程（含安装上传）** | **~1.2–1.5 h** | |
| 费用 | **¥1.4–3.9 / 次会话**（4090 ¥1.65/h） | fold-0 会话 ¥3.88（含调试）；fold-4 会话 ¥1.43（纯执行） |

## 选型建议

- **首选：RTX 4090 24GB（¥1.65/h）**——已两次验证，余量充足（峰值 13.6/24GB）。
- **省钱替代：RTX 3090 24GB（¥1.1/h）**——同 cu124 wheel 原生支持，慢约 20–40%，单 fold 全程约 ¥1.5–2。
- **不可用**：RTX 3080（10GB < 13.6GB 峰值）；V100-16GB（16GB 余量过薄 + sm_70 与新 CUDA 的兼容性坑，见 D-096 选型记录）。
- **避免**：RTX 5090（sm_120 需更新 torch，未验证）；A100/A800（性能过剩，贵）。
- 配置底线：显存 ≥16GB（建议 24GB）、内存 ≥32GB、磁盘 ≥20GB 空闲、CPU ≥8 核。

## 操作要点（细节见 `gpu_fold4_20260905/README.md`）

1. 远端必须：`PYTHONPATH=/home/huyudi/013_spatial`、`R04_STORAGE_RESERVE_BYTES=100000000000`、`CUBLAS_WORKSPACE_CONFIG=:4096:8`。
2. manifest 用绝对路径 → 远端必须是 `/home/huyudi/013_spatial`（新实例先 `mkdir -p`）。
3. rsync 上传用 `-azR`（漏 R 会把目录扁平化到根）；回传方向远端→本地，用 `/./` 锚点。
4. pip 用 `--isolated` + tuna 镜像（pypi.nvidia.cn 外链超时）。
5. 后台启动 ssh 命令要 `< /dev/null`；ssh 输出过管道会吞退出码。
6. 退租走 CLI：`ai-galaxy-compute release-plan <name>` → `release --preview-file --approval-token`（token 取同次 preview 文件）；MCP release 工具拒放非 MCP 创建实例。实例名用 `ai-galaxy-compute instances` 查（Container_name）。
7. 同 pin 同卡型的不同实例 runtime_environment_hash 逐位相同（已验证两台均为 fdb3b0ec…）。
