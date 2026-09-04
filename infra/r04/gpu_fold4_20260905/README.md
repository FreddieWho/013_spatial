# R-04 fold-4 GPU 作业手册（2026-09-05 预置）

目的：在租用 GPU 实例上按冻结协议完成 **fold 4 的 K=0/K=3 fit 8400 + 推断 2400×2**，与昨晚 fold-0 结果配对，闭合 fold 异质性裁决的证据缺口（roadmap 2026-09-04 条目）。结果仍为代表性诊断，`selected_k` 保持 `null`。

授权：D-096（GPU 路线，fold-0）+ D-098（fold-4 复跑范围）。协议与 D-080/D-084 完全同参，仅 `--fold-values 4`。

## 执行步骤（拿到实例地址后）

设 `IP=<实例IP>` `PORT=<SSH端口>` `SSHPASS=<密码>`（不写进仓库）。

```bash
SSH="sshpass -e ssh -o StrictHostKeyChecking=accept-new -p $PORT root@$IP"

# 1) 上传代码 + 数据（~651MB，约 6 分钟 @1.8MB/s）
#    manifest 用绝对路径，远端路径必须是 /home/huyudi/013_spatial（需先 mkdir -p）
#    数据为四个 panel_cache_training* 目录（47 段全部指向 part1/2/3 下的 npz）
$SSH 'mkdir -p /home/huyudi/013_spatial'
cd /home/huyudi/013_spatial
rsync -azR -e "sshpass -p $SSHPASS ssh -o StrictHostKeyChecking=accept-new -p $PORT" \
  r04 scripts infra/r04/panel_cache_training infra/r04/panel_cache_training_part1 \
  infra/r04/panel_cache_training_part2 infra/r04/panel_cache_training_part3 \
  infra/r04/gene_universe.txt infra/r04/environment.lock.json infra/r04/gpu_fold4_20260905 \
  root@$IP:/home/huyudi/013_spatial/

# 2) 装环境（幂等，~5 分钟；torch 2.5.1 PyPI 版 = cu124，避开 pypi.nvidia.cn 超时）
$SSH 'bash /home/huyudi/013_spatial/infra/r04/gpu_fold4_20260905/setup_remote.sh'

# 3) 后台启动主运行（预检 fail-closed → 配对跑 ~1-1.5h → 端点审计）
$SSH 'nohup bash /home/huyudi/013_spatial/infra/r04/gpu_fold4_20260905/run_fold4.sh > /tmp/run_fold4.log 2>&1 & echo PID=$!'

# 4) 本地监视（单次终态通知，不轮询）：
#    bg_run: ssh 循环 3 分钟心跳检测 PREFLIGHT/RUN_ALL_DONE/进程消失
#    注意：pgrep/pkill -f 会匹配到自己这条 ssh 命令，用 PID 文件判断进程死活。

# 5) 完成后回传（rsync -azR 方向：远端 → 本地，注意 /./ 相对锚点）
rsync -azR -e "sshpass -p $SSHPASS ssh -p $PORT" \
  "root@$IP:/home/huyudi/013_spatial/./infra/r04/torch_representative_fold4_20260905_k03" \
  "root@$IP:/home/huyudi/013_spatial/./infra/r04/torch_k0_fold4_checkpoint_audit_*.json" \
  "root@$IP:/home/huyudi/013_spatial/./infra/r04/torch_k3_fold4_checkpoint_audit_*.json" \
  "root@$IP:/home/huyudi/013_spatial/./infra/r04/environment.lock.gpu_*.json" \
  /home/huyudi/013_spatial/

# 6) 退租（MCP release 拒放非 MCP 创建实例；用 CLI）
~/.local/bin/ai-galaxy-compute release-plan <instance_name> > /tmp/rp.json
TOKEN=$(python3 -c "import json;print(json.load(open('/tmp/rp.json'))['approval_token'])")
~/.local/bin/ai-galaxy-compute release --preview-file /tmp/rp.json --approval-token "$TOKEN"
# instance_name 从 ~/.local/bin/ai-galaxy-compute instances 查（Container_name 字段）
```

## 预期与判读

- 预估：全程 ~2h（安装+预检 ~15min，主跑 ~1-1.5h，回传 ~5min），费用 ~¥3.5；VRAM 峰值预期 ≤11GB。
- fold 4 规模：训练 37 段 / 90,812 spot / 2 患者；留出 10 段 / 22,532 spot / 2 患者。
- 判读锚点：TF 五 seed fold-4 K3−K0 全正（均值 +20.05），fold-0 全负（均值 −30.17；Torch 复现 −28.96）。Torch fold-4 若显著为正 → fold 异质性在 Torch 复现，结论仍是 `K_UNDECIDED`；若偏负 → 与 TF 冲突，记 ISSUES 并停下报告。
- minibatch fit 门禁预期仍警告（正常）；`inference_platform` 应收敛；端点审计应为 `DIAGNOSTIC_ONLY` 且 hash 全 VERIFIED、重复差 sd=0。
- 产物：`infra/r04/torch_representative_fold4_20260905_k03/`（cell.json ×2 + run/vram 日志）+ 两份端点审计 json + 新 env lock。

## 已知坑（昨晚踩过，已固化进脚本）

- 必须环境变量：`PYTHONPATH=/home/huyudi/013_spatial`、`R04_STORAGE_RESERVE_BYTES=100000000000`、`CUBLAS_WORKSPACE_CONFIG=:4096:8`。
- pip 用 `--isolated` + tuna 镜像；download.pytorch.org 的 nvidia 外链在该网络超时。
- k_search 不接受 k=3-only（聚合需 K=0 cell）→ 必须 `--k-values 0,3` 同跑。
- rsync 回传方向与 `/./` 锚点不要写反。
