# 受控线程探针交接（review N-01 / required_fix.minimal_verification）

范围：只对 `source_0` 的 EMT（target index 3）做受控探针；未改正式启动器、未全量重拟合、未 evaluate/summarize、未放宽容差、未换求解器、未做其它扫描。
脚本与结果（本目录）：`probe_source0_emt_threads.py`（可复用小脚本）、`probe_results.tsv`、`probe_results.json`、`probe.log`、`sampled_rows.npz`、`contract.json`（机器合同副本，仅供脚本读取）。

## 精确命令

```bash
cd /home/huyudi/013_spatial
PYTHONNOUSERSITE=1 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 \
  CUDA_VISIBLE_DEVICES=-1 /opt/anaconda3/bin/python \
  infra/field_blas_probe_20261019/probe_source0_emt_threads.py
# 返回码 0；墙钟 172 s
```
线程由 `threadpool_limits(limits=k)` 设定，并用 `threadpoolctl.threadpool_info()` 记录**实际生效**线程数（非仅环境变量）：`{1: ['1'], 8: ['8']}`。

## 设计身份核对（不预先假定一致；与 fit_diagnostics 逐项比对）

| 量 | 探针实测 | `splits/source_0/fit_diagnostics.tsv`（point_ridge/EMT） | 一致 |
|---|---|---|---|
| 训练行 | **25803** | 25803 | ✓ |
| 原列数 | **1500** | 1500 | ✓ |
| 训练常量列剔除 | **2** | 2 | ✓ |
| 缩放列 | 0 | 0 | ✓ |
| 上下文列 | 6 | 6 | ✓ |
| dtype | float32 | — | — |
| solver / alpha | auto / 10.0 | auto / 10.0 | ✓ |
| 训练 RMSE | 0.0611007661（t=1）/ 0.0611007623（t=8） | 0.06110076233744621 | ✓（t=8 位级一致） |

同时记录设计指纹：`design_sha256 / y_sha256 / weights_sha256 / keep_mask_sha256`（见 `probe_results.json`）；抽样行已存 `sampled_rows.npz`；`W/splits.json` 与冻结 `PRED/splits.json` **哈希相同**（无划分路径差异）；基因矩阵来自 `infra/field_followup_20261018/genes/NCBI642_t3.npy`（哈希已记录）。
**结论：设计与保存路径不存在可解释差异**（本项按审核要求单独核对，未归因环境）。

## 受控结果（同一设计，仅线程不同；每个设置重复 2 次）

单位：预测最大绝对差。`rep0`=同设置重复差；`1↔8`=线程对照；`frozen`=冻结 `gene_ridge`/`gene_ridge_smooth_half`；`stored`=当前运行存盘包。

| 设置 | 片/读出 | rep0（同设置噪声） | 1↔8 | vs frozen | vs stored |
|---|---|---|---|---|---|
| t=1 | NCBI642/point | 0.0 | **2.08e-4** | 2.09e-5 | **2.08e-4** |
| t=1 | NCBI642/smooth | 0.0 | **2.07e-4** | 1.83e-5 | **2.07e-4** |
| t=1 | NCBI643/point | 1.19e-7 | **2.08e-4** | 2.40e-5 | **2.08e-4** |
| t=1 | NCBI643/smooth | 5.96e-8 | **2.07e-4** | 2.16e-5 | **2.07e-4** |
| t=8 | NCBI642/point | 0.0 | **2.08e-4** | 2.28e-4 | **2.38e-7** |
| t=8 | NCBI642/smooth | 0.0 | **2.07e-4** | 2.24e-4 | **1.19e-7** |
| t=8 | NCBI643/point | 5.96e-8 | **2.08e-4** | 2.29e-4 | **4.77e-7** |
| t=8 | NCBI643/smooth | 0.0 | **2.07e-4** | 2.26e-4 | **2.38e-7** |

设计条件数（本次脚本重算，补上审核指出缺失的配套脚本）：缩放后 Gram 特征值 **8.77e-12 / 4.51e6 → 原始 rcond 1.9e-18**；ridge(α=10) 系统 rcond **2.2e-6**；系数 L2 范数 7.35e-2。

## 判读（供 Grok 判断，我不宣布结论）

1. **线程设置能复现约 2e-4 的差距**：同一设计、同抽样、同权重、同 dtype、同 alpha/solver、同常量列 mask，仅把生效线程从 1 改为 8，point 与平滑预测同时变化 **≈2.07–2.08e-4**；而同设置重复噪声仅 ≤1.2e-7。
2. **存盘包 ≡ 8 线程实现**：t=8 预测与存盘包差 **1.2e-7–4.8e-7**；t=1 与存盘差 **2.07e-4**。即全量 fit 实际以 8 线程 BLAS 生效（启动器里 `OPENBLAS/MKL=2` 被 `threadpool_limits(8)` 覆盖），这不是仅凭同设置换求解器的小差异推断。
3. **冻结基线接近 1 线程行为但尚不完全一致**：t=1 与冻结包仍差 **1.8–2.4e-5**（比 t=8 的 2.2–2.3e-4 小一个量级）。冻结轮次合同是 `threads=2, workers=3`；t=2 未被本次探针覆盖，**这是尚未解释的残差**，但量级已在 1e-4 界内、1e-5 界外（每片 RMSE 差待同法核对）。
4. 因此：`2e-4` 量级差异可归因于 BLAS 生效线程；**剩余 2e-5 量级的残差尚未定因**（最可能是冻结轮的 2 线程设置或该轮 BLAS 版本/库差异）。按审核 `if_not_confirmed` 要求，不换求解器压界、不静默扩大计算。

## 下一步（等复核批准，本次未执行）

- 最小追加探针：同一脚本增加 `threads=2`（冻结轮次合同的设置）与一次 `threadpool_info` 库清单比对，检验 t=2 是否把与冻结包的差异压到 ~1e-6 量级；以及是否需要在启动器层把并发线程固定为与冻结轮次一致（`threadpool_limits(2)`），再按 `recompute_if_confirmed` 决定是否重算 44 个拟合/500 包并重跑审计。
- 未做：改启动器、全量重拟合、evaluate/summarize。

## 保留与未提交

- 本目录脚本与结果为新增运行记录；旧失败记录（`field_raw_gene_spatial_20261019/archive/*.log`）与当前 500 个预测包原样保留，未被覆盖。
- 500 个预测包仍为未提交大文件（逐包 sha256 见 `field_raw_gene_spatial_20261019/fit_package_inventory.tsv`）。
