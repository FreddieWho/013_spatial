# 受控线程探针交接（review N-01 / required_fix.minimal_verification）

**追加轮（threads=1/2/8）**：在 1/8 原结果保持不变的前提下（原文件另存为 `probe_results_threads_1_8_original.{json,tsv}`、`sampled_rows_threads_1_8_original.npz`），用同一脚本补测旧冻结合同规定的 **2 线程**，输出到 `probe_results_t128.{json,tsv}`、`probe_t128.log`。命令：

```bash
cd /home/huyudi/013_spatial
PYTHONNOUSERSITE=1 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 \
  CUDA_VISIBLE_DEVICES=-1 /opt/anaconda3/bin/python \
  infra/field_blas_probe_20261019/probe_source0_emt_threads.py --threads 1 2 8 --out-suffix _t128
# 返回码 0
```

生效线程由 `threadpool_limits` + `threadpool_info()` 逐库记录：t=1 → mkl/openmp 均 1；t=2 → 均 2；t=8 → 均 8。每个设置重复 2 次，同设置重复噪声 ≤1.8e-7。

## 追加轮结果（同设计，仅线程不同）

两两对照（max_abs / 每片 RMSE 绝对差，NCBI642 与 NCBI643 的 point 与 smooth 一致）：

| 对照 | max_abs | 每片 RMSE 绝对差 |
|---|---:|---:|
| **t1 ↔ t2** | 1.04e-6 – 1.61e-6 | 5.1e-7 – 6.0e-7 |
| **t1 ↔ t8** | 2.07e-4 – 2.08e-4 | 1.94e-4 – 2.00e-4 |
| **t2 ↔ t8** | 2.07e-4 | 1.93e-4 – 2.00e-4 |

按原定界限逐项报告（每个设置取该读出跨两片的最差值）：

| 设置 vs 参照 | point 最差 max_abs（限 1e-4） | point 最差每片RMSE差（限 1e-5） | smooth 最差 max_abs | smooth 最差每片RMSE差 |
|---|---:|---:|---:|---:|
| t=1 vs 冻结 | 2.40e-5 **通过** | 1.25e-5 **未通过** | 2.16e-5 通过 | 1.25e-5 **未通过** |
| t=2 vs 冻结 | 2.44e-5 **通过** | 1.31e-5 **未通过** | 2.22e-5 通过 | 1.31e-5 **未通过** |
| t=8 vs 冻结 | 2.29e-4 **未通过** | 2.08e-4 **未通过** | 2.27e-4 未通过 | 2.11e-4 **未通过** |
| t=1 vs 存盘 | 2.08e-4 **未通过** | 1.98e-4 **未通过** | 2.07e-4 未通过 | 2.00e-4 **未通过** |
| t=2 vs 存盘 | 2.07e-4 **未通过** | 1.98e-4 **未通过** | 2.07e-4 未通过 | 2.00e-4 **未通过** |
| t=8 vs 存盘 | 4.8e-7 **通过** | 5.9e-9 **通过** | 2.4e-7 通过 | 7.1e-9 通过 |

设计身份（追加轮复核，均与 `fit_diagnostics` 一致）：**25803 行 / 1500 原列 / 剔除 2 常量列 / scaled=0 / ctx=6 / float32 / alpha=10 / solver=auto**；运行 `splits.json` 与冻结哈希相同；Gram 原始 rcond 1.9e-18、ridge rcond 2.2e-6。

## 残差仍未定因（关键结论）

- **已确认（线程解释超界部分）**：>1e-4 的差距完全由生效线程数解释——存盘包 ≡ 8 线程（t8↔存盘 ≤4.8e-7），而 1/2 线程与 8 线程相差 2.07e-4；旧冻结轮次按合同用 2 线程，两者落在同一侧。
- **仍未确认（小残差）**：所有设置（1、2、8 线程）与冻结包都稳定相差 **≈2.2–2.4e-5（max_abs，通过 1e-4）与 ≈1.25–1.31e-5（每片 RMSE，未通过 1e-5）**；1↔2 线程只有 1.6e-6，说明这不是线程噪声，而是我的重拟合族与冻结产物之间一个**稳定的 ~1e-5 偏移**。**逐项结论：即使采用冻结轮次自身的 2 线程设置，每片 RMSE 界限（1e-5）仍未满足（最好 1.25e-5，超出约 25%）。**
- **库版本已排除**：冻结轮次记录的运行环境与当前**完全相同**（python 3.11.5 / numpy 1.26.4 / scipy 1.13.1 / sklearn 1.5.1 / h5py 3.11.0；见 `field_source_holdout_review_20261018/environment.json` 与 `field_discovery_20261018/environment.json`），故残差也不能归因库版本。
- 因此按审核 `if_not_confirmed` 分支：**不换求解器压界、不全量重拟合、不改正式启动器**；需要 Grok 判断残差的下一步（例如核对冻结轮次设计数组是否可复原/是否有未保存的设计差异，或该 1e-5 界限在本设计上是否可达）。

以下为第一轮 1/8 探针原文（保留）。

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
