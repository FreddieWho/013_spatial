# 同原始基因输入空间对照：全量运行交接（在冻结基线复现关卡处停止）

授权：`review/reviewed_implementation_round2.json` = `ACCEPTED_FOR_FROZEN_DESCRIPTIVE_RUN`，`allow_full_run=true`，R-01 已关闭，未关闭结果问题为空；只允许按 Git `0e6e5f5` 及其冻结合同运行。
结果：**fit 全量完成（44/44 case，500 个预测包，全部收敛、全部有限）；随后按任务约定执行逐 case 冻结基线复现核对——11/44 case 超出预定界，链条在 evaluate/summarize 之前停止。未放宽任何容差、未更换求解器或参数、未做广泛重跑。本轮无科学结论，交 Grok 判断。**

## 1. 版本与范围

- 运行前核实：`git diff 0e6e5f5` 对已审路径（脚本、启动器、机器合同、交接、科学合同）**为空**；HEAD = `0e6e5f5`。
- 本提交（仅我拥有的新成果与运行记录）：`infra/field_raw_gene_spatial_20261019/{audit_frozen_agreement.py,audit_frozen_agreement.tsv,audit_frozen_agreement.json,audit_numerical_probe.json,fit_package_inventory.tsv,fit_full.log,fit_time.txt,fit_diagnostics(逐划分),full_run_stages.log}` 与本文件。
- 未提交（工作区）：500 个预测包（40.3 MB，逐包 sha256 见 `fit_package_inventory.tsv`）、`splits/*/split_receipt.json` 之外的中间产物（无）；未修改主代理合同/根文档。

## 2. 实际命令、返回码、阶段计时

```bash
bash scripts/run_field_raw_gene_spatial_20261019.sh fit      # 见 fit_full.log / fit_time.txt
/opt/anaconda3/bin/python infra/field_raw_gene_spatial_20261019/audit_frozen_agreement.py
# 链条随后在 audit 返回码 1 处按约定停止，未执行 evaluate / summarize
```
| 阶段 | 返回码 | 墙钟 | 进程 User | 进程 Sys | CPU% | 最大 RSS |
|---|---:|---:|---:|---:|---:|---:|
| fit | 0 | **1048 s (17:27.98)** | **5130.39 s** | **812.89 s** | 567% | **5.19 GB** |
| audit | **1** | 1 s | — | — | — | — |
| evaluate / summarize | 未执行 | — | — | — | — | — |

说明：进程 CPU 秒取自 `/usr/bin/time -v`（含被回收子进程），**不是**把 8 线程墙钟当成 CPU；`fit` 的实际并行度 ≈5.7 核。磁盘：运行前 844.5 GB 可用，运行后 **827.4 GB**（净增约 17 GB，仍 ≫ 200 GB 下限）。

## 3. 覆盖与规模（真实行数）

- 44/44（11 划分 × 4 目标）case 均有 fit 诊断；**176** 条诊断行 = 44 × 4 个拟合模型；**500** 个预测包 = 125 次切片访问 × 4 目标（每包 5 模型）；全部 `converged=True`、`solver=auto`、**0 个非有限数组**。
- 因此实现交接中旧成本口径的 **5000 行已按 2500 行纠正**（full15 单视图：125×4×5=2500）；本轮 evaluate 未运行，所以逐片记录尚未产生。
- 训练行常量列剔除（仅训练行决定）：44 个 case 中 20 个剔除 0 列、12 个剔除 1 列、8 个剔除 2 列、4 个剔除 3 列——被剔除者恒为某些折里恒定的器官/平台指示列；按合同该列在训练中系数为 0，剔除与保留在预测上等价。

## 4. 冻结基线复现核对（停止原因）

容差在看结果前写定（本轮未改）：逐点 `1e-4`、每片 RMSE 绝对差 `1e-5`（point 与固定平滑同界）。逐 case 结果：`audit_frozen_agreement.tsv`；汇总：`audit_frozen_agreement.json`。

**11/44 case 超出；通过者也在界边**（通过 case 最差：逐点 7.27e-5、每片 RMSE 8.86e-6）：

| 划分 | 目标 | 片数 | point max\|Δ\| | point RMSE差 | smooth max\|Δ\| | smooth RMSE差 |
|---|---|---:|---:|---:|---:|---:|
| source_0 | EMT | 2 | **2.29e-4** | **2.08e-4** | 2.26e-4 | 2.11e-4 |
| source_4 | B_AXIS | 34 | 1.10e-4 | 5.06e-6 | 5.74e-5 | 5.21e-6 |
| source_5 | OXPHOS | 1 | 1.06e-4 | 1.56e-5 | 7.70e-5 | 1.51e-5 |
| source_5 | B_AXIS | 1 | 1.05e-4 | 4.88e-5 | 1.00e-4 | 4.92e-5 |
| internal_4 | OXPHOS | 12 | 5.89e-5 | 1.03e-5 | 5.46e-5 | 9.86e-6 |
| source_0 | B_AXIS | 2 | 5.75e-5 | 4.47e-5 | 5.50e-5 | 4.56e-5 |
| source_0 | OXPHOS | 2 | 5.61e-5 | 3.05e-5 | 5.08e-5 | 3.07e-5 |
| internal_1 | OXPHOS | 15 | 5.22e-5 | 1.55e-5 | 4.18e-5 | 1.53e-5 |
| internal_1 | B_AXIS | 15 | 5.13e-5 | 1.94e-5 | 4.47e-5 | 2.18e-5 |
| source_2 | EMT | 2 | 4.16e-5 | 1.23e-5 | 2.91e-5 | 1.25e-5 |
| source_0 | IFN-γ | 2 | 3.89e-5 | 2.10e-5 | 3.65e-5 | 2.28e-5 |

未超界的 case 中，逐点最差 7.27e-5 与界 1e-4 只差 1.4 倍，每片 RMSE 最差 8.86e-6 与界 1e-5 只差 1.1 倍——即整个语料都贴着这两条界，不是个别离群。

## 5. 数值影响证据（定位到环境敏感的近似零方向，而非参数/管线差异）

运行日志保留 **36 条 `LinAlgWarning: Ill-conditioned matrix (rcond=3.15e-8 … 4.51e-8)`**（来自三块设计；point 设计本身未触发）。

对最差 case（source_0 / EMT）做了**一次有界、只读的重拟合探针**（不写任何冻结产物，`audit_numerical_probe.json`）：

| 量 | NCBI642 | NCBI643 |
|---|---:|---:|
| 同设计换求解器（cholesky↔lsqr） | 5.45e-6 | 6.35e-6 |
| 新重拟合 vs 冻结 `gene_ridge` | 2.09e-5 | 2.40e-5 |
| **运行存盘 vs 冻结 `gene_ridge`** | **2.28e-4** | **2.29e-4** |
| **运行存盘 vs 新重拟合（同代码同种子）** | **2.08e-4** | **2.08e-4** |

- point 设计的原始 Gram 最小/最大特征值 = **6.23e-11 / 4.51e6 → 原始 rcond ≈ 1.4e-17（数值上奇异）**；ridge（α=10）把系统条件数抬到 ≈2e-6，但**近零方向完全由舍入决定**，其上分到的系数权重 ≈O(1)；系数 L2 范数 7.35e-2，cholesky 与 lsqr 的训练 RMSE 相同到 1e-10。
- 关键判读依据：同一次重拟合内部两种求解器只差 ~6e-6，**而运行存盘与新重拟合（同 Git 代码、同种子、同数据、同权重、同 alpha、同列剔除）差 ~2e-4** ⇒ 差异不是参数或管线差别，而是**该病态设计下近零方向舍入导致的环境（BLAS/线程）敏感性**；冻结产物恰好落在另一个舍入实现上。
- 该敏感性对**候选/背景/置换**更需留意：三块设计的 rcond（3e-8）比 point 更差，但同一次运行内四种模型共享同一环境，配对差值内部一致；受影响的是「相对冻结基线/跨环境复现」这一层面。
- 上表所有量都远大于数值噪声下限，也不满足预定的 1e-4 / 1e-5 复现界限；按任务约定**停在此处交 Grok 判断**，不自行放宽容差、不换求解器/参数、不重跑。

## 6. 保留物 / 复核所需文件

- 关卡证据：`audit_frozen_agreement.py`、`audit_frozen_agreement.tsv`、`audit_frozen_agreement.json`、`audit_frozen_agreement.log`、`audit_numerical_probe.json`。
- 运行记录：`fit_full.log`（含 36 条警告原文）、`fit_time.txt`（墙钟/CPU/RSS）、`full_run_stages.log`（阶段返回码）、`splits/*/fit_diagnostics.tsv`（176 行）、`fit_package_inventory.tsv`（500 包逐包 sha256/size）。
- 授权与历史：`review/reviewed_implementation_round2.json`、`review/reviewed_implementation_round1.json`、`review/preflight.md`、`implementation_handoff.md`。
- 未提交：500 个预测包（40.3 MB，清单见 inventory）、`smoke/` 与大中间产物；未提交状态如实记录，不构成本轮实现版本。

## 7. NOT_RUN / 未产生

- **evaluate 与 summarize 未运行**（按停止约定）；因此无 `prediction_metrics.tsv`/`macro_*/unit_*/*paired*` 等 2500 行逐片记录与配对表，也无该轮的 freeze manifest。
- 未运行：全量配对比较、readout_claim_summary（本轮不存在）、任何显著性、独立 marker、人工 TLS、物理尺度、视野外推断。
- 未做：容差放宽、求解器更换、参数/划分/缩放调整、广泛重跑、另开重复训练（仅对最差 case 做了一次只读重拟合探针）。

## 8. 交付后冻结

本文件交付后停止修改实现/合同/检查/运行记录，等待 Grok 就"冻结基线复现超界 + 近零方向环境敏感性"给出处置判断（例如：是否需要在合同层固定 BLAS 线程/求解器以实现跨环境可复现、或如何界定哪些读出在此分辨率下可解释）；在得到指令前不继续 evaluate/summarize、不宣布科学结论。
