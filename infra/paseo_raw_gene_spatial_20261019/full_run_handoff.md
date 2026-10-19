# 同原始基因输入空间对照：全量描述性结果交接

**状态：PENDING_INDEPENDENT_REVIEW。** 本轮数字均为**未审核描述性结果**，不构成科学结论，不支持或否定 HF-01/HF-04；主代理负责解释，Grok 复核关闭前不得引用为结论。
授权：D-275（旧逐位复现失败**保留**为已披露失败，不再作为读取 500 包与计算的先决条件）；`review/result_comparison_review.json`（`allow_scientific_claim=false`）。实现冻结版：Git `0e6e5f5`（R-01 修复，已通过第二轮复核）；本轮只做必要最小代码修改，未重拟合、未换线程/求解器、未新增模型。

## 1. 本轮最重要修复：评价产物被旧探针污染（结果有效性问题）

- **现象**：首轮 evaluate 后内部臂只有 57 片/51 单位、`section_rows` 2200（应 2500）；`evaluation/internal_0/prediction_metrics.tsv` 时间戳 **16:27:52**，其余划分均为 18:33。
- **根因**：`check` 命令的 tiny 探针调用 `run_split_evaluate` 时只把预测来源目录重定向到 `checks/tiny/`，**指标仍写入 `WORK/evaluation/internal_0`**（test_limit=1），并留下 `COMPLETE` 标记；随后全量 evaluate 命中该标记而跳过 internal_0（16 片只算了 1 片）。
- **修复（最小）**：`run_split_evaluate` 新增 `metrics_root` 参数（探针写入自身目录，永不落生产目录）；新增 `evaluation_receipt.json` 记录 `n_sections/rows/paired/expected_*`，跳过前逐项校验，**陈旧/不完整目录一律重算而不被接受**。
- **处置**：陈旧目录**归档未删除** → `archive/evaluation_stale_partial_1843/`；随后重跑 evaluate+summarize。
- **重算后核对**（`verify_results.py`，17 项 / 0 失败）：`section_rows` **2500**（125 访问 × 4 目标 × 5 模型）、`paired_rows` **2000**、`unit_rows` **2240**、`unit_paired` **1792**；内部 **72 片/65 单位**、来源 **53 片/47 单位**；确认身份 **14 片/8 单位**、未知 **58 片**；来源身份不足 **19 片**仅留内部；器官外推划分 = {source_0, source_5}。
- 该污染也解释了首轮验证脚本的 `local-detail Pearson` 重算差 0.233 与 q75 差 0.184：它们落在 internal_0 的**陈旧行**上（那些行来自 tiny 3 片训练模型）；重算后两项重算误差 ≤1e-9（通过的 17 项含这两项）。

## 2. 冻结基线复现：保留的失败与已知数值边界（不改写成通过）

- 逐 case 核对（`audit_frozen_agreement.{json,tsv}`）：44/44 覆盖，**11 个 case 超出旧逐点 1e-4 / 每片 RMSE 1e-5 界**，最差 source_0/EMT **2.29e-4 / 2.08e-4**；容差**未放宽、未改写成通过**。按 D-275，审计脚本新增 `--no-fail`：**仍输出 `status: FAIL`**，仅不再作为描述性计算的硬阻塞（`disclosure` 字段写明）。
- 单 case 受控线程证据（`infra/field_blas_probe_20261019/`）：仅 source_0/EMT 的 point 及其派生 smooth。1↔2 线程 1.0–1.6e-6；1/2 ↔ 8 线程 **2.07e-4**；存盘包 ≡ 8 线程（≤4.8e-7）。**边界**：该证据不覆盖 candidate/background/permuted，也不覆盖另外 43 个 case；不假定同批数值误差相消。冻结包在任意线程设置下与重拟合仍有 ~1.25–1.31e-5 每片 RMSE 的稳定残差（未定因，非线程、非库版本）。

## 3. 结果表（同一批 500 包；内部/来源分开；65/47 单位口径）

完整表：**`RESULTS_TABLES.md`**（覆盖与分母／主次读出／候选减四基线／来源分层／方向合取／限制）与 **`RESULTS_SOURCE_TABLE.tsv`**（四目标 × 6 来源，候选减同点）。表内方向标注遵循 C-01：**有符号半径差只作描述、不标方向**，范围一律用“逐片绝对值再按完整单位聚合”的读数；C-02：每行给出可评单位数，来源分层单独列出，器官外推单独标注。

描述性要点（**不构成结论**）：
- 局部细节 Pearson（内部臂，单位等权）：EMT local 0.6485 > point 0.6252、background 0.6465；IFN-γ local 0.5390 > point 0.4904、background 0.5299；OXPHOS local 0.6933 > point 0.6773 但 < background 0.6985；B_AXIS local 0.3283 < point 0.3449 < background 0.3635。
- q75 Jaccard：IFN-γ/EMT/OXPHOS 内部臂 local 略高于 point（+0.005～+0.027）；B_AXIS 略低（−0.004）。
- 置换对照（permuted）在全部四目标、两臂均**差于**候选，说明置换确实破坏了信息；但**整片背景对照在 B_AXIS/EMT/OXPHOS 与候选同量级或更好**。
- 来源臂：候选相对同点的局部细节差异在 source_2/4（肠）为负、source_0/5（肝/乳腺，器官外推）为正负参差；**未与内部臂合并**，按来源逐项见 `RESULTS_SOURCE_TABLE.tsv`。
- 方向合取（`readout_claim_summary.tsv`，仅描述）：内部臂 IFN-γ 与 EMT 的 local-detail/q75/质心/半径绝对偏差四项“同时优于四对照”为 True（EMT 的 RMSE 否）；B_AXIS 与 OXPHOS 为 False；来源臂逐目标可见。**这不是阈值、不是显著性、不含 PASS 判定**。
- 数值参考：主要读出配对差为 0.002–0.09 量级，远大于单 case 已知 2e-4 波动参考；但该参考只覆盖一个 point case，故仅把 |Δ|≤2e-4 的行标 `~num`，其余仍**未定**而非“通过”。

## 4. 有效 N、缺失与限制

- 内部臂：72 片 / **65 单位**（含 `singleton::TENX13` 两片合并为一单位、`source_key` 冲突标 `mixed`）；来源臂：53 片 / **47 单位** / 6 来源。
- 分母（`denominators.tsv`）：五个模型的主读出均为**完整 260（内部）／188（来源）单位组，0 部分、0 缺失**；`train_mean` 类常量基线不在本轮模型集内（本轮无该模型）。
- 身份：**确认 14 片 / 8 单位**，未知 **58 片 / 57 单位**；切片数不冒称患者数，未知不升级。
- 来源身份不足 19 片（含全部 6 片肺）仅留内部，不进来源分母；乳腺（source_5）与肝脏（source_0）为器官外推，单独标注。
- 受限读出：`rec_calibration_slope`/`rec_sd_ratio` 仅描述（C-01 方向不设）；有符号半径差仅描述；接口固定的 `rec_decay_half_pitch_delta` 不属本轮指标集。

## 5. 命令、计时、资源

```bash
cd /home/huyudi/013_spatial
bash scripts/run_field_raw_gene_spatial_20261019.sh evaluate     # exit 0
bash scripts/run_field_raw_gene_spatial_20261019.sh summarize    # exit 0
PYTHONNOUSERSITE=1 OMP_NUM_THREADS=2 CUDA_VISIBLE_DEVICES=-1 \
  /opt/anaconda3/bin/python infra/field_raw_gene_spatial_20261019/verify_results.py    # exit 0
/opt/anaconda3/bin/python infra/field_raw_gene_spatial_20261019/build_results_tables.py
/opt/anaconda3/bin/python infra/field_raw_gene_spatial_20261019/audit_frozen_agreement.py --no-fail
```
| 阶段 | 退出码 | 墙钟 | 进程 User | 进程 Sys | CPU% | 最大 RSS |
|---|---:|---:|---:|---:|---:|---:|
| fit（含 500 包，先前完成） | 0 | 1048 s | 5130.4 s | 812.9 s | 567% | 5.19 GB |
| audit（--no-fail，披露） | 0 | 6 s | — | — | — | — |
| evaluate | 0 | 25 s | 68.5 s | 2.5 s | — | 281 MB |
| summarize | 0 | 4 s | 3.4 s | 0.2 s | — | 196 MB |
| verify_results | 0 | 4 s | — | — | — | — |

说明：进程 CPU 取自 `/usr/bin/time -v`（含被回收子进程），**不是**把并行墙钟当成 CPU 分钟。磁盘：运行前 844.5 GB → 结束后 **827 GB 可用**（≫200 GB）。全程 CPU、本地数据；人工 TLS/annotations/incoming/HookNet 未读取。

## 6. 交付物与复核清单

- 实现（本轮最小修改）：`scripts/field_raw_gene_spatial_20261019.py`（metrics_root 隔离、evaluation_receipt 校验、C-01 方向集、claim 行方向与可评单位列）。
- 脚本：`infra/field_raw_gene_spatial_20261019/{audit_frozen_agreement.py,verify_results.py,build_results_tables.py}`。
- 结果：`RESULTS_TABLES.md`、`RESULTS_SOURCE_TABLE.tsv`、`prediction_metrics.tsv`(2500)、`paired_metrics.tsv`(2000)、`unit_metrics.tsv`(2240)、`unit_paired.tsv`(1792)、`macro_summary.tsv`、`macro_paired.tsv`、`readout_claim_summary.tsv`、`denominators.tsv`、`summary.json`、`REPORT.md`、`freeze_manifest.json`、`evaluation/*/evaluation_receipt.json`。
- 核对与披露：`audit_frozen_agreement.{json,tsv,log}`、`verify_results.{json,log}`、`full_run_stages.log`、`*_time.txt`、`fit_package_inventory.tsv`、`archive/evaluation_stale_partial_1843/`（污染证据，保留）。
- 未提交（工作区大文件）：500 个预测包（40 MB，清单在 `fit_package_inventory.tsv`）；其余本轮产物见 Git。

## 7. NOT_RUN

人工 TLS 与结构验证、独立 marker 隔离与重训、未知患者身份补齐、正式空间显著性、原场物理尺度与视野外推断、输入规模搜索、shuffle 参考（本轮 `do_shuffle=False`）。数值边界之外的一切“改善/失败”解读均属未定，等待独立复核。
