# FA-13 交付修正说明（2026-10-23）

状态：**DELIVERED_PENDING_REVIEW / BLOCKED_REVIEW**。独立 Grok `xhigh` 新代理 `82d66e60-3f71-4310-8840-14777ad3151f` 在开始后遭 xAI 403（额度/订阅不足），**未完成任何独立复核**。以下四项为主代理对交付表述的核对要求，不是“Grok 审核发现”；本人的核对也不替代独立审核。所有科学结果仍不得用于结论。

本轮不做真实重拟合以外的扩样；第 1 项为纯措辞与命名修正（数值未变），第 2 项为解析错误纠正并重生成资产表（主拟合未重跑，但为保持 `summary.json` 与代码一致做了同参数重跑，数值逐位一致），第 3 项为边界内补查，第 4 项为可比性/边界表述修正。

## 1 删除“≥6/8 为预先固定通过标准”的不实来源

- **问题**：报告称 `>=6/8` 是预先固定判定；实际冻结合同 `3a07f06` 只写“报告效应大小与单位间一致性；无稳定增量则冻结该分支”，**没有任何数值门槛或 6/8 规则**。
- **修正**：
  - `field_analyzer/conditional_relations.py`：删除 `STABLE_UNITS`/`stable`/`decision_rule`，改为 `DESCRIPTIVE_MAJORITY_UNITS`；`summarise()` 输出 `descriptive_majority_positive`、`descriptive_note`（明示“executor-added、未列入冻结合同、非预注册、非 AF 通过标准”）与 `majority_units`。
  - `scripts/field_analyzer_conditional_20261023.py`：`settings.majority_units` 取代 `stable_units`。
  - `summary.json`：重跑后不含 `stable`/`decision_rule` 键；`n_positive`、`median_delta_r2`、`median_r2_m0/m1`、`unit_delta_r2` 等原始数字**逐位未变**（重跑前后自动比对：0 处差异）。
  - `REPORT_FA13.md`：判定口径段改写为“合同无数值标准，只报告”；表头“稳定”改为“多数为正（描述性，非合同）”，逐行计数按 `summary.json` 修正（2/8、2/8、5/8、4/8、5/8、3/8、7/8、7/8）。
  - `run_manifest.json`：limitations 明示同一条款。
- 未新增任何生物学或 RMSE 百分比门槛。

## 2 Xenium 基因轴量级错误的纠正（解析错误）

- **问题**：上一版报告称 Xenium 基因轴 10 756–16 915、EMT≈100/IFNG≈105/OXPHOS≈45。实际提取器把 `gene_panel.json` 中所有 `name` 字段（含公司名、面板名等非基因字段）都当基因名，导致数量虚高；且**面板设计目标数 ≠ run 实际检测轴**。
- **交叉核对（只读少量 metadata）**：
  - `data/derived/xenium_colon_base/cell_feature_matrix.h5`：只读 `matrix/features/name` → **541 个输出特征名条目**（基因/控制探针分类待核验）；`gene_panel.json` 设计目标 **325**；`experiment.xenium` 声明 predesigned 325 / custom 0。
  - `data/derived/xenium_colon_addon`：检测轴同为 **541**；设计目标 425（325+100 custom）。
  - 检测轴上的面板覆盖：base EMT 10、IFNG 6、OXPHOS 1；addon EMT 28、IFNG 10、OXPHOS 2（分母 200）。
- **修正**：
  - 提取器重写为 `_read_zip_panel_metadata`：仅读 `experiment.xenium`、`gene_panel.json`（按 `type.descriptor == "gene"` 取 `type.data.name`）、`metrics_summary.csv`、`cells.csv.gz`；归档内**不再声称检测轴**，统一 `UNVERIFIED_DETECTED_AXIS`。
  - 新增已解包交验行（`platform=xenium_unpacked`，`detected_axis=VERIFIED_H5_FEATURES`），证据即上面的 541/325/425 与覆盖数。
  - `local_asset_fit.tsv` 以同一有界范围重生成（94 行：归档 92 + 交验 2）：`detected_axis` 分布 = `VERIFIED_H5_FEATURES` 2 / `UNVERIFIED_DETECTED_AXIS` 35 / `NOT_CHECKED` 57（后者为 Visium HD 51 + 6 个 post-Xenium CytAssist 的 tar.gz）。
  - `REPORT_FA13.md` §5 重写；主代理进一步收窄解释：只核实了两份已解包输出，541 是特征名条目数，不直接等同基因数；较少的程序基因命中限制原定义的分子支持，但不能仅凭覆盖比例判定不可分析，也不能外推为全部本地 Xenium 不可用。其余归档轴未核验，当前清单不用于选数据。
- 可追踪证据：`local_asset_fit.tsv` 的两行 `xenium_unpacked`（含 `n_target_entries`、`panel_design_gene_targets_n`、`panel_declared_targets_n`、`detected_axis_n`、`detected_panel_hits`、`experiment_*` 字段）与本节数字一一对应。

## 3 Visium HD 基因轴：NOT_CHECKED（不解释成“不覆盖”）

- 在 `data/derived/` 与 `data/other_sources/10x_genomics/` 内做有界搜索（`features.tsv*`、`filtered_feature_bc_matrix.h5`、`gene_panel.json`，深度 4）：命中的 `filtered_feature_bc_matrix.h5` 全部属于 `data/derived/r04/st_crc_cms/`（Visium ST-CRC，不是 Visium HD）；未发现任何已解包的 Visium HD features。
- HD 的基因轴只在 6–26 GB 的 `_binned_outputs.tar.gz` 内；`_spatial.tar.gz` 仅含图像/scalefactors。
- 因此 `local_asset_fit.tsv` 中 51 行 HD 一律记 `detected_axis=NOT_CHECKED`，报告明示“未读 ≠ 不覆盖”。未解包、未下载、未扩大范围。

## 4 与 FA-12 的可比性及 AF-03 边界

- `REPORT_FA13.md` 与 `equivalence_audit.md` 新增明示：FA-12 的 `same_neighbour` 用的是**目标 witness（核验半边）在 O 的值**；本轮为基因内部核验改用**特征半边**、核验半边只作监督/离线评价。两者不是同一输入预算，故本轮结果**不能**表述为“超过/未超过 FA-12 同预算最强参照”。
- 描述目标改为：**在当前限制的基因输入预算内**，检查目标自身邻域之外是否还有跨场条件增量。
- 新增明示：本轮是**随机整点 10% 留出**，不是连续区域缺失；即使出现增量也不能直接升级 AF-03 的连续区域恢复，只能作为后续设计线索。

## 未关闭范围

- 独立 Grok `xhigh` 复核：**BLOCKED_REVIEW**（xAI 403），未进行；本文件与代码实现均未经独立复核。
- Visium HD 基因轴：`NOT_CHECKED`（需解包/读取 binned outputs，超出本轮边界）。
- Xenium 归档检测轴：`UNVERIFIED_DETECTED_AXIS`（仅两份已解包输出可核实）。
- 主拟合结果本身未变：主关系 EMT←IFNG（A→B、条件化）M0 0.4989 / M1 0.4893、ΔR² 中位 −0.0019、2/8 单位为正；换半边 +0.0015、5/8；对照 IFN←OXPHOS 条件化 +0.0021、5/8，未条件化 +0.0102、7/8（仅敏感性）。

## 主代理结果解释补充

上述范围修正属于交付表述核对，不是独立代码审核。主拟合与两个基因方向的数值未改；指定 Grok 的完成复核仍被 403 阻塞。
