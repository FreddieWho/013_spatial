# FA-13 交付：目标自身邻域之外的跨场增量（DELIVERED_PENDING_REVIEW）

日期：2026-10-23。合同：`docs/FA13_CONDITIONAL_RELATION_20261023.md`（D-284）。执行：Pi `opencode-go/deepseek-v4.1-flash`（thinking `max`）。复核：待独立 Grok `xhigh`。
**本报告只作描述性工程与统计交付；独立复核关闭前不用于任何 AF 判断。**

## 1 先读与等价性核对

按合同要求先核对 FA-12 / FA-05 / 旧 HF / 旧 FC 是否存在“目标自身邻域 + 另一场邻域”的等价嵌套比较。
逐项结论见 `equivalence_audit.md`，摘要：

| 来源 | 关键差距 | 判定 |
|---|---|---|
| FA-12 `spatial_lag` | 三个预测器各自单独评价，无 M0/M1 嵌套；且 `_prepare` 在掩膜前 residualise/zscore | NOT_EQUIVALENT |
| FA-05 `multifield_ridge` | 只用其它场，无目标自身邻域项 | NOT_EQUIVALENT |
| 旧 X/WX/VX（原始基因空间对照） | 同场同点+同场邻域，无跨场块 | NOT_EQUIVALENT |
| 旧 HF M0/M1/M2 | 不同队列，且审计记录隐藏表达泄漏 | NOT_EQUIVALENT，数值不用 |
| 旧 FC01/02/03 | 输入阶梯/字典探针，无嵌套自身邻域设计 | NOT_EQUIVALENT |

**没有完全等价的有效结果**，因此按合同实现并运行。

## 2 实现与固定口径

新增 `field_analyzer/conditional_relations.py`、`scripts/field_analyzer_conditional_20261023.py`、`tests/test_fa_conditional_relations.py`。
**未复用** `real_relations.section_scores` / `pipeline.build_inputs` 的全片路径（那会在掩膜前拟合标准化与深度系数）；只复用 `spatial_lag.holdout_mask` 生成与 FA-12 相同的 Q，其余全部从原始计数与观测掩膜重建。

- **数据**：FA-00 的 14 片 / 8 个研究内单位（`identity_review.tsv`），只读原始 counts、`in_tissue`、`obsm/spatial`；不读人工 TLS、不下载。
- **查询点 Q**：直接调用 FA-12 的 `spatial_lag.holdout_mask(n, 0.1, seed 0, section_id)`（按片确定），并在 `in_tissue` 内取；14 片共 2 198 个 Q 点（每片 32–273）。
- **场评分**：逐基因 `log1p(1e4 × count / 背景文库)`，背景文库为**排除 EMT+IFNG+OXPHOS 全部面板基因**后的 36 011 个基因计数；每场分数为其半边基因的均值。
- **基因半边**：`roles.build_field_splits(seed 0)`；特征半边再扣掉“该方向的核验并集”（共享基因按并集规则排除）。实测每片特征半边 EMT 133 / IFNG 131 / OXPHOS 133 个基因，核验半边 67 / 65 / 67。
- **只在 O 拟合**：深度残差系数（分片，O 行）、z 标准化（分片，O 行）、Ridge（训练单位 O 行）全部只用 O；Q 的表达/深度/QC 不进入任何特征或拟合。离线目标变换用同一 O 拟合系数，对 Q 应用（仅用于评价，不返回模型）。
- **特征**：k=8 个最近**观测**邻居的均值，**自身永远排除**；无有效邻居的点计数并报告（本轮 Q 侧为 0）。
- **模型**：M0 = 目标场特征半边的邻居均值；M1 = M0 + 预测场特征半边的邻居均值；带截距岭回归 `alpha=1.0`，同折同训练行、同 Q 评价。主口径 = 深度条件化；未条件化版本为敏感性，不参与判定。
- **判定**：预先固定为描述性一致性——主关系 EMT←IFNG、特征/核验半边 A/B、条件化版，要求 ≥6/8 单位 ΔR²>0 且同一要求在换半边后也成立；不设 RMSE 百分比门槛、不报确认性显著性。

## 3 主结果（单位等权，留一单位，同 Q）

| 关系 | 半边 | 口径 | M0 R²(中位) | M1 R²(中位) | ΔR²中位 | ΔR²>0 单位 | 稳定 |
|---|---|---|---|---|---|---|---|
| **EMT←IFNG（主）** | A→B | 条件化 | 0.4989 | 0.4893 | **−0.0019** | **2/8** | 否 |
| EMT←IFNG | A→B | 未条件化 | 0.5200 | 0.5128 | −0.0006 | 2/8 | 否 |
| EMT←IFNG | B→A | 条件化 | 0.5236 | 0.5266 | +0.0015 | 5/8 | 否 |
| EMT←IFNG | B→A | 未条件化 | 0.5685 | 0.5684 | −0.0000 | 4/8 | 否 |
| IFN←OXPHOS（深度敏感对照） | A→B | 条件化 | 0.3723 | 0.3713 | +0.0021 | 5/8 | 否 |
| IFN←OXPHOS | B→A | 条件化 | 0.3136 | 0.3175 | −0.0006 | 3/8 | 否 |
| IFN←OXPHOS | A→B | 未条件化 | 0.3612 | 0.3760 | +0.0102 | 7/8 | 是（仅敏感性口径） |
| IFN←OXPHOS | B→A | 未条件化 | 0.3641 | 0.3789 | +0.0078 | 7/8 | 是（仅敏感性口径） |

单特征参照（说明“单独弱于自身邻域、联合仍无增量”）：主关系 A→B 条件化下，same-only 中位 R² **0.499**、cross-only **0.165**、联合 M1 **0.489**；换半边 same-only 0.524、cross-only 0.167、联合 0.527。IFN←OXPHOS 的 cross-only 仅 0.034（条件化 A→B）。

逐单位（主关系 A→B 条件化）：ΔR² = +0.0007 / −0.0046 / −0.0075 / −0.0198 / −0.0020 / −0.0017 / +0.0003 / −0.0002；M1 的跨场系数在多数单位为负（−0.021 / −0.037 / −0.064 / +0.067 / −0.024 / −0.021 / −0.020 / −0.012），即加入跨场邻居没有提供可用信息。

**判读（按合同边界）**：在目标自身邻域已进入模型的条件下，另一场的可见邻域没有带来稳定条件增量；两个基因方向一致地接近 0（|ΔR²|≤0.002，均未达 6/8）。这不支持“跨场邻居在自身邻域之外提供额外预测”，也不否定“单独跨场弱于同场”之外的任何机制；无匹配空间零分布与独立结构证据，不报机制、不当确认性显著。

**深度敏感对照**：只有**未做深度条件化**时 IFN←OXPHOS 才有稳定表面增量（+0.010，7/8）；加入 O 上拟合的深度条件化后消失（+0.002，5/8）。这与“该关系历史深度敏感”一致，因此不计为生物证据，只作技术敏感性记录。

## 4 必要校验（先跑，全部通过）

`smoke_checks.json`：

| 检查 | 结果 |
|---|---|
| 小型正例：cross 单独 0.465 < same 单独 0.531 < 联合 0.975 | 通过（增量 +0.444，6/6 单位） |
| 反例（目标只依赖自身邻域）：ΔR² 中位 1.0e-4，4/6 为正 | 通过（无条件增量） |
| 毒化隐藏原始计数/深度（合成小片）：特征、校准参数、系数、Q 预测 **逐位相同**；离线目标改变 | 通过（`features_identical=True`、`predictions_identical=True`、`truth_changed=True`） |
| 同一毒化在**真实片** NCBI643（2 261 点 / 240 Q / 240 隐藏）上重做 | 通过（特征、校准、系数、预测全部逐位相同；truth_changed=True） |
| Q 掩膜与 FA-12 一致（同 section/seed/10%） | 单元测试通过 |
| 邻域排除自身、k=8、单位不跨折、基因半边互斥与共享基因排除、数组/条码对应 | 单元测试 7 项全部通过 |

测试命令与结果：`python -m pytest tests/test_fa_conditional_relations.py -q` → **7 passed in 0.57s**。
真实运行命令：`python scripts/field_analyzer_conditional_20261023.py smoke` 与 `... run`（后者墙钟 11.6 s，8 线程 CPU，`CUDA_VISIBLE_DEVICES=-1`）。

## 5 本地 Xenium / Visium HD 只读核对

产物 `local_asset_fit.tsv`（92 个归档行：Xenium 41、Visium HD 51），扫描为**有界只读**：zip 顺序本地头扫描（不依赖中央目录）、tar.gz 流式读取小成员、最多读 2×10⁵ 行细胞表；**不解包、不复制、不载入表达矩阵**。

- Xenium：41 个归档中 16 个有中央目录、19 个缺失（归档不完整/流式写），6 个未判定；35 个可读出基因轴（`gene_panel.json`），34 个可读坐标范围（细胞表 `x_centroid/y_centroid`）。基因轴规模 10 756–16 915；面板覆盖 EMT ≈100/200、IFNG ≈105/200、OXPHOS ≈45/200（如 Brain GBM 面板 OXPHOS 为 0）。
- Visium HD：51 个归档（binned outputs 6–26 GB）基因轴未读（`NOT_READ_BIG_ARCHIVE`，按不载入大矩阵的边界）；`_spatial.tar.gz` 只含图像与缩放信息，无坐标表，坐标范围记为空。
- 能支持：Xenium 面板内 EMT/IFN 的部分基因可做同样的“自身邻域 vs +跨场”检验；OXPHOS 覆盖面普遍不足。
- 不能支持：由这些归档确认患者/病例身份（metadata 只有 run/sample/cell-seg 名称）；把 ROI 当患者；把“Base/Addon 同病例连续切片”当作事实（未知保持未知）；Visium HD 的基因轴覆盖（未读）。
- 单位/ROI/连续切片三者分开记录，未用任何人工 TLS 标签。

## 6 资源与边界

- CPU 8 线程（`OMP_NUM_THREADS=8`），未用 GPU；真实运行墙钟 11.6 s；磁盘余量实测 1 486.1 GB（≥200 GB）。
- 未做（记录为 NOT_RUN）：匹配空间零分布/置换零模型、独立结构证据、跨平台或来源扩样、参数搜索（k、alpha、方向、场组合）、B 路线精修、任何机制或因果表述。
- 限制：8 个研究内单位少、片内 Q 点少（最小 32）；片不是患者；主关系是看过历史结果后的探索性跟进；Q 的预测场同点表达本轮不可读，因此不能声称已控制该未观测变量；深度条件化本身是片内可见表达的无监督校准（已报告），它会带走部分真实生物变化。
- 结论限制在**本次固定模型/观测条件**；不升级 AF-01/AF-02，不称“桥”“非局部连接”“因果”。

## 7 对假设的影响

无正向支持：本次固定表示（k=8、alpha=1、线性、10% 整点缺失、8 单位）下，目标自身邻域之外的跨场条件增量未出现。该负结果限制的是**这个模型与观测条件**，不是“所有空间关系不存在”；AF-01（跨场关系辨识）不因本节点升级，AF-03 的条件恢复能力也不因本节点获得支持。

## 8 产物清单

- `equivalence_audit.md`、`REPORT_FA13.md`（本文件）、`run_manifest.json`、`local_asset_fit.tsv`
- `unit_results.tsv`（64 行 = 2 关系 × 2 半边 × 2 口径 × 8 单位）、`section_results.tsv`（112 行）、`section_inputs.tsv`（14 片 Q/O/基因半边与背景基因数）、`summary.json`（含单特征参照）、`smoke_checks.json`
- 代码：`field_analyzer/conditional_relations.py`、`scripts/field_analyzer_conditional_20261023.py`、`tests/test_fa_conditional_relations.py`
