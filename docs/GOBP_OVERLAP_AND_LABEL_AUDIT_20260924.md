# GO BP 候选重叠与本地结构标注审计（2026-09-24）

## 结论

全量 GO BP 阶段留下的45个分子预测候选不是45项独立机制发现。它们的基因集合整体较大，部分程序高度重叠；但没有一个输入基因出现在全部45项中。这个清点说明需要控制候选间重复，不能替代单基因基线，因此还不能据此排除某些程序由少数强预测基因驱动。

本地结构标签比此前待办的概述更完整：R-02登记的80条可审计标签来源均能在本地找到并通过 SHA-256 核对，且都映射到冻结外层切分角色。GSE175540 有18位患者对应的外部 TLS 标签，可按冻结角色评估患者级外部验证；18条记录都没有组织块编号，不能支撑块级独立性主张。

## 45个候选的基因重叠

按上一阶段 `candidate_context.tsv` 中冻结的 `molecular_robust_gt20=true` 选择45项，并从冻结代表定义读取 input/readout 基因。每项覆盖100–1023个检测基因（中位数499）；input 基因数中位数250，readout 中位数249。45项共有11,746个 input 基因出现次数、3,835个不同基因；readout 有11,717次、3,867个不同基因。按每项 input 与 readout 的并集计，共有4,818个不同基因。

| 对比集合 | 990对程序的 Jaccard 中位数 |
| --- | ---: |
| input 基因 | 0.0287 |
| readout 基因 | 0.0283 |
| 每项 input∪readout | 0.0601 |

对每项 input∪readout 集合，121/990对的 Jaccard≥0.25，42/990对≥0.50。

input 中最常见的基因是 TGFB1（21/45项）；IL1B、TNF 各18项，HMGB1、IGF1 各16项。没有单个基因覆盖全部候选。多数程序对之间的重叠较低，但42对的 input∪readout Jaccard 至少为0.50，说明候选中存在明显重复。重叠计数不等于通路独立性检验，也没有按 GO 父子关系合并术语。

所以，现阶段可以说“有45项达到上一阶段分子预测筛选标准”，不能说“发现了45条独立生物机制”，也不能仅凭这次集合重叠断言预测不由单个基因驱动。使用相同冻结切分和评分的单基因基线仍未运行。

## 本地结构标签与资格

| 冻结谱系 | 可审计标签来源 | 结构类别 | 物理样本 / 患者数 | R-02角色 | 文件哈希核验 | 块级资格 |
| --- | ---: | --- | ---: | --- | ---: | --- |
| HTAN Vanderbilt CRC | 42 | TLS、肿瘤—基质边界 | 36 / 24 | training | 42/42 | 42/42 eligible |
| GSE175540 | 18 | TLS | 18 / 18 | external_validation | 18/18 | 0/18 eligible |
| ST_CRC_CMS | 12 | 肿瘤—基质边界 | 12 / 6 | internal_validation | 12/12 | 12/12 eligible |
| TLS_VISIUM_USZ | 8 | TLS | 8 / 8 | external_validation | 8/8 | 8/8 eligible |
| **合计** | **80** | — | — | — | **80/80** | — |

来源数、物理样本数和患者数是不同计数单位。GSE175540 的18条来源都有冻结患者链接和 `external_validation` 角色；`block_id` 为空，`block_level_eligible=no`。这保留了 R-02 定义的患者级外部验证用途，但不允许把18个样本升级成18个组织块的独立重复。

GSE175540 本地有23个逐点 TLS 标注 CSV.gz、24个表达矩阵 H5 和24个 tissue-position 文件。R-02只将其中18条来源审计为可重放标签；另外3个标注文件只有汇总型 T_agg 信息、2个零 TLS 样本不产生实例行，另有1个样本没有沉积标注。进入新的外部评分前，应先核对这18个样本的基因覆盖、表达矩阵与标签的样本键对应及必要 QC。该队列为肾癌，与主 CRC 队列不同；结果应按跨癌种患者级背景验证解释。

## 对假设和下一步的影响

本次清点没有新增模型证据，不改变 `docs/plan.md` 中任何假设的支持状态，也不解除真实空间检验仍为 `NOT_CALIBRATED` 的状态。它收紧了45项候选的解释方式，并确认 GSE175540 可在冻结 R-02 患者级角色下进行资格检查。

下一步先用相同冻结切分和评分运行45项候选的训练侧单基因基线，量化其相对最佳单基因预测的增量；另对GSE175540完成18个可审计样本的基因覆盖、样本键对应及必要 QC，再决定是否进入外部评分。两项工作均不需要新下载或 GPU。本轮没有做模型拟合、没有添加外部数据；剩余磁盘约1.39 TB。

## 可复核文件

- 重算程序：[audit_gobp_overlap_labels.py](../scripts/audit_gobp_overlap_labels.py)
- 输入与配对明细：[go45_programs.tsv](../infra/gobp_followup_20260924/go45_programs.tsv)、[go45_pairwise_overlap.tsv](../infra/gobp_followup_20260924/go45_pairwise_overlap.tsv)
- 基因频次：[go45_gene_frequency.tsv](../infra/gobp_followup_20260924/go45_gene_frequency.tsv)
- 标签来源及队列汇总：[gt_sources.tsv](../infra/gobp_followup_20260924/gt_sources.tsv)、[gt_cohort_summary.tsv](../infra/gobp_followup_20260924/gt_cohort_summary.tsv)、[gt_lineage_summary.tsv](../infra/gobp_followup_20260924/gt_lineage_summary.tsv)
- GSE175540 本地资源清单：[gse175540_local_inventory.tsv](../infra/gobp_followup_20260924/gse175540_local_inventory.tsv)
- 输入/输出 SHA-256 与计数收据：[receipt.json](../infra/gobp_followup_20260924/receipt.json)

重现命令：`python scripts/audit_gobp_overlap_labels.py`。它只读取仓库内冻结 GO 定义、候选表与 R-02 注册文件，并核验标签来源文件哈希。

## 2026-09-25 范围更正

本报告把 GSE175540 的后续用途概括为“跨癌种”背景，表述过宽。45项候选使用 ST-CRC 与 USZ 数据筛选；USZ 来源记录说明其中包括3个肾癌和5个肺癌样本。GSE175540 是肾细胞癌队列，因此它提供独立患者队列的外部标签，但不是未见癌种验证，也不是跨平台验证。更完整的样本键、标签和基因覆盖限制见[资格核查报告](GOBP_GSE175540_QUALIFICATION_20260925.md)。
