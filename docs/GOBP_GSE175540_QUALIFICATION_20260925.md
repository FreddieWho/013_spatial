# GSE175540 患者级外部评分资格核查（2026-09-25）

## 结论

GSE175540 的标签来源可追溯到18位不同病人，18/18标签文件均与R-02登记哈希一致，表达矩阵、坐标文件和标签文件的样本编号也逐一对应。严格要求标签和表达矩阵完整对应、且病人同时有明确 TLS 与 NO_TLS 标签后，可用于两类比较的样本为16位病人。

45项候选的全部输入和读出基因符号都能在18份表达矩阵中找到；但14项候选至少包含一个重复映射到多个基因编号的符号。按不自行挑选或合并feature的严格规则，只有31项候选在18份矩阵里都能无歧义地对应到唯一feature。因此目前可确认的评分资格是 **16位病人 × 31项候选**。这只是数据资格，不是外部生物结果；本阶段没有计算任何表达分数、效应量或显著性。

## 标签与样本对应

| 检查 | 结果 | 处理边界 |
| --- | ---: | --- |
| R-02冻结外部角色 | 18/18位不同患者 | 全部患者级；0条具块级资格 |
| 标签文件哈希 | 18/18匹配登记值 | 原始标签文件未改动 |
| 样本编号对应 | 18/18标签、矩阵、坐标文件编号相同 | 三类文件按同一 GSM 编号配对 |
| 矩阵结构与计数检查 | 18/18通过 | 矩阵形状、稀疏索引一致；保存的UMI为正整数 |
| 明确 TLS 与 NO_TLS 两类 | 17份文件含两类 | 空白标签均记为未知，不当作 NO_TLS |
| 严格完整两类患者 | 16/18位 | 两个例外见下文 |

有两处必须保留的缺口：GSM5924043（`frozen_a_3`）有1个 TLS 点和8个 NO_TLS 点没有进入过滤后的表达矩阵，虽然这9个条码有组织坐标；它的明确标签表达覆盖率为1310/1319（99.32%），但按完整对应规则不纳入严格两类集合。GSM5924049（`frozen_b_18`）有78个 TLS 点、没有明确的 NO_TLS 点，另有1108个空白标签；空白不是阴性，因此该病人不参加 TLS/NO_TLS 两类比较。

## 候选基因覆盖

45/45项的冻结输入与读出符号均在所有18份矩阵的 `Gene Expression`、`GRCh38` feature 名称中出现，没有使用同义词或外部别名映射。问题在于部分名称对应多个 feature ID：`TBCE` 影响5项候选，`MATR3` 影响8项，`TMSB15B` 影响1项，合计14项候选。31项候选在18份矩阵中都有无歧义的一对一映射。

这一步没有挑选重复名称中的某一个编号，也没有把多个编号的计数相加。如此处理是为了避免静默改变冻结基因集合；后续若希望评分全部45项，必须先明确冻结按 exact symbol 聚合多个 feature ID 的规则，再重新核验实现。否则应把受影响候选标为不可测，不能让程序默认选中其中一个编号。

## 对外部验证范围的影响

GSE175540 的本地 GEO 系列元数据将其描述为肾细胞癌空间数据。45项候选的筛选使用了 ST-CRC 与 USZ；USZ 的本地来源记录明确包含3个肾癌和5个肺癌样本。因此 GSE175540 是新的患者队列，但不是未见过的癌种确认，也不是跨平台验证。它未来能提供的证据应限定为独立患者来源中的 TLS 标签关联；不能据此声称跨癌种或跨平台迁移。

这次只完成数据资格清点，H-01、H-03、H-05均无新生物证据，`docs/plan.md` 保持不变；真实空间 p/T2P2 仍为 `NOT_CALIBRATED`/未运行。本轮没有下载新数据或使用 GPU，剩余可用磁盘约1.3 TB，高于项目下限。

## 可复核文件

- [资格核查脚本](../scripts/qualify_gse175540_gobp.py)
- [运行合同](../infra/gobp_external_qualification_20260925/run_contract.json)
- [运行收据与哈希](../infra/gobp_external_qualification_20260925/receipt.json)
- [逐患者标签和矩阵 QC](../infra/gobp_external_qualification_20260925/sample_qc.tsv)
- [逐候选基因覆盖表](../infra/gobp_external_qualification_20260925/program_coverage.tsv)
- [候选级资格汇总](../infra/gobp_external_qualification_20260925/program_qualification_summary.tsv)
- [重复符号与 feature ID 明细](../infra/gobp_external_qualification_20260925/feature_symbol_collisions.tsv)
- [USZ 样本来源记录](../infra/structure-registry/provenance/USZ_ZENODO_14620362_TLS_VISIUM.md)

重现命令：`python scripts/qualify_gse175540_gobp.py freeze`（首次冻结）后运行 `python scripts/qualify_gse175540_gobp.py run`。完成后再次运行 `run` 会校验现有输入与输出哈希，不会重算或改写数据。
