# P-stromal：成纤维/基质重塑活性

- **解释标签**：identity+state｜**最强替代解释**：基质丰度｜**保留限制**：quad敏感：部分残余是线性欠拟合
- **定义**：input 9基因（COL1A1, COL1A2, COL3A1, COL6A3, DCN, SPARC, IGFBP7, VIM, C1R）；readout 6基因（COL6A2, COL5A1, FBN1, BGN, LUM, FN1）；input/readout零重叠；冻结2026-09-18T16:33:33Z hash=675df82fde39b306
- **发现集**（Vanderbilt 30病人）：QC残余Moran中位0.221，30/30病人残余>0（全表见program_patient_effects.tsv）
- **外部三层**（分子corr / dM1 / dM2中位）：
| 队列 | 分子再现 | dM1 | dM2 |
|---|---|---|---|
| ST-CRC | 0.892 | 68.8% | 2.2% |
| USZ | 0.814 | 49.5% | 2.1% |
- **T5B TLS隐藏核心**：见masked_structure_results.tsv（mhc2在5/6可评片≥B轴；其余程序未超B轴）
- **原始/调整图**：cards_data 4切片（raw/resid_Q/resid_QC三层分数已存npz，作图见triage_log）
- **失败样本**：USZ肺片分子corr 0.26–0.34（epi）、stress USZ 3/8病人corr<0.25；ST-CRC个别病人dM1为负（M0已很好或input外推翻车，见incremental_prediction.tsv）
