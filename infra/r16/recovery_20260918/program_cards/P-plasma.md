# P-plasma：浆细胞丰度

- **解释标签**：identity｜**最强替代解释**：Plasma轴（不可分）｜**保留限制**：只验分子再现，不争增量
- **定义**：input 11基因（IGHA1, IGKC, IGLC2, IGLC3, JCHAIN, IGHV3-30, IGHV3-35, IGKV1D-12, IGKV3OR2-268, DEPRECATED_ENSG00000211685, DEPRECATED_ENSG00000211890）；readout 6基因（IGHM, IGHA2, MZB1, SDC1, TNFRSF17, CD38）；input/readout零重叠；冻结2026-09-18T16:33:33Z hash=675df82fde39b306
- **发现集**（Vanderbilt 30病人）：QC残余Moran中位0.246，30/30病人残余>0（全表见program_patient_effects.tsv）
- **外部三层**（分子corr / dM1 / dM2中位）：
| 队列 | 分子再现 | dM1 | dM2 |
|---|---|---|---|
| ST-CRC | 0.596 | NA | NA |
| USZ | 0.659 | NA | NA |
- **T5B TLS隐藏核心**：见masked_structure_results.tsv（mhc2在5/6可评片≥B轴；其余程序未超B轴）
- **原始/调整图**：cards_data 4切片（raw/resid_Q/resid_QC三层分数已存npz，作图见triage_log）
- **失败样本**：USZ肺片分子corr 0.26–0.34（epi）、stress USZ 3/8病人corr<0.25；ST-CRC个别病人dM1为负（M0已很好或input外推翻车，见incremental_prediction.tsv）
