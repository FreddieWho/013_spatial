# P-stress：即早应激反应

- **解释标签**：state/technical｜**最强替代解释**：操作梯度未排除｜**保留限制**：USZ分化；技术标签未摘
- **定义**：input 2基因（FOS, DEPRECATED_ENSG00000170345）；readout 6基因（JUN, JUNB, EGR1, IER2, DUSP1, ATF3）；input/readout零重叠；冻结2026-09-18T16:33:33Z hash=675df82fde39b306
- **发现集**（Vanderbilt 30病人）：QC残余Moran中位0.291，30/30病人残余>0（全表见program_patient_effects.tsv）
- **外部三层**（分子corr / dM1 / dM2中位）：
| 队列 | 分子再现 | dM1 | dM2 |
|---|---|---|---|
| ST-CRC | 0.488 | 59.8% | 15.8% |
| USZ | 0.330 | 20.1% | 14.5% |
- **T5B TLS隐藏核心**：见masked_structure_results.tsv（mhc2在5/6可评片≥B轴；其余程序未超B轴）
- **原始/调整图**：cards_data 4切片（raw/resid_Q/resid_QC三层分数已存npz，作图见triage_log）
- **失败样本**：USZ肺片分子corr 0.26–0.34（epi）、stress USZ 3/8病人corr<0.25；ST-CRC个别病人dM1为负（M0已很好或input外推翻车，见incremental_prediction.tsv）
