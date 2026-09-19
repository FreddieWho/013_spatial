# P-epi：肠上皮/杯状分化态

- **解释标签**：identity+state｜**最强替代解释**：Epi丰度本身｜**保留限制**：分化梯度 vs 纯丰度不可分（Visium分辨率限制）
- **定义**：input 16基因（EPCAM, CEACAM5, CEACAM6, CLDN4, CLDN7, KRT8, KRT18, MUC13, TFF3, ELF3, FXYD3, LGALS4, ANXA2, PHGR1, LGALS3BP, GPX2）；readout 6基因（MUC2, MUC4, TFF1, CLDN3, KRT19, CDH1）；input/readout零重叠；冻结2026-09-18T16:33:33Z hash=675df82fde39b306
- **发现集**（Vanderbilt 30病人）：QC残余Moran中位0.340，30/30病人残余>0（全表见program_patient_effects.tsv）
- **外部三层**（分子corr / dM1 / dM2中位）：
| 队列 | 分子再现 | dM1 | dM2 |
|---|---|---|---|
| ST-CRC | 0.831 | 44.0% | -0.1% |
| USZ | 0.582 | 72.2% | 7.5% |
- **T5B TLS隐藏核心**：见masked_structure_results.tsv（mhc2在5/6可评片≥B轴；其余程序未超B轴）
- **原始/调整图**：cards_data 4切片（raw/resid_Q/resid_QC三层分数已存npz，作图见triage_log）
- **失败样本**：USZ肺片分子corr 0.26–0.34（epi）、stress USZ 3/8病人corr<0.25；ST-CRC个别病人dM1为负（M0已很好或input外推翻车，见incremental_prediction.tsv）
