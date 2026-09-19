# P-smmhc：平滑肌/肌成纤维收缩态

- **解释标签**：identity+state｜**最强替代解释**：平滑肌丰度｜**保留限制**：小而干净，但M2=0
- **定义**：input 3基因（ACTA2, FLNA, MYL9）；readout 5基因（TAGLN, CNN1, DES, MYH11, TPM2）；input/readout零重叠；冻结2026-09-18T16:33:33Z hash=675df82fde39b306
- **发现集**（Vanderbilt 30病人）：QC残余Moran中位0.151，30/30病人残余>0（全表见program_patient_effects.tsv）
- **外部三层**（分子corr / dM1 / dM2中位）：
| 队列 | 分子再现 | dM1 | dM2 |
|---|---|---|---|
| ST-CRC | 0.782 | 52.6% | 0.1% |
| USZ | 0.590 | 15.0% | 3.4% |
- **T5B TLS隐藏核心**：见masked_structure_results.tsv（mhc2在5/6可评片≥B轴；其余程序未超B轴）
- **原始/调整图**：cards_data 4切片（raw/resid_Q/resid_QC三层分数已存npz，作图见triage_log）
- **失败样本**：USZ肺片分子corr 0.26–0.34（epi）、stress USZ 3/8病人corr<0.25；ST-CRC个别病人dM1为负（M0已很好或input外推翻车，见incremental_prediction.tsv）
