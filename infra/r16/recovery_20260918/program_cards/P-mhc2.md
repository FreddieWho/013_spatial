# P-mhc2：抗原呈递活性（APC）

- **解释标签**：state｜**最强替代解释**：B/Mye丰度｜**保留限制**：唯一非丰度候选；TLS超B轴线索
- **定义**：input 2基因（CD74, HLA-DRA）；readout 5基因（HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DQA1, CD86）；input/readout零重叠；冻结2026-09-18T16:33:33Z hash=675df82fde39b306
- **发现集**（Vanderbilt 30病人）：QC残余Moran中位0.139，30/30病人残余>0（全表见program_patient_effects.tsv）
- **外部三层**（分子corr / dM1 / dM2中位）：
| 队列 | 分子再现 | dM1 | dM2 |
|---|---|---|---|
| ST-CRC | 0.798 | 54.3% | 9.4% |
| USZ | 0.723 | 21.4% | 8.2% |
- **T5B TLS隐藏核心**：见masked_structure_results.tsv（mhc2在5/6可评片≥B轴；其余程序未超B轴）
- **原始/调整图**：cards_data 4切片（raw/resid_Q/resid_QC三层分数已存npz，作图见triage_log）
- **失败样本**：USZ肺片分子corr 0.26–0.34（epi）、stress USZ 3/8病人corr<0.25；ST-CRC个别病人dM1为负（M0已很好或input外推翻车，见incremental_prediction.tsv）

- **L-009修订（2026-09-19，D-134）**：细窗963配对head-to-head，mhc2 hit 0.606 vs B轴0.672（McNemar p=3.3e-08）；T5B的超B轴结论不成立。本程序退回场活性描述，不做TLS定位器主张。
