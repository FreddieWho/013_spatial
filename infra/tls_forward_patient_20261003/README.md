# 全量前向患者单位汇总（D-186）

同一患者、同条件的可测切片取均值，再患者间取中位数。USZ保留肿瘤样本单位；不混长度、半宽、路径、mask或端点类型。无候选/缺值不是零。

三组结果：bins（晕/边分箱）、contrasts（晕/边对比）、bridges（桥）。每组都有`*_unit_index.tsv`与`*_unit_values.npz`、`*_summary_index.tsv`与`*_summary.npz`。数组的set_ids固定全部1691 GO及51对照，行由对应索引给出。unit_index.source_rows指向旧integrated表的零起始行，所有旧标签/校准/参考支持可追溯。

unit_values每个metric为患者内可测切片均值，metric_n_sections为实际贡献数；summary里metric_median是患者中位数，metric_old_section_median为旧切片中位数，metric_n_units为可测患者数，metric_n_complete_units要求该患者全部登记切片均贡献。n_positive_units/n_negative_units仅为相对数值零的方向计数，只对有符号效应适合作生物描述；均值、标准差、参考分位的符号计数不是发现数。参考z/分位的汇总不是联合检验或p。

`all_program_bidirectional_patient_screen.tsv`含全部1742项×两前向mask，与D-185来源留出的患者均值对齐；前向核心mask与反向圆形查询mask不是同一个几何定义。`fixed_readout_examples.tsv`只展示已冻结B轴及两个既有GO例子，不代替全量。`contract.json`绑定输入与身份SHA；`complete.json`及`tests.txt`给出完整性/聚合检查。原始源资格见旧`infra/tls_multiscale_gobp_20261002/integrated/halo_edge_source_completion.tsv`与bridge calibration目录；没有几何行的源不因未列入condition索引而变成阴性。
