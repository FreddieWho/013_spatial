# 共同查询的尺度比较 D-188

[报告](../../docs/TLS_PAIRED_MASK_COMPARISON_20261003.md)。all_program_paired_summary.tsv含1742联合模型+geometry_depth/train_selected_GO，三个验证/来源组合共5232行；径向/角向单臂不在此比较内。all_unit_paired_metrics.tsv.gz与all_section_paired_metrics.tsv.gz分别保留患者/样本及切片指标。共同正负查询、同源spot_index和标签一致性由common_query_identity.tsv追溯，common_query_support.tsv保留缺类和未匹配支持。GSE在来源留出与患者留出各出现一次，不能把记录行数当独立查询/患者数。

contract.json与complete.json绑定原预测SHA及完整性。原完整人群指标并列在汇总中，共同子集不是新的盲验或纯距离因果效应；不重训、不翻转方向。c_2 Frozen因500μm无正例不能进入共同可评子集，其完整1000μm失败保留。所有正式p留空。
