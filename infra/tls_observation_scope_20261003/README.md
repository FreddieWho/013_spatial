# 冻结设计的范围核查 D-187

[报告](../../docs/TLS_OBSERVATION_SCOPE_20261003.md)。cohort_scope.tsv是入口；all_component_scope.tsv保留每个切片TLS连通域×三mask，all_section_scope.tsv与all_unit_scope.tsv保留切片与患者/肿瘤样本单位。所有结构数都是切片标注实例，不是独立患者或跨片去重后的真实TLS数。component_exclusion_reasons.tsv与query_exclusion_reasons.tsv的分母不同，不可混加。visible_support_by_section_class.tsv是每片每类查询的组织支持，非患者等权统计。geometry_vs_model_support.tsv区分几何双类支持和来源留出训练资格；cross_mask_component_overlap.tsv按同片同连通域匹配，未配准跨片。必要直径上限限制核心遮蔽，不限制桥长度。

contract.json记录身份/几何/查询表SHA，complete.json记录计数核对与资格上限检查。没有读取表达、重选程序、改变模型或新造标签。
