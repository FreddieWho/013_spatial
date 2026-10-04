# c_2固定范围诊断

[结果报告](../../docs/TLS_FORWARD_PATIENT_AND_C2_AUDIT_20261003.md)。contract.json绑定两源原始文件、评分缓存和固定预测SHA。complete.json记录原始counts全1742项复算、源标签/坐标/顺序核对、三mask查询/固定特征及1000μm固定读数预测重构的一致性。

source_comparison.tsv与query_support.tsv描述样本与资格。raw_score_by_annotation.tsv只检查标注区域读数，不是细胞比例或独立病理验证。R1000_query_diagnostics.tsv中的可见TLS统计仅用于审计，未作为模型输入；barcode_order是源矩阵行号。feature_contributions.tsv（实际文件名R1000_feature_contributions.tsv）是固定模型各维对正负平均分差的贡献，不是因果贡献。5个Frozen正查询来自同一TLS区域，不是5个独立结构。未翻转方向、改标签或改mask，成对制备差异原因仍未确定。
