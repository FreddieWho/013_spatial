# 患者分组验证

[结果报告](../../docs/TLS_PATIENT_VALIDATION_20261003.md)。`contract.json`固定D-185及57折；`complete.json`记录训练/测试分组和预测/指标SHA。3个c_2联合留出模型在models/；54个精确等价折引用旧路径。`all_program_summary.tsv`为全部程序入口，`all_unit_metrics.tsv.gz`与`all_section_metrics.tsv.gz`保留患者/切片层，`c2_pair_all_program_comparison.tsv.gz`保留成对失败。固定读数无事后新筛选。GSE患者等权，USZ肿瘤样本等权；训练仍原切片/类别权重。内部验证不是桥、空间校准或定位确认。
