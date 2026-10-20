# FA-04 逐节点复核（xai/grok-4.7，xhigh）

只审 FA-04：**真实 14 片的 A/B 分子核验**。重点是"有没有会让结果解释变错的错误"。

## 范围

- 代码：`field_analyzer/real_validation.py`、`field_analyzer/roles.py`、`field_analyzer/pipeline.py`（外部文库与场评分部分）
- 产物：`infra/field_analyzer_20261020/real_validation/{map_agreement.tsv,cross_field_validation.tsv,scale_agreement.tsv,summary.json}`
- 报告：`infra/field_analyzer_20261020/REPORT_FA04.md`

## 必答

1. **A/B 隔离是否对称且无后门**：`half_specs` 的两套 spec，外部文库是否都排除了 A∪B？B 是否有任何路径通过 witness、分母、score_genes 进入 A？`feature_mask` 在两次调用里是否一致？
2. **两半是否真的"独立"**：A 半与 B 半只是同一批文献集合的不相交子集，报告有没有把"两半一致"说成独立验证（应只是分子内部一致性）？
3. **图一致性 0.776 / 0.578 的口径**：Pearson 是否只在可见点算？`depth_corr` 用的是哪个文库（应是与目标无关的那个）？0.79 的 OXPHOS 深度相关是否被正确解释为混杂而不是生物信号？
4. **汇总单位**：`summarise()` 先单位内再单位间取中位吗？14 片 / 8 单位的不均衡是否影响报告里的中位数？
5. **尺度一致性**：`field_range` 与 FA-01 的拟合是同一实现吗？两半 range 差 0.72–4.44 pitch 的解读是否过度（应只说明不稳）？
6. **符号一致性的解释**：报告说"6 对符号 8/8 单位一致"，但共享文库/深度会诱导正相关；这句是否被写成了超出证据的重复性主张？
7. **数字一致性**：抽 6 个数与 TSV 对齐（0.776、0.926/0.793/0.738/0.578、0.049/0.109/0.161/0.220/0.040/0.022、−0.06/0.10/0.44/0.79、0.72–4.44）。

## 交付

写 `infra/field_analyzer_20261020/review/fa04_node_review.json`，`scope.read`/`scope.ran` 必须非空（未跑就写 `[]` 并在 `not_checked` 说明）。字段同 FA-01 模板（`node` 填 `FA-04`，findings 前缀用 `J`）。

不要修改仓库里除该 JSON 以外的任何文件。
