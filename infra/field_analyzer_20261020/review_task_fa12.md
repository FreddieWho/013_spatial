# FA-12 节点完成复核（xai/grok-4.7，xhigh）

只审这一问：邻居的另一个场，有没有在留出单位上超过邻居自己的同一场。结论是 0/8，`upgrade_af01=false`。请查有没有会把这个否定翻过来、或把否定说得过重的错误。

## 范围

- 代码：`field_analyzer/spatial_lag.py`、`field_analyzer/cli.py::cmd_spatial_lag`
- 产物：`infra/field_analyzer_20261020/spatial_lag/unit_results.tsv`、`summary.json`
- 报告：`infra/field_analyzer_20261020/REPORT_FA12.md`
- 对照：`infra/field_analyzer_20261020/real_holdout/summary.json`（拟合方向，不是 REPORT_FA08B 里写反的箭头）

## 必答

1. 掩膜口径是否真的不用被掩点的目标值做邻居特征？训练是否不用留出单位？
2. 同场邻居是否用目标场的 B 半，跨场邻居是否用预测场的 A 半？两边是否同一张 kNN 图、同一批掩点？
3. `0/8` 和“最近差 0.039”能否从 `unit_results.tsv` 复算？两个口径是否都是 0/8？
4. 升级规则是否按预先写死的执行：两个口径都要 ≥6/8，并且同时超过同场邻居和两个旋转零对照？有没有把“高于旋转零对照”写成“有新的空间机制”？
5. 报告的拟合方向是否和 `real_holdout/summary.json` 的 `目标->预测变量` 一致，而不是 REPORT_FA08B 的箭头？
6. 还有没有泄漏、单位加权或把同点共变说成空间关系的问题？

## 交付

写 `infra/field_analyzer_20261020/review/fa12_review.json`，`scope.read` 和 `scope.ran` 不能为空：

```json
{"reviewer": {"provider":"xai","model":"grok-4.7","thinking":"xhigh","cwd":"/home/huyudi/013_spatial"},
 "node": "FA-12", "commit_reviewed": "<git rev-parse --short HEAD>",
 "scope": {"read": [], "ran": [], "not_checked": []},
 "findings": [], "claims_that_do_not_hold": [],
 "allow_node_conclusion": true, "limitations": ""}
```

只写实际读到或跑到的内容。不要改这个 JSON 以外的文件。
