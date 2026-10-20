# FA-00 节点完成复核（xai/grok-4.7，xhigh）

你是本项目的独立复核代理。本轮只做一件事：判断 **FA-00 已交付的工作里有没有会影响结果正确性的错误**。不是泛化审计，不是逐项审批；按实际风险做必要检查。

## 复核范围

- 代码：`field_analyzer/`（`observation.py`、`adapters.py`、`geometry.py`、`roles.py`、`pipeline.py`、`poisoning.py`、`cli.py`），提交 `cb7256c`
- 测试：`tests/test_fa_observation.py`、`test_fa_adapters.py`、`test_fa_geometry.py`、`test_fa_roles.py`、`test_fa_poisoning.py`
- 产物：`infra/field_analyzer_20261020/{REPORT_FA00.md, reuse_resolution.tsv, section_asset_check.tsv, unit_role_assertions.tsv, poisoning_report_ZEN36.json, poisoning_report_NCBI642.json, geometry_smoke/, run_manifest.json}`
- 合同：`docs/FA_CONTRACT_ADDENDUM_20261020.md`（权威，优先于施工包正文）；`docs/field_analyzer_killtest_20261020/05_COMMON_PROTOCOL.md`

## 必须回答的问题（逐条给结论）

1. **隐藏信息是否真的进不去**：`pipeline.build_inputs` 里有没有任何路径让隐藏行参与计算（文库、深度、检出、整片均值/方差、场评分、邻接权重、查询资格、以及任何归一化分母）？请指出具体行；若没有，说明你如何确认。
2. **毒化证据是否有效**：`poisoning.run_poisoning` 的比较对象是否覆盖了这些入口？"评价必须随动"的反向断言是否真的能发现"评价没读隐藏区"？负对照（`test_leaky_slide_statistic_would_be_caught`）是否有检出力？
3. **witness（B 组）隔离**：从 `roles.split_field_genes` 到 `InputSpec.feature_mask` 到外部文库分母，B 是否可能以任何代数方式进入 A 或进入评分？并集处理是否正确？
4. **坐标/轴/单位**：`adapters.read_h5ad_counts` 是否杜绝了静默转置？xy、barcode、`source_rows` 的对应关系是否被强制？方阵或形状巧合时会发生什么？单位是否被硬编码？
5. **几何**：孤立点是否确实为 NaN（无 self fallback）？最大边距是否会因数据而误判为"组织断开"？空洞与"被遮但可能有组织"是否严格区分？局部方向的不可测判定是否会漏判？
6. **单位与分折**：`unit_role_assertions.tsv` 的 8 单位/14 片、同单位同折、重复 mask 不计额外患者，是否与 `identity_review.tsv` 一致？有没有把 unknown 身份升级成独立患者？
7. **口径是否会误导后续 FA-01/03**：场评分用"外部文库排除 A∪B"是否引入了新的偏倚（例如分母与目标无关却与组织类型相关）？`ridge_predict` 作为不变性参照是否会被误读成模型能力？`unit_role_assertions` 的 5 折轮转是否会被误读成正式交叉验证？
8. **报告与产物是否一致**：`REPORT_FA00.md` 里的数字（14/14、22 数组、max|Δ|=0、Δ0.513、1.29–1.40 GB、438/4 测试）能否从产物复算？有没有夸大或漏披露？
9. **是否还有会让 FA-01 起步就错的遗留项**。

## 你可以做的检查

可以读代码、跑测试（`/opt/anaconda3/bin/python -m pytest tests/test_fa_*.py -q`）、跑 `python -m field_analyzer.cli poison --section ZEN36`、复算产物。不要修改仓库任何文件（除输出复核 JSON），不要读人工 TLS 标签，不要下载、不要 GPU。

## 交付

写 `infra/field_analyzer_20261020/review/fa00_completion.json`：

```json
{
  "reviewer": {"provider": "xai", "model": "grok-4.7", "thinking": "xhigh", "cwd": "/home/huyudi/013_spatial"},
  "commit_reviewed": "<git rev-parse --short HEAD>",
  "scope": {"read": [...], "ran": [...], "not_checked": [...]},
  "findings": [{"id": "F1", "severity": "BLOCKER|MAJOR|MINOR", "topic": "...", "evidence": "...", "required_fix": "..."}],
  "answers": {"1": "...", "2": "...", "...": "..."},
  "claims_that_do_not_hold": [],
  "allow_next_node": true,
  "blocking_items": [],
  "limitations": "..."
}
```

`allow_next_node` 为 true 的条件：没有 BLOCKER，且 MAJOR 都能在不改变结论范围的条件下于下一节点前修完。若发现会改变结果解释的问题（泄漏、坐标错位、单位跨折、汇总口径错误），置 false 并写清最小修复。

**只报你实际读到和实际跑出来的内容。** 不承诺"绝对无缺陷"；无法从代码判断的实现行为写 `NOT_ASSESSABLE`。
