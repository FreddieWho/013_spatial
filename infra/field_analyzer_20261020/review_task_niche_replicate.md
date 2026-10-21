# FA-11 增补复核：跨片重复 + P1 关闭确认（xai/grok-4.7，xhigh）

你上一轮审到 `fdf8eac`，判 `allow_node_conclusion: true`，给出 P1（旋转零对照偏低不能写成已证明的贴边伪影），并明确**跨片 0.998 不在该轮通过范围内**。本轮只审这一块 + P1 是否关闭。

## 范围

- 代码：`field_analyzer/niche_graph.py::argmax_niches`、`cli.py::cmd_niche_replicate`
- 产物：`infra/field_analyzer_20261020/niche_replicate/{section_stats.tsv,enrichment_by_section.tsv,summary.json}`
- 报告：`REPORT_FA11.md` §4b 与 §3 纠正 2（本次改动）
- 数据：`data/derived/xenium_colon_base`、`data/derived/xenium_colon_addon`（各 541 基因 panel，两套 panel 有 100 个基因不同）

## 必答

1. **P1 是否关闭**：纠正 2 现在是否只写"观测低于旋转零对照、原因未单独验证、不能读成真实结构更弱"？
2. **跨片比较是否同口径**：两片是否只用**共有基因**评分？生态位是否用不依赖拟合的规则（argmax + MIXED 阈值）从而可跨片对齐？`MIXED` 类是否会造成两片富集矩阵的机械相关？
3. **子采样**：两片各取 10 万细胞是否独立抽样、是否影响富集值？两片统计量 4.022/4.022 几乎相同，这是巧合还是可能的实现问题（例如两次都用了同一 seed 导致相同子集？）？
4. **零对照**：两片各自做了置换零对照（z=+385/+287），零对照是否每片单独生成、没有被复用？
5. **0.998 与 42/42**：是否能从 `enrichment_by_section.tsv` 复算？42 对是 (7×7−7)/1 还是包含对角线（应为 7×7−7=42）？
6. **限定**：报告是否已写明这两片是**同一病例**的连续切片、不是跨病人？有没有把 0.998 写成"跨病人可重复"？
7. 还有没有别的过度表述（尤其把"共定位"说成"连接"）？

## 交付

写 `infra/field_analyzer_20261020/review/niche_replicate_review.json`，`scope.read`/`scope.ran` 非空：

```json
{"reviewer": {"provider":"xai","model":"grok-4.7","thinking":"xhigh","cwd":"/home/huyudi/013_spatial"},
 "node": "FA-11", "commit_reviewed": "<git rev-parse --short HEAD>",
 "scope": {"read": [...], "ran": [...], "not_checked": [...]},
 "p1_closed": true,
 "findings": [], "claims_that_do_not_hold": [],
 "allow_cross_section_claim": true, "limitations": "..."}
```

只报实际读到/跑到的内容；可跑 `niche-replicate` 复算。不要改仓库里除该 JSON 以外的文件。
