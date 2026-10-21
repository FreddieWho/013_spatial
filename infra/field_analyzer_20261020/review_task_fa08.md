# FA-08A / FA-08B 节点完成复核（xai/grok-4.7，xhigh）

两个真实材料跟进节点。重点：**有没有会让结论变错的错误**。特别盯 FA-08A 的零对照是否公平、FA-08B 的留出是否泄漏。

## 范围

| 节点 | 代码 | 产物 | 报告 |
|---|---|---|---|
| FA-08A 真实位移 | `real_relations.py`（`displacement_curve`/`displacement_rows`/`summarise_displacement`）、`cli.py::cmd_real_displacement` | `real_displacement/{displacement.tsv,summary.json}` | `REPORT_FA08A.md` |
| FA-08B 留一单位 | `real_relations.py`（`louo_holdout`/`zscore`/`_fit_linear`/`_r2`）、`cli.py::cmd_real_holdout` | `real_holdout/{louo.tsv,summary.json}` | `REPORT_FA08B.md` |

## 必答

**FA-08A**
1. **零对照是否公平**：旋转目标场 137° 后，匹配容差/最小配对/搜索范围是否与真实工况一致？旋转是否无意中改变了点密度或重叠（使零对照更容易/更难出大峰）？
2. 真实峰"贴边界"的计数是否与 TSV 一致（7–8/14）？边界判定阈值 11.9 是否合理？
3. **A/B 两半不一致**的结论是否给出可复算证据（哪几片、什么方向）？
4. 报告说"与 OXPHOS 无关的三对峰在零位移、contrast=0.000"——这与 TSV 一致吗？这是否意味着估计器有鉴别力、只是对宽梯度场失效？
5. 方向约定（target = source 平移 d）是否与代码一致，报告数值是否已按该约定？
6. "真实位移不可用"是否被证据充分支持，还是只能说"幅值不可信、contrast 仍有信息"？

**FA-08B**
7. **有没有泄漏**：拟合用 A 半基因、评估用 B 半基因；片内 z-score 是否在训练/测试之间共享了测试片的均值/方差（若有，是什么后果）？训练时是否用了留出单位的任何信息？
8. 留一单位是否按 `unit_id` 分组（同单位多片不跨折）？
9. `r2_within`（片内 oracle）是否只在留出单位内拟合？`r_half_agreement` 是否已从错误的 `2r−1` 改为 Pearson `r`（上一版曾错）？
10. 三对 R² 0.19–0.27、7–8/8 单位为正，是否与 `louo.tsv` 一致？"达到片内上限 78%–99%"的计算是否正确？
11. **结论强度**：报告说"首次获得单位级迁移的有限正面证据"，是否越界（8 个单位、同平台、片内 z-score 去掉绝对亮度）？有没有把"同点耦合可迁移"说成"空间关系成立"？

## 交付

写 `infra/field_analyzer_20261020/review/fa08_review.json`，`scope.read`/`scope.ran` 非空：

```json
{"reviewer": {"provider":"xai","model":"grok-4.7","thinking":"xhigh","cwd":"/home/huyudi/013_spatial"},
 "nodes": ["FA-08A","FA-08B"], "commit_reviewed": "<git rev-parse --short HEAD>",
 "scope": {"read": ["..."], "ran": ["..."], "not_checked": ["..."]},
 "findings": [{"id":"N1","node":"FA-08A","severity":"BLOCKER|MAJOR|MINOR","topic":"...","evidence":"...","required_fix":"..."}],
 "answers": {"1":"...", "...":"..."},
 "claims_that_do_not_hold": [], "allow_node_conclusion": true, "limitations": "..."}
```

只报实际读到/跑到的内容；可跑命令复算。不要改仓库里除该 JSON 以外的文件。
