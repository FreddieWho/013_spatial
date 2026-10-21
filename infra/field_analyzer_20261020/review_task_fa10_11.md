# FA-10 / FA-11 节点完成复核（xai/grok-4.7，xhigh）

两条新架构（A 计数似然因子模型、C 生态位共定位）。重点是**有没有会让结论变错的错误**，以及报告里的两个"纠正"是否成立。

## 范围

| 节点 | 代码 | 产物 | 报告 |
|---|---|---|---|
| FA-10 A | `field_analyzer/factor_model.py`、`cli.py::cmd_factor_model` | `factor_model/{spot_holdout.tsv,louo_spot_holdout.tsv,loadings.npz,summary.json}` | `REPORT_FA10.md` |
| FA-11 C | `field_analyzer/niche_graph.py`、`cli.py::cmd_niche` | `niche/{contact_enrichment.tsv,labels.npz,summary.json}` | `REPORT_FA11.md` |

## 必答

**FA-10**
1. `depth` 是否真的只用**非建模基因**估计（无循环）？各变体是否用同一基因集、同一留出 spot、同一 depth？
2. `smooth_0` 与 `depth_only` 数值完全相同——这是设计使然（留出点因子为 0）还是 bug？
3. `gene_knn` 基线是否只用**训练** spot 的值、并且是逐基因（不共享信息）？它的 k=6 是否与因子模型的图 k=6 一致（公平）？
4. 留一单位（LOUO）是否用训练单位的载荷、留出单位的因子在训练 spot 上拟合后再传播？有没有用到留出点计数？
5. "因子模型比 depth-only 好 15%、但逐基因 kNN 好 27%"这两个数是否与 `spot_holdout.tsv` 可复算一致？"12/14 片最优"是否与 per-section best 一致？
6. 拟合是 log-rate 上的最小二乘近似，报告是否已明确这不是完整 Poisson MLE，且不足以解释 14% 的差距？

**FA-11**
7. **纠正 1 是否成立**：对称 kNN 图是否必然给出对称的计数矩阵？`max|C−Cᵀ|=0` 是否被代码与产物同时支持？据此说"方向不可识别"是否过头（是否存在其他可行定义）？
8. **纠正 2 是否成立**：旋转零对照的边界伪影（把标签场旋转后回贴到同一轮廓）是否真的是观测低于零对照（z=−2.42）的原因？还有没有别的原因？
9. 置换零对照 z=+220 是否只反映空间自相关（即平凡结果）？报告有没有把它当成"连接证据"（应只当共定位下限）？
10. 半区一致性 r=0.986、符号 20/20 是否与产物一致？半区划分是否会造成两侧样本量/组成差异从而夸大一致性？
11. 生态位定义（k-means on 6 programs，未用空间信息）是否可能产生"标签本身依赖空间"的循环？报告是否已说明 panel 缺基质/血管 marker 的限制？

## 交付

写 `infra/field_analyzer_20261020/review/fa10_11_review.json`，`scope.read`/`scope.ran` 非空：

```json
{"reviewer": {"provider":"xai","model":"grok-4.7","thinking":"xhigh","cwd":"/home/huyudi/013_spatial"},
 "nodes": ["FA-10","FA-11"], "commit_reviewed": "<git rev-parse --short HEAD>",
 "scope": {"read": [...], "ran": [...], "not_checked": [...]},
 "findings": [{"id":"P1","node":"FA-10","severity":"BLOCKER|MAJOR|MINOR","topic":"...","evidence":"...","required_fix":"..."}],
 "answers": {"1":"...", "...":"..."},
 "claims_that_do_not_hold": [], "allow_node_conclusion": true, "limitations": "..."}
```

只报实际读到/跑到的内容。不要改仓库里除该 JSON 以外的文件。
