# FA-07A / FA-07B 节点完成复核（xai/grok-4.7，xhigh）

两个并发跟进节点。重点只有一件事：**有没有会让结论变错的错误**。特别注意 FA-07B 声称推翻了 FA-03B 的归因——这是本轮最需要盯的地方。

## 范围

| 节点 | 代码 | 产物 | 报告 |
|---|---|---|---|
| FA-07A 深度条件化 | `field_analyzer/depth_conditioning.py`、`cli.py::cmd_depth_conditioned` | `depth_conditioning/{map_agreement_before_after.tsv,relations_before_after.tsv,summary.json}` | `REPORT_FA07A.md` |
| FA-07B 观测信号 | `field_analyzer/signal_power.py`、`route_b.py::cross_correlation_peak`、`synthetic.py`（markers_per_field / sigma） | `signal_sweep/{sweep.tsv,summary.json}` | `REPORT_FA07B.md` |
| 更正 | — | `REPORT_FA03B.md` 顶部更正框、`REPORT_FA06.md` AF-02 行 | — |

## 必答

**FA-07A**
1. `residualise` 是否只用了可见点、是否与评分共用同一个外部文库？残差是否真的只扣了深度（`1+d+d²`），有没有顺手扣掉结构？
2. "三对存活 / IFN~OXPHOS 翻转 / EMT~OXPHOS 抑制"是否与 `relations_before_after.tsv` 按单位中位复算一致？翻转对在条件化后符号一致单位数是几个？
3. OXPHOS 的 A/B 一致性 0.790→0.369 是否在同样的点上算、用同样的两半文库？这个下降能否支持"其空间图大量由深度提供"？
4. 报告有没有把描述性条件化说成因果校正？

**FA-07B**
5. `score_fields` 是否只用计数与基因名，没有用 latent truth 或 structures？三种观测（markers/depth/σ）之间除这些参数外是否真的只变了观测？
6. **互相关峰 6–8 pitch 是否可信**：符号约定与 ridge 的 (0,8) 是否只是平移约定差异？互相关的峰值是否落在搜索边界（±12）从而可能是假峰？`min_overlap=500` 会不会让边界平移获得不公平优势？
7. **"ridge 失败于共线性"这个归因是否被证据支持**，还是可能由实现细节（零填充、α=1、步长 4）造成？报告是否已把这一点写成"未穷尽改进"？
8. "12/12 组合"的计数是否与 `sweep.tsv` 一致？same_point 是否 12/12 为 0？
9. "marker 12→96、深度 ×3 不改变任何结论"是否成立（有没有任何一项指标确实变好/变差）？
10. **FA-03B 的更正是否过头**：把"观测抹掉信号"改成"估计器无法定位"，是否削弱了原本正确的观察（宽场下对比度很弱）？AF-02 的否定范围收窄是否被产物支持？

## 交付

写 `infra/field_analyzer_20261020/review/fa07_review.json`，`scope.read`/`scope.ran` 必须非空：

```json
{"reviewer": {"provider":"xai","model":"grok-4.7","thinking":"xhigh","cwd":"/home/huyudi/013_spatial"},
 "nodes": ["FA-07A","FA-07B"], "commit_reviewed": "<git rev-parse --short HEAD>",
 "scope": {"read": [...], "ran": [...], "not_checked": [...]},
 "findings": [{"id":"M1","node":"FA-07B","severity":"BLOCKER|MAJOR|MINOR","topic":"...","evidence":"...","required_fix":"..."}],
 "answers": {"1":"...", "...": "..."},
 "claims_that_do_not_hold": [], "allow_node_conclusion": true, "limitations": "..."}
```

只报实际读到/跑到的内容；可跑测试与命令复算。不要改仓库里除该 JSON 以外的文件。
