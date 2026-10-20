# FA-01 逐节点复核（xai/grok-4.7，xhigh）

只审 FA-01：**尺度/方向/跨场关系可测性**。重点是"有没有会让结果解释变错的错误"。

## 范围

- 代码：`field_analyzer/scales.py`，以及 `field_analyzer/cli.py` 的 `cmd_scales` / `_load_synth_sample` 相关部分（如涉及）
- 产物：`infra/field_analyzer_20261020/scales/{scale_curves.tsv,scale_fits.tsv,direction_curves.tsv,direction_anisotropy.tsv,cross_field.tsv,summary.json}`
- 报告：`infra/field_analyzer_20261020/REPORT_FA01.md`
- 合同：`docs/FA_CONTRACT_ADDENDUM_20261020.md`、`docs/field_analyzer_killtest_20261020/05_COMMON_PROTOCOL.md` §7、`08_REAL_DATA_AND_SCALE.md` §3

## 必答

1. **原场 vs 残差尺度是否真的分开**：`quadratic_trend` 扣的是趋势，`fit_variogram` 对两者都用同一模型；报告有没有把残差尺度当成原场尺度使用？
2. **NOT_IDENTIFIED_WITHIN_FOV 判定是否正确**：`reached` 用 `max(gamma) >= nugget+0.8*sill`，`range >= 0.8*span` 也判未识别。这两个条件是否会**错误地**把真实有平台的场判成未识别，或反过来？
3. **拟合是否有可识别性陷阱**：指数模型的 nugget/sill/range 是否可能互相补偿（例如 sill→0 时 range 无意义）？报告里的 30/56 OK 是否受此影响？
4. **方向箱**：4 个 45° 箱 + ≥30 配对的可测阈值是否会造成方向性结论的系统偏差？`direction_anisotropy.tsv` 的 2–6 pitch 带宽选择是否会挑出有利结果？
5. **跨场指标**：`cross_field_curves` 报的是标准化交叉协方差（不是增量相关），h→0 即同点相关；报告有没有把它误读成"位移"或"因果"？（FA-01 报告里"不声明位移"这句是否与产物一致）
6. **统计单位**：14 片来自 8 个单位，`summary.json` 是否按片汇总就当成独立？报告的中位数是否被片数不均衡带偏？
7. **数字一致性**：抽 6 个数与 TSV/JSON 对齐（30/56、26 例、range 中位 2.52/5.01/4.18/5.49、93.75%、0.995、跨场 r）。

## 交付

写 `infra/field_analyzer_20261020/review/fa01_node_review.json`，**`scope.read` 必须列出你实际打开的文件，`scope.ran` 列出你实际执行的命令**（这两项不得为空；若某项未跑就写 `[]` 并在 `not_checked` 说明原因）：

```json
{
  "reviewer": {"provider": "xai", "model": "grok-4.7", "thinking": "xhigh", "cwd": "/home/huyudi/013_spatial"},
  "node": "FA-01",
  "commit_reviewed": "<git rev-parse --short HEAD>",
  "scope": {"read": ["..."], "ran": ["..."], "not_checked": ["..."]},
  "findings": [{"id": "H1", "severity": "BLOCKER|MAJOR|MINOR", "topic": "...", "evidence": "...", "required_fix": "..."}],
  "answers": {"1": "...", "...": "..."},
  "claims_that_do_not_hold": [],
  "allow_node_conclusion": true,
  "limitations": "..."
}
```

不要修改仓库里除该 JSON 以外的任何文件。只报实际读到/跑到的内容；不能确认的写 `NOT_ASSESSABLE`。
