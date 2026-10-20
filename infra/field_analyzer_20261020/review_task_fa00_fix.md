# FA-00 复核意见 F1–F4 的修复复核（xai/grok-4.7，xhigh）

上一轮你审 `cb7256c` 给出 `allow_next_node: true` 与 4 项发现（F1/F2 MAJOR，F3/F4 MINOR）。本轮只复核这 4 项是否已正确关闭，不重新审全节点。

## 修复内容（工作区，未提交）

- **F1** `field_analyzer/poisoning.py`：`truth_max_abs_delta` 改为"无有限差值时写 `null`"，新增 `n_finite_hidden_truth`；`metric_changed` 只接受**有限且非零**的差；全缺失判失败。回归测试 `tests/test_fa_poisoning.py::test_all_missing_hidden_truth_cannot_satisfy_the_reverse_assertion`。
- **F2** `field_analyzer/adapters.py`：`AxisSpec` 增加 `allow_square=False`、`square_reason`；`n_obs == n_genes` 时默认 `ContractError`（"square matrix: obs and gene axes are indistinguishable..."），显式声明才读取并写入 `report.notes`。回归测试 `tests/test_fa_adapters.py::test_square_matrix_is_rejected_without_explicit_declaration`。
- **F3** `infra/field_analyzer_20261020/REPORT_FA00.md`：数字改为从 JSON 复算——NCBI642 隐藏点 130、ZEN36 B_AXIS 真值差 0.547、EMT 1.091、场基因并集（存在）634、文库 35,967（定义并集 639 已注明）。
- **F4** `field_analyzer/pipeline.py`：`detection` 与 `depth` 使用同一个目标无关分母（`library_mask`），`meta.detection_scope` 记录；回归测试 `tests/test_fa_adapters.py::test_detection_denominator_excludes_field_genes`。

已重跑：`poison --section ZEN36` 与 `--section NCBI642`，两者 `passed=true`、三个场景 `max|Δ|=0`、`metric_changed=true`；FA 定向测试 45 passed。

## 请回答

1. F1 的判定现在能否拒绝"全缺失"这种空测试？给出你检查到的代码路径。
2. F2 的方阵拒绝是否覆盖了稠密与稀疏两条读入路径？显式声明时是否把假设写进产物？
3. F3 现在报告里的四个数字是否与 `poisoning_report_*.json`、`reuse_resolution.tsv` 一致？还有没有残留的不一致？
4. F4 的 detection 分母是否已经与目标无关？是否还有别的量仍用含 A 的分母？
5. `tests/test_fa_*.py` 是否真的覆盖了你在 F2/F1 指出的空路径（不只是文件存在）。

## 交付

覆盖写 `infra/field_analyzer_20261020/review/fa00_fix_review.json`：

```json
{
  "reviewer": {"provider": "xai", "model": "grok-4.7", "thinking": "xhigh", "cwd": "/home/huyudi/013_spatial"},
  "commit_reviewed": "<git rev-parse --short HEAD>",
  "items": [{"id": "F1", "closed": true, "evidence": "...", "residual_risk": "..."}],
  "new_findings": [],
  "allow_next_node": true,
  "limitations": "..."
}
```

只报你实际读到、实际跑出来的内容；不能确认的写 `NOT_ASSESSABLE`。不要修改仓库里除该 JSON 以外的文件。
