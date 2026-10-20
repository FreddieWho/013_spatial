# FA 轮次调度记录（2026-10-20 起）

本文件记录每次派发的 agent ID、实际 provider／model／推理档位、工作目录、任务文件与范围，供审核与复现追溯。凭据不写入。

## 生效分工（AGENTS.md + D-280 + D-282）

| 角色 | Paseo agent | Provider | Model | Thinking | Cwd |
|---|---|---|---|---|---|
| 实现 | `41f406bd-bd3a-44e6-8dd8-2fa5c0ea93a6` | pi | `opencode-go/deepseek-v4.1-flash` | `max` | `/home/huyudi/013_spatial` |
| 独立审核 | `a5bcb012-6618-488e-8f88-f700773fd4c4` | pi | `xai/grok-4.7` | `high` | `/home/huyudi/013_spatial` |

实际核对时间：2026-10-20，`paseo inspect` 返回 Provider=pi、Model=xai/grok-4.7、Thinking=high、Status=running；实现代理同一轮 `paseo ls` 显示 Thinking=max、Status=idle。若后续运行时模型/档位不可用，按项目约束报阻塞，不静默降级。

## 派发记录

### 2026-10-20｜预审（Grok high）

- **任务文件：** `infra/field_analyzer_20261020/preflight_task.md`
- **派发命令：** `paseo send a5bcb012 --prompt-file infra/field_analyzer_20261020/preflight_task.md --no-wait`
- **范围：** 只读预审 `05_COMMON_PROTOCOL.md` 观测合同、`reuse_register.tsv` 复用裁定、FA-00 实现范围（对应 I-055）。不改仓库、不跑模型、不读人工 TLS。
- **期望产出：** `infra/field_analyzer_20261020/review/preflight_fa00.json`，含每项 PASS/GAP/BLOCKER、`allow_fa00_implementation`、必要修复与最小范围建议。
- **等待方式：** 单个持久 `paseo wait`（后台任务 `Grok preflight wait`），不轮询、不循环读日志。
- **结果（2026-10-20）：** `allow_fa00_implementation = true`，无 BLOCKER，3 PASS / 7 GAP，审查提交 `8b174e2`。GAP 已由 `docs/FA_CONTRACT_ADDENDUM_20261020.md` 关闭；其中第 8 项指出包内精度式 `Q = τI + ΣLᵀWL` 漏掉观测项 `MᵀR⁻¹M`，属 FA-03B 前置修正。预审自述未覆盖实现、未运行毒化反例、未核验本地文件。

### 2026-10-20｜FA-00 实现（DeepSeek max，已派发）

- **任务文件：** `infra/field_analyzer_20261020/implementation_task_fa00.md`（已按补充条款更新为权威口径）
- **派发命令：** `paseo send 41f406bd --prompt-file infra/field_analyzer_20261020/implementation_task_fa00.md --no-wait`
- **开工前实测：** Provider=pi，Model=`opencode-go/deepseek-v4.1-flash`，Thinking=`max`，Status=idle，Cwd=`~/013_spatial`。
- **范围：** `field_analyzer/` 最小接口（observation/adapters/geometry）、`reuse_resolution.tsv`（五态资产）、`unit_role_assertions.tsv`、五入口毒化反例 + 评价随动反向断言、含被掩膜场值的小几何样本、`run_manifest.json`、`REPORT_FA00.md`。
- **禁止项：** 不读人工 TLS、不改 `scripts/`、不开 GPU、不下载、不进入 B/C、不跑真实拟合。
- **等待方式：** 单个持久 `paseo wait`（后台任务 `FA-00 impl wait`，上限 3 小时），不轮询。
- **完成后动作：** 交 `a5bcb012`（Grok high）复核实现与毒化证据，审核未关闭前不进入 FA-02 之外的任何结果解释。

### 待派发

- **FA-02 合成生成器（DeepSeek max）：** 与 FA-00 可部分并行，但生成器与 S0/S1/S3/S5 反例的接口合同需在 FA-00 接口事实明确后再派，避免两条线各自造数据接口。

## 记录规则

- 每次派发写：时间、agent、任务文件、范围、前置条件、期望产出、等待方式。
- 不记录凭据，不记录模型内部推理。
- 审核结论与必修项存 `infra/field_analyzer_20261020/review/`，范围用 Git 提交与必要 diff 标明，不另建手工哈希清单。
