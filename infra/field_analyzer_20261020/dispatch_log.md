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

### 待派发（视预审结论）

- **FA-00 实现（DeepSeek max）：** 任务文件 `infra/field_analyzer_20261020/implementation_task_fa00.md`，前置条件为预审 `allow_fa00_implementation: true`。
- **FA-02 合成生成器：** 与 FA-00 可并行，但其生成器与 S0/S1/S3/S5 反例的具体合同需要在 FA-00 接口冻结后再派，避免两条线各自造数据接口。

## 记录规则

- 每次派发写：时间、agent、任务文件、范围、前置条件、期望产出、等待方式。
- 不记录凭据，不记录模型内部推理。
- 审核结论与必修项存 `infra/field_analyzer_20261020/review/`，范围用 Git 提交与必要 diff 标明，不另建手工哈希清单。
