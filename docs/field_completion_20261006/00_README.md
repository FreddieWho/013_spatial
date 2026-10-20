# 013_spatial｜空间分子场条件恢复：目标与路线更新包

版本：2026-10-06。用途：将本次讨论形成的科学目标写入现有项目，指导后续执行；不是实验完成报告，也不是新的平行项目。

> **主线：学习信号场之间的条件空间关系，从可见场的任意可用子集，预测未观测场的分子身份、空间分布及不确定性。固定字典先验证恢复任务，分子组成发现随后升级；TLS 是有独立结构参照的正案例，easy / harder 是核心，extreme 是探索。**

## 文件与使用

| 文件 | 职责 |
|---|---|
| [01_SCIENCE.md](01_SCIENCE.md) | 科学目标、场的抽象、模型输入输出、假设和结论边界；不规定代码、算法实现或执行顺序。 |
| [02_ROADMAP.md](02_ROADMAP.md) | 分阶段路线、资产复用、证据要求与 agent 的执行裁量。 |
| [03_UPDATE_PROJECT_PROMPT.md](03_UPDATE_PROJECT_PROMPT.md) | 交给仓库 agent，实际更新 plan / roadmap / decisions / STATUS 等上层文件；本轮仅迁移文档。 |

将解压后的 `field_completion_20261006/` 放入仓库 `docs/`，把第三份文件全文交给 agent。无需先手工替换现有文件。迁移完成后，项目的唯一现行入口仍是 `STATUS.md`、`docs/plan.md`、`docs/roadmap.md`；本包保留为本次更新的来源快照，不另建一套长期维护的上层系统。

## 依据与仓库起点

科学设计依据本次对话中用户逐步确认的目标；证据边界和路线组织是对这些约定的整理，不代表假设已经成立。仓库核对基准为 `main@5a2d97aa8d1ef3315b32773dabf89f6d717b3ce7`。本包没有重新运行实验、核验全部底层产物，或独立复核此前对话中的外部论文和工具性能，因此不承诺“首个方法”、特定架构必胜或自动标注必能扩样 10–100 倍。

已核对的关键事实：现有入口仍以 TLS 外围预测和桥校准为主；GSE 的 18 片已按 17 位患者分组，USZ 的 8 个肿瘤样本仍有患者身份限制；既有预测属于内部证据，不能改名为本次新模型的结果。[R1–R4]

既有 GO、Hallmark、B-axis、TLS 标注、坐标和遮蔽审计可复用，但旧缓存是否符合新任务的可见信息约束，必须按任务检查。旧桥检验失败保留，不再充当全部恢复任务的前置门槛。当前 GPU、外部获取和存储约束不因本包自动解除。[R1, R2]

### 可追溯来源（固定到核对版本）

- [R1｜AGENTS.md：文档职责与资源约束](https://github.com/FreddieWho/013_spatial/blob/5a2d97aa8d1ef3315b32773dabf89f6d717b3ce7/AGENTS.md)
- [R2｜STATUS.md：当前桥校准工作面](https://github.com/FreddieWho/013_spatial/blob/5a2d97aa8d1ef3315b32773dabf89f6d717b3ce7/STATUS.md)
- [R3｜docs/plan.md：原科学目标和 H-01～H-08](https://github.com/FreddieWho/013_spatial/blob/5a2d97aa8d1ef3315b32773dabf89f6d717b3ce7/docs/plan.md)
- [R4｜TLS_PATIENT_VALIDATION_20261003.md：患者分组、内部预测与限制](https://github.com/FreddieWho/013_spatial/blob/5a2d97aa8d1ef3315b32773dabf89f6d717b3ce7/docs/TLS_PATIENT_VALIDATION_20261003.md)
- [R5｜docs/roadmap.md：原节点、既有工作及历史状态](https://github.com/FreddieWho/013_spatial/blob/5a2d97aa8d1ef3315b32773dabf89f6d717b3ce7/docs/roadmap.md)
- [R6｜docs/decisions.md：仅追加的决策记录](https://github.com/FreddieWho/013_spatial/blob/5a2d97aa8d1ef3315b32773dabf89f6d717b3ce7/docs/decisions.md)

第三份 prompt 要求 agent 再读执行时的实际 HEAD；本包基准不是要求回退代码或覆盖更新进展。
