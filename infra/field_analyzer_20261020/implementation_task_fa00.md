# FA-00 实现任务（交指定执行代理 opencode-go/deepseek-v4.1-flash，max）

开始前必须读：`AGENTS.md`、`docs/plan.md`、`docs/roadmap.md`（顶部 FA 节点表）、`docs/field_analyzer_killtest_20261020/05_COMMON_PROTOCOL.md`、`.../07_ALGORITHM_ROUTES.md` §0、`.../08_REAL_DATA_AND_SCALE.md` §1/§6、`.../reuse_register.tsv`、`docs/ISSUES.md` I-053/I-055，以及 `infra/field_analyzer_20261020/review/preflight_fa00.json`（Grok 预审结论，必须 `allow_fa00_implementation: true` 才能开始）。

若预审为 false 或含未关闭 BLOCKER：只做被允许的子项，其余报阻塞，不自行改合同。

## 项目边界（不得越界）

- 仅本地已有数据、仅 CPU（默认 ≤16 线程，`CUDA_VISIBLE_DEVICES=-1`），磁盘余量必须 ≥200,000,000,000 bytes，未测峰值内存前不扩大规模。
- 不下载任何外部数据/权重/依赖；不联系作者；不使用 GPU。
- 不读取、不派生、不引用人工 TLS 标签或由标签产生的端点/候选；不读 HookNet 输出与自动弱标签。
- `scripts/` 下既有脚本是冻结证据：**不修改、不重排、不改字节**。新代码只放 `field_analyzer/`，产物只放 `infra/field_analyzer_20261020/`。
- 不改 `docs/plan.md`、`docs/roadmap.md`、`docs/decisions.md`、`STATUS.md`、`TODO.md`；如需改科学表述，先在交付报告里提出，由主代理处理。
- 不把"接口写完/单测通过"表述为 AF 假设获得支持。

## 交付物

### 1. `field_analyzer/` 最小包

- `observation.py`：`ObservationBundle` 数据类，字段至少 `sample_id/unit_id/source_id/patient_status`、`barcodes/source_rows/xy/coordinate_unit/scale_provenance`、`tissue_support/visible_mask`、`field_ids/gene_sets/score_definition`、`values_visible`、`auxiliary_visible`、`input_provenance`。掩膜必须把不可见值置 `NaN` 并物理分离到其他文件，禁止 `0` 冒充缺失；不得在 bundle 上暴露隐藏区任何派生量（QC、深度、检出率、全片 PCA、预平滑）。
- `adapters.py`：只读允许的原始 counts 与 metadata。**显式** barcode 轴与单位映射：读入时断言矩阵形状与 barcode 数一致，出现方向不确定就报错而不是转置。重复 gene symbol、缺失 gene、零深度样本的处理必须写明并留痕。参考 `scripts/field_discovery_20261018.py` 的 `strings/read_matrix/assert_grouped`，但必须去掉形状启发式转置。
- `geometry.py`：单位、邻接、距离对枚举、局部方向（结构张量或可见场梯度），区分**已知组织空洞**与**被遮但可能有组织**的区域；不得把不可见区域边当成真实空洞，也不得用 ConvexHull 把缺口补成实测组织。`tissue_graph` 类实现禁止硬编码 130.01 μm。
- 测试：`tests/test_field_analyzer_*.py`，覆盖非方阵/方阵、CSR/CSC、barcode 重排、重复 symbol、缺基因、带洞组织、近接曲线、零长度段。

### 2. 资产复用裁定 `infra/field_analyzer_20261020/reuse_resolution.tsv`

列：`path, what, availability, exists, bytes, mtime, manifest_version, reuse_decision, provenance, notes`。

- `availability` 取 `VERIFIED_LOCAL / REFERENCED_LOCAL_NOT_VERIFIED / NOT_PRESENT / TRACKED_SOURCE_ONLY`；`reuse_decision` 取 `REUSE / ADAPT / REBUILD / HISTORICAL_ONLY / SEALED`，与 `reuse_register.tsv` 对齐，不一致必须写理由。
- 至少覆盖：`data/other_sources/hest1k/st/{section_id}.h5ad`、`infra/field_discovery_review_20261018/identity_review.tsv`、`infra/field_source_holdout_20261018/`、`infra/tls_multiscale_gobp_20261002/panel_controls.json`、`scripts/` 中四个待适配脚本、`infra/r16/t2p2_20260928/`。
- 清点方式：相对路径 + 存在性 + 字节数 + 已有 manifest 版本；**不做全量逐字节 hash**，必要时只用头尾摘要做快速比对并注明它不能证明物理独立。
- 真实片数与单位数按实际 manifest 统计，不得为凑"8 单位/14 片"而放行任何片；不符就写实际值 + 原因。

### 3. 信息入口毒化反例（本节点最重要）

最小样本上建立断言并做成测试或可复跑脚本，输出 `infra/field_analyzer_20261020/poisoning_report.json`：

- (a) 把隐藏区 H 内表达替换为极大值 / 随机置换，(b) 置换核验组 B，(c) 改动一个含人工 TLS 的历史文件路径的存在性，三者在**不改变 O 与允许 metadata** 的前提下重建全部模型输入，要求：输入指纹、资格集合、方向/端点（若已实现）、拟合参数、预测数组全部逐位不变；只有最终评价可以变。
- 断言必须覆盖"资格/查询点选择不由 H 表达决定"和"图乘不把隐藏 0 当真实 0、不用孤立点 self fallback 复制 H 值"。
- 反例应能在小数据上秒级到分钟级跑完，不需要真实全量。

### 4. 小几何样本

用合成 64×64 规则格 + 一个带洞/窄颈的变体，验证：坐标单位来源、邻接构建、距离对枚举、局部方向可测与不可测判定。输出 `infra/field_analyzer_20261020/geometry_smoke/*.tsv`（真实数值列，不要只放图），并在报告中说明哪些量在有洞组织上不可评。

### 5. `infra/field_analyzer_20261020/run_manifest.json`

按 `docs/field_analyzer_killtest_20261020/schemas/run_manifest.schema.json` 填写；`status` 用 `PARTIAL` 或 `COMPLETED`，`results_scope` 在无科学结果时用 `ENGINEERING_ONLY`，`observation_contract.hidden_expression_available=false`、`witness_available_to_model=false`、`human_tls_used=false`。

### 6. `infra/field_analyzer_20261020/REPORT_FA00.md`

人可读，必须逐条回答：实际读了哪些路径（及字节数）、哪些资产缺失或不可用、实现了哪些接口与断言、毒化反例的实测结果（通过/未通过及数值）、几何样本的可评/不可评、资源消耗（墙钟、峰值内存、线程、磁盘余量）、明确失败与未做项、以及"本节点对哪条 AF 有影响"（预期答案：无直接影响，属基础设施）。

## 验收（必须自证，不得只写声明）

1. `python -c "import field_analyzer"` 及全部新增测试通过；给出测试命令与实际输出摘要。
2. 毒化反例的逐位不变性有数值证据（数组哈希或 max|Δ|），不是"没有报错"。
3. `reuse_resolution.tsv` 每行的 `exists` 都与同一命令的实际结果一致。
4. 无 GPU、无下载、线程数、磁盘余量有实测记录。
5. 不改任何 `scripts/` 文件（用 `git status` 证明）。

## 交付后的动作

交付后停在本节点：不进入 B/C 实现、不跑真实拟合、不写科学结论。由主代理交独立审核代理复核，审核未关闭前不把结果用于任何 AF 判断。若发现预审合同与本地实际冲突，写清冲突并停下，不自行降低要求。

提交信息用中文、说明改了什么与为什么；提交前更新 `STATUS.md`/`TODO.md`/`docs/roadmap.md` 的节点状态由主代理负责，你只提供事实清单。
