# 文档包自身检查（不是实验结果）

日期：2026-10-20。以下检查仅针对文档、接口规格和交付完整性。没有运行013_spatial的新模型、真实数据或公共方法。

| 检查 | 结果 |
|---|---|
| `yaml_documents_parse` | PASS |
| `proposed_not_executed` | PASS |
| `synthetic_seed_roles_disjoint` | PASS |
| `core_within_smallest_mask` | PASS |
| `resource_boundaries_consistent` | PASS |
| `no_hidden_or_human_inputs` | PASS |
| `json_schema_valid` | PASS |
| `all_source_ids_resolved` | PASS |
| `no_stale_chat_citation_tags` | PASS |
| `balanced_code_fences` | PASS |
| `all_main_documents_present` | PASS |
| `six_worker_prompts_present` | PASS |
| `reuse_sources_resolved` | PASS |
| `no_computational_project_code_bundled` | PASS |

## 人工语义复核重点
已在文档中区分：旧含同点X的空间模型与新块遮挡；原场尺度与残差尺度；实测连接与推断连接；几何已知与坐标生成；A/B分子核验与人工病理；GitHub源码与未验本地缓存；经典核/图算法与尚待证明的新贡献。

已补入：B联合条件求解中的截距/背景右端项；旧候选缓存的行与barcode/xy再对齐；C首版可执行的共同支持评分，以及它可能被简单聚合解释的停止条件。

未验证：服务器文件存在性、环境与资源、源码未读分支、旧结果独立重训、所有近邻工具的运行性能。最终使用前仍需FA-00与指定独立审核。
