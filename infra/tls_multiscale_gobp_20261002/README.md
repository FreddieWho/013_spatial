# 全量GO BP、多尺度TLS场与桥

本地执行完整；正式双向关系与桥确认未成立。先读[结果报告](../../docs/TLS_MULTISCALE_GOBP_RESULTS_20261003.md)和[执行验收](../../docs/TLS_MULTISCALE_ACCEPTANCE_20261003.md)。

- `contract.json`：冻结范围、输入与尺度；1691项项目压缩GO BP，不是官方GO-slim。
- `summary/`：前向汇总、长度/半宽可测性、真实标注支持。
- `reverse/summary/`：全部程序的84任务留出结果、资格、配对基线及训练内选集。
- `controls/`、`bridge_calibration/`：全量原位桥参考、失败的必要校准、功效和不可测状态。
- `halo_edge/`：26片×两任务的原位背景及必要校准。
- `integrated/program_evidence_index.tsv`：1742项共同入口。
- `integrated/all_program_bidirectional_screen.tsv`：3484行双向证据；`all_program_exploratory_ranking.tsv`保留全部mask/信息臂及未排名项。
- `integrated/all_bins.tsv`与`halo_edge_bins.npz`：完整曲线；`all_contrasts.tsv`与`halo_edge_contrasts.npz`：对比；长表`all_source_program_contrasts.tsv.gz`。
- `integrated/bridge_groups.tsv`与`bridge_all_programs.npz`：1559组×1742项桥视图；端点尺度见`all_endpoint_pair_scales.tsv`。
- [完整图表入口](figures/index.md)：1742页曲线、完整桥矩阵、全部来源代表图及反向总览。
- `acceptance_audit.json`：输入/缓存/预测SHA、覆盖/形状/完整性验收；不等于科学确认。
- `missing_evidence/patient_mapping_request.tsv`：患者对应模板，未知字段保持空白。

正式p留空；参考秩、分位和z均是诊断量。切片不是已确认独立患者。代表通路为探索性选择，代表连线只按几何资格选取。核心外组成信息保留；不作组成之外机制、解剖连接、位置概率或因果主张。
