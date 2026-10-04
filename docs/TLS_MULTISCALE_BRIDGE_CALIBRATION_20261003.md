# 多尺度桥的原位参考与必要校准（2026-10-03）

当前结果：旋转参考方法未通过必要校准，不能据此确认分子桥，也不能用未检出证明没有桥。全量参考与诊断计算已经完成；这项失败不会被更换半径、只挑最强GO或降低阈值覆盖。

## 全量执行范围

26片、全部1691个压缩GO BP及51个固定对照；5048对端点形成80768种路径/mask/宽度组合，41058种可评。原位旋转共得到404049个合格几何参考，全部程序均已评分。长度含不足500μm至8mm以上，半宽100/250/500/1000μm；宽片与狭长连接分开解释。

6449种几何满足至少19个参考且四个方向区间各至少2个的描述性参考门槛；其中6300种同时满足长度至少为走廊全宽。参考不足的其他组合全部保留得分与排除原因，不按表达挑选或抹成阴性。TUM/UNKNOWN比例匹配仅为次级敏感性，主对照保留组成信息。

## 为什么不能直接把参考分位当显著性

严格诊断使用真实方向及全部39个旋转方向都可测的同一候选集合，每个方向都执行相同的路径、mask、宽度、长度及端点类型搜索。只有8片、1242种几何满足这一条件，约占全部可评几何的3.0%。其余几何不在此必要校准的覆盖范围内。

在实际组织点和mask上，六类无TLS定向连接的已知场各模拟400次。48个设置中35个未通过误报率的预定上界门槛。8片中7片至少一项失败；另一片仅通过单场必要条件，不等于真实多GO联合场或患者级确认。偏差也出现在平稳高斯场中，不能只归咎于最复杂的压力场。

0.5与1.0可见域SD的注入共1776个设置，每个200次；0.5SD的888个设置中767个未达到预定检出能力下界。阴性不能排除有用桥效应。不同设置共享背景，不把设置数量当独立生物重复。

| 切片 | 共同几何数 | 误报未过 / 已测设置 | 0.5SD功效未过 | 状态 |
|---|---:|---:|---:|---|
| usz-KC1 | 511 | 4 / 6 | 199 | ROTATION_EXCHANGEABILITY_GATE_FAILED |
| usz-KC2 | 16 | 4 / 6 | 15 | ROTATION_EXCHANGEABILITY_GATE_FAILED |
| usz-KC3 | 157 | 5 / 6 | 148 | ROTATION_EXCHANGEABILITY_GATE_FAILED |
| usz-LC1 | 151 | 6 / 6 | 183 | ROTATION_EXCHANGEABILITY_GATE_FAILED |
| usz-LC2 | 136 | 6 / 6 | 112 | ROTATION_EXCHANGEABILITY_GATE_FAILED |
| usz-LC3 | 248 | 6 / 6 | 55 | ROTATION_EXCHANGEABILITY_GATE_FAILED |
| usz-LC4 | 22 | 4 / 6 | 52 | ROTATION_EXCHANGEABILITY_GATE_FAILED |
| usz-LC5 | 0 | — / — | — | NO_COMMON_39_DIRECTION_GEOMETRY |
| gse175540-GSM5924030_ffpe_c_2 | 0 | — / — | — | NO_COMMON_39_DIRECTION_GEOMETRY |
| gse175540-GSM5924031_ffpe_c_3 | 1 | 0 / 6 | 3 | NECESSARY_UNIVARIATE_GATE_PASSED_ONLY |
| gse175540-GSM5924032_ffpe_c_4 | 0 | — / — | — | NO_COMMON_39_DIRECTION_GEOMETRY |
| gse175540-GSM5924033_ffpe_c_7 | 0 | — / — | — | NO_COMMON_39_DIRECTION_GEOMETRY |
| gse175540-GSM5924035_ffpe_c_20 | 0 | — / — | — | NO_EVALUABLE_GEOMETRY |
| gse175540-GSM5924037_ffpe_c_34 | 0 | — / — | — | NO_EVALUABLE_GEOMETRY |
| gse175540-GSM5924038_ffpe_c_36 | 0 | — / — | — | NO_EVALUABLE_GEOMETRY |
| gse175540-GSM5924039_ffpe_c_39 | 0 | — / — | — | NO_EVALUABLE_GEOMETRY |
| gse175540-GSM5924040_ffpe_c_45 | 0 | — / — | — | NO_EVALUABLE_GEOMETRY |
| gse175540-GSM5924041_ffpe_c_51 | 0 | — / — | — | NO_COMMON_39_DIRECTION_GEOMETRY |
| gse175540-GSM5924043_frozen_a_3 | 0 | — / — | — | NO_EVALUABLE_GEOMETRY |
| gse175540-GSM5924044_frozen_a_15 | 0 | — / — | — | NO_EVALUABLE_GEOMETRY |
| gse175540-GSM5924046_frozen_b_1 | 0 | — / — | — | NO_COMMON_39_DIRECTION_GEOMETRY |
| gse175540-GSM5924049_frozen_b_18 | 0 | — / — | — | NO_EVALUABLE_GEOMETRY |
| gse175540-GSM5924050_frozen_c_2 | 0 | — / — | — | NO_COMMON_39_DIRECTION_GEOMETRY |
| gse175540-GSM5924051_frozen_c_5 | 0 | — / — | — | NO_EVALUABLE_GEOMETRY |
| gse175540-GSM5924052_frozen_c_23 | 0 | — / — | — | NO_COMMON_39_DIRECTION_GEOMETRY |
| gse175540-GSM5924053_frozen_c_57 | 0 | — / — | — | NO_EVALUABLE_GEOMETRY |

## 对假设的影响

H-08仍未获确认。这次完整审计说明当前旋转参考不能提供可靠的普遍桥检验，不说明桥不存在。正向图上连续、端点调整后仍高或旋转参考中排名靠前，均只能列为描述性线索。严格患者映射仍受I-033限制。

H-01与H-05的反向内部预测另见[全量反向结果](TLS_MULTISCALE_REVERSE_RESULTS_20261003.md)，不能用那里的AUC替代桥的空间证据，也不能用桥校准失败抹去已观测的预测信息。

## 完整结果

- [全部源与全部程序诊断](../infra/tls_multiscale_gobp_20261002/bridge_calibration/all_source_program_diagnostics.tsv)：45292行；正式p字段全部留空，秩只作诊断。
- [误报校准全部设置](../infra/tls_multiscale_gobp_20261002/bridge_calibration/null_calibration.tsv)
- [功效全部设置](../infra/tls_multiscale_gobp_20261002/bridge_calibration/power_calibration.tsv)
- [校准汇总](../infra/tls_multiscale_gobp_20261002/bridge_calibration/calibration_summary.json)
- [原位参考全源汇总](../infra/tls_multiscale_gobp_20261002/controls/source_summary.tsv)
- [多尺度可测性](../infra/tls_multiscale_gobp_20261002/summary/length_width_support.tsv)

R-18的晕/边背景、全GO前向/反向整合和完整图表交付现已完成，见[完整结果](TLS_MULTISCALE_GOBP_RESULTS_20261003.md)。原目标的正式确认仍未达标，本桥校准失败不会因交付完整而解除。

**D-189更新：** 桥误报/功效失效条件已细化，见[诊断报告](TLS_BRIDGE_FAILURE_AUDIT_20261003.md)。原旋转方法仍不支持正式确认，不由成功子集或无噪声示例解除。
