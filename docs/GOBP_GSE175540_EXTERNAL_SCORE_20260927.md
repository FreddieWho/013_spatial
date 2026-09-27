# GSE175540 患者级外部 TLS 标签关联（2026-09-27）

## 问题与结论

45 个 GO BP 分子预测候选的 input 分数，在 16 位新患者的 TLS 点与 NO_TLS 点之间有没有差异？**整体没有**：31 项严格程序的中位 SMD 为 raw +0.045、归一化 +0.051，中位 AUC 0.51–0.52；raw 评分下没有符号检验 p<0.05 的程序。唯一的例外是 `GOBP_ADAPTIVE_IMMUNE_RESPONSE`（归一化中位 SMD +1.37，16/16 同向，AUC 0.84，Holm 后显著）， acompañada de 淋巴/炎症相关程序的弱同向——这与“TLS 是免疫聚集体”一致，应作阳性对照一致性读，不作新发现。

分子可预测性（45 项的入选依据）不等于 TLS 标签关联。本轮是两者在独立患者中的第一次分离描述。

## 方法（D-148 冻结，评分前未看数）

- 两臂并行：Arm1 严格 16×31；Arm2 同 16 患者×45 程序，两聚合路线 AGG-SUM（重复符号先加总）与 AGG-DROP（剔除歧义符号）。
- 评分 raw 主＋深度归一化敏（lib_excl 两臂一致）；端点主 SMD＋敏 AUC；跨患者中位＋符号计数＋双侧符号 p，每臂每评分 SMD 家族内 Holm。
- 空白标签未知；患者为单位；无块级、无跨癌种主张（USZ 含肾癌样本）。

## 结果

| 臂 | 评分 | 程序数 | 中位SMD | 中位AUC | 符号p<0.05 | Holm显著 |
|---|---|---|---|---|---|---|
| strict | raw | 31 | +0.045 | 0.510 | 0 | 0 |
| strict | 归一化 | 31 | +0.051 | 0.518 | 4 | 1（ADAPTIVE_IMMUNE_RESPONSE） |
| agg_sum | raw | 45 | +0.050 | 0.518 | 0 | 0 |
| agg_sum | 归一化 | 45 | +0.088 | 0.524 | 7 | 1（同上） |
| agg_drop | raw | 45 | +0.050 | 0.518 | 0 | 0 |
| agg_drop | 归一化 | 45 | +0.088 | 0.524 | 7 | 1（同上） |

- SUM 与 DROP 的程序中位 SMD 最大差 0.005：聚合选择不影响结论， frozen 敏感性问题已回答。
- raw 下 ADAPTIVE_IMMUNE_RESPONSE 仅 11/16（SMD +0.66，Holm 不显著）；归一化把它推到 16/16——评分处理敏感，与 45 项筛选时的 629→45 收缩同类。
- 其余 7 个归一化 raw-p<0.05 均为淋巴/炎症/黏附相关，Holm 后全不显著；只作描述。

## 边界

1. 45 项名单用 ST-CRC/USZ 全队列结果后选；本轮患者虽新，候选不是，仍属探索性。
2. 本轮只测标签关联，不测 readout 预测——对 H-03（复合生态预测）无新增支持；H-01/H-05 不变，plan.md 不变。
3. ADAPTIVE_IMMUNE_RESPONSE 的显著是免疫基因在 TLS 富集的预期内结果，不能升级为定位能力或机制发现。
4. 真实空间 p/T2P2 仍未做（TODO 未勾选项保持）。

## 文件

- 合同：`infra/gobp_external_score_20260927/run_contract.json`
- 逐患者端点：`patient_endpoints.tsv`（3872 行）；臂汇总：`arm_summary.tsv`（242 行）；收据：`receipt.json`（已复核）
- 脚本：`scripts/external_score_gse175540.py`（freeze/run/verify；重跑只校验不重算）
- 同伴报告：`docs/GOBP_GSE175540_SCORE_20260928.md`（同一批numbers的主终点口径复述：严格raw无家族级显著；归一化适应性免疫作敏感性读；结论一致）
- 磁盘余量>1.2TB；CPU 本地；无 GPU、新数据、旧模型重训。
