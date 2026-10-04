# 完整图表与源数据索引

全部图表为描述性或内部验证证据；未确认分子桥或正式双向关系。完整PDF曲线图集不以代表通路替代全量。

## 全部前向曲线

[逐程序图集索引](forward_atlas/index.tsv)：15册PDF，共1742页（1691 GO及51对照）；每页为两来源晕、TUM侧及明确非TUM侧四面板，两种核心隔离带。灰/空缺是不可测，曲线不跨空箱连接。[逐页内容审计](forward_atlas/pdf_content_audit.json)检查全部页程序名与面板；[B轴预览](forward_atlas/B_AXIS_preview.pdf)。数值见上一级`integrated/all_bins.tsv`和`halo_edge_bins.npz`。

## 完整桥矩阵

每张包含全部1742项×1559组；行与列分别见[程序索引](program_row_index.tsv)和[几何组索引](../integrated/bridge_groups.tsv)。

- [原始中心与侧带差](bridge_all_raw_median.pdf)
- [端点晕调整后差](bridge_all_adjusted_median.pdf)
- [相对原位参考](bridge_all_native_z_median.pdf)：诊断z，不是p。
- [正向沿程一致性](bridge_all_positive_continuity.pdf)
- [负向沿程一致性](bridge_all_negative_continuity.pdf)
- [两侧不对称](bridge_all_side_asymmetry.pdf)

数值见`../integrated/bridge_all_programs.npz`。含宽片区及全部已测尺度，不将所有合格几何命名为桥。

## 全来源代表空间图

三项读数均展示全部26片，共78面板；按D-180选择，属于探索性展示。已知TLS核心及130μm隔离带隐藏，灰色表示资格不足；线只表示按metadata几何选择的候选，不是已发现桥。

- B轴固定对照：[USZ](representative_01_USZ.pdf) / [GSE175540](representative_01_GSE175540.pdf)。
- γδ T细胞活化GO：[USZ](representative_02_USZ.pdf) / [GSE175540](representative_02_GSE175540.pdf)。
- B细胞受体信号GO：[USZ](representative_03_USZ.pdf) / [GSE175540](representative_03_GSE175540.pdf)。

[候选几何选择记录](representative_candidate_geometry.tsv)。GO名称不证明具体细胞类型或细胞内状态。

## 反向预测与尺度支持

- [多尺度可测性](multiscale_support.pdf)：包含短距离；点数不足不是没有桥。
- [逐切片留出](reverse_leave_section.pdf)与[跨来源留出](reverse_leave_cohort.pdf)：全部三mask和两来源，不声称患者独立。
- [全部1691 GO](inverse_all_GO.pdf)与[51固定对照](inverse_all_controls.pdf)：径向、角向、联合三臂，灰色不可测。

反向源矩阵为`inverse_heatmap_source.npz`，完整表位于`../reverse/summary/`。总览及代表图提供同名SVG、PDF、600dpi TIFF和PNG；全量曲线按PDF分册交付。图表规则见`figure_contract.json`和`forward_figure_contract.json`，核查见`visual_qa.json`、`forward_visual_qa.json`及各完成凭证。目视检查按布局抽查，完整曲线另外逐页做文本内容审计。
