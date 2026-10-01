# 无效产物：CRC 参考的 NNLS 组成

这些 `*.npz` 来自 v9.0/v9.1-nnls：用 R-04 冻结的 **GSE236581 CRC 单细胞** 参考
对 USZ（肾/肺）与 GSE175540（ccRCC）切片做 per-spot NNLS。结果退化
（usz-LC3 均值 B=0.63、Epi=0.00），判为**领域不匹配、无效**，不承担任何结论。
有效版本（v9.1 marker-proxy 协变量）不产生 per-spot npz。
详见 docs/COMPOSITION_CONTROL_FIELD_20261001.md 与 D-166。
