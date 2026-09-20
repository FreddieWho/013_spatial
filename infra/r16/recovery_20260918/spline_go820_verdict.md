# Spline × 820 双队列通路：形态由切片主导，通路身份被抹平

6560 fits（820 通路 × 8 USZ 片），口径同 D-137（ns df=3 ~ TLS 距离，TUM spots，Spearman+15% 地板）。

## 主结论

**形态几乎完全由切片决定，不由通路决定。** 没有任何一条通路在 ≥3 片上稳定 high-to-low 或 low-to-high。

| 切片 | 主导形态 | 纯度 |
|---|---|---|
| LC3 | nonmonotonic | 100% |
| LC4 | nonmonotonic | 99.4% |
| LC5 | nonmonotonic | 85% |
| KC3 | flat | 99.9% |
| LC2 | flat | 99.9% |
| KC1 | high-to-low | 71% |
| LC1 | high-to-low | 57% |
| KC2 | low-to-high | 63% |

通路多数票：nonmonotonic 725 / flat 95 / 单调 0。LC3 上随机通路对的曲线 Spearman ρ≈1.0——820 条曲线是同一条形状乘以不同振幅。

## 对 L-011 晕环的修订

LC3 上从 actin 到 Wnt 到 wound healing，归一化曲线都是 `1.00 → 0.00 → 0.21 → 0.45 → 0.32`。L-011 的 peri-TLS 晕环是 **LC3/LC4 切片几何**，不是 mhc2 特异。mhc2 振幅仍大于 B 轴（D-137 数字仍在），但形状不是它的身份。

## 原因（方法，不是生物学）

Spline 拟合的是 **未残差的 log1p 均值**，没做 Q/C 回归。深度、细胞密度、肿瘤-基质混沿 TLS 距离的共变被所有通路共享。论文用 AUCell（秩）+ tumor-annotated spots；我们限了 TUM，但没去掉组成/深度。GO 基因高度重叠加重了共变。

## 不支持的主张

- 「820 通路里普遍存在 TLS 距离梯度」——不支持（切片主导）。
- 「晕环是 mhc2 的组织化信息」——降级为切片描述（L-011 修订）。
- 灵感论文的 8 癌种保守 high-to-low：在本 8 片 USZ 上未复现为通路特异。

## 若继续

必须先对程序分做 Q+C 残差再拟合 spline，否则再铺 6870 条也是同一张切片形状。残差后若通路形状仍共线，才是基因重叠/组成不可分。
