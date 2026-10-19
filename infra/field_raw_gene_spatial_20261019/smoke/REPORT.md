# 同原始基因输入的空间对照（描述性）

固定候选与控制，无调参、无通过阈值、无 p 值；只读冻结输入。
范围绝对偏差先在片内取绝对值，再按完整单位合同聚合。

## 臂级宏平均（完整单位等权）

| 模型 | 目标 | 臂 | 单位数 | 局部细节Pearson | q75 Jaccard | 质心位移 | 半径绝对偏差 | sd比 | RMSE |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| point_ridge | B_AXIS | internal | 2 | 0.2525 | 0.3161 | 1.7195 | 10.0629 | 0.9025 | 0.1037 |
| point_ridge_smooth_half | B_AXIS | internal | 2 | 0.2585 | 0.3606 | 1.6533 | 3.9967 | 0.6578 | 0.0891 |
| local_ridge | B_AXIS | internal | 2 | 0.1498 | 0.2620 | 1.8040 | 7.0673 | 1.0463 | 0.1316 |
| background_ridge | B_AXIS | internal | 2 | 0.2552 | 0.3085 | 1.8663 | 6.1802 | 0.8626 | 0.1107 |
| permuted_ridge | B_AXIS | internal | 2 | 0.1400 | 0.2439 | 1.7708 | 12.0682 | 1.2797 | 0.1780 |
| point_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | 2 | 0.7893 | 0.6752 | 0.5514 | 1.5512 | 0.9671 | 0.0675 |
| point_ridge_smooth_half | HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | 2 | 0.7669 | 0.6704 | 0.2922 | 3.7273 | 0.8522 | 0.0743 |
| local_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | 2 | 0.6886 | 0.5877 | 0.5148 | 2.7561 | 1.0434 | 0.0887 |
| background_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | 2 | 0.7991 | 0.6866 | 0.3438 | 2.1282 | 0.9323 | 0.0863 |
| permuted_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | 2 | 0.6182 | 0.5373 | 0.8216 | 2.8651 | 1.0082 | 0.1228 |
| point_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | 2 | 0.5889 | 0.5052 | 0.5837 | 3.3127 | 0.9638 | 0.0661 |
| point_ridge_smooth_half | HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | 2 | 0.5839 | 0.5578 | 0.7510 | 7.6739 | 0.8250 | 0.0612 |
| local_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | 2 | 0.4763 | 0.4510 | 0.6197 | 6.0030 | 1.0234 | 0.0795 |
| background_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | 2 | 0.6078 | 0.5071 | 0.5516 | 7.3186 | 0.9096 | 0.0831 |
| permuted_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | 2 | 0.4118 | 0.4008 | 0.9822 | 8.5465 | 1.0939 | 0.0882 |
| point_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | 2 | 0.6299 | 0.6105 | 1.2049 | 3.7635 | 1.0451 | 0.0959 |
| point_ridge_smooth_half | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | 2 | 0.6299 | 0.6417 | 1.2926 | 2.7321 | 0.8848 | 0.0857 |
| local_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | 2 | 0.5153 | 0.5434 | 1.2392 | 6.2459 | 1.1220 | 0.1132 |
| background_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | 2 | 0.6343 | 0.6172 | 1.1489 | 3.8223 | 1.0209 | 0.0938 |
| permuted_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | 2 | 0.4413 | 0.4786 | 1.2896 | 5.4300 | 1.1744 | 0.1372 |
| point_ridge | B_AXIS | source | 2 | 0.1501 | 0.3183 | 0.6062 | 2.7524 | 1.7072 | 0.0676 |
| point_ridge_smooth_half | B_AXIS | source | 2 | 0.1546 | 0.3480 | 0.9347 | 1.7083 | 1.2308 | 0.0505 |
| local_ridge | B_AXIS | source | 2 | 0.0960 | 0.2003 | 1.4275 | 5.5902 | 2.3238 | 0.2878 |
| background_ridge | B_AXIS | source | 2 | 0.1496 | 0.3249 | 0.7462 | 3.0142 | 1.7610 | 0.0996 |
| permuted_ridge | B_AXIS | source | 2 | 0.0295 | 0.2070 | 0.4310 | 10.2719 | 2.8646 | 0.4390 |
| point_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | 2 | 0.3456 | 0.2564 | 0.5752 | 9.8112 | 0.8165 | 0.2729 |
| point_ridge_smooth_half | HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | 2 | 0.3516 | 0.2891 | 0.7016 | 6.8365 | 0.5880 | 0.2714 |
| local_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | 2 | 0.2697 | 0.1952 | 0.8590 | 12.3294 | 1.1359 | 0.6831 |
| background_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | 2 | 0.3513 | 0.2684 | 0.5325 | 13.0057 | 0.7877 | 0.7703 |
| permuted_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | 2 | 0.1835 | 0.1965 | 0.6212 | 14.7181 | 1.2911 | 0.1013 |
| point_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | source | 2 | 0.2833 | 0.2765 | 0.3265 | 1.0137 | 1.6020 | 0.1208 |
| point_ridge_smooth_half | HALLMARK_INTERFERON_GAMMA_RESPONSE | source | 2 | 0.2864 | 0.3085 | 0.6015 | 4.3232 | 1.1154 | 0.1172 |
| local_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | source | 2 | 0.1887 | 0.1803 | 0.7895 | 3.0852 | 2.2623 | 0.0972 |
| background_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | source | 2 | 0.2945 | 0.2801 | 0.3120 | 4.5620 | 1.5298 | 0.4958 |
| permuted_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | source | 2 | 0.1365 | 0.2175 | 0.1727 | 5.4892 | 2.6644 | 0.3715 |
| point_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | 2 | 0.3495 | 0.4082 | 0.4627 | 1.0187 | 1.1338 | 0.1991 |
| point_ridge_smooth_half | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | 2 | 0.3685 | 0.4597 | 0.8533 | 0.5768 | 0.8240 | 0.1958 |
| local_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | 2 | 0.2304 | 0.2768 | 0.9174 | 5.6048 | 1.5119 | 0.3169 |
| background_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | 2 | 0.3458 | 0.4106 | 0.4364 | 1.1167 | 1.0882 | 0.3598 |
| permuted_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | 2 | 0.1689 | 0.2787 | 0.5950 | 7.5931 | 1.8600 | 0.1526 |

## 候选减基线（单位等权配对差值）

| 目标 | 臂 | 基线 | 局部细节Δ | q75 Jaccard Δ | 质心位移 Δ | 半径绝对偏差 Δ |
|---|---|---|---:|---:|---:|---:|
| B_AXIS | internal | point_ridge | -0.1027 | -0.0541 | +0.0845 | -2.9956 |
| B_AXIS | internal | point_ridge_smooth_half | -0.1086 | -0.0987 | +0.1507 | +3.0706 |
| B_AXIS | internal | background_ridge | -0.1054 | -0.0465 | -0.0624 | +0.8871 |
| B_AXIS | internal | permuted_ridge | +0.0098 | +0.0181 | +0.0332 | -5.0009 |
| B_AXIS | source | point_ridge | -0.0540 | -0.1180 | +0.8213 | +2.8377 |
| B_AXIS | source | point_ridge_smooth_half | -0.0586 | -0.1476 | +0.4928 | +3.8819 |
| B_AXIS | source | background_ridge | -0.0536 | -0.1246 | +0.6813 | +2.5760 |
| B_AXIS | source | permuted_ridge | +0.0666 | -0.0067 | +0.9965 | -4.6818 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | point_ridge | -0.1007 | -0.0874 | -0.0366 | +1.2050 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | point_ridge_smooth_half | -0.0784 | -0.0827 | +0.2225 | -0.9711 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | background_ridge | -0.1106 | -0.0988 | +0.1710 | +0.6279 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | permuted_ridge | +0.0704 | +0.0504 | -0.3068 | -0.1089 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | point_ridge | -0.0759 | -0.0613 | +0.2838 | +2.5182 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | point_ridge_smooth_half | -0.0819 | -0.0939 | +0.1574 | +5.4929 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | background_ridge | -0.0816 | -0.0732 | +0.3265 | -0.6764 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | permuted_ridge | +0.0862 | -0.0013 | +0.2378 | -2.3887 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | point_ridge | -0.1127 | -0.0542 | +0.0360 | +2.6903 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | point_ridge_smooth_half | -0.1077 | -0.1068 | -0.1313 | -1.6709 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | background_ridge | -0.1315 | -0.0560 | +0.0682 | -1.3156 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | permuted_ridge | +0.0644 | +0.0502 | -0.3624 | -2.5435 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | point_ridge | -0.0946 | -0.0962 | +0.4630 | +2.0715 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | point_ridge_smooth_half | -0.0977 | -0.1282 | +0.1880 | -1.2380 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | background_ridge | -0.1058 | -0.0998 | +0.4774 | -1.4769 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | permuted_ridge | +0.0522 | -0.0372 | +0.6168 | -2.4041 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | point_ridge | -0.1146 | -0.0670 | +0.0342 | +2.4824 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | point_ridge_smooth_half | -0.1146 | -0.0983 | -0.0535 | +3.5138 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | background_ridge | -0.1190 | -0.0737 | +0.0902 | +2.4236 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | permuted_ridge | +0.0740 | +0.0649 | -0.0504 | +0.8159 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | point_ridge | -0.1191 | -0.1314 | +0.4547 | +4.5861 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | point_ridge_smooth_half | -0.1381 | -0.1829 | +0.0641 | +5.0280 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | background_ridge | -0.1154 | -0.1338 | +0.4809 | +4.4880 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | permuted_ridge | +0.0615 | -0.0019 | +0.3224 | -1.9883 |

## 仅描述：候选是否同时优于四个对照（不构成阈值或显著性）

- internal / B_AXIS / rec_local_pattern_pearson：同时有利 = False
- internal / B_AXIS / rec_jaccard_q75：同时有利 = False
- internal / HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION / rec_local_pattern_pearson：同时有利 = False
- internal / HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION / rec_jaccard_q75：同时有利 = False
- internal / HALLMARK_INTERFERON_GAMMA_RESPONSE / rec_local_pattern_pearson：同时有利 = False
- internal / HALLMARK_INTERFERON_GAMMA_RESPONSE / rec_jaccard_q75：同时有利 = False
- internal / HALLMARK_OXIDATIVE_PHOSPHORYLATION / rec_local_pattern_pearson：同时有利 = False
- internal / HALLMARK_OXIDATIVE_PHOSPHORYLATION / rec_jaccard_q75：同时有利 = False
- source / B_AXIS / rec_local_pattern_pearson：同时有利 = False
- source / B_AXIS / rec_jaccard_q75：同时有利 = False
- source / HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION / rec_local_pattern_pearson：同时有利 = False
- source / HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION / rec_jaccard_q75：同时有利 = False
- source / HALLMARK_INTERFERON_GAMMA_RESPONSE / rec_local_pattern_pearson：同时有利 = False
- source / HALLMARK_INTERFERON_GAMMA_RESPONSE / rec_jaccard_q75：同时有利 = False
- source / HALLMARK_OXIDATIVE_PHOSPHORYLATION / rec_local_pattern_pearson：同时有利 = False
- source / HALLMARK_OXIDATIVE_PHOSPHORYLATION / rec_jaccard_q75：同时有利 = False

证据角色：既有开发材料上的探索性机制对照；预测与输入共享可见基因，含输入回声；
打乱邻域为单次机制破坏，不作 p 值或正式零模型；未同时优于四个对照的读出不得宣称依赖真实邻域排列。
NOT_RUN：人工 TLS、独立 marker 隔离、未知患者身份、正式空间显著性、物理尺度、视野外推断、
打乱参考描述（本轮关闭 shuffle 计算）。

