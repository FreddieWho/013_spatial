# 同原始基因输入的空间对照（描述性）

固定候选与控制，无调参、无通过阈值、无 p 值；只读冻结输入。
范围绝对偏差先在片内取绝对值，再按完整单位合同聚合。

## 臂级宏平均（完整单位等权）

| 模型 | 目标 | 臂 | 单位数 | 局部细节Pearson | q75 Jaccard | 质心位移 | 半径绝对偏差 | sd比 | RMSE |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| point_ridge | B_AXIS | internal | 65 | 0.3449 | 0.4223 | 1.1726 | 3.4572 | 0.9256 | 0.1039 |
| point_ridge_smooth_half | B_AXIS | internal | 65 | 0.3506 | 0.4513 | 1.1505 | 3.8602 | 0.7293 | 0.0976 |
| local_ridge | B_AXIS | internal | 65 | 0.3283 | 0.4260 | 1.1863 | 4.0149 | 0.9505 | 0.1046 |
| background_ridge | B_AXIS | internal | 65 | 0.3635 | 0.4336 | 1.1235 | 3.4882 | 0.8146 | 0.1044 |
| permuted_ridge | B_AXIS | internal | 65 | 0.3108 | 0.4026 | 1.2197 | 4.0888 | 0.8993 | 0.1078 |
| point_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | 65 | 0.6773 | 0.6470 | 0.3147 | 2.9162 | 0.9690 | 0.0601 |
| point_ridge_smooth_half | HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | 65 | 0.6493 | 0.6350 | 0.5197 | 3.3819 | 0.8769 | 0.0623 |
| local_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | 65 | 0.6933 | 0.6427 | 0.3257 | 3.2155 | 0.9642 | 0.0562 |
| background_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | 65 | 0.6985 | 0.6554 | 0.3143 | 3.0557 | 0.9456 | 0.0585 |
| permuted_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | 65 | 0.6565 | 0.6288 | 0.3407 | 2.8712 | 0.9575 | 0.0645 |
| point_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | 65 | 0.4904 | 0.4976 | 0.6564 | 3.6052 | 0.9607 | 0.0648 |
| point_ridge_smooth_half | HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | 65 | 0.4857 | 0.5236 | 0.6966 | 4.4710 | 0.8185 | 0.0622 |
| local_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | 65 | 0.5390 | 0.5245 | 0.5550 | 3.5699 | 0.9507 | 0.0593 |
| background_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | 65 | 0.5299 | 0.5129 | 0.6911 | 3.8394 | 0.9011 | 0.0629 |
| permuted_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | 65 | 0.4781 | 0.4838 | 0.7238 | 3.5820 | 0.9348 | 0.0688 |
| point_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | 65 | 0.6252 | 0.6084 | 0.6522 | 2.4576 | 0.9900 | 0.0805 |
| point_ridge_smooth_half | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | 65 | 0.6117 | 0.6190 | 0.8419 | 2.8316 | 0.8263 | 0.0783 |
| local_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | 65 | 0.6485 | 0.6239 | 0.5824 | 2.2633 | 0.9710 | 0.0777 |
| background_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | 65 | 0.6465 | 0.6191 | 0.6313 | 2.2885 | 0.9364 | 0.0749 |
| permuted_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | 65 | 0.6095 | 0.5980 | 0.6802 | 2.3512 | 0.9629 | 0.0857 |
| point_ridge | B_AXIS | source | 47 | 0.2550 | 0.3645 | 1.3126 | 3.5877 | 0.9196 | 0.1135 |
| point_ridge_smooth_half | B_AXIS | source | 47 | 0.2670 | 0.4044 | 1.1945 | 4.2325 | 0.6942 | 0.1049 |
| local_ridge | B_AXIS | source | 47 | 0.2254 | 0.3366 | 1.6256 | 4.6701 | 1.1807 | 0.1493 |
| background_ridge | B_AXIS | source | 47 | 0.2704 | 0.3748 | 1.2750 | 3.9414 | 0.8659 | 0.1142 |
| permuted_ridge | B_AXIS | source | 47 | 0.1925 | 0.3185 | 1.4448 | 4.9969 | 1.1159 | 0.1905 |
| point_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | 47 | 0.6279 | 0.6298 | 0.4309 | 3.1939 | 0.9412 | 0.0692 |
| point_ridge_smooth_half | HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | 47 | 0.6144 | 0.6280 | 0.6333 | 3.5313 | 0.8532 | 0.0696 |
| local_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | 47 | 0.5938 | 0.5914 | 0.5103 | 3.6054 | 0.9984 | 0.0765 |
| background_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | 47 | 0.6379 | 0.6277 | 0.4521 | 3.6614 | 0.8877 | 0.0716 |
| permuted_ridge | HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | 47 | 0.5244 | 0.5649 | 0.5187 | 3.2065 | 0.9267 | 0.1023 |
| point_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | source | 47 | 0.3855 | 0.4575 | 0.7741 | 3.8992 | 0.8801 | 0.0703 |
| point_ridge_smooth_half | HALLMARK_INTERFERON_GAMMA_RESPONSE | source | 47 | 0.3931 | 0.5071 | 0.8217 | 3.4292 | 0.7265 | 0.0664 |
| local_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | source | 47 | 0.4185 | 0.4351 | 0.9011 | 3.4508 | 1.0363 | 0.0927 |
| background_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | source | 47 | 0.4378 | 0.4882 | 0.6821 | 3.6418 | 0.8738 | 0.0689 |
| permuted_ridge | HALLMARK_INTERFERON_GAMMA_RESPONSE | source | 47 | 0.3309 | 0.4162 | 0.7839 | 4.7888 | 1.0079 | 0.0819 |
| point_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | 47 | 0.5766 | 0.5975 | 0.7924 | 2.9156 | 0.9875 | 0.0877 |
| point_ridge_smooth_half | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | 47 | 0.5755 | 0.6222 | 0.8993 | 2.9916 | 0.8225 | 0.0835 |
| local_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | 47 | 0.5530 | 0.5625 | 1.0154 | 2.5596 | 1.0645 | 0.1263 |
| background_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | 47 | 0.5771 | 0.5989 | 0.7925 | 2.9078 | 0.9664 | 0.0760 |
| permuted_ridge | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | 47 | 0.4796 | 0.5381 | 0.9049 | 2.9491 | 1.0822 | 0.1047 |

## 候选减基线（单位等权配对差值）

| 目标 | 臂 | 基线 | 局部细节Δ | q75 Jaccard Δ | 质心位移 Δ | 半径绝对偏差 Δ |
|---|---|---|---:|---:|---:|---:|
| B_AXIS | internal | point_ridge | -0.0167 | +0.0037 | +0.0138 | +0.5577 |
| B_AXIS | internal | point_ridge_smooth_half | -0.0223 | -0.0254 | +0.0358 | +0.1547 |
| B_AXIS | internal | background_ridge | -0.0353 | -0.0076 | +0.0629 | +0.5267 |
| B_AXIS | internal | permuted_ridge | +0.0174 | +0.0234 | -0.0334 | -0.0739 |
| B_AXIS | source | point_ridge | -0.0297 | -0.0279 | +0.3129 | +1.0824 |
| B_AXIS | source | point_ridge_smooth_half | -0.0416 | -0.0678 | +0.4311 | +0.4375 |
| B_AXIS | source | background_ridge | -0.0450 | -0.0381 | +0.3505 | +0.7287 |
| B_AXIS | source | permuted_ridge | +0.0329 | +0.0181 | +0.1807 | -0.3268 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | point_ridge | +0.0160 | -0.0043 | +0.0111 | +0.2993 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | point_ridge_smooth_half | +0.0440 | +0.0077 | -0.1940 | -0.1664 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | background_ridge | -0.0052 | -0.0127 | +0.0114 | +0.1598 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | permuted_ridge | +0.0368 | +0.0139 | -0.0150 | +0.3443 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | point_ridge | -0.0341 | -0.0384 | +0.0794 | +0.4115 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | point_ridge_smooth_half | -0.0206 | -0.0365 | -0.1230 | +0.0741 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | background_ridge | -0.0441 | -0.0363 | +0.0582 | -0.0560 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | permuted_ridge | +0.0694 | +0.0265 | -0.0084 | +0.3989 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | point_ridge | +0.0485 | +0.0269 | -0.1014 | -0.0352 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | point_ridge_smooth_half | +0.0533 | +0.0009 | -0.1415 | -0.9011 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | background_ridge | +0.0091 | +0.0116 | -0.1361 | -0.2695 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | permuted_ridge | +0.0609 | +0.0407 | -0.1688 | -0.0121 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | point_ridge | +0.0331 | -0.0223 | +0.1270 | -0.4484 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | point_ridge_smooth_half | +0.0254 | -0.0719 | +0.0793 | +0.0216 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | background_ridge | -0.0192 | -0.0530 | +0.2190 | -0.1910 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | permuted_ridge | +0.0877 | +0.0190 | +0.1171 | -1.3380 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | point_ridge | +0.0233 | +0.0155 | -0.0698 | -0.1943 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | point_ridge_smooth_half | +0.0368 | +0.0049 | -0.2595 | -0.5683 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | background_ridge | +0.0020 | +0.0047 | -0.0489 | -0.0252 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | permuted_ridge | +0.0390 | +0.0258 | -0.0978 | -0.0879 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | point_ridge | -0.0236 | -0.0350 | +0.2230 | -0.3560 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | point_ridge_smooth_half | -0.0226 | -0.0597 | +0.1161 | -0.4320 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | background_ridge | -0.0241 | -0.0364 | +0.2230 | -0.3482 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | permuted_ridge | +0.0734 | +0.0244 | +0.1105 | -0.3895 |

## 仅描述：候选是否同时优于四个对照（不构成阈值或显著性）

- internal / B_AXIS / rec_local_pattern_pearson：同时有利 = False
- internal / B_AXIS / rec_jaccard_q75：同时有利 = False
- internal / HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION / rec_local_pattern_pearson：同时有利 = True
- internal / HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION / rec_jaccard_q75：同时有利 = True
- internal / HALLMARK_INTERFERON_GAMMA_RESPONSE / rec_local_pattern_pearson：同时有利 = True
- internal / HALLMARK_INTERFERON_GAMMA_RESPONSE / rec_jaccard_q75：同时有利 = True
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

