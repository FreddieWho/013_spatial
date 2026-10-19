# 同原始基因输入空间对照：描述性结果表（REVIEWED：只允许描述性报告，不得声称真实邻域排列带来恢复；见 review/result_comparison_review.json）

来源：`infra/field_raw_gene_spatial_20261019/` 的 evaluate/summarize 产物；同一批 500 个预测包、5 个模型；内部 65 单位 / 来源 47 单位；无显著性、无 PASS 阈值。
已知数值边界：single-case（source_0/EMT point）1/2 线程 vs 8 线程的逐点差约 2.1e-4、每片 RMSE 差约 1.9e-4；该证据只覆盖这一个 case 的 point 模型，未覆盖 candidate/background/permuted，也不假定同批误差相消。下表中 |Δ| ≤ 2e-4 的配对差标 `~num`（处在已知数值波动内，结论未确定）。

## 1. 覆盖与分母（完整 / 部分 / 缺失单位组）

| 臂 | 模型 | 单位组 | 主读出可评 | 部分 | 缺失 |
|---|---|---:|---:|---:|---:|
| internal | point_ridge | 260 | 260 | 0 | 0 |
| internal | point_ridge_smooth_half | 260 | 260 | 0 | 0 |
| internal | local_ridge | 260 | 260 | 0 | 0 |
| internal | background_ridge | 260 | 260 | 0 | 0 |
| internal | permuted_ridge | 260 | 260 | 0 | 0 |
| source | point_ridge | 188 | 188 | 0 | 0 |
| source | point_ridge_smooth_half | 188 | 188 | 0 | 0 |
| source | local_ridge | 188 | 188 | 0 | 0 |
| source | background_ridge | 188 | 188 | 0 | 0 |
| source | permuted_ridge | 188 | 188 | 0 | 0 |

## 2. 主读出与次读出（臂级单位等权宏平均；内部/来源分开）

| 目标 | 臂 | 模型 | local-detail Pearson | q75 Jaccard | q80 / q90 Jaccard | 质心位移(pitch) | 半径绝对偏差(pitch) | sd比 | RMSE | 去均值RMSE |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| B_AXIS | internal | point_ridge | 0.3449 | 0.4223 | 0.3974 / 0.3435 | 1.1726 | 3.4572 | 0.9256 | 0.1039 | 0.0882 |
| B_AXIS | internal | point_ridge_smooth_half | 0.3506 | 0.4513 | 0.4270 / 0.3720 | 1.1505 | 3.8602 | 0.7293 | 0.0976 | 0.0805 |
| B_AXIS | internal | local_ridge | 0.3283 | 0.4260 | 0.4029 / 0.3466 | 1.1863 | 4.0149 | 0.9505 | 0.1046 | 0.0880 |
| B_AXIS | internal | background_ridge | 0.3635 | 0.4336 | 0.4078 / 0.3523 | 1.1235 | 3.4882 | 0.8146 | 0.1044 | 0.0835 |
| B_AXIS | internal | permuted_ridge | 0.3108 | 0.4026 | 0.3778 / 0.3239 | 1.2197 | 4.0888 | 0.8993 | 0.1078 | 0.0907 |
| B_AXIS | source | point_ridge | 0.2550 | 0.3645 | 0.3457 / 0.3020 | 1.3126 | 3.5877 | 0.9196 | 0.1135 | 0.0867 |
| B_AXIS | source | point_ridge_smooth_half | 0.2670 | 0.4044 | 0.3843 / 0.3378 | 1.1945 | 4.2325 | 0.6942 | 0.1049 | 0.0760 |
| B_AXIS | source | local_ridge | 0.2254 | 0.3366 | 0.3158 / 0.2676 | 1.6256 | 4.6701 | 1.1807 | 0.1493 | 0.1039 |
| B_AXIS | source | background_ridge | 0.2704 | 0.3748 | 0.3543 / 0.3109 | 1.2750 | 3.9414 | 0.8659 | 0.1142 | 0.0835 |
| B_AXIS | source | permuted_ridge | 0.1925 | 0.3185 | 0.2942 / 0.2452 | 1.4448 | 4.9969 | 1.1159 | 0.1905 | 0.1041 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | point_ridge | 0.6252 | 0.6084 | 0.5829 / 0.5192 | 0.6522 | 2.4576 | 0.9900 | 0.0805 | 0.0717 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | point_ridge_smooth_half | 0.6117 | 0.6190 | 0.5935 / 0.5241 | 0.8419 | 2.8316 | 0.8263 | 0.0783 | 0.0693 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | local_ridge | 0.6485 | 0.6239 | 0.5952 / 0.5300 | 0.5824 | 2.2633 | 0.9710 | 0.0777 | 0.0681 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | background_ridge | 0.6465 | 0.6191 | 0.5951 / 0.5289 | 0.6313 | 2.2885 | 0.9364 | 0.0749 | 0.0684 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | permuted_ridge | 0.6095 | 0.5980 | 0.5720 / 0.5066 | 0.6802 | 2.3512 | 0.9629 | 0.0857 | 0.0733 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | point_ridge | 0.5766 | 0.5975 | 0.5752 / 0.5113 | 0.7924 | 2.9156 | 0.9875 | 0.0877 | 0.0692 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | point_ridge_smooth_half | 0.5755 | 0.6222 | 0.5995 / 0.5294 | 0.8993 | 2.9916 | 0.8225 | 0.0835 | 0.0634 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | local_ridge | 0.5530 | 0.5625 | 0.5384 / 0.4761 | 1.0154 | 2.5596 | 1.0645 | 0.1263 | 0.0807 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | background_ridge | 0.5771 | 0.5989 | 0.5778 / 0.5127 | 0.7925 | 2.9078 | 0.9664 | 0.0760 | 0.0685 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | permuted_ridge | 0.4796 | 0.5381 | 0.5153 / 0.4482 | 0.9049 | 2.9491 | 1.0822 | 0.1047 | 0.0862 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | point_ridge | 0.4904 | 0.4976 | 0.4677 / 0.4013 | 0.6564 | 3.6052 | 0.9607 | 0.0648 | 0.0530 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | point_ridge_smooth_half | 0.4857 | 0.5236 | 0.4937 / 0.4255 | 0.6966 | 4.4710 | 0.8185 | 0.0622 | 0.0498 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | local_ridge | 0.5390 | 0.5245 | 0.4959 / 0.4311 | 0.5550 | 3.5699 | 0.9507 | 0.0593 | 0.0498 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | background_ridge | 0.5299 | 0.5129 | 0.4825 / 0.4107 | 0.6911 | 3.8394 | 0.9011 | 0.0629 | 0.0503 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | permuted_ridge | 0.4781 | 0.4838 | 0.4506 / 0.3800 | 0.7238 | 3.5820 | 0.9348 | 0.0688 | 0.0540 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | point_ridge | 0.3855 | 0.4575 | 0.4292 / 0.3484 | 0.7741 | 3.8992 | 0.8801 | 0.0703 | 0.0526 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | point_ridge_smooth_half | 0.3931 | 0.5071 | 0.4781 / 0.3957 | 0.8217 | 3.4292 | 0.7265 | 0.0664 | 0.0475 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | local_ridge | 0.4185 | 0.4351 | 0.3996 / 0.3130 | 0.9011 | 3.4508 | 1.0363 | 0.0927 | 0.0580 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | background_ridge | 0.4378 | 0.4882 | 0.4565 / 0.3739 | 0.6821 | 3.6418 | 0.8738 | 0.0689 | 0.0480 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | permuted_ridge | 0.3309 | 0.4162 | 0.3790 / 0.2934 | 0.7839 | 4.7888 | 1.0079 | 0.0819 | 0.0598 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | point_ridge | 0.6773 | 0.6470 | 0.6042 / 0.4848 | 0.3147 | 2.9162 | 0.9690 | 0.0601 | 0.0495 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | point_ridge_smooth_half | 0.6493 | 0.6350 | 0.5923 / 0.4702 | 0.5197 | 3.3819 | 0.8769 | 0.0623 | 0.0519 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | local_ridge | 0.6933 | 0.6427 | 0.5983 / 0.4780 | 0.3257 | 3.2155 | 0.9642 | 0.0562 | 0.0489 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | background_ridge | 0.6985 | 0.6554 | 0.6112 / 0.4892 | 0.3143 | 3.0557 | 0.9456 | 0.0585 | 0.0480 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | permuted_ridge | 0.6565 | 0.6288 | 0.5818 / 0.4566 | 0.3407 | 2.8712 | 0.9575 | 0.0645 | 0.0514 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | point_ridge | 0.6279 | 0.6298 | 0.5912 / 0.4802 | 0.4309 | 3.1939 | 0.9412 | 0.0692 | 0.0543 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | point_ridge_smooth_half | 0.6144 | 0.6280 | 0.5923 / 0.4745 | 0.6333 | 3.5313 | 0.8532 | 0.0696 | 0.0549 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | local_ridge | 0.5938 | 0.5914 | 0.5487 / 0.4373 | 0.5103 | 3.6054 | 0.9984 | 0.0765 | 0.0630 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | background_ridge | 0.6379 | 0.6277 | 0.5875 / 0.4753 | 0.4521 | 3.6614 | 0.8877 | 0.0716 | 0.0550 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | permuted_ridge | 0.5244 | 0.5649 | 0.5181 / 0.3974 | 0.5187 | 3.2065 | 0.9267 | 0.1023 | 0.0651 |

注：有符号半径差（`rec_patch75/90_largest_radius_pitch_delta`）只作描述，负值不代表范围改善（C-01）；范围结论一律以“每片绝对值再按完整单位聚合”的读数为准。

## 3. 候选减四基线（臂级单位等权配对差；C-01 方向标注）

| 目标 | 臂 | 基线 | 判据 | local-detail Δ | q75 Jaccard Δ | 质心位移 Δ | 半径绝对偏差 Δ | RMSE Δ |
|---|---|---|---|---:|---:|---:|---:|---:|
| B_AXIS | internal | point_ridge | beyond-2e-4（n_est=65） | -0.0167 | +0.0037 | +0.0138 | +0.5577 | +0.0007 |
| B_AXIS | internal | point_ridge_smooth_half | beyond-2e-4（n_est=65） | -0.0223 | -0.0254 | +0.0358 | +0.1547 | +0.0069 |
| B_AXIS | internal | background_ridge | ~num（n_est=65） | -0.0353 | -0.0076 | +0.0629 | +0.5267 | +0.0002 ~num |
| B_AXIS | internal | permuted_ridge | beyond-2e-4（n_est=65） | +0.0174 | +0.0234 | -0.0334 | -0.0739 | -0.0033 |
| B_AXIS | source | point_ridge | beyond-2e-4（n_est=47） | -0.0297 | -0.0279 | +0.3129 | +1.0824 | +0.0358 |
| B_AXIS | source | point_ridge_smooth_half | beyond-2e-4（n_est=47） | -0.0416 | -0.0678 | +0.4311 | +0.4375 | +0.0444 |
| B_AXIS | source | background_ridge | beyond-2e-4（n_est=47） | -0.0450 | -0.0381 | +0.3505 | +0.7287 | +0.0350 |
| B_AXIS | source | permuted_ridge | beyond-2e-4（n_est=47） | +0.0329 | +0.0181 | +0.1807 | -0.3268 | -0.0412 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | point_ridge | beyond-2e-4（n_est=65） | +0.0233 | +0.0155 | -0.0698 | -0.1943 | -0.0028 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | point_ridge_smooth_half | beyond-2e-4（n_est=65） | +0.0368 | +0.0049 | -0.2595 | -0.5683 | -0.0006 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | background_ridge | beyond-2e-4（n_est=65） | +0.0020 | +0.0047 | -0.0489 | -0.0252 | +0.0029 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | internal | permuted_ridge | beyond-2e-4（n_est=65） | +0.0390 | +0.0258 | -0.0978 | -0.0879 | -0.0080 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | point_ridge | beyond-2e-4（n_est=47） | -0.0236 | -0.0350 | +0.2230 | -0.3560 | +0.0386 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | point_ridge_smooth_half | beyond-2e-4（n_est=47） | -0.0226 | -0.0597 | +0.1161 | -0.4320 | +0.0428 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | background_ridge | beyond-2e-4（n_est=47） | -0.0241 | -0.0364 | +0.2230 | -0.3482 | +0.0503 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | source | permuted_ridge | beyond-2e-4（n_est=47） | +0.0734 | +0.0244 | +0.1105 | -0.3895 | +0.0216 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | point_ridge | beyond-2e-4（n_est=65） | +0.0485 | +0.0269 | -0.1014 | -0.0352 | -0.0056 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | point_ridge_smooth_half | beyond-2e-4（n_est=65） | +0.0533 | +0.0009 | -0.1415 | -0.9011 | -0.0029 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | background_ridge | beyond-2e-4（n_est=65） | +0.0091 | +0.0116 | -0.1361 | -0.2695 | -0.0036 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | internal | permuted_ridge | beyond-2e-4（n_est=65） | +0.0609 | +0.0407 | -0.1688 | -0.0121 | -0.0095 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | point_ridge | beyond-2e-4（n_est=47） | +0.0331 | -0.0223 | +0.1270 | -0.4484 | +0.0224 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | point_ridge_smooth_half | beyond-2e-4（n_est=47） | +0.0254 | -0.0719 | +0.0793 | +0.0216 | +0.0263 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | background_ridge | beyond-2e-4（n_est=47） | -0.0192 | -0.0530 | +0.2190 | -0.1910 | +0.0238 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | source | permuted_ridge | beyond-2e-4（n_est=47） | +0.0877 | +0.0190 | +0.1171 | -1.3380 | +0.0108 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | point_ridge | beyond-2e-4（n_est=65） | +0.0160 | -0.0043 | +0.0111 | +0.2993 | -0.0039 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | point_ridge_smooth_half | beyond-2e-4（n_est=65） | +0.0440 | +0.0077 | -0.1940 | -0.1664 | -0.0061 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | background_ridge | beyond-2e-4（n_est=65） | -0.0052 | -0.0127 | +0.0114 | +0.1598 | -0.0023 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | internal | permuted_ridge | beyond-2e-4（n_est=65） | +0.0368 | +0.0139 | -0.0150 | +0.3443 | -0.0082 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | point_ridge | beyond-2e-4（n_est=47） | -0.0341 | -0.0384 | +0.0794 | +0.4115 | +0.0073 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | point_ridge_smooth_half | beyond-2e-4（n_est=47） | -0.0206 | -0.0365 | -0.1230 | +0.0741 | +0.0069 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | background_ridge | beyond-2e-4（n_est=47） | -0.0441 | -0.0363 | +0.0582 | -0.0560 | +0.0050 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | source | permuted_ridge | beyond-2e-4（n_est=47） | +0.0694 | +0.0265 | -0.0084 | +0.3989 | -0.0257 |

`~num`：该差落在单一 case 已知数值波动（2e-4）之内，按审核要求标为未确定；`beyond-2e-4` 仅表示数值上超出该单 case 波动参考，不是显著性。

## 4. 来源分层（来源臂，按 held_source；乳腺/肝脏为器官外推）

派生自 `unit_paired.tsv`（候选减同点，单位等权；每个来源内单位平均），
四目标 × 6 来源的完整行见 `RESULTS_SOURCE_TABLE.tsv`。

| 目标 | 来源 | 单位组(可评) | local-detail Δ | q75 Jaccard Δ | 质心位移 Δ | 半径绝对偏差 Δ | 器官外推 |
|---|---|---:|---:|---:|---:|---:|---|
| B_AXIS | https://www.biorxiv.org/content/10.1101/2022 | 2(2) | 0.00490 | -0.07014 | 0.37572 | 6.55691 | 是 |
| B_AXIS | https://www.biorxiv.org/content/10.1101/2023 | 6(6) | -0.00375 | -0.04064 | 0.38492 | 1.53752 | 否 |
| B_AXIS | https://www.biorxiv.org/content/10.1101/2024 | 2(2) | -0.10956 | -0.03522 | 1.18451 | 2.60525 | 否 |
| B_AXIS | https://www.cell.com/cell/fulltext/S0092-867 | 2(2) | -0.02133 | 0.01188 | 0.05703 | 4.60672 | 否 |
| B_AXIS | https://www.cell.com/cell/fulltext/S0092-867 | 34(34) | -0.03253 | -0.02256 | 0.23468 | 0.39296 | 否 |
| B_AXIS | https://www.nature.com/articles/s41467-023-4 | 1(1) | -0.01446 | -0.11234 | 1.18456 | 0.75051 | 是 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | https://www.biorxiv.org/content/10.1101/2022 | 2(2) | 0.00852 | -0.02137 | 0.20111 | 0.92830 | 是 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | https://www.biorxiv.org/content/10.1101/2023 | 6(6) | 0.00270 | 0.00210 | 0.02757 | -0.45802 | 否 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | https://www.biorxiv.org/content/10.1101/2024 | 2(2) | -0.02653 | -0.00302 | 0.34802 | 0.92517 | 否 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | https://www.cell.com/cell/fulltext/S0092-867 | 2(2) | 0.01148 | 0.01814 | -0.11773 | 0.51828 | 否 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | https://www.cell.com/cell/fulltext/S0092-867 | 34(34) | -0.03286 | -0.04970 | 0.27174 | -0.61024 | 否 |
| HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | https://www.nature.com/articles/s41467-023-4 | 1(1) | 0.00457 | 0.04547 | 0.21605 | 2.01930 | 是 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | https://www.biorxiv.org/content/10.1101/2022 | 2(2) | 0.03850 | 0.00161 | -0.04139 | -3.19059 | 是 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | https://www.biorxiv.org/content/10.1101/2023 | 6(6) | 0.00150 | 0.00030 | 0.09151 | 0.81356 | 否 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | https://www.biorxiv.org/content/10.1101/2024 | 2(2) | -0.00958 | 0.02485 | -0.50582 | 3.62678 | 否 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | https://www.cell.com/cell/fulltext/S0092-867 | 2(2) | 0.01499 | 0.01598 | 0.14055 | 2.53260 | 否 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | https://www.cell.com/cell/fulltext/S0092-867 | 34(34) | 0.04256 | -0.03180 | 0.17991 | -1.15177 | 否 |
| HALLMARK_INTERFERON_GAMMA_RESPONSE | https://www.nature.com/articles/s41467-023-4 | 1(1) | 0.01081 | -0.05511 | 0.11665 | 7.26696 | 是 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | https://www.biorxiv.org/content/10.1101/2022 | 2(2) | 0.00599 | 0.02245 | -0.08257 | 1.52021 | 是 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | https://www.biorxiv.org/content/10.1101/2023 | 6(6) | 0.01631 | -0.01823 | -0.00109 | -0.29867 | 否 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | https://www.biorxiv.org/content/10.1101/2024 | 2(2) | 0.00652 | -0.03511 | 0.61472 | -1.26887 | 否 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | https://www.cell.com/cell/fulltext/S0092-867 | 2(2) | 0.00497 | 0.02571 | -0.20346 | 1.51534 | 否 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | https://www.cell.com/cell/fulltext/S0092-867 | 34(34) | -0.05249 | -0.05017 | 0.09154 | 0.52266 | 否 |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | https://www.nature.com/articles/s41467-023-4 | 1(1) | 0.04944 | -0.01417 | -0.03174 | -0.16914 | 是 |

注：分来源必须在多目标、多基线下逐一核对（完整分来源汇总可再按需派生），不得把这个分层当成显著性或成功率（C-02）。

## 5. 描述性方向合取（不是阈值、不是显著性）

| 臂 | 目标 | 读出 | 方向 | 同时优于四个对照 | 最小可评单位数 |
|---|---|---|---|---|---:|
| internal | B_AXIS | rec_local_pattern_pearson | higher_is_better | False | 65 |
| internal | B_AXIS | rec_jaccard_q75 | higher_is_better | False | 65 |
| internal | B_AXIS | rec_centroid_displacement_pitch | lower_is_better | False | 65 |
| internal | B_AXIS | rec_patch75_largest_radius_abs_delta | lower_is_better | False | 65 |
| internal | B_AXIS | rec_rmse | lower_is_better | False | 65 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | rec_local_pattern_pearson | higher_is_better | True | 65 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | rec_jaccard_q75 | higher_is_better | True | 65 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | rec_centroid_displacement_pitch | lower_is_better | True | 65 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | rec_patch75_largest_radius_abs_delta | lower_is_better | True | 65 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | rec_rmse | lower_is_better | False | 65 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | rec_local_pattern_pearson | higher_is_better | True | 65 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | rec_jaccard_q75 | higher_is_better | True | 65 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | rec_centroid_displacement_pitch | lower_is_better | True | 65 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | rec_patch75_largest_radius_abs_delta | lower_is_better | True | 65 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | rec_rmse | lower_is_better | True | 65 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | rec_local_pattern_pearson | higher_is_better | False | 65 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | rec_jaccard_q75 | higher_is_better | False | 65 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | rec_centroid_displacement_pitch | lower_is_better | False | 65 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | rec_patch75_largest_radius_abs_delta | lower_is_better | False | 65 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | rec_rmse | lower_is_better | True | 65 |
| source | B_AXIS | rec_local_pattern_pearson | higher_is_better | False | 47 |
| source | B_AXIS | rec_jaccard_q75 | higher_is_better | False | 47 |
| source | B_AXIS | rec_centroid_displacement_pitch | lower_is_better | False | 47 |
| source | B_AXIS | rec_patch75_largest_radius_abs_delta | lower_is_better | False | 47 |
| source | B_AXIS | rec_rmse | lower_is_better | False | 47 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | rec_local_pattern_pearson | higher_is_better | False | 47 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | rec_jaccard_q75 | higher_is_better | False | 47 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | rec_centroid_displacement_pitch | lower_is_better | False | 47 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | rec_patch75_largest_radius_abs_delta | lower_is_better | True | 47 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | rec_rmse | lower_is_better | False | 47 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | rec_local_pattern_pearson | higher_is_better | False | 47 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | rec_jaccard_q75 | higher_is_better | False | 47 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | rec_centroid_displacement_pitch | lower_is_better | False | 47 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | rec_patch75_largest_radius_abs_delta | lower_is_better | False | 47 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | rec_rmse | lower_is_better | False | 47 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | rec_local_pattern_pearson | higher_is_better | False | 47 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | rec_jaccard_q75 | higher_is_better | False | 47 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | rec_centroid_displacement_pitch | lower_is_better | False | 47 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | rec_patch75_largest_radius_abs_delta | lower_is_better | False | 47 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | rec_rmse | lower_is_better | False | 47 |

## 6. 限制与未确定项

- 未知患者：内部 72 片中 14 片（8 单位）为研究内确认身份，其余 58 片（57 单位）未知；切片数不冒称患者数。
- 来源身份不足的 19 片只保留内部结果，不进入来源臂；来源臂 47 单位覆盖 6 个来源。
- 乳腺（source_5）与肝脏（source_0）为器官外推，未与同器官来源混合解释。
- 冻结基线复现：11/44 case 超出旧逐点 1e-4 / 每片 RMSE 1e-5 界，按 D-275 保留为已披露的失败，不改写成通过、不放宽界限。
- 数值波动证据只覆盖 source_0/EMT 的 point 模型（单 case）；不推广到 candidate/background/permuted，也不假定同批误差相消。
- 本轮 NOT_RUN：人工 TLS、独立 marker、正式显著性、物理尺度、视野外推断、shuffle 参考（本轮关闭）。
