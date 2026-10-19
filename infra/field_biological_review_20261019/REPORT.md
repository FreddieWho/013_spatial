# 实测场资格与冻结预测恢复：描述性清单

只描述，不判定成功/失败，不设生物学通过线，不给出 p 值。
长度单位一律为该切片经验点间距（pitch）；分位阈值是操作阈值，不是生物真值。
三本账分开：测量资格（coverage/支持）、输入回声关系（与可见模块场）、
已开发预测的探索性再分析（实测 vs 预测）。三者的列名与 evidence_role 不同，不共用“成功”一词。

## 0. 语料与身份

- 主清单：来源留出 `infra/field_source_holdout_20261018` 的 1000 份冻结预测包（内部 72 片/65 单位、来源 53 片/47 单位/6 来源）；`train_indices.npz` 不进入评分。
- 实测清单覆盖全部 72 片 × 4 目标；来源身份不足的 19 片保留在内部描述，不进入来源分母。
- 研究内确认身份 14 片 / 8 个单位；其余为未知单位，不升级为独立患者。
- 状态计数（实测）：{"DESCRIBED": 2, "DESCRIBED_WITH_NOT_ESTIMABLE_METRICS": 286}
- 状态计数（预测）：{"DESCRIBED": 6, "DESCRIBED_WITH_NOT_ESTIMABLE_METRICS": 3994, "NOT_ESTIMABLE_CONSTANT_PREDICTION": 1000}

## 1. 分子检测覆盖（连续报告，不分档删除）

- B_AXIS：72 片；检出基因数中位 47；任意目标基因检出点占比 中位 0.999（0.612–1.000）；有效基因数中位 2.87（最小 1.19）；首位基因占比中位 0.535（最大 0.918）
- HALLMARK_OXIDATIVE_PHOSPHORYLATION：72 片；检出基因数中位 191；任意目标基因检出点占比 中位 1.000（0.919–1.000）；有效基因数中位 75.21（最小 4.20）；首位基因占比中位 0.055（最大 0.482）
- HALLMARK_INTERFERON_GAMMA_RESPONSE：72 片；检出基因数中位 195；任意目标基因检出点占比 中位 1.000（0.950–1.000）；有效基因数中位 24.66（最小 10.83）；首位基因占比中位 0.120（最大 0.269）
- HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION：72 片；检出基因数中位 199；任意目标基因检出点占比 中位 1.000（0.962–1.000）；有效基因数中位 32.79（最小 7.66）；首位基因占比中位 0.100（最大 0.299）

完整 72 片 × 4 目标逐片数值见 `measured_fields.tsv`；分层中位/极值见 `measured_summary.tsv`。

## 2. 恢复描述（单位宏平均，整臂；不跨信息预算判胜负）

| 臂 | 目标 | 视图 | 模型 | 模型类别 | 单位 | RMSE | 去均值RMSE | 逐点相关 | 局部图样相关 | q75 Jaccard | q75 打乱中位 | 质心位移(pitch) | 自相关差 |
|---|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| internal | B_AXIS | full15 | fixed_local_hgb | module_local_hgb | 65 | 0.1320 | 0.1009 | 0.5107 | 0.2597 | 0.3353 | 0.1427 | 1.7350 | 0.0388 |
| internal | B_AXIS | full15 | fixed_point_hgb | module_point_hgb | 65 | 0.1362 | 0.1026 | 0.5028 | 0.2583 | 0.3335 | 0.1427 | 1.5505 | -0.0170 |
| internal | B_AXIS | full15 | gene_ridge | raw_gene_ridge | 65 | 0.1039 | 0.0882 | 0.6306 | 0.3449 | 0.4223 | 0.1429 | 1.1726 | -0.0218 |
| internal | B_AXIS | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 65 | 0.0976 | 0.0805 | 0.6805 | 0.3506 | 0.4513 | 0.1427 | 1.1505 | 0.2170 |
| internal | B_AXIS | full15 | train_mean | constant_prior | 65 | 0.1692 | 0.1166 | NA | NA | NA | NA | NA | NA |
| internal | B_AXIS | small5 | fixed_local_hgb | module_local_hgb | 65 | 0.1512 | 0.1100 | 0.4100 | 0.2048 | 0.2821 | 0.1427 | 1.9823 | 0.0392 |
| internal | B_AXIS | small5 | fixed_point_hgb | module_point_hgb | 65 | 0.1498 | 0.1103 | 0.4210 | 0.2227 | 0.2880 | 0.1430 | 1.8772 | -0.0159 |
| internal | B_AXIS | small5 | gene_ridge | raw_gene_ridge | 65 | 0.1190 | 0.0981 | 0.5531 | 0.2753 | 0.3726 | 0.1433 | 1.3561 | -0.0347 |
| internal | B_AXIS | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 65 | 0.1111 | 0.0879 | 0.6174 | 0.2867 | 0.4081 | 0.1429 | 1.2903 | 0.1955 |
| internal | B_AXIS | small5 | train_mean | constant_prior | 65 | 0.1692 | 0.1166 | NA | NA | NA | NA | NA | NA |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | fixed_local_hgb | module_local_hgb | 65 | 0.1441 | 0.1090 | 0.6885 | 0.4489 | 0.4381 | 0.1431 | 1.3866 | -0.0018 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | fixed_point_hgb | module_point_hgb | 65 | 0.1495 | 0.1131 | 0.6828 | 0.4541 | 0.4321 | 0.1424 | 1.3911 | -0.0409 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | gene_ridge | raw_gene_ridge | 65 | 0.0805 | 0.0717 | 0.8499 | 0.6252 | 0.6084 | 0.1427 | 0.6522 | -0.0265 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 65 | 0.0783 | 0.0693 | 0.8669 | 0.6117 | 0.6190 | 0.1424 | 0.8419 | 0.2189 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | train_mean | constant_prior | 65 | 0.2340 | 0.1529 | NA | NA | NA | NA | NA | NA |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | fixed_local_hgb | module_local_hgb | 65 | 0.1615 | 0.1170 | 0.6633 | 0.4214 | 0.4129 | 0.1426 | 1.5876 | 0.0021 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | fixed_point_hgb | module_point_hgb | 65 | 0.1632 | 0.1210 | 0.6492 | 0.4321 | 0.3985 | 0.1432 | 1.5315 | -0.0429 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | gene_ridge | raw_gene_ridge | 65 | 0.1083 | 0.0901 | 0.7800 | 0.5176 | 0.5283 | 0.1427 | 0.8626 | -0.0406 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 65 | 0.0996 | 0.0801 | 0.8214 | 0.5181 | 0.5609 | 0.1430 | 1.0278 | 0.1960 |
| internal | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | train_mean | constant_prior | 65 | 0.2340 | 0.1529 | NA | NA | NA | NA | NA | NA |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | fixed_local_hgb | module_local_hgb | 65 | 0.0731 | 0.0561 | 0.7636 | 0.5063 | 0.4575 | 0.1426 | 0.7720 | 0.0411 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | fixed_point_hgb | module_point_hgb | 65 | 0.0784 | 0.0585 | 0.7516 | 0.4858 | 0.4462 | 0.1428 | 0.7767 | -0.0125 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | gene_ridge | raw_gene_ridge | 65 | 0.0648 | 0.0530 | 0.7972 | 0.4904 | 0.4976 | 0.1423 | 0.6564 | -0.0141 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 65 | 0.0622 | 0.0498 | 0.8204 | 0.4857 | 0.5236 | 0.1423 | 0.6966 | 0.2392 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | train_mean | constant_prior | 65 | 0.1452 | 0.0915 | NA | NA | NA | NA | NA | NA |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | fixed_local_hgb | module_local_hgb | 65 | 0.0918 | 0.0640 | 0.6775 | 0.4735 | 0.3879 | 0.1428 | 1.0281 | 0.0372 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | fixed_point_hgb | module_point_hgb | 65 | 0.0925 | 0.0642 | 0.6860 | 0.4586 | 0.3873 | 0.1424 | 0.9904 | -0.0087 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | gene_ridge | raw_gene_ridge | 65 | 0.0751 | 0.0604 | 0.7441 | 0.4290 | 0.4425 | 0.1427 | 0.7308 | -0.0285 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 65 | 0.0706 | 0.0548 | 0.7814 | 0.4312 | 0.4808 | 0.1430 | 0.7798 | 0.2148 |
| internal | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | train_mean | constant_prior | 65 | 0.1452 | 0.0915 | NA | NA | NA | NA | NA | NA |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | fixed_local_hgb | module_local_hgb | 65 | 0.0790 | 0.0563 | 0.9072 | 0.6755 | 0.5792 | 0.1429 | 0.4803 | 0.0446 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | fixed_point_hgb | module_point_hgb | 65 | 0.0789 | 0.0572 | 0.9052 | 0.6590 | 0.5818 | 0.1431 | 0.4629 | 0.0091 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | gene_ridge | raw_gene_ridge | 65 | 0.0601 | 0.0495 | 0.9258 | 0.6773 | 0.6469 | 0.1429 | 0.3147 | 0.0116 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 65 | 0.0623 | 0.0519 | 0.9206 | 0.6493 | 0.6350 | 0.1423 | 0.5197 | 0.2884 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | train_mean | constant_prior | 65 | 0.2227 | 0.1406 | NA | NA | NA | NA | NA | NA |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | fixed_local_hgb | module_local_hgb | 65 | 0.0854 | 0.0621 | 0.8881 | 0.6263 | 0.5531 | 0.1428 | 0.5130 | 0.0492 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | fixed_point_hgb | module_point_hgb | 65 | 0.0877 | 0.0643 | 0.8806 | 0.5975 | 0.5521 | 0.1428 | 0.4756 | -0.0038 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | gene_ridge | raw_gene_ridge | 65 | 0.0674 | 0.0554 | 0.9085 | 0.6265 | 0.6082 | 0.1426 | 0.3435 | -0.0003 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 65 | 0.0677 | 0.0556 | 0.9093 | 0.6087 | 0.6116 | 0.1424 | 0.5421 | 0.2710 |
| internal | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | train_mean | constant_prior | 65 | 0.2227 | 0.1406 | NA | NA | NA | NA | NA | NA |
| source | B_AXIS | full15 | fixed_local_hgb | module_local_hgb | 47 | 0.1347 | 0.0965 | 0.4125 | 0.2073 | 0.2939 | 0.1428 | 1.8746 | 0.0349 |
| source | B_AXIS | full15 | fixed_point_hgb | module_point_hgb | 47 | 0.1375 | 0.0968 | 0.4115 | 0.2043 | 0.2917 | 0.1427 | 1.7848 | -0.0166 |
| source | B_AXIS | full15 | gene_ridge | raw_gene_ridge | 47 | 0.1135 | 0.0867 | 0.5402 | 0.2550 | 0.3645 | 0.1427 | 1.3126 | -0.0280 |
| source | B_AXIS | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.1049 | 0.0760 | 0.6128 | 0.2670 | 0.4044 | 0.1427 | 1.1945 | 0.1709 |
| source | B_AXIS | full15 | train_mean | constant_prior | 47 | 0.1611 | 0.1013 | NA | NA | NA | NA | NA | NA |
| source | B_AXIS | small5 | fixed_local_hgb | module_local_hgb | 47 | 0.1435 | 0.1000 | 0.3337 | 0.1470 | 0.2549 | 0.1427 | 2.0591 | 0.0267 |
| source | B_AXIS | small5 | fixed_point_hgb | module_point_hgb | 47 | 0.1432 | 0.0996 | 0.3632 | 0.1845 | 0.2633 | 0.1432 | 1.9598 | -0.0134 |
| source | B_AXIS | small5 | gene_ridge | raw_gene_ridge | 47 | 0.1212 | 0.0909 | 0.4848 | 0.2131 | 0.3350 | 0.1433 | 1.3756 | -0.0321 |
| source | B_AXIS | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.1138 | 0.0808 | 0.5622 | 0.2268 | 0.3744 | 0.1425 | 1.3191 | 0.1626 |
| source | B_AXIS | small5 | train_mean | constant_prior | 47 | 0.1611 | 0.1013 | NA | NA | NA | NA | NA | NA |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | fixed_local_hgb | module_local_hgb | 47 | 0.1661 | 0.1119 | 0.6672 | 0.4087 | 0.4307 | 0.1429 | 1.6178 | -0.0025 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | fixed_point_hgb | module_point_hgb | 47 | 0.1455 | 0.1095 | 0.6538 | 0.4026 | 0.4153 | 0.1425 | 1.5164 | -0.0478 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | gene_ridge | raw_gene_ridge | 47 | 0.0876 | 0.0692 | 0.8304 | 0.5766 | 0.5975 | 0.1428 | 0.7924 | -0.0303 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.0835 | 0.0634 | 0.8604 | 0.5755 | 0.6222 | 0.1424 | 0.8994 | 0.1865 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | full15 | train_mean | constant_prior | 47 | 0.2329 | 0.1378 | NA | NA | NA | NA | NA | NA |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | fixed_local_hgb | module_local_hgb | 47 | 0.1905 | 0.1113 | 0.6528 | 0.3686 | 0.4191 | 0.1429 | 1.5962 | 0.0098 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | fixed_point_hgb | module_point_hgb | 47 | 0.1708 | 0.1143 | 0.6221 | 0.3846 | 0.3876 | 0.1431 | 1.5929 | -0.0520 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | gene_ridge | raw_gene_ridge | 47 | 0.1145 | 0.0823 | 0.7693 | 0.4948 | 0.5310 | 0.1428 | 0.9205 | -0.0408 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.1079 | 0.0723 | 0.8178 | 0.5033 | 0.5748 | 0.1431 | 1.0521 | 0.1684 |
| source | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION | small5 | train_mean | constant_prior | 47 | 0.2329 | 0.1378 | NA | NA | NA | NA | NA | NA |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | fixed_local_hgb | module_local_hgb | 47 | 0.0864 | 0.0606 | 0.6223 | 0.4546 | 0.3777 | 0.1427 | 1.1186 | 0.0354 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | fixed_point_hgb | module_point_hgb | 47 | 0.0884 | 0.0584 | 0.6665 | 0.3772 | 0.4063 | 0.1428 | 0.9369 | -0.0173 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | gene_ridge | raw_gene_ridge | 47 | 0.0703 | 0.0526 | 0.7354 | 0.3855 | 0.4575 | 0.1428 | 0.7741 | -0.0192 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.0664 | 0.0475 | 0.7876 | 0.3931 | 0.5071 | 0.1425 | 0.8217 | 0.1933 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | full15 | train_mean | constant_prior | 47 | 0.1633 | 0.0792 | NA | NA | NA | NA | NA | NA |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | fixed_local_hgb | module_local_hgb | 47 | 0.1175 | 0.0631 | 0.6043 | 0.4270 | 0.3637 | 0.1424 | 1.2572 | 0.0377 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | fixed_point_hgb | module_point_hgb | 47 | 0.1145 | 0.0623 | 0.6366 | 0.3591 | 0.3771 | 0.1426 | 1.0148 | -0.0132 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | gene_ridge | raw_gene_ridge | 47 | 0.0857 | 0.0595 | 0.6694 | 0.3355 | 0.3974 | 0.1424 | 0.8552 | -0.0290 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.0803 | 0.0519 | 0.7314 | 0.3458 | 0.4451 | 0.1425 | 0.9343 | 0.1748 |
| source | HALLMARK_INTERFERON_GAMMA_RESPONSE | small5 | train_mean | constant_prior | 47 | 0.1633 | 0.0792 | NA | NA | NA | NA | NA | NA |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | fixed_local_hgb | module_local_hgb | 47 | 0.0830 | 0.0653 | 0.8796 | 0.6352 | 0.5473 | 0.1426 | 0.6520 | 0.0402 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | fixed_point_hgb | module_point_hgb | 47 | 0.0826 | 0.0640 | 0.8839 | 0.6104 | 0.5694 | 0.1432 | 0.5948 | 0.0066 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | gene_ridge | raw_gene_ridge | 47 | 0.0692 | 0.0543 | 0.9107 | 0.6279 | 0.6298 | 0.1426 | 0.4308 | 0.0032 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.0696 | 0.0549 | 0.9122 | 0.6144 | 0.6279 | 0.1420 | 0.6333 | 0.2512 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | full15 | train_mean | constant_prior | 47 | 0.2065 | 0.1400 | NA | NA | NA | NA | NA | NA |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | fixed_local_hgb | module_local_hgb | 47 | 0.0915 | 0.0673 | 0.8705 | 0.6068 | 0.5543 | 0.1431 | 0.6452 | 0.0418 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | fixed_point_hgb | module_point_hgb | 47 | 0.0933 | 0.0683 | 0.8675 | 0.5771 | 0.5665 | 0.1428 | 0.5837 | -0.0014 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | gene_ridge | raw_gene_ridge | 47 | 0.0717 | 0.0576 | 0.9003 | 0.5990 | 0.6087 | 0.1427 | 0.4152 | -0.0021 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | gene_ridge_smooth_half | raw_gene_ridge_posthoc_neighbourhood_smoothing | 47 | 0.0713 | 0.0570 | 0.9048 | 0.5899 | 0.6161 | 0.1426 | 0.6306 | 0.2418 |
| source | HALLMARK_OXIDATIVE_PHOSPHORYLATION | small5 | train_mean | constant_prior | 47 | 0.2065 | 0.1400 | NA | NA | NA | NA | NA | NA |

打乱中位是保持取值边际、打乱空间排列的描述性扰动参考；不保持空间支撑与密度，未校准，
不构成零模型、不产生 p 值、不用于筛选关系。
`gene_ridge_smooth_half` 的自相关上升是构造结果（后处理平滑），不得作为空间关系或恢复证据；
本阶段不做同点/邻域基线对照，也不在信息预算之间判胜负。

## 3. 边界与未运行

- 分位阈值（q75/q80/q90）是操作阈值，取严格大于分位值；逐场报告 area_frac 与 at_threshold_frac，
  不同场的点集大小可以不同。
- 预测场是可见模块/基因的函数；预测与可见模块的相关部分由输入结构决定，不能读作独立生物学证据。
- 相关、自相关与尺度在任一方方差低于合同 sd_floor 时为 NOT_ESTIMABLE，不填 0、不填 1；RMSE 仍保留。
- 数学上可算不等于有足够生物学参照：检出覆盖、单基因主导、有效基因数、边界/孤立点支撑与未确认的物理尺度
  仍会限制解释；这些证据逐项并列，不合成单一分数，也不用统一百分比删除切片/目标/行。
- 不得把"可计算"格子升级为已具生物学评价资格；NOT_TESTABLE/NOT_ESTIMABLE 也不读作场不存在或恢复失败。
- 生物学可用性由主代理结合上述连续证据判断，本清单只提供描述性事实与不可评价状态。
- 未训练/未重训、未做独立 marker 隔离、未读取人工 TLS、未做正式显著性、未做物理尺度换算。
- 已知患者仅 8 个研究内单位；切片数不冒称患者数。

