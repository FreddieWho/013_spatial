#!/usr/bin/env bash
# Deep-stack priority subset: 100 Synapse files (88 Visium + 8 CODEX + 4 scRNA)
# Generated from portal manifest x deep3d biospecimen IDs. Status: TRIGGERED 2026-09-11
set -euo pipefail
synapse get syn53215828 --downloadLocation synapse/syn53215828  # V10X_auxiliary_BRCA23/20201103_ST_HT206B1-S1Fc1U2.tif
synapse get syn53215848 --downloadLocation synapse/syn53215848  # V10X_auxiliary_BRCA23/20201103_ST_HT206B1-S1Fc1U3.tif
synapse get syn53215842 --downloadLocation synapse/syn53215842  # V10X_auxiliary_BRCA23/20201103_ST_HT206B1-S1Fc1U4.tif
synapse get syn61465442 --downloadLocation synapse/syn61465442  # V10X_auxiliary-spatial_clone/20201103_ST_HT206B1-S1Fc1U5.tif
synapse get syn61465364 --downloadLocation synapse/syn61465364  # V10X_auxiliary-spatial_clone/A1_HT397B1-S1H3A1U1.tif
synapse get syn61465339 --downloadLocation synapse/syn61465339  # V10X_auxiliary-spatial_clone/B1_HT397B1-S1H3A1U21.tif
synapse get syn61465475 --downloadLocation synapse/syn61465475  # V10X_auxiliary-spatial_clone/BRCA_206B1.tif
synapse get syn61465393 --downloadLocation synapse/syn61465393  # V10X_auxiliary-spatial_clone/C1_HT397B1-S1H3A1.tif
synapse get syn61465369 --downloadLocation synapse/syn61465369  # V10X_auxiliary-spatial_clone/C1_HT397B1-S1H3A1U41.tif
synapse get syn61465449 --downloadLocation synapse/syn61465449  # V10X_auxiliary-spatial_clone/D1_HT397B1-S1H2A4.tif
synapse get syn61465420 --downloadLocation synapse/syn61465420  # V10X_auxiliary-spatial_clone/D1_HT397B1-S1H3A1U61.tif
synapse get syn53215135 --downloadLocation synapse/syn53215135  # V10X_level_3_BRCA23/HT206B1-S1Fc1U2Z1B1-barcodes.tsv.gz
synapse get syn53215213 --downloadLocation synapse/syn53215213  # V10X_level_3_BRCA23/HT206B1-S1Fc1U2Z1B1-features.tsv.gz
synapse get syn53215292 --downloadLocation synapse/syn53215292  # V10X_level_3_BRCA23/HT206B1-S1Fc1U2Z1B1-matrix.mtx.gz
synapse get syn53215250 --downloadLocation synapse/syn53215250  # V10X_level_3_BRCA23/HT206B1-S1Fc1U2Z1B1-scalefactors_json.json
synapse get syn53215278 --downloadLocation synapse/syn53215278  # V10X_level_3_BRCA23/HT206B1-S1Fc1U2Z1B1-tissue_lowres_image.png
synapse get syn53215161 --downloadLocation synapse/syn53215161  # V10X_level_3_BRCA23/HT206B1-S1Fc1U2Z1B1-tissue_positions_list.csv
synapse get syn53215195 --downloadLocation synapse/syn53215195  # V10X_level_3_BRCA23/HT206B1-S1Fc1U3Z1B1-barcodes.tsv.gz
synapse get syn53215117 --downloadLocation synapse/syn53215117  # V10X_level_3_BRCA23/HT206B1-S1Fc1U3Z1B1-features.tsv.gz
synapse get syn53215308 --downloadLocation synapse/syn53215308  # V10X_level_3_BRCA23/HT206B1-S1Fc1U3Z1B1-matrix.mtx.gz
synapse get syn53215097 --downloadLocation synapse/syn53215097  # V10X_level_3_BRCA23/HT206B1-S1Fc1U3Z1B1-scalefactors_json.json
synapse get syn53215276 --downloadLocation synapse/syn53215276  # V10X_level_3_BRCA23/HT206B1-S1Fc1U3Z1B1-tissue_lowres_image.png
synapse get syn53215104 --downloadLocation synapse/syn53215104  # V10X_level_3_BRCA23/HT206B1-S1Fc1U3Z1B1-tissue_positions_list.csv
synapse get syn53215269 --downloadLocation synapse/syn53215269  # V10X_level_3_BRCA23/HT206B1-S1Fc1U4Z1B1-barcodes.tsv.gz
synapse get syn53215160 --downloadLocation synapse/syn53215160  # V10X_level_3_BRCA23/HT206B1-S1Fc1U4Z1B1-features.tsv.gz
synapse get syn53215307 --downloadLocation synapse/syn53215307  # V10X_level_3_BRCA23/HT206B1-S1Fc1U4Z1B1-matrix.mtx.gz
synapse get syn53215259 --downloadLocation synapse/syn53215259  # V10X_level_3_BRCA23/HT206B1-S1Fc1U4Z1B1-scalefactors_json.json
synapse get syn53215089 --downloadLocation synapse/syn53215089  # V10X_level_3_BRCA23/HT206B1-S1Fc1U4Z1B1-tissue_lowres_image.png
synapse get syn53215176 --downloadLocation synapse/syn53215176  # V10X_level_3_BRCA23/HT206B1-S1Fc1U4Z1B1-tissue_positions_list.csv
synapse get syn61462356 --downloadLocation synapse/syn61462356  # V10X_level_3-spatial_clone/HT206B1-S1Fc1U5Z1B1-barcodes.tsv.gz
synapse get syn61462085 --downloadLocation synapse/syn61462085  # V10X_level_3-spatial_clone/HT206B1-S1Fc1U5Z1B1-features.tsv.gz
synapse get syn61462518 --downloadLocation synapse/syn61462518  # V10X_level_3-spatial_clone/HT206B1-S1Fc1U5Z1B1-matrix.mtx.gz
synapse get syn61462118 --downloadLocation synapse/syn61462118  # V10X_level_3-spatial_clone/HT206B1-S1Fc1U5Z1B1-scalefactors_json.json
synapse get syn61462309 --downloadLocation synapse/syn61462309  # V10X_level_3-spatial_clone/HT206B1-S1Fc1U5Z1B1-tissue_lowres_image.png
synapse get syn61462047 --downloadLocation synapse/syn61462047  # V10X_level_3-spatial_clone/HT206B1-S1Fc1U5Z1B1-tissue_positions_list.csv
synapse get syn61462061 --downloadLocation synapse/syn61462061  # V10X_level_3-spatial_clone/HT206B1-U1_ST_Bn1-barcodes.tsv.gz
synapse get syn61462218 --downloadLocation synapse/syn61462218  # V10X_level_3-spatial_clone/HT206B1-U1_ST_Bn1-features.tsv.gz
synapse get syn61462504 --downloadLocation synapse/syn61462504  # V10X_level_3-spatial_clone/HT206B1-U1_ST_Bn1-matrix.mtx.gz
synapse get syn61462415 --downloadLocation synapse/syn61462415  # V10X_level_3-spatial_clone/HT206B1-U1_ST_Bn1-scalefactors_json.json
synapse get syn61462430 --downloadLocation synapse/syn61462430  # V10X_level_3-spatial_clone/HT206B1-U1_ST_Bn1-tissue_lowres_image.png
synapse get syn61462237 --downloadLocation synapse/syn61462237  # V10X_level_3-spatial_clone/HT206B1-U1_ST_Bn1-tissue_positions_list.csv
synapse get syn61462145 --downloadLocation synapse/syn61462145  # V10X_level_3-spatial_clone/HT397B1-S1H2Fs4U1Bp1-barcodes.tsv.gz
synapse get syn61462397 --downloadLocation synapse/syn61462397  # V10X_level_3-spatial_clone/HT397B1-S1H2Fs4U1Bp1-features.tsv.gz
synapse get syn61462475 --downloadLocation synapse/syn61462475  # V10X_level_3-spatial_clone/HT397B1-S1H2Fs4U1Bp1-matrix.mtx.gz
synapse get syn61462301 --downloadLocation synapse/syn61462301  # V10X_level_3-spatial_clone/HT397B1-S1H2Fs4U1Bp1-scalefactors_json.json
synapse get syn61462191 --downloadLocation synapse/syn61462191  # V10X_level_3-spatial_clone/HT397B1-S1H2Fs4U1Bp1-tissue_lowres_image.png
synapse get syn61462255 --downloadLocation synapse/syn61462255  # V10X_level_3-spatial_clone/HT397B1-S1H2Fs4U1Bp1-tissue_positions_list.csv
synapse get syn61462019 --downloadLocation synapse/syn61462019  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U1Bp1-barcodes.tsv.gz
synapse get syn61462180 --downloadLocation synapse/syn61462180  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U1Bp1-features.tsv.gz
synapse get syn61462467 --downloadLocation synapse/syn61462467  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U1Bp1-matrix.mtx.gz
synapse get syn61462278 --downloadLocation synapse/syn61462278  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U1Bp1-scalefactors_json.json
synapse get syn61462216 --downloadLocation synapse/syn61462216  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U1Bp1-tissue_lowres_image.png
synapse get syn61462133 --downloadLocation synapse/syn61462133  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U1Bp1-tissue_positions.csv
synapse get syn61462105 --downloadLocation synapse/syn61462105  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U21Bp1-barcodes.tsv.gz
synapse get syn61462306 --downloadLocation synapse/syn61462306  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U21Bp1-features.tsv.gz
synapse get syn61462443 --downloadLocation synapse/syn61462443  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U21Bp1-matrix.mtx.gz
synapse get syn61462243 --downloadLocation synapse/syn61462243  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U21Bp1-scalefactors_json.json
synapse get syn61462189 --downloadLocation synapse/syn61462189  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U21Bp1-tissue_lowres_image.png
synapse get syn61462246 --downloadLocation synapse/syn61462246  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U21Bp1-tissue_positions.csv
synapse get syn61462109 --downloadLocation synapse/syn61462109  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U41Bp1-barcodes.tsv.gz
synapse get syn61462113 --downloadLocation synapse/syn61462113  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U41Bp1-features.tsv.gz
synapse get syn61462258 --downloadLocation synapse/syn61462258  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U41Bp1-matrix.mtx.gz
synapse get syn61462199 --downloadLocation synapse/syn61462199  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U41Bp1-scalefactors_json.json
synapse get syn61462355 --downloadLocation synapse/syn61462355  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U41Bp1-tissue_lowres_image.png
synapse get syn61462029 --downloadLocation synapse/syn61462029  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U41Bp1-tissue_positions.csv
synapse get syn61462150 --downloadLocation synapse/syn61462150  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U61Bp1-barcodes.tsv.gz
synapse get syn61462028 --downloadLocation synapse/syn61462028  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U61Bp1-features.tsv.gz
synapse get syn61462245 --downloadLocation synapse/syn61462245  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U61Bp1-matrix.mtx.gz
synapse get syn61462161 --downloadLocation synapse/syn61462161  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U61Bp1-scalefactors_json.json
synapse get syn61462417 --downloadLocation synapse/syn61462417  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U61Bp1-tissue_lowres_image.png
synapse get syn61462196 --downloadLocation synapse/syn61462196  # V10X_level_3-spatial_clone/HT397B1-S1H3A1U61Bp1-tissue_positions.csv
synapse get syn61462328 --downloadLocation synapse/syn61462328  # V10X_level_3-spatial_clone/HT397B1-S1H3Fs1U1Bp1-barcodes.tsv.gz
synapse get syn61462091 --downloadLocation synapse/syn61462091  # V10X_level_3-spatial_clone/HT397B1-S1H3Fs1U1Bp1-features.tsv.gz
synapse get syn61462514 --downloadLocation synapse/syn61462514  # V10X_level_3-spatial_clone/HT397B1-S1H3Fs1U1Bp1-matrix.mtx.gz
synapse get syn61462406 --downloadLocation synapse/syn61462406  # V10X_level_3-spatial_clone/HT397B1-S1H3Fs1U1Bp1-scalefactors_json.json
synapse get syn61462170 --downloadLocation synapse/syn61462170  # V10X_level_3-spatial_clone/HT397B1-S1H3Fs1U1Bp1-tissue_lowres_image.png
synapse get syn61462184 --downloadLocation synapse/syn61462184  # V10X_level_3-spatial_clone/HT397B1-S1H3Fs1U1Bp1-tissue_positions_list.csv
synapse get syn61462684 --downloadLocation synapse/syn61462684  # V10X_level_4-spatial_clone/HT206B1-S1Fc1U2Z1B1-SeuratObj.rds
synapse get syn61462766 --downloadLocation synapse/syn61462766  # V10X_level_4-spatial_clone/HT206B1-S1Fc1U3Z1B1-SeuratObj.rds
synapse get syn61462656 --downloadLocation synapse/syn61462656  # V10X_level_4-spatial_clone/HT206B1-S1Fc1U4Z1B1-SeuratObj.rds
synapse get syn61462685 --downloadLocation synapse/syn61462685  # V10X_level_4-spatial_clone/HT206B1-S1Fc1U5Z1B1-SeuratObj.rds
synapse get syn61462616 --downloadLocation synapse/syn61462616  # V10X_level_4-spatial_clone/HT206B1-U1_ST_Bn1-SeuratObj.rds
synapse get syn61462732 --downloadLocation synapse/syn61462732  # V10X_level_4-spatial_clone/HT397B1-S1H2Fs4U1Bp1-SeuratObj.rds
synapse get syn61462714 --downloadLocation synapse/syn61462714  # V10X_level_4-spatial_clone/HT397B1-S1H3A1U1Bp1-SeuratObj.rds
synapse get syn61462543 --downloadLocation synapse/syn61462543  # V10X_level_4-spatial_clone/HT397B1-S1H3A1U21Bp1-SeuratObj.rds
synapse get syn61462559 --downloadLocation synapse/syn61462559  # V10X_level_4-spatial_clone/HT397B1-S1H3A1U41Bp1-SeuratObj.rds
synapse get syn61462536 --downloadLocation synapse/syn61462536  # V10X_level_4-spatial_clone/HT397B1-S1H3A1U61Bp1-SeuratObj.rds
synapse get syn61462648 --downloadLocation synapse/syn61462648  # V10X_level_4-spatial_clone/HT397B1-S1H3Fs1U1Bp1-SeuratObj.rds
synapse get syn64720589 --downloadLocation synapse/syn64720589  # HTAN WUSTL/CODEX_Imaging_level_3_Segmentation-spatial_clone/HT397B1-S1H3A1-U12_c
synapse get syn64720591 --downloadLocation synapse/syn64720591  # HTAN WUSTL/CODEX_Imaging_level_3_Segmentation-spatial_clone/HT397B1-S1H3A1-U2_ce
synapse get syn64720588 --downloadLocation synapse/syn64720588  # HTAN WUSTL/CODEX_Imaging_level_3_Segmentation-spatial_clone/HT397B1-S1H3A1-U22_c
synapse get syn64720587 --downloadLocation synapse/syn64720587  # HTAN WUSTL/CODEX_Imaging_level_3_Segmentation-spatial_clone/HT397B1-S1H3A1-U31_c
synapse get syn64720594 --downloadLocation synapse/syn64720594  # HTAN WUSTL/CODEX_Imaging_level_4-spatial_clone/HT397B1-U12.txt
synapse get syn64720593 --downloadLocation synapse/syn64720593  # HTAN WUSTL/CODEX_Imaging_level_4-spatial_clone/HT397B1-U2.txt
synapse get syn64720590 --downloadLocation synapse/syn64720590  # HTAN WUSTL/CODEX_Imaging_level_4-spatial_clone/HT397B1-U22.txt
synapse get syn64720592 --downloadLocation synapse/syn64720592  # HTAN WUSTL/CODEX_Imaging_level_4-spatial_clone/HT397B1-U31.txt
synapse get syn61463921 --downloadLocation synapse/syn61463921  # snRNA-Seq_level_3-spatial_clone/HT397B1-S1H4A4Y1N1Z1_1Bmn1_1-barcodes.tsv.gz
synapse get syn61463882 --downloadLocation synapse/syn61463882  # snRNA-Seq_level_3-spatial_clone/HT397B1-S1H4A4Y1N1Z1_1Bmn1_1-features.tsv.gz
synapse get syn61463962 --downloadLocation synapse/syn61463962  # snRNA-Seq_level_3-spatial_clone/HT397B1-S1H4A4Y1N1Z1_1Bmn1_1-matrix.mtx.gz
synapse get syn61465425 --downloadLocation synapse/syn61465425  # snRNA-Seq_level_4-spatial_clone/snRNA_L4__HT397B1-S1H4_combo.rds
