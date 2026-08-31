# Provenance: Heiser et al. 2023 (Cell) HTAN Vanderbilt CRC spot-level pathology annotations

- **Source publication:** Heiser CN, Simmons AJ, Revetta F, et al. "Molecular cartography uncovers evolutionary and microenvironmental dynamics in sporadic colorectal tumors." *Cell* 2023 Dec 7;186(25):5620–5637.e16. doi:10.1016/j.cell.2023.11.006. PMID 38065082; PMCID PMC10756562 (author manuscript, NIHMS1944717).
- **Annotation assets:** GitHub repository `Ken-Lau-Lab/spatial_CRC_atlas`, `resources/ST/{capture}_pathology_annotation.csv` (42 files, two columns `Barcode,pathology_annotation`), `resources/ST/visium_sample_key.csv`, `resources/WES/*_ROIs_*.csv`.
- **Content addressing:** repository cloned 2026-07-31 at commit `64e585453514801a3b730f86aec90bbc0f0595da` into `data/other_sources/htan/spatial_CRC_atlas_repo/`; per-file SHA-256 recorded in `gt_source_audit.tsv`. Full-text XML snapshot: `infra/structure-registry/provenance/PMC10756562_Heiser2023_Cell_fulltext.xml` (SHA-256 `410bb1e34fe08a9db889f3103c8dd41f3ac712c7218b92ad9682c8ae811ffb0d`, NCBI efetch db=pmc id=10756562, retrieved 2026-07-31).
- **R-01 lineage:** `HTAN_VANDERBILT_CRC` (training role, R01-v1 freeze).

## Annotation modality and generation method (verbatim quotes)

From the STAR Methods and Results text of PMC10756562:

> "Regions of interest (ROIs) for ST were chosen based on histological annotation of FFPE blocks, targeting tumor areas with morphology indicative of various stages of malignancy and transition points between them."

> "Samples with concurrent pre-malignant, malignant, and invasive regions were identified by a pathologist (STAR Methods: Experimental model and study participant details)."

> "In ST, we grouped Visium microwells by patient and used stromal regions and adjacent normal epithelium, manually annotated in the 10X Genomics Loupe Browser, to provide a normal background for calling CNVs in tumor regions."

> "Registration of LCM-WES with ST was performed manually by creating masks for each LCM ROI in the 10X Genomics Loupe Browser that were used to subset Visium microwells to LCM ROIs for downstream analysis."

> "Loupe Browser version 6.4.0, 10X Genomics" (Key Resources Table).

> "Data and code availability — All data have been deposited to the HTAN Data Coordinating Center Data Portal at the National Cancer Institute: https://data.humantumoratlas.org/ (under the HTAN Vanderbilt Atlas). Code used to perform analyses in this manuscript is available at https://github.com/Ken-Lau-Lab/spatial_CRC_atlas (for WES, scRNA, and ST data) and https://github.com/Ken-Lau-Lab/spatial_CRC_atlas_imaging (for MxIF and spatial registration). Annotated, pre-processed data compatible with the codebase can be downloaded from OSF Storage: https://osf.io/hftq2/."

Author contributions list data curation by, among others, F.R. (Frank Revetta, pathology informatics) and M.K.W. (M. Kay Washington, GI pathologist, Dept. of Pathology, Microbiology and Immunology, Vanderbilt University Medical Center).

## Audit assessment (2026-07-31)

- **Modality:** manual region annotation on H&E morphology in 10X Loupe Browser at Visium spot (microwell) resolution; independent of the expression assay. Not expression clustering, not deconvolution, not image-model inference.
- **Coordinate system:** Visium spot barcode (16 nt + `-1`); replayable to array coordinates via the OSF h5ad `obs/_index, array_row, array_col` (identity/coordinate reads only, per D-033).
- **Label vocabulary (observed):** `carcinoma`, `carcinoma_border`, `carcinoma_edge`, `adenoma`, `adenoma_border`, `normal_mucosa`, `smooth_muscle`, `lymphoid_follicle`; empty cells = unannotated (unknown, never negative).
- **Crosswalk:** capture area → h5ad → block → patient via `visium_sample_key.csv` (41 keyed capture areas), identical on the keyed subset to the verified R-01 evidence asset `repo/data_meta/ST_CRC_cohort_meta2.csv` (48 rows incl. 12 `duplicated=Yes`).
- **Residual provenance gap (honest limitation):** the per-file annotator identity (which person drew which CSV) is not stated file-by-file in the paper or repository; the method text and team composition (GI pathologist + pathology informatics co-authors, "identified by a pathologist", "manually annotated in the 10X Genomics Loupe Browser") support manual pathology annotation as the modality, but individual annotator attribution is not auditable. This is registered as boundary/provenance uncertainty, not as a fail condition.
- **Exclusions observed at registration:** `7794_4` annotation CSV has no sample-key or meta2 row (physical link unverifiable, fail-closed); `7003_7` is keyed but its OSF h5ad is not in the local inventory (barcode/coordinate replay impossible without one additional OSF fetch; deferred).
