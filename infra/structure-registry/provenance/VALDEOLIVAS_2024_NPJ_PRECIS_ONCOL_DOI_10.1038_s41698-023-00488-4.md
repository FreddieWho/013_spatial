# Provenance: Valdeolivas et al. 2024 (npj Precision Oncology) ST_CRC_CMS pathologist spot categorization

- **Source publication:** Valdeolivas A, Amberg B, Giroud N, et al. "Profiling the heterogeneity of colorectal cancer consensus molecular subtypes using spatial transcriptomics." *npj Precision Oncology* 2024;8:10. doi:10.1038/s41698-023-00488-4.
- **Source record:** Zenodo 7760264, "Profiling the Heterogeneity of Colorectal Cancer Consensus Molecular Subtypes using Spatial Transcriptomics: datasets", creator Valdeolivas, Alberto; published 2023-01-19. doi:10.5281/zenodo.7760264. Record metadata snapshot: `data/other_sources/zenodo_st_crc_cms/zenodo_record_7760264.json`; companion GitHub README snapshot: `data/other_sources/zenodo_st_crc_cms/ST_CRC_CMS_README.md` (both retrieved 2026-08-03).
- **Annotation assets:** `Pathology_SpotAnnotations.zip` → 14 `Pathologist_Annotations_{sample}.csv` files in `data/other_sources/zenodo_st_crc_cms/Pathology_SpotAnnotations/`; per-sample Space Ranger `tissue_positions_list.csv` members extracted verbatim from the 14 deposit sample zips into `data/other_sources/zenodo_st_crc_cms/tissue_positions/`.
- **Content addressing:** per-file SHA-256 recorded in `gt_source_audit.tsv`.
- **R-01 lineage:** `ST_CRC_CMS` (internal validation role, R01-v1 freeze); 14 physical units `ST_CRC_CMS::{sample}` across 7 patients, two serial sections per patient.

## Annotation modality and generation method (verbatim quotes)

From the Zenodo record description:

> "This contents the raw Spatial Transcriptomics data, spot categorization made by pathologist, the results of the deconvolution and intermediary files required to run the analysis described in our manuscript"

From the companion GitHub README:

> "We processed fresh-frozen resection samples obtained from seven CRC patients for ST using 10x Genomics VISIUM aiming at exploring spatial molecular heterogeneity in CRC. We considered two serial sections per patient to generate technical replicates."

## Audit assessment (2026-08-03)

- **Modality:** pathologist per-spot categorization on H&E (QuPath/Loupe tooling per the publication methods), deposited as per-barcode category CSVs; independent of the expression assay. Not expression clustering, not deconvolution, not image-model inference. Per D-036 the human pathologist annotators are registered as high-confidence at this stage; per-file annotator identity is attributed only for A938797 (`Pathologist_KH` header), registered as provenance uncertainty, not a fail condition.
- **Coordinate system:** Visium spot barcode (16 nt + `-1`); replayable to array coordinates via the deposit `tissue_positions_list.csv` files (in-tissue filter applied). Annotated spot sets equal the in-tissue spot sets exactly for every sample (verified at registration).
- **Label vocabulary (observed, fail-closed):** 39 verbatim deposit labels (misspelling variants included) frozen in `scripts/r02_validation_gt.py` (`STCRC_LABEL_VOCABULARY`); empty cells = unknown, never negative.
- **Scoped structure:** only the explicit `tumor&stroma` label family (`tumor&stroma`, `tumor&stroma IC med to high`, `tumor&stroma_IC low`, `tumor&stroma_IC med to high`) is registered as TUMOR_STROMA_BOUNDARY. The `IC aggregate*` labels are immune-cell aggregates, not verified TLS, and are fail-closed excluded from the TLS scope; every other label is frozen deposit vocabulary but unscoped.
- **Scoped exclusions:** both A798015 serial sections (`Rep1`, `Rep2`) contain zero in-tissue `tumor&stroma`-family spots and yield no audit rows (audited negatives, Heiser-lineage convention).
- **Result:** 12 auditable GT sources, 551 tumor-stroma boundary connected-component instances across 6 patients (12 serial sections; each patient's two serial sections form one explicit tumor block).
