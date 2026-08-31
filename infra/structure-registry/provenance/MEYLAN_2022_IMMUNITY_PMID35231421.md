# Provenance: Meylan et al. 2022 (Immunity) GSE175540 ccRCC per-spot TLS annotations

- **Source publication:** Meylan F, et al. "Tertiary lymphoid structures generate and propagate anti-tumor antibody-producing plasma cells in renal cell cancer." *Immunity* 2022;55(3):527–541.e5. doi:10.1016/j.immuni.2022.02.001. PMID 35231421.
- **Annotation assets:** GEO series GSE175540 raw data (`GSE175540_RAW.tar`, extracted to `data/GEO/GSE175540/raw/`): 23 `{GSM}_{ffpe|frozen}_{alias}_TLS_annotation.csv.gz` per-spot annotation files plus the matching Space Ranger `{GSM}_{alias}_tissue_positions_list.csv.gz` coordinate files for all 24 Visium samples (GSM5924030–GSM5924053, 12 FFPE + 12 fresh-frozen, BioProject PRJNA732692).
- **Content addressing:** series metadata snapshot `data/GEO/GSE175540/GSE175540_family.soft.gz` (retrieved 2026-08-03); per-file SHA-256 recorded in `gt_source_audit.tsv`.
- **R-01 lineage:** `GEO::GSE175540` (external validation role, R01-v1 freeze); 24 physical units `external_geo::GSM59…`, one sample per patient (`GSE175540::{alias}`).

## Annotation modality and generation method (evidence)

From the local GEO series metadata snapshot:

> "!Series_title = Tertiary lymphoid structures generate and propagate anti-tumor antibody-producing plasma cells in renal cell cancer"

> "!Series_overall_design = Spatial transcriptomics of fresh frozen ccRCC human tumors"

The publication's methods (D-037 literature census) identify TLS on CD3/CD20 immunostaining and the authors deposited the resulting per-spot TLS calls as the raw-data annotation CSVs audited here. Annotators are human domain experts of the study; per D-036 any human annotator is registered as high-confidence at this stage.

## Audit assessment (2026-08-03)

- **Modality:** human expert TLS annotation guided by CD3/CD20 immunostaining, deposited as per-barcode labels; independent of the expression assay. Not expression clustering, not deconvolution, not image-model inference.
- **Coordinate system:** Visium spot barcode (16 nt + `-1`); replayable to array coordinates via the same deposit's `tissue_positions_list.csv.gz` (in-tissue filter applied).
- **Label vocabulary (observed, fail-closed):** header `Barcode,TLS_2_cat` with values {`TLS`, `NO_TLS`, empty} (20 files); header `Barcode,TLS` with values {`T_agg`, empty} (3 files). Empty cells = unknown, never negative.
- **Scoped exclusions:** the three `T_agg`-only files (GSM5924042 `a_1`, GSM5924045 `a_17`, GSM5924047 `b_7`) carry T-cell aggregate labels; T_agg is not TLS and those rows are registered `NOT_AUDITABLE_GT`. GSM5924034 (`ffpe_c_10`) has no deposited annotation file. Two files (`ffpe_c_21`, `frozen_b_13`) contain zero in-tissue TLS spots and yield no audit rows (audited negatives, Heiser-lineage convention).
- **Result:** 18 auditable sample-level GT sources, 35 TLS connected-component instances across 18 patients.
- **Atlas count reconciliation:** the pan-cancer atlas Table S4 registers 30 KIRC TLS IDs (R_P1–R_P24 sample axis); our replayed component count is 35. The difference is definition granularity (atlas histological structure counts vs. hex-connected spot components); the registry uses the recomputed, replayable counts, and the GSM↔R_P per-sample crosswalk is not published, so the duplicate link between the two deposits is study-level only.
