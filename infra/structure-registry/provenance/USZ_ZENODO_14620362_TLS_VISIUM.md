# Provenance: USZ TLS Visium dataset (Zenodo 14620362) h5ad obs ground-truth labels

- **Source record:** Zenodo 14620362, "10x Visium Spatial Transcriptomics Dataset: Kidney (3) and Lung (5) Cancer with Tertiary Lymphoid Structures", creators Dawo, Sebastian; Nonchev, Kalin; Silina, Karina; published 2025-01-09. doi:10.5281/zenodo.14620362. Record metadata snapshot: `data/other_sources/zenodo_usz_tls_visium/zenodo_record_14620362.json` (retrieved 2026-08-03).
- **Annotation assets:** `TLS_VISIUM_USZ.zip` (2,086,265,101 bytes; selectively extracted 2026-08-03, high-resolution tif images skipped) → `data/other_sources/zenodo_usz_tls_visium/TLS_VISIUM_USZ/h5ad_preprocessed/{KC1,KC2,KC3,LC1,LC2,LC3,LC4,LC5}.h5ad`, with per-sample sidecar JSON `{"SAMPLE": "KC1", …}`.
- **Content addressing:** per-file SHA-256 recorded in `gt_source_audit.tsv`.
- **R-01 lineage:** `TLS_VISIUM_USZ` (external validation role, R01-v1 freeze); 8 physical units `TLS_VISIUM_USZ::{sample}`, one sample per patient.

## Annotation modality and generation method (verbatim quotes)

From the Zenodo record description:

> "The kidney and lung cancer dataset with tertiary lymphoid structures (TLS) consists of 5 μm thick FFPE sections from kidney (3) and lung (5) tumors obtained from the Institute of Pathology at the University Hospital of Zurich, mounted onto Visium slides with the Human Probe Set v1. The samples were stained with hematoxylin and eosin (H&E) and subsequently processed for sequencing following the manufacturer's recommendations."

> "The Visium spots were annotated in the corresponding H&E images by expert researchers (K.S. and S.D.) and included manual annotations with the following labels: TLS, Immune, Tumor, Normal, and Unassigned."

> "K.N performed the manual fiducial alignment and tissue detection on the high-resolution images using Loupe Browser 8.0.0 and preprocesseed the AnnData objects."

## Audit assessment (2026-08-03)

- **Modality:** manual expert annotation of Visium spots on H&E images, stored inside the deposit h5ad `obs/ground_truth` categorical column; independent of the expression assay. Not expression clustering, not deconvolution, not image-model inference. Per D-036 the human expert annotators are registered as high-confidence at this stage.
- **Scoped reads (D-038):** only `obs/_index`, `obs/x_array`, `obs/y_array` and `obs/ground_truth` (categories + codes) are read from the h5ad files; no expression values (`X`, `layers`, `obsm`) and no image pixels are read.
- **Coordinate system:** Visium spot barcode replayed to `x_array`/`y_array` array coordinates from the same files (all spots in-tissue, `in_tissue == 1`).
- **Label vocabulary (observed, fail-closed):** {`TLS`, `INFL`, `TUM`, `NOR`, `UNASSIGNED`, `LN`}; `LN` (33 spots, LC3 only) extends the five classes named in the record description and is outside the TLS scope. `UNASSIGNED` cells are unknown, never negative.
- **Result:** 8 auditable GT sources (one per sample), 108 TLS connected-component instances across 8 patients (KC1–3 kidney, LC1–5 lung).
