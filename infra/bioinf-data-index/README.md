# External bioinformatics metadata index

This directory records the approved metadata-only external inputs used for
R-01 identity resolution. `index.tsv` is the checksum inventory,
`manifest.json` records provenance and scope incidents, and `summary.json`
provides machine-readable retention totals.

Only official metadata and an official patient-to-GSM identity crosswalk are
retained here. No expression matrix, image, or biological result file is
retained. Local GEO raw archives outside this directory are used only as
header-level GSM locators; their members are not opened by the R-01 extractor.
