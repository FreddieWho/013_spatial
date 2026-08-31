# R-02 structure registry

This is a metadata/schema-only, fail-closed control plane. It does not read
expression matrices or image pixels, does not use the network, and does not
use a GPU. The only h5ad access is the per-spot barcode/index coordinate
columns (`obs/_index`, `obs/array_row`, `obs/array_col`) needed to replay
annotation geometry (D-033) plus the per-spot human `ground_truth` label
column of the eight USZ Zenodo 14620362 h5ad files (D-038).

The scoped TLS inventory contains 87 source-reported TLS IDs: HTAN
Vanderbilt CRC 44, GSE226997 4, GSE274103 8, GSE274557 1 (Atlas Table S4)
and GSE175540 30 (Meylan 2022 Table S4 R-IDs). These are cross-checked
against Table S2 counts rather than expanded from count ordinals and remain
non-confirmatory inventory rows.

Auditable instance-level GT now exists for the training-lineage HTAN
Vanderbilt CRC unit and for three public-deposit validation lineages
(D-036/D-037/D-038):

- HTAN Vanderbilt CRC: Heiser et al. 2023 (PMC10756562) commit-pinned
  per-spot pathology annotations, crosswalked piece-by-piece through the
  GitHub sample key and the R-01 meta2 evidence asset, replayed to array
  coordinates through the local OSF h5ad spot indexes. 42 piece/rule GT
  sources are auditable; one piece is replay-blocked (missing local h5ad)
  and one annotation capture is unlinked (its R-01 unit is excluded).
- GSE175540 (Meylan et al. 2022, Immunity, PMID 35231421): author-deposited
  per-spot TLS annotation CSVs in the GEO raw data, replayed through the
  local Space Ranger tissue positions. 18 sample-level sources are
  auditable (35 TLS components, 18 patients); three T_agg-only files are
  registered NOT_AUDITABLE, one sample has no deposited annotation file,
  and two zero-TLS samples yield no audit rows.
- TLS_VISIUM_USZ (Zenodo 14620362): per-spot human `ground_truth` labels
  inside the eight deposited h5ad files. 8 sources are auditable (108 TLS
  components, 8 patients).
- ST_CRC_CMS (Zenodo 7760264, Valdeolivas et al. 2024): pathologist
  per-spot category CSVs for 14 sections from 7 CRC patients. 12 sources
  are auditable (551 tumor-stroma boundary components); two replicate
  sections without boundary labels yield no audit rows.

Together these give 889 confirmatory structure instances: 147 TLS (HTAN 4,
GSE175540 35, USZ 108) and 742 tumor-stroma boundary components (HTAN 191,
ST_CRC_CMS 551). Instances in preneoplastic or normal-mucosa context remain
registered as non-confirmatory context only. GSE226997, GSE274103 and
GSE274557 still have no auditable GT; vasculature and necrosis have no
public GT in any examined source.

Nineteen local STOmicsDB TLS annotation files are inventoried by path,
checksum, and compressed CSV header only. They are outside the R-01 frozen
core and are excluded from the 87 scoped summaries. No STOmics row is
accepted as confirmatory GT.

Outer splits use a conservative patient-wide envelope, including all explicit
blocks from the same patient. Explicit blocks remain separately registered for
block-level counting. R-01 patient-linked physical-specimen exceptions stay
block-empty and are not block-level evidence. The one
explicit TENX section pair is registered as an identity link, not a confirmed
cross-section structure correspondence.

The stored gate is expected to be `PARTIAL_GT_READY`, with TLS and
TUMOR_STROMA_BOUNDARY frozen as claim-bearing structures for the training
lineage and the three public-deposit validation lineages. The
`SECOND_CLAIM_BEARING_STRUCTURE_NOT_FROZEN` blocker is cleared; the
remaining validation-lineage GT gaps are tracked in `docs/ISSUES.md`.
