# R-01 sample registry

This directory is the versioned control surface for roadmap node R-01. It
answers who a physical sample is, where each asset came from, which assets may
share a leakage group, and which analysis role a logical data unit may carry.
It does not decide whether a structure annotation is valid ground truth
(R-02), perform registration (R-03), inspect expression matrices, or select a
model.

## Evidence boundary

Only metadata bundled with the local data, plus user-approved official
metadata/crosswalks indexed under `infra/bioinf-data-index/`, may support
identity facts. Structure facts remain limited to data-carried metadata and
are audited separately in R-02. Every normalized patient, block, specimen,
section, serial relationship, or
structure field must retain its raw value, metadata location and evidence
grade:

- `E3_explicit`: bundled metadata explicitly states the relationship.
- `E2_corroborated`: at least two consistent bundled metadata fields support it.
- `E1_weak`: filename, directory, isolated alias or incomplete mapping only.
- `E0_unknown`: no auditable evidence.
- `EC_conflict`: bundled metadata sources conflict.

E1/E0/EC records cannot enter claim-bearing validation. Different checksums do
not prove physical independence. A possible cross-source match remains in one
conservative leakage group or outside confirmation until resolved.

R-01 also permits one fail-closed `patient_linked_physical_specimen` identity
granularity at E2: official patient mapping, explicit physical-tissue
description, and a stable specimen locator must corroborate one another.
`block_id` remains empty. GSM, BioSample, sample, slide, capture area, section,
directory and filename values cannot be promoted to block. Discovery of a
multi-block specimen, repeated block sections, mapping conflict or cross-study
reuse revokes the equivalence and requires a gate rebuild.

The gate does not trust either explicit patient/block rows or specimen
equivalences by themselves. It checks patient/block or
patient/specimen/basis evidence values, raw keys and values, active checksummed
source assets, actual small-file SHA-256 values, canonical source-lineage
duplicate-group membership, and role-to-lineage leakage-group agreement.
Missing or mismatched references remove the row from eligibility and are listed under
`identity_integrity_errors`, `source_asset_errors`, or `role_freeze_errors`.

## Registry tables

The frozen registry will use five normalized TSV tables:

1. `source_assets.tsv`: source, accession, local path, processing lineage,
   bytes and optional checksum.
2. `physical_units.tsv`: study, patient, block, section and serial-section
   fields, each linked to identity evidence.
3. `identity_evidence.tsv`: raw metadata key/value, source location, evidence
   grade, interpretation and conflict status.
4. `duplicate_groups.tsv`: exact asset, reprocessed section, serial block,
   same-patient and possible-match relationships.
5. `role_freeze.tsv`: 6–10 logical units, conservative leakage group, one
   primary role, allowed/forbidden uses and freeze version.

The metadata inventory generated here is only a source-level preflight. It
flags field names that may encode identity, section geometry, or circular GT;
it does not accept those fields as truth.

## Leakage and contamination model

| Risk | Why it contaminates evidence | R-01 disposition |
| --- | --- | --- |
| One physical sample appears in HEST and GEO/10x/HTAN or another aggregator | Database names differ but the patient, block or section is not independent | Place all confirmed or possible matches in one leakage group until resolved |
| Serial sections or reprocessed versions cross outer folds | Nearby sections share patient, block and often the same physical structure | Split at the conservative patient/block/serial/duplicate group |
| A metadata label, TLS score, segmentation or embedding was derived from an allowed model input | The target contains information from the predictor, creating circular validation | Record provenance only; R-02 must reject, mask or downgrade it |
| A discovery dataset is later renamed external validation | Analysis history already influenced representation, threshold or model choice | Freeze one primary role and append any later role change |
| Similar filenames or aliases are merged without sufficient bundled evidence | False merging reduces effective sample size and can invent leakage groups | Preserve raw values; use E1/E0/EC and quarantine rather than force a match |
| Different checksums are treated as proof of independence | Reformatting, cropping and reprocessing change bytes without changing the physical sample | Use checksums only for exact assets; resolve physical identity from metadata |

## Reproducible build

Run the stages in this order. All stages are CPU-only, inspect only allowlisted
small bundled metadata, and write versioned TSV/JSON control artifacts:

```bash
python scripts/r01_inventory_metadata.py
python scripts/r01_extract_atlas.py
python scripts/r01_extract_hest.py
python scripts/r01_extract_htan.py
python scripts/r01_cross_source_duplicates.py
python scripts/r01_summarize_units.py
python scripts/r01_extract_tenx_explicit.py
python scripts/r01_extract_external_geo.py
python scripts/r01_build_registry.py
python scripts/r01_freeze_roles.py --freeze-date 2026-07-31
python scripts/r01_validate_gate.py
python scripts/r01_prepare_metadata_request.py --approval-status APPROVED_METADATA_ONLY
python scripts/update_bioinf_data_index.py --generated-at 2026-07-31T06:16:28+08:00
python -m unittest discover -s tests -p 'test_r01_*.py'
env PYTHONPATH=. pytest -q tests
```

The current machine gate is `COMPLETE_WITH_EXCLUSIONS`: six canonical logical
units pass identity eligibility and have one frozen role each. Four GEO units
use the documented E2 specimen equivalence; HTAN CRC and 10x Breast Block A
retain explicit patient-block identity. GSE211956 is the locked independent
external lineage. This completes identity infrastructure only and does not
validate any structure ground truth. Before any
later unpacking or derived-data task, recheck that at least 2 TB remains
available. The validator returns zero when it successfully writes a gate
artifact, including a blocked artifact; automation must inspect
`r01_gate.json::status` rather than treating process exit alone as scientific
success.

`r01_metadata_request.tsv` is a proposed metadata-only acquisition queue, not
a role freeze. Its ranking uses bundled identity coverage, provenance
locators, and conflict counts; it excludes TLS/structure outcomes and model
results. Shared accession, PMID, or DOI aliases are grouped conservatively,
and every proposed group contributes zero eligible units until an acquired
patient/block crosswalk passes duplicate and independence audit.
