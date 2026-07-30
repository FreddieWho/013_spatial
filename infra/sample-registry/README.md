# R-01 sample registry

This directory is the versioned control surface for roadmap node R-01. It
answers who a physical sample is, where each asset came from, which assets may
share a leakage group, and which analysis role a logical data unit may carry.
It does not decide whether a structure annotation is valid ground truth
(R-02), perform registration (R-03), inspect expression matrices, or select a
model.

## Evidence boundary

Only metadata bundled with the local data may support identity or structure
facts. Every normalized patient, block, section, serial relationship, or
structure field must retain its raw value, metadata location and evidence
grade:

- `E3_explicit`: bundled metadata explicitly states the relationship.
- `E2_corborated`: at least two consistent bundled metadata fields support it.
- `E1_weak`: filename, directory, isolated alias or incomplete mapping only.
- `E0_unknown`: no auditable evidence.
- `EC_conflict`: bundled metadata sources conflict.

E1/E0/EC records cannot enter claim-bearing validation. Different checksums do
not prove physical independence. A possible cross-source match remains in one
conservative leakage group or outside confirmation until resolved.

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

## Commands

```bash
python -m unittest discover -s tests -p 'test_*.py'
python scripts/r01_inventory_metadata.py
```

Both commands are CPU-only and write only a small TSV. Before any later
unpacking or derived-data task, recheck that at least 2 TB remains available.
