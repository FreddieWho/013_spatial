# R-02 structure registry

This is a metadata/schema-only, fail-closed control plane. It does not read
expression matrices, image pixels, use the network, or use a GPU.

The scoped TLS inventory contains 57 source-reported Table S4 TLS IDs:
HTAN Vanderbilt CRC 44, GSE226997 4, GSE274103 8, and GSE274557 1. These are
cross-checked against Table S2 counts rather than expanded from count ordinals.
All 57 structure rows remain
non-confirmatory because replayable instance geometry and the full
physical/provenance/overlap boundary are not closed.

Nineteen local STOmicsDB TLS annotation files are inventoried by path,
checksum, and compressed CSV header only. They are outside the R-01 frozen
core and are excluded from the 57 scoped summaries. No STOmics row is accepted
as confirmatory GT.

Outer splits use a conservative patient-wide envelope, including all explicit
blocks from the same patient. Explicit blocks remain separately registered for
block-level counting. R-01 patient-linked physical-specimen exceptions stay
block-empty and are not block-level evidence. The one
explicit TENX section pair is registered as an identity link, not a confirmed
cross-section structure correspondence.

The stored gate is expected to be `HARD_BLOCKED_NO_AUDITABLE_GT`, with zero
confirmatory instances and no second claim-bearing structure frozen.
