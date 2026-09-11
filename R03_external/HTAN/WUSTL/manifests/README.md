# WUSTL download batches (NOT_TRIGGERED)
Batches W1_visium, W2_visium_hd, W3_xenium, W4_cosmx, W5_codex,
W6_histology_imaging, W7_other_spatial, W8_auxiliary_processed.
Each batch will store: Portal filter description, Synapse-ID download script
(W*_synapse_download.sh), file list, and status TSV.
Blocked on: Portal manifests + Synapse PAT (NEED_USER_SYNAPSE_TOKEN).
Rule: synapse/<synapse_id>/<filename> layout; no silent overwrite; retries kept.
