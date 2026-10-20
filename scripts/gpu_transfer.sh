#!/bin/bash
# Transfer GPU probe bundle to the compute box. Run from repo root.
set -e
: "${SSHPASS:?SSHPASS env var required (never store in repo)}"
SSH="sshpass -e ssh -o StrictHostKeyChecking=no -p 21312"
DST="root@180.127.11.169:/root/probe/"
cd /home/huyudi/013_spatial
for f in external/istar infra/gpu_probe_20261007 infra/istar_probe_20261007 \
         infra/tls_weaklabel_probe_20261006/weights repo/data/TLS_signature; do
  rsync -az --exclude='__pycache__' -e "$SSH" "$f" "$DST"
done
rsync -az -e "$SSH" external/HIPT/HIPT_4K/Checkpoints/ "$DST/external/istar/checkpoints/"
rsync -az -e "$SSH" \
  data/other_sources/hest1k/wsis/TENX13.tif \
  data/other_sources/hest1k/wsis/NCBI642.tif \
  data/other_sources/hest1k/st/TENX13.h5ad \
  data/other_sources/hest1k/st/NCBI642.h5ad \
  "$DST/data/"
for g in GSM5924053_frozen_c_57; do
  rsync -az -e "$SSH" \
    data/GEO/GSE175540/raw/${g}_tissue_hires_image.png.gz \
    data/GEO/GSE175540/raw/${g}_tissue_positions_list.csv.gz \
    data/GEO/GSE175540/raw/${g}_filtered_feature_bc_matrix.h5 \
    data/GEO/GSE175540/raw/${g}_TLS_annotation.csv.gz \
    data/GEO/GSE175540/raw/${g}_scalefactors_json.json.gz \
    "$DST/data/"
done
echo TRANSFER_DONE
$SSH root@180.127.11.169 'du -sh /root/probe; ls /root/probe/external/istar/checkpoints/'
