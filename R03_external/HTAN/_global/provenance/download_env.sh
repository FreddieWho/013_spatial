# Source before any synapse download work: API via proxy (required),
# bulk S3 bytes direct (fast path). No secrets in this file.
# API (repo-prod) direct == 403 region block; via 127.0.0.x proxy == 200.
# S3 direct reachable (~0.6s handshake); AWS has no domestic mirror,
# so bulk bytes cross the Pacific either way — direct avoids proxy bottleneck.
export NO_PROXY="$NO_PROXY,amazonaws.com,s3.amazonaws.com,googleapis.com,storage.googleapis.com"
export no_proxy="$no_proxy,amazonaws.com,s3.amazonaws.com,googleapis.com,storage.googleapis.com"
