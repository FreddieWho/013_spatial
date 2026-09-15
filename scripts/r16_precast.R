#!/usr/bin/env Rscript
# R-16 joint-embedding Arm P: PRECAST (intrinsic CAR spatial model).
# Discovery role ONLY (Vanderbilt bridge input). Exploratory; no claims.
# Pipeline: Seurat list (row/col metadata from coords.csv) ->
# CreatePRECASTObject(customGenelist=HVG-10k) ->
# AddAdjList(type=fixed_number, k=6; scale-free, no pixel-unit assumption) ->
# PRECAST seqK + MBIC SelectModel -> save RDS + cluster labels CSV.
# Usage: Rscript scripts/r16_precast.R <bridge_dir> <out_dir> [Kmin] [Kmax]

suppressPackageStartupMessages({
  library(Seurat)
  library(PRECAST)
  library(jsonlite)
})

SEED <- 20260914
set.seed(SEED)
KNN_FIXED <- 6L

args <- commandArgs(trailingOnly = TRUE)
bridge_dir <- args[1]
out_dir <- args[2]
kmin <- ifelse(length(args) >= 3, as.integer(args[3]), 4L)
kmax <- ifelse(length(args) >= 4, as.integer(args[4]), 10L)
# ngene>0: keep top-ngene genes by pooled variance (smoke/cost control).
# Full runs use all common genes; any cap is disclosed in run.log + registry.
ngene <- ifelse(length(args) >= 5, as.integer(args[5]), 0L)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logf <- file(file.path(out_dir, "run.log"), open = "wt")
logmsg <- function(...) { cat(..., "\n", file = logf); flush(logf) }
logmsg("PRECAST ", as.character(packageVersion("PRECAST")),
       " Seurat ", as.character(packageVersion("Seurat")), " seed ", SEED)

# NOTE: simplifyVector=FALSE keeps manifest as a list of row-lists;
# default simplification collapses it to a data.frame (length() = ncols!).
manifest <- fromJSON(file.path(bridge_dir, "manifest.json"), simplifyVector = FALSE)
message("sections: ", length(manifest))

seuList <- lapply(seq_along(manifest), function(i) {
  m <- manifest[[i]]
  d <- file.path(bridge_dir, m$stem)
  counts <- Read10X(data.dir = d, gene.column = 1)
  so <- CreateSeuratObject(counts = counts, project = m$stem,
                           min.cells = 3, min.features = 200)
  co <- read.csv(file.path(d, "coords.csv"))
  rownames(co) <- co$barcode
  cb <- colnames(so)
  so$row <- co[cb, "y"]
  so$col <- co[cb, "x"]
  so$section <- m$stem
  so$patient <- m$patient_id
  so
})
logmsg("loaded: ", length(seuList), " sections")

common_genes <- Reduce(intersect, lapply(seuList, rownames))
logmsg("common genes across sections: ", length(common_genes))
if (ngene > 0 && length(common_genes) > ngene) {
  v <- rowMeans(do.call(cbind, lapply(seuList, function(s) {
    m <- as.matrix(GetAssayData(s, slot = "counts")[common_genes, ])
    (m - rowMeans(m))^2
  })))
  common_genes <- names(sort(v, decreasing = TRUE))[seq_len(ngene)]
  logmsg("capped to top-", ngene, " variance genes")
}

po <- CreatePRECASTObject(seuList, customGenelist = common_genes,
                          rawData.preserve = FALSE, verbose = TRUE)
po <- AddAdjList(po, type = "fixed_number", number = KNN_FIXED)
po <- AddParSetting(po, coreNum = 8, seed = SEED,
                    maxIter = 20, verbose = TRUE)
KGRID <- seq(kmin, kmax)
logmsg("fitting seqK: ", paste(KGRID, collapse = ","))
po <- PRECAST(po, K = KGRID)
po <- SelectModel(po)
logmsg("model selection done")

saveRDS(po, file.path(out_dir, "precast_obj.rds"))

# ---- cluster extraction (primary path + structural sidecar fallback) ------
sidecar <- list(resList_names = names(po@resList),
                resList_class = class(po@resList))
write_json(sidecar, file.path(out_dir, "resList_sidecar.json"),
           auto_unbox = TRUE, pretty = TRUE)
# resList$cluster is an UNNAMED per-section list of n×1 matrices aligned
# with po@seulist (filtered cells); barcodes come from colnames(seulist).
# (Smoke 2026-09-16 caught the old names()-based branch writing empty output.)
tryCatch({
  clu <- po@resList$cluster
  stopifnot(is.list(clu), length(clu) == length(manifest))
  tbl <- do.call(rbind, lapply(seq_along(clu), function(i) {
    bcs <- colnames(po@seulist[[i]])
    stopifnot(length(bcs) == nrow(clu[[i]]))
    data.frame(section = manifest[[i]]$stem, barcode = bcs,
               orig_barcode = bcs, cluster = as.character(clu[[i]][, 1]),
               stringsAsFactors = FALSE)
  }))
  write.csv(tbl, file.path(out_dir, "labels_armP.csv"), row.names = FALSE)
  logmsg("labels extracted: ", nrow(tbl), " spots, K groups: ",
         paste(sort(unique(tbl$cluster)), collapse = ","))
}, error = function(e) {
  logmsg("EXTRACTION_FALLBACK_NEEDED: ", conditionMessage(e))
})
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
close(logf)
message("OK")
