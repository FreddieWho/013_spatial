#!/usr/bin/env Rscript
# R-16 joint-embedding R track (D-119/D-120): dedicated single-cell methods.
# Arm S: Seurat v4 SCTransform + anchor integration (reference-based).
# Arm H: merged SCTransform? No: NormalizeData/ScaleData/PCA + Harmony.
# Discovery role ONLY (Vanderbilt bridge input). Exploratory; no claims.
# Usage: Rscript scripts/r16_seurat_integration.R <bridge_dir> <out_dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(sctransform)
  library(harmony)
  library(jsonlite)
})

SEED <- 20260914
set.seed(SEED)
# PrepSCTIntegration ships ~22GB of globals to future workers (whole object
# list); default 500MB cap aborts. RAM is ample (1TB); raise the cap.
options(future.globals.maxSize = 40 * 1024^3)
# NOTE (Seurat behavior, not a bug in this script): feature names containing
# underscores are rewritten with dashes on load. Gene-level joins between R
# outputs and Python artifacts must normalize '_' -> '-' first.
RESOLUTIONS <- c(0.25, 0.5, 1.0)
N_SCT_FEATURES <- 3000
N_ANCHOR_FEATURES <- 3000
N_REFERENCE <- 5

args <- commandArgs(trailingOnly = TRUE)
bridge_dir <- args[1]
out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
logf <- file(file.path(out_dir, "run.log"), open = "wt")
logmsg <- function(...) { cat(..., "\n", file = logf, append = FALSE); flush(logf) }
logmsg("Seurat ", as.character(packageVersion("Seurat")),
       " sctransform ", as.character(packageVersion("sctransform")),
       " harmony ", as.character(packageVersion("harmony")),
       " seed ", SEED)

# NOTE: simplifyVector=FALSE keeps manifest as a list of row-lists;
# default simplification collapses it to a data.frame (length() = ncols!).
manifest <- fromJSON(file.path(bridge_dir, "manifest.json"), simplifyVector = FALSE)
message("sections: ", length(manifest))

# ---- load + SCTransform per section (checkpointed) -------------------------
ckpt_dir <- file.path(out_dir, "checkpoints")
dir.create(ckpt_dir, showWarnings = FALSE, recursive = TRUE)
objs_ckpt <- file.path(ckpt_dir, "objs_sct.rds")
if (file.exists(objs_ckpt)) {
  objs <- readRDS(objs_ckpt)
  logmsg("resumed SCT objects from checkpoint: ", length(objs), " sections")
} else {
objs <- lapply(seq_along(manifest), function(i) {
  m <- manifest[[i]]
  counts <- Read10X(data.dir = file.path(bridge_dir, m$stem), gene.column = 1)
  so <- CreateSeuratObject(counts = counts, project = m$stem,
                           min.cells = 3, min.features = 200)
  so$section <- m$stem
  so$patient <- m$patient_id
  so <- SCTransform(so, variable.features.n = N_SCT_FEATURES,
                    verbose = FALSE, seed.use = SEED)
  gc()
  so
})
logmsg("SCTransform done: ", length(objs), " sections")
saveRDS(objs, objs_ckpt)
logmsg("checkpoint saved: objs_sct")
}

strip_prefix <- function(bc, sec) sub(paste0("^", sec, "_"), "", bc)

# ---- Arm S: SKIPPED_PER_D119 (2026-09-16) ----------------------------------
# IntegrateData fails at 47-section scale (kNN weight normalize length
# mismatch over ~5.4M anchors); downgraded to method-negative-control
# per D-119 failover. Anchors checkpoint retained. Wrap in if(FALSE).
if (FALSE) {
n_spots <- vapply(objs, ncol, integer(1))
ref_idx <- order(n_spots, decreasing = TRUE)[seq_len(min(N_REFERENCE, length(objs)))]
features <- SelectIntegrationFeatures(object.list = objs, nfeatures = N_ANCHOR_FEATURES)
objs <- PrepSCTIntegration(object.list = objs, anchor.features = features, verbose = FALSE)
anchors_ckpt <- file.path(ckpt_dir, "anchors.rds")
if (file.exists(anchors_ckpt)) {
  anchors <- readRDS(anchors_ckpt)
  logmsg("resumed anchors from checkpoint")
} else {
anchors <- FindIntegrationAnchors(object.list = objs, normalization.method = "SCT",
                                 anchor.features = features, reference = ref_idx,
                                 verbose = FALSE)
saveRDS(anchors, anchors_ckpt)
logmsg("checkpoint saved: anchors")
}
logmsg("anchors found: ", nrow(anchors@anchors))
# IntegrateData collects ~50GB of globals (whole object list + closure env)
# at future() CREATION time, regardless of backend plan — plan(sequential)
# does NOT help (verified: same globals error under sequential, 2026-09-16).
# RAM is ~1TB, so lift the cap entirely instead of raising it piecemeal.
options(future.globals.maxSize = Inf)
int_s <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose = FALSE)
int_s <- RunPCA(int_s, npcs = 30, verbose = FALSE, seed.use = SEED)
int_s <- FindNeighbors(int_s, dims = 1:30, verbose = FALSE)
for (res in RESOLUTIONS) {
  int_s <- FindClusters(int_s, resolution = res, verbose = FALSE,
                        random.seed = SEED,
                        cluster.name = paste0("SCT_", gsub("\\.", "", as.character(res))))
}
s_labels <- data.frame(
  barcode = colnames(int_s),
  section = int_s$section,
  orig_barcode = mapply(strip_prefix, colnames(int_s), int_s$section, USE.NAMES = FALSE),
  SCT_025 = int_s$SCT_025, SCT_05 = int_s$SCT_05, SCT_10 = int_s$SCT_10
)
write.csv(s_labels, file.path(out_dir, "labels_armS.csv"), row.names = FALSE)
logmsg("Arm S done")
} # end SKIPPED_PER_D119

# ---- Arm H: merged log-normalize + Harmony --------------------------------
merged <- merge(objs[[1]], objs[-1], add.cell.ids = vapply(objs, function(o) o$section[1], character(1)))
DefaultAssay(merged) <- "RNA"
merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, nfeatures = N_SCT_FEATURES, verbose = FALSE)
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE, seed.use = SEED)
merged <- RunHarmony(merged, group.by.vars = "section", verbose = FALSE)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30, verbose = FALSE)
# Seurat v4 FindClusters has no cluster.name arg (v5 API; silently ignored).
# Store each resolution's Idents into its HARM_* column explicitly.
# NOTE: suffixes are positional, NOT derived from as.character(res):
# as.character(1.0)=="1" would give HARM_1, breaking the _025/_05/_10
# convention used by the scorer and the GraphST arm (fixed 2026-09-16).
HARM_SUFFIX <- c("025", "05", "10")
for (k in seq_along(RESOLUTIONS)) {
  res <- RESOLUTIONS[k]
  merged <- FindClusters(merged, resolution = res, verbose = FALSE, random.seed = SEED)
  cname <- paste0("HARM_", HARM_SUFFIX[k])
  merged[[cname]] <- Idents(merged)
}
h_labels <- data.frame(
  barcode = colnames(merged),
  section = merged$section,
  orig_barcode = mapply(strip_prefix, colnames(merged), merged$section, USE.NAMES = FALSE),
  HARM_025 = merged$HARM_025, HARM_05 = merged$HARM_05, HARM_10 = merged$HARM_10
)
write.csv(h_labels, file.path(out_dir, "labels_armH.csv"), row.names = FALSE)
logmsg("Arm H done")

params <- list(seed = SEED, resolutions = RESOLUTIONS,
               sct_features = N_SCT_FEATURES, anchor_features = N_ANCHOR_FEATURES,
               arm_S = "SKIPPED_PER_D119", arm_H = "merged_logNormalize_Harmony",
               n_sections = length(objs))
write_json(params, file.path(out_dir, "params.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
close(logf)
message("OK")
