#!/usr/bin/env Rscript
# 13_whitepaper_embeddings.R — export UMAP / integration embeddings for whitepaper
# Group B, from the committed _05 objects. Subsampled for plotting (embeddings are a
# visualization; full-data structure is preserved in the subsample). No re-clustering.
#
#   Rscript R/13_whitepaper_embeddings.R

suppressPackageStartupMessages({ library(Seurat) })
ROOT <- Sys.getenv("XENIUM_ROOT", unset = getwd())
OBJ  <- file.path(ROOT, "outputs", "objects")
set.seed(1); NSUB <- 60000L

DS <- list(
  list(label = "kidney_RCC_protein",  lab_col = "final_ref_cell_type"),
  list(label = "kidney_preview_PRCC", lab_col = "final_ref_cell_type"),
  list(label = "cln_cosmx",           lab_col = "insitutype"))

# ---- per-dataset UMAP + label ----------------------------------------------
for (d in DS) {
  obj <- readRDS(file.path(OBJ, paste0(d$label, "_05_refann_kidney.rds")))
  emb <- Embeddings(obj, "umap"); lab <- as.character(obj[[d$lab_col]][, 1])
  keep <- if (ncol(obj) > NSUB) sort(sample(ncol(obj), NSUB)) else seq_len(ncol(obj))
  out <- data.frame(umap1 = emb[keep, 1], umap2 = emb[keep, 2], label = lab[keep])
  write.csv(out, file.path(OBJ, paste0("wp_umap_", d$label, ".csv")), row.names = FALSE)
  message("wrote wp_umap_", d$label, ".csv (", nrow(out), " cells, ",
          length(unique(out$label)), " labels)")
  rm(obj); gc(verbose = FALSE)
}

# ---- cLN Harmony before/after (slides mixing without erasing biology) -------
obj <- readRDS(file.path(OBJ, "cln_cosmx_05_refann_kidney.rds"))
set.seed(1); keep <- sort(sample(ncol(obj), NSUB))
sub <- subset(obj, cells = colnames(obj)[keep]); rm(obj); gc(verbose = FALSE)
# BEFORE: UMAP on raw PCA (no integration). AFTER: UMAP on Harmony embedding.
sub <- RunUMAP(sub, reduction = "pca",     dims = 1:30, reduction.name = "umap_pca",
               verbose = FALSE)
sub <- RunUMAP(sub, reduction = "harmony", dims = 1:30, reduction.name = "umap_harm",
               verbose = FALSE)
ep <- Embeddings(sub, "umap_pca"); eh <- Embeddings(sub, "umap_harm")
out <- data.frame(pca_umap1 = ep[, 1], pca_umap2 = ep[, 2],
                  harm_umap1 = eh[, 1], harm_umap2 = eh[, 2],
                  slide = as.character(sub$slidename),
                  label = as.character(sub$insitutype))
write.csv(out, file.path(OBJ, "wp_harmony_cln.csv"), row.names = FALSE)
message("wrote wp_harmony_cln.csv (", nrow(out), " cells, ",
        length(unique(out$slide)), " slides)")
message("== 13 done ==")
