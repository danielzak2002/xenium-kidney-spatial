#!/usr/bin/env Rscript
# 01_load_qc.R — Phase A step 1: load Xenium, imaging-ST QC, save object.
#
#   XENIUM_DATASET=preview Rscript R/01_load_qc.R     # default
#   XENIUM_DATASET=big     Rscript R/01_load_qc.R
#
# Xenium QC philosophy (NOT scRNA-seq): targeted low-plex => low per-cell counts.
# We do NOT hard-filter on min counts/genes. We FLAG (not drop) elevated
# negative-control / blank-codeword cells and candidate segmentation merges.
# The only removal is zero-RNA-count cells, which cannot be normalized.

suppressPackageStartupMessages({
  library(Seurat); library(arrow); library(ggplot2); library(patchwork)
})
source(file.path(Sys.getenv("XENIUM_ROOT", unset = getwd()), "R", "config.R"))

cfg <- ensure_dirs(get_config())
set.seed(cfg$seed)
message("== 01_load_qc :: ", cfg$label, " ==")

# ---- Load (lean: no per-molecule coords, no polygon segmentations) ----------
# Scale-up note (BIG): LoadXenium 5.5.0 slot.map has no entry for the 15
# 'Deprecated Codeword' features in the RCC panel and may error; if so, load the
# matrix manually and drop that feature class before CreateSeuratObject.
obj <- tryCatch(
  LoadXenium(cfg$data_dir, assay = cfg$assay,
             cell.centroids = TRUE, molecule.coordinates = FALSE),
  error = function(e) stop("LoadXenium failed for ", cfg$label, ": ",
                           conditionMessage(e),
                           "\n(If dataset='big', see the Deprecated Codeword note above.)"))
message("  loaded ", ncol(obj), " cells; assays: ", paste(Assays(obj), collapse = ", "))

# ---- Authoritative per-cell QC metrics from cells.parquet (lean) ------------
cells <- as.data.frame(arrow::read_parquet(
  cfg$cells_parquet,
  col_select = c("cell_id","x_centroid","y_centroid","transcript_counts",
                 "control_probe_counts","control_codeword_counts",
                 "unassigned_codeword_counts","cell_area","nucleus_area")))
idx <- match(colnames(obj), cells$cell_id)
if (anyNA(idx)) stop(sum(is.na(idx)), " cells in matrix not found in cells.parquet")
cells <- cells[idx, ]

den <- with(cells, transcript_counts + control_probe_counts +
              control_codeword_counts + unassigned_codeword_counts)
den[den == 0] <- NA_real_
obj$x_centroid   <- cells$x_centroid
obj$y_centroid   <- cells$y_centroid
obj$cell_area    <- cells$cell_area
obj$nucleus_area <- cells$nucleus_area
obj$nuc_cyto_ratio <- cells$nucleus_area / cells$cell_area
obj$neg_frac     <- (cells$control_probe_counts + cells$control_codeword_counts) / den
obj$blank_frac   <- cells$unassigned_codeword_counts / den
obj$neg_frac[is.na(obj$neg_frac)]     <- 0
obj$blank_frac[is.na(obj$blank_frac)] <- 0

# ---- QC flags (recorded, not applied except zero-count removal) -------------
nc <- obj$nCount_Xenium
count_thr <- quantile(nc,           cfg$seg_merge_count_quantile, na.rm = TRUE)
area_thr  <- quantile(obj$cell_area, cfg$seg_merge_area_quantile,  na.rm = TRUE)
obj$flag_zero  <- nc < cfg$min_counts_keep
obj$flag_neg   <- obj$neg_frac   > cfg$neg_frac_flag
obj$flag_blank <- obj$blank_frac > cfg$blank_frac_flag
obj$flag_seg_merge <- (nc > count_thr) & (obj$cell_area > area_thr)

# ---- Apply minimal removal: zero-count cells only ---------------------------
n0 <- ncol(obj)
keep <- which(!obj$flag_zero)
if (cfg$drop_flagged)
  keep <- intersect(keep, which(!(obj$flag_neg | obj$flag_blank | obj$flag_seg_merge)))
obj <- subset(obj, cells = colnames(obj)[keep])
message("  removed ", n0 - ncol(obj), " cell(s) (zero-count",
        if (cfg$drop_flagged) " + flagged" else "", "); kept ", ncol(obj))

# ---- QC summary table -------------------------------------------------------
qc <- data.frame(
  dataset            = cfg$label,
  n_cells_loaded     = n0,
  n_cells_kept       = ncol(obj),
  median_counts      = median(obj$nCount_Xenium),
  median_genes       = median(obj$nFeature_Xenium),
  mean_neg_frac      = mean(obj$neg_frac),
  mean_blank_frac    = mean(obj$blank_frac),
  median_cell_area   = median(obj$cell_area),
  median_nuc_cyto    = median(obj$nuc_cyto_ratio, na.rm = TRUE),
  n_flag_neg         = sum(obj$flag_neg),
  n_flag_blank       = sum(obj$flag_blank),
  n_flag_seg_merge   = sum(obj$flag_seg_merge),
  seg_merge_count_thr = as.numeric(count_thr),
  seg_merge_area_thr  = as.numeric(area_thr))
write.csv(qc, tab_path(cfg, "qc_summary.csv"), row.names = FALSE)
message("  wrote ", tab_path(cfg, "qc_summary.csv"))
print(t(qc))

# ---- Figures ----------------------------------------------------------------
md <- obj[[]]

p_vln <- VlnPlot(obj, c("nCount_Xenium","nFeature_Xenium"), pt.size = 0,
                 ncol = 2) & theme(legend.position = "none")
save_png(p_vln, fig_path(cfg, "qc_counts_violin.png"), width = 8, height = 4)

p_ctrl <- (ggplot(md, aes(neg_frac)) +
             geom_histogram(bins = 60) +
             geom_vline(xintercept = cfg$neg_frac_flag, colour = "red", linetype = 2) +
             labs(x = "neg-control fraction / cell", y = "cells") + theme_bw()) |
  (ggplot(md, aes(blank_frac)) +
     geom_histogram(bins = 60) +
     geom_vline(xintercept = cfg$blank_frac_flag, colour = "red", linetype = 2) +
     labs(x = "blank-codeword fraction / cell", y = "cells") + theme_bw())
save_png(p_ctrl, fig_path(cfg, "qc_control_fractions.png"), width = 9, height = 4)

p_seg <- ggplot(md, aes(cell_area, nCount_Xenium, colour = flag_seg_merge)) +
  geom_point(size = 0.3, alpha = 0.4) +
  geom_hline(yintercept = count_thr, linetype = 2) +
  geom_vline(xintercept = area_thr,  linetype = 2) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "red"),
                      name = "candidate\nseg merge") +
  labs(x = "cell area (um^2)", y = "RNA counts") + theme_bw()
save_png(p_seg, fig_path(cfg, "qc_segmentation_merges.png"), width = 7, height = 5)

p_sp <- ggplot(md, aes(x_centroid, y_centroid, colour = log1p(nCount_Xenium))) +
  geom_point(size = 0.15) + coord_fixed() +
  scale_colour_viridis_c(name = "log1p(counts)") +
  labs(x = NULL, y = NULL, title = paste0(cfg$label, " — counts in situ")) +
  theme_void() + theme(legend.position = "right")
save_png(p_sp, fig_path(cfg, "qc_spatial_counts.png"), width = 8, height = 7)

# ---- Persist ----------------------------------------------------------------
saveRDS(obj, cfg$obj_01)
message("  wrote ", cfg$obj_01)
write_session_info(cfg, "01")
message("== 01 done ==")
