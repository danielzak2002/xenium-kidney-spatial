#!/usr/bin/env Rscript
# 01_load_qc.R — Phase A step 1: load (Xenium or CosMx), spatial-ST QC, save object.
#
#   XENIUM_DATASET=preview Rscript R/01_load_qc.R          # Xenium (default)
#   XENIUM_DATASET=big     Rscript R/01_load_qc.R
#   XENIUM_DATASET=cln_cosmx Rscript R/01_load_qc.R        # CosMx (Danaher cLN)
#   XENIUM_SUBSAMPLE=30000 XENIUM_DATASET=cln_cosmx Rscript R/01_load_qc.R  # validation subset
#
# QC philosophy (NOT scRNA-seq): targeted panels => low per-cell counts. FLAG (not
# drop) elevated negative-control cells and candidate segmentation merges; the only
# removal is zero-count cells. CosMx neg-probes are not in the matrix, so the per-
# cell negative rate comes from annot$negmean (mapped onto the same flag scheme).

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork)
})
source(file.path(Sys.getenv("XENIUM_ROOT", unset = getwd()), "R", "config.R"))

cfg <- ensure_dirs(get_config())
set.seed(cfg$seed)
message("== 01_load_qc :: ", cfg$label, " (", cfg$platform, ") ==")
assay  <- cfg$assay
nc_col <- paste0("nCount_", assay); nf_col <- paste0("nFeature_", assay)

# ---- Load + per-cell QC source metrics (platform-specific) ------------------
if (cfg$platform == "cosmx") {
  obj <- load_cosmx(cfg)
  obj$cell_area    <- obj$Area                     # CosMx morphology area
  obj$nucleus_area <- NA_real_; obj$nuc_cyto_ratio <- NA_real_
  obj$blank_frac   <- NA_real_                      # no blank codewords in this matrix
  ng <- nrow(obj); nc0 <- obj[[nc_col]][, 1]
  # neg metric (adjustment: no assumed probe count): per-probe background (negmean)
  # relative to mean per-gene signal (nCount/n_genes). >1 => background dominates.
  obj$neg_frac <- obj$negmean * ng / pmax(nc0, 1)
  neg_flag_thr <- cfg$cosmx_neg_ratio_flag
} else {
  library(arrow)
  obj <- load_xenium_lean(cfg)
  cells <- as.data.frame(arrow::read_parquet(cfg$cells_parquet,
    col_select = c("cell_id","x_centroid","y_centroid","transcript_counts",
                   "control_probe_counts","control_codeword_counts",
                   "unassigned_codeword_counts","cell_area","nucleus_area")))
  idx <- match(colnames(obj), cells$cell_id)
  if (anyNA(idx)) stop(sum(is.na(idx)), " cells not found in cells.parquet")
  cells <- cells[idx, ]
  den <- with(cells, transcript_counts + control_probe_counts +
                control_codeword_counts + unassigned_codeword_counts)
  den[den == 0] <- NA_real_
  obj$x_centroid <- cells$x_centroid; obj$y_centroid <- cells$y_centroid
  obj$cell_area  <- cells$cell_area;  obj$nucleus_area <- cells$nucleus_area
  obj$nuc_cyto_ratio <- cells$nucleus_area / cells$cell_area
  obj$neg_frac   <- (cells$control_probe_counts + cells$control_codeword_counts) / den
  obj$blank_frac <- cells$unassigned_codeword_counts / den
  obj$neg_frac[is.na(obj$neg_frac)] <- 0; obj$blank_frac[is.na(obj$blank_frac)] <- 0
  neg_flag_thr <- cfg$neg_frac_flag
}
message("  loaded ", ncol(obj), " cells; assay ", assay, " (", nrow(obj), " features)")

# Optional validation subsetting (applies to any platform).
# XENIUM_SAMPLES: comma-separated sample IDs (cfg$sample_col) to keep.
# XENIUM_SUBSAMPLE: random N cells to keep.
samp_keep <- Sys.getenv("XENIUM_SAMPLES", "")
if (nzchar(samp_keep) && cfg$multi_sample) {
  want <- trimws(strsplit(samp_keep, ",")[[1]]); sv <- obj[[cfg$sample_col]][, 1]
  obj <- subset(obj, cells = colnames(obj)[sv %in% want])
  message("  >> SUBSET to samples {", paste(want, collapse = ", "), "}: ", ncol(obj), " cells")
}
ss <- as.integer(Sys.getenv("XENIUM_SUBSAMPLE", "0"))
if (ss > 0 && ss < ncol(obj)) {
  set.seed(cfg$seed); obj <- subset(obj, cells = sample(colnames(obj), ss))
  message("  >> SUBSAMPLED to ", ncol(obj), " cells (validation run)")
}

# ---- QC flags (recorded, not applied except zero-count removal) -------------
nc <- obj[[nc_col]][, 1]; nf <- obj[[nf_col]][, 1]
count_thr <- quantile(nc,           cfg$seg_merge_count_quantile, na.rm = TRUE)
area_thr  <- quantile(obj$cell_area, cfg$seg_merge_area_quantile,  na.rm = TRUE)
obj$flag_zero  <- nc < cfg$min_counts_keep
obj$flag_neg   <- obj$neg_frac   > neg_flag_thr
obj$flag_blank <- !is.na(obj$blank_frac) & obj$blank_frac > cfg$blank_frac_flag
obj$flag_seg_merge <- (nc > count_thr) & (obj$cell_area > area_thr)

# Low-quality signal: MTRNR2L mito-pseudogene fraction (self-skips if absent, as
# on the CosMx panel) + complexity (nFeature/nCount).
cmat <- GetAssayData(obj, assay = assay, layer = "counts")
mt_rows <- grep("^MTRNR2L", rownames(cmat), value = TRUE)
mt_counts <- (if (length(mt_rows)) Matrix::colSums(cmat[mt_rows, , drop = FALSE])
              else rep(0, ncol(obj)))
obj$mtrnr2l_frac <- as.numeric(mt_counts) / pmax(nc, 1)
obj$complexity   <- nf / pmax(nc, 1)
obj$flag_lowq    <- obj$mtrnr2l_frac > cfg$lowq_mtrnr2l_frac | obj$complexity < cfg$lowq_complexity
message("  MTRNR2L features: ", length(mt_rows), "; flag_lowq cells: ", sum(obj$flag_lowq))

# ---- Apply minimal removal: zero-count cells only ---------------------------
n0 <- ncol(obj)
keep <- which(!obj$flag_zero)
if (cfg$drop_flagged)
  keep <- intersect(keep, which(!(obj$flag_neg | obj$flag_blank | obj$flag_seg_merge)))
obj <- subset(obj, cells = colnames(obj)[keep])
message("  removed ", n0 - ncol(obj), " zero-count cell(s); kept ", ncol(obj))

# ---- QC summary -------------------------------------------------------------
nck <- obj[[nc_col]][, 1]
qc <- data.frame(
  dataset = cfg$label, platform = cfg$platform,
  n_cells_loaded = n0, n_cells_kept = ncol(obj),
  median_counts = median(nck), median_genes = median(obj[[nf_col]][, 1]),
  mean_neg_frac = mean(obj$neg_frac),
  mean_blank_frac = if (all(is.na(obj$blank_frac))) NA else mean(obj$blank_frac),
  median_cell_area = median(obj$cell_area),
  n_flag_neg = sum(obj$flag_neg), n_flag_blank = sum(obj$flag_blank),
  n_flag_seg_merge = sum(obj$flag_seg_merge), n_flag_lowq = sum(obj$flag_lowq),
  n_samples = if (cfg$multi_sample) length(unique(obj[[cfg$sample_col]][, 1])) else 1L)
write.csv(qc, tab_path(cfg, "qc_summary.csv"), row.names = FALSE)
print(t(qc))
if (cfg$multi_sample) {
  st <- as.data.frame.matrix(table(obj[[cfg$sample_col]][, 1], obj[[cfg$condition_col]][, 1]))
  message("  samples x condition:"); print(st)
}

# ---- Figures ----------------------------------------------------------------
md <- obj[[]]; xl <- if (cfg$coord_units == "mm") "x (mm)" else NULL
p_vln <- VlnPlot(obj, c(nc_col, nf_col), pt.size = 0, ncol = 2) &
  theme(legend.position = "none")
save_png(p_vln, fig_path(cfg, "qc_counts_violin.png"), width = 8, height = 4)

negp <- ggplot(md, aes(neg_frac)) + geom_histogram(bins = 60) +
  geom_vline(xintercept = neg_flag_thr, colour = "red", linetype = 2) +
  labs(x = if (cfg$platform == "cosmx") "neg ratio (negmean*nGenes/nCount)" else
       "neg-control fraction / cell", y = "cells") + theme_bw()
save_png(negp, fig_path(cfg, "qc_control_fractions.png"), width = 6, height = 4)

p_seg <- ggplot(md, aes(cell_area, .data[[nc_col]], colour = flag_seg_merge)) +
  geom_point(size = 0.3, alpha = 0.4) +
  geom_hline(yintercept = count_thr, linetype = 2) +
  geom_vline(xintercept = area_thr, linetype = 2) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "red"),
                      name = "candidate\nseg merge") +
  labs(x = "cell area", y = "RNA counts") + theme_bw()
save_png(p_seg, fig_path(cfg, "qc_segmentation_merges.png"), width = 7, height = 5)

p_sp <- ggplot(md, aes(x_centroid, y_centroid, colour = log1p(.data[[nc_col]]))) +
  geom_point(size = 0.15) + coord_fixed() +
  scale_colour_viridis_c(name = "log1p(counts)") +
  labs(x = xl, y = NULL, title = paste0(cfg$label, " — counts in situ")) +
  theme_void() + theme(legend.position = "right")
save_png(p_sp, fig_path(cfg, "qc_spatial_counts.png"), width = 8, height = 7)

# ---- Persist ----------------------------------------------------------------
saveRDS(obj, cfg$obj_01)
message("  wrote ", cfg$obj_01)
write_session_info(cfg, "01")
message("== 01 done ==")
