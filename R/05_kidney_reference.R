#!/usr/bin/env Rscript
# 05_kidney_reference.R — Phase A step 5: KIDNEY/EPITHELIAL reference layer,
# layered on top of the immune labels from 04. Uses the KPMP/HuBMAP Azimuth
# human-kidney reference (Zenodo 10694842, CC BY 4.0; cached git-ignored) at
# annotation.l1 (16 nephron / stroma / endothelial / immune classes).
#
#   XENIUM_DATASET=preview Rscript R/05_kidney_reference.R
#
# NOTE on method: this reference's ref.Rds ships only the supervised-PCA reduction
# (the gene matrix is emptied), so it cannot be run correlation-style via SingleR.
# We therefore use Azimuth's own mapping (SCTransform + projection onto the reference
# SPCA + label transfer) and gate by its prediction score. The immune layer (04)
# stays SingleR/Monaco. The kidney label is adopted ONLY for non-immune, non-(tumour/
# LowQ) clusters; immune clusters keep their Monaco label; Tumor_RCC / LowQ_MTRNR2L /
# Mast are never overridden. Output: a single combined final_ref_cell_type.

suppressPackageStartupMessages({
  library(Seurat); library(Azimuth); library(ggplot2); library(dplyr)
})
source(file.path(Sys.getenv("XENIUM_ROOT", unset = getwd()), "R", "config.R"))

# Azimuth's label transfer routes the (already in-memory) reference + anchor weights
# through future. Pin a sequential plan (no worker copies) and lift the global-export
# guard — at 465k cells the closure exceeds any fixed cap, but nothing is duplicated.
future::plan("sequential")
options(future.globals.maxSize = 1e12)

cfg <- ensure_dirs(get_config())
set.seed(cfg$seed)
message("== 05_kidney_reference (Azimuth KPMP) :: ", cfg$label, " ==")

obj  <- readRDS(cfg$obj_04)          # carries immune ref_cell_type + marker cell_type
DefaultAssay(obj) <- cfg$assay
clus <- obj$seurat_clusters

# ---- Azimuth mapping onto the KPMP kidney reference -------------------------
if (!file.exists(file.path(cfg$ref_kidney_dir, "ref.Rds")))
  stop("Kidney reference not found in ", cfg$ref_kidney_dir,
       " — download per DATA.md (Zenodo 10694842: ref.Rds + idx.annoy).")
pred_col  <- paste0("predicted.", cfg$ref_kidney_level)
score_col <- paste0(pred_col, ".score")
azi <- RunAzimuth(obj, reference = cfg$ref_kidney_dir,
                  annotation.levels = cfg$ref_kidney_level,
                  assay = cfg$assay, verbose = FALSE)
am <- azi[[]]
obj$kidney_l1            <- am[colnames(obj), pred_col]
obj$kidney_score         <- am[colnames(obj), score_col]
obj$kidney_mapping_score <- am[colnames(obj), "mapping.score"]
obj$kidney_lowconf       <- obj$kidney_score < cfg$ref_kidney_score_min
obj$kidney_pruned        <- ifelse(obj$kidney_lowconf, NA, obj$kidney_l1)
message("  Azimuth done on ", ncol(obj), " cells; ",
        length(unique(na.omit(obj$kidney_l1))), " ", cfg$ref_kidney_level,
        " classes; low-confidence (score<", cfg$ref_kidney_score_min, "): ",
        sum(obj$kidney_lowconf))

# ---- Per-cluster consensus + layered reconciliation -------------------------
consensus <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(c(NA, 0))
  tb <- sort(table(x), decreasing = TRUE); c(names(tb)[1], tb[1] / length(x)) }
cl_levels   <- levels(clus)
marker_type <- vapply(cl_levels, function(k)
  as.character(obj$cell_type[which(clus == k)[1]]), character(1))
immune_lab  <- vapply(cl_levels, function(k)              # from 04 (immune layer)
  as.character(obj$ref_cell_type[which(clus == k)[1]]), character(1))
kcons <- t(vapply(cl_levels, function(k)
  consensus(obj$kidney_pruned[clus == k]), character(2)))

# Adopt kidney l1 for non-immune, non-keep clusters with adequate agreement. The
# reference's coarse "Immune" bucket is NOT an adoptable refinement (immune is the
# Monaco layer's job); a non-immune marker cluster that Azimuth confidently calls
# "Immune" is instead FLAGGED as a conflict and keeps its marker label for review.
kid_lab    <- kcons[, 1]
kid_agree  <- as.numeric(kcons[, 2])
non_immune <- !(marker_type %in% IMMUNE_TYPES)
eligible   <- non_immune & !(marker_type %in% cfg$ref_keep_marker) &
              kid_agree >= cfg$ref_consensus_min & !is.na(kid_lab)
adopt_kidney   <- eligible & kid_lab != "Immune"
kidney_immune_conflict <- eligible & kid_lab == "Immune"
final <- immune_lab                              # immune labels + tumour/LowQ/Mast/marker
final[adopt_kidney] <- kid_lab[adopt_kidney]     # override non-immune with nephron label

ann <- data.frame(
  cluster          = cl_levels,
  n_cells          = as.integer(table(clus)),
  marker_type      = marker_type,
  immune_label     = immune_lab,
  kidney_l1        = kid_lab,
  kidney_agreement = round(kid_agree, 3),
  mean_kidney_score = round(as.numeric(tapply(obj$kidney_score, clus, mean)[cl_levels]), 3),
  adopted_kidney   = adopt_kidney,
  kidney_immune_conflict = kidney_immune_conflict,
  final_ref_cell_type = final,
  row.names = NULL)
write.csv(ann, tab_path(cfg, "final_reference_labels.csv"), row.names = FALSE)
message("  clusters adopting a kidney label: ", sum(adopt_kidney), " / ", sum(non_immune),
        " non-immune; kidney-vs-marker 'Immune' conflicts (kept marker, review): ",
        paste(ann$cluster[kidney_immune_conflict], collapse = ", "))

cl2final <- setNames(ann$final_ref_cell_type, ann$cluster)
obj$final_ref_cell_type <- factor(unname(cl2final[as.character(clus)]))
ctc <- as.data.frame(table(final_ref_cell_type = obj$final_ref_cell_type))
write.csv(ctc, tab_path(cfg, "final_celltype_counts.csv"), row.names = FALSE)

# ---- Figures ----------------------------------------------------------------
save_png(DimPlot(obj, reduction = "umap", group.by = "final_ref_cell_type",
                 label = TRUE, repel = TRUE) +
           ggtitle(paste0(cfg$label, " — final layered cell types (immune + kidney)")),
         fig_path(cfg, "umap_final_celltypes.png"), width = 12, height = 8)

md <- obj[[]]
p_sp <- ggplot(md, aes(x_centroid, y_centroid, colour = final_ref_cell_type)) +
  geom_point(size = 0.2) + coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 1)) +
  labs(x = NULL, y = NULL,
       title = paste0(cfg$label, " — final cell types in situ")) +
  theme_void()
save_png(p_sp, fig_path(cfg, "spatial_final_celltypes.png"), width = 12, height = 9)

# ---- Persist ----------------------------------------------------------------
saveRDS(obj, cfg$obj_05)
message("  wrote ", cfg$obj_05)
write_session_info(cfg, "05")
message("== 05 done ==")
