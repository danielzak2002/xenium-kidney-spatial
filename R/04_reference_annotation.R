#!/usr/bin/env Rscript
# 04_reference_annotation.R — Phase A step 4: reference-based annotation, LAYERED
# onto the marker-based clusters from 02. IMMUNE LAYER (this script): SingleR with
# the celldex Monaco immune reference (blood-derived, fine subsets).
#
#   XENIUM_DATASET=preview Rscript R/04_reference_annotation.R
#
# Targeted 377-gene panel -> SingleR (correlation over shared genes) is robust
# where anchor/integration methods are not. We label per-cell, then take a
# per-cluster majority CONSENSUS against the Leiden clusters and gate by within-
# cluster agreement. Reference labels are applied ONLY to immune clusters; tumour,
# LowQ_MTRNR2L and Mast keep their marker labels (no analog in a blood reference).
# Kidney/epithelial subtypes are a separate (later) reference layer.

suppressPackageStartupMessages({
  library(Seurat); library(SingleR); library(celldex)
  library(SummarizedExperiment); library(ggplot2); library(dplyr)
})
source(file.path(Sys.getenv("XENIUM_ROOT", unset = getwd()), "R", "config.R"))

cfg <- ensure_dirs(get_config())
set.seed(cfg$seed)
message("== 04_reference_annotation (immune/SingleR) :: ", cfg$label, " ==")

obj <- readRDS(cfg$obj_02)
DefaultAssay(obj) <- cfg$assay
clus <- obj$seurat_clusters

# ---- Reference (cached by celldex/gypsum in the user cache, outside the repo) -
ref <- tryCatch(celldex::fetchReference("monaco_immune", cfg$ref_monaco_version),
                error = function(e) celldex::MonacoImmuneData())
assay(ref, "logcounts") <- as.matrix(assay(ref, "logcounts"))  # realize (drop HDF5 backing)
fine2main <- tapply(as.character(ref$label.main), ref$label.fine,
                    function(x) names(sort(table(x), decreasing = TRUE))[1])
message("  reference: Monaco immune (", ncol(ref), " samples, ",
        length(unique(ref$label.fine)), " fine labels)")

# ---- Per-cell SingleR on log-normalized Xenium counts -----------------------
query <- GetAssayData(obj, assay = cfg$assay, layer = "data")
pred <- SingleR(test = query, ref = ref, labels = ref$label.fine,
                de.method = "classic")
obj$singler_fine   <- pred$labels
obj$singler_main   <- unname(fine2main[pred$labels])
obj$singler_pruned <- pred$pruned.labels            # NA = low-confidence call
obj$ref_lowconf    <- is.na(pred$pruned.labels)
n_shared <- length(intersect(rownames(query), rownames(ref)))
message("  SingleR done on ", ncol(obj), " cells (", n_shared,
        " shared genes); low-confidence (pruned NA): ", sum(obj$ref_lowconf))

# ---- Per-cluster consensus (majority of non-NA pruned fine labels) ----------
consensus <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(c(NA, 0))
  tb <- sort(table(x), decreasing = TRUE); c(names(tb)[1], tb[1] / length(x)) }
cl_levels <- levels(clus)
marker_type <- vapply(cl_levels, function(k)
  as.character(obj$cell_type[which(clus == k)[1]]), character(1))
fine_cons <- t(vapply(cl_levels, function(k)
  consensus(obj$singler_pruned[clus == k]), character(2)))
main_cons <- t(vapply(cl_levels, function(k)
  consensus(obj$singler_main[clus == k & !obj$ref_lowconf]), character(2)))

# ---- Layered reconciliation -------------------------------------------------
# Per-cell fine labels fragment within a cluster (subtypes), so we adopt the
# DOMINANT fine label when it is lineage-concordant with the marker type (via
# fine->main mapping) — this refines (Plasma->Plasmablasts, T->CD8/CD4 subsets,
# NK->NK) without relabelling across lineages. Non-concordant immune clusters
# keep the marker label and are flagged for review; ref_keep_marker types
# (tumour/LowQ/Mast) are never overridden; non-immune clusters await the kidney layer.
allowed_main <- list(
  T_cell = c("CD4+ T cells", "CD8+ T cells", "T cells"),
  NK_cell = "NK cells", B_cell = "B cells", Plasma = "B cells",
  Myeloid = c("Monocytes", "Dendritic cells"), DC = "Dendritic cells",
  pDC = "Dendritic cells", Neutrophil = "Neutrophils")
dom_fine      <- fine_cons[, 1]
dom_fine_main <- unname(fine2main[dom_fine])
is_immune  <- marker_type %in% IMMUNE_TYPES
concordant <- mapply(function(mt, fm) mt %in% names(allowed_main) &&
                       !is.na(fm) && fm %in% allowed_main[[mt]],
                     marker_type, dom_fine_main)
adopt <- is_immune & !(marker_type %in% cfg$ref_keep_marker) & concordant

ann <- data.frame(
  cluster        = cl_levels,
  n_cells        = as.integer(table(clus)),
  marker_type    = marker_type,
  singler_main   = main_cons[, 1],
  main_agreement = round(as.numeric(main_cons[, 2]), 3),
  singler_fine   = dom_fine,
  fine_agreement = round(as.numeric(fine_cons[, 2]), 3),
  concordant     = concordant,
  adopted        = adopt,
  low_consensus  = round(as.numeric(fine_cons[, 2]), 3) < cfg$ref_consensus_min,
  frac_lowconf   = round(as.numeric(tapply(obj$ref_lowconf, clus, mean)[cl_levels]), 3),
  row.names = NULL)
ann$ref_cell_type <- ifelse(adopt, dom_fine, marker_type)
write.csv(ann, tab_path(cfg, "cluster_reference_labels.csv"), row.names = FALSE)
message("  clusters adopting a SingleR label: ", sum(adopt), " / ", length(adopt),
        " (immune: ", sum(is_immune), "); non-concordant immune (kept marker): ",
        paste(ann$cluster[is_immune & !concordant & !(marker_type %in% cfg$ref_keep_marker)],
              collapse = ", "))

cl2ref <- setNames(ann$ref_cell_type, ann$cluster)
obj$ref_cell_type <- factor(unname(cl2ref[as.character(clus)]))
ctc <- as.data.frame(table(ref_cell_type = obj$ref_cell_type))
write.csv(ctc, tab_path(cfg, "reference_celltype_counts.csv"), row.names = FALSE)

# ---- Figures ----------------------------------------------------------------
save_png(DimPlot(obj, reduction = "umap", group.by = "ref_cell_type",
                 label = TRUE, repel = TRUE) +
           ggtitle(paste0(cfg$label, " — layered reference cell types")),
         fig_path(cfg, "umap_reference_celltypes.png"), width = 11, height = 7)

md <- obj[[]]
p_sp <- ggplot(md, aes(x_centroid, y_centroid, colour = ref_cell_type)) +
  geom_point(size = 0.2) + coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 1)) +
  labs(x = NULL, y = NULL, title = paste0(cfg$label, " — reference cell types in situ")) +
  theme_void()
save_png(p_sp, fig_path(cfg, "spatial_reference_celltypes.png"), width = 11, height = 8)

# marker-cluster x SingleR-main agreement among immune cells (tile heatmap)
imm <- md[md$cell_type %in% IMMUNE_TYPES & !md$ref_lowconf, ]
if (nrow(imm)) {
  tab <- as.data.frame(prop.table(table(marker = droplevels(factor(imm$cell_type)),
                                        singler_main = imm$singler_main), margin = 1))
  p_h <- ggplot(tab, aes(singler_main, marker, fill = Freq)) +
    geom_tile() + scale_fill_viridis_c(name = "prop") +
    theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "SingleR main label", y = "marker-based type",
         title = paste0(cfg$label, " — marker vs SingleR (immune cells)"))
  save_png(p_h, fig_path(cfg, "marker_vs_singler_heatmap.png"), width = 9, height = 6)
}

# ---- Persist ----------------------------------------------------------------
saveRDS(obj, cfg$obj_04)
message("  wrote ", cfg$obj_04)
write_session_info(cfg, "04")
message("== 04 done ==")
