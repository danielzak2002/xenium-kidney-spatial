#!/usr/bin/env Rscript
# 02_cluster_annotate.R — Phase A step 2: normalize -> cluster -> annotate ->
# spatial smoke-test. Consumes the object written by 01_load_qc.R.
#
#   XENIUM_DATASET=preview Rscript R/02_cluster_annotate.R
#
# Normalization = LogNormalize (operator choice; SCT switch in config).
# Clustering = native igraph Leiden on the SNN graph (no python/reticulate).
# Annotation is marker-based and PROVISIONAL: a cluster x cell-type score matrix
# + argmax assignment, left for expert review. B-cell / plasma emphasis.

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
  library(RANN); library(Matrix)
})
source(file.path(Sys.getenv("XENIUM_ROOT", unset = getwd()), "R", "config.R"))

cfg <- ensure_dirs(get_config())
set.seed(cfg$seed)
message("== 02_cluster_annotate :: ", cfg$label, " ==")

obj <- readRDS(cfg$obj_01)
DefaultAssay(obj) <- cfg$assay

# ---- Normalize / reduce -----------------------------------------------------
if (cfg$norm_method == "SCT") {
  obj <- SCTransform(obj, assay = cfg$assay, verbose = FALSE)
} else {
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = cfg$norm_scale_factor, verbose = FALSE)
  # Targeted panel: use ALL genes as features (HVG selection is inappropriate).
  VariableFeatures(obj) <- if (cfg$use_all_features) rownames(obj)
    else VariableFeatures(FindVariableFeatures(obj, nfeatures = cfg$n_variable_features))
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
}
obj <- RunPCA(obj, npcs = cfg$npcs_compute, verbose = FALSE,
              features = VariableFeatures(obj))
obj <- RunUMAP(obj, dims = cfg$dims, verbose = FALSE)
obj <- FindNeighbors(obj, dims = cfg$dims, verbose = FALSE)

# ---- Cluster (native Leiden) ------------------------------------------------
graph_name <- paste0(DefaultAssay(obj), "_snn")
clusters <- run_leiden(obj, graph_name, cfg$cluster_resolution, cfg$seed,
                       cfg$leiden_objective, cfg$min_cluster_size,
                       emb = Embeddings(obj, "pca")[, cfg$dims])
obj$seurat_clusters <- clusters
Idents(obj) <- clusters
message("  Leiden (res=", cfg$cluster_resolution, ") -> ",
        nlevels(clusters), " clusters; sizes: ",
        paste(as.integer(table(clusters)), collapse = ", "))

# ---- Marker DE --------------------------------------------------------------
markers <- FindAllMarkers(obj, only.pos = cfg$de_only_pos, min.pct = cfg$de_min_pct,
                          logfc.threshold = cfg$de_logfc, verbose = FALSE)
write.csv(markers, tab_path(cfg, "cluster_markers_all.csv"), row.names = FALSE)
top <- markers %>% group_by(cluster) %>%
  slice_max(avg_log2FC, n = cfg$de_top_n, with_ties = FALSE) %>% ungroup()
write.csv(top, tab_path(cfg, "cluster_markers_top.csv"), row.names = FALSE)
message("  wrote marker tables (", nrow(markers), " rows)")

# ---- Provisional annotation: cluster x cell-type score ----------------------
# Score = mean z-scored (scale.data) expression of each marker set per cell,
# averaged per cluster; argmax => provisional label. Transparent, no binning.
sets <- lapply(marker_sets(), intersect, rownames(obj))
sets <- sets[lengths(sets) > 0]
scaled <- GetAssayData(obj, layer = "scale.data")
present <- lapply(sets, intersect, rownames(scaled))
present <- present[lengths(present) > 0]
cell_scores <- sapply(present, function(g)
  if (length(g) == 1) scaled[g, ] else Matrix::colMeans(scaled[g, , drop = FALSE]))
clus <- obj$seurat_clusters
clus_scores <- t(apply(cell_scores, 2, function(v) tapply(v, clus, mean)))  # type x cluster
clus_scores <- t(clus_scores)                                              # cluster x type
assign_type <- colnames(clus_scores)[max.col(clus_scores, ties.method = "first")]
names(assign_type) <- rownames(clus_scores)

ann <- data.frame(cluster = rownames(clus_scores), round(clus_scores, 3),
                  assigned = assign_type, n_cells = as.integer(table(clus)),
                  check.names = FALSE)
write.csv(ann, tab_path(cfg, "cluster_celltype_scores.csv"), row.names = FALSE)
obj$cell_type <- factor(unname(assign_type[as.character(clus)]))
ct_counts <- as.data.frame(table(cell_type = obj$cell_type))
write.csv(ct_counts, tab_path(cfg, "celltype_counts.csv"), row.names = FALSE)
message("  provisional cell types: ",
        paste(sprintf("%s=%d", ct_counts$cell_type, ct_counts$Freq), collapse = "  "))

# ---- Figures ----------------------------------------------------------------
save_png(DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
                 label = TRUE) + ggtitle(paste0(cfg$label, " — clusters")),
         fig_path(cfg, "umap_clusters.png"), width = 8, height = 6)
save_png(DimPlot(obj, reduction = "umap", group.by = "cell_type", label = TRUE,
                 repel = TRUE) + ggtitle(paste0(cfg$label, " — provisional cell types")),
         fig_path(cfg, "umap_celltypes.png"), width = 9, height = 6)

md <- obj[[]]
p_spc <- ggplot(md, aes(x_centroid, y_centroid, colour = cell_type)) +
  geom_point(size = 0.2) + coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  labs(x = NULL, y = NULL, title = paste0(cfg$label, " — cell types in situ")) +
  theme_void()
save_png(p_spc, fig_path(cfg, "spatial_celltypes.png"), width = 9, height = 7)

dot_feats <- unique(unlist(lapply(present, head, 3)))
save_png(DotPlot(obj, features = dot_feats, group.by = "seurat_clusters") +
           RotatedAxis() + ggtitle(paste0(cfg$label, " — markers x cluster")),
         fig_path(cfg, "dotplot_markers.png"),
         width = max(8, length(dot_feats) * 0.28), height = 7)

# B/plasma highlight in situ (project goal)
b_genes <- intersect(c("MS4A1","CD79A","MZB1"), rownames(obj))
if (length(b_genes)) {
  bp <- FetchData(obj, vars = b_genes, layer = "data")
  md$bcell_signal <- rowMeans(bp)
  p_b <- ggplot(md[order(md$bcell_signal), ],
                aes(x_centroid, y_centroid, colour = bcell_signal)) +
    geom_point(size = 0.25) + coord_fixed() +
    scale_colour_viridis_c(option = "magma",
                           name = paste(b_genes, collapse = "/")) +
    labs(x = NULL, y = NULL, title = paste0(cfg$label, " — B/plasma signal")) +
    theme_void()
  save_png(p_b, fig_path(cfg, "spatial_bcell_signal.png"), width = 8, height = 7)
}

# ---- Spatial smoke-test: immune-neighbor fraction (Phase B does the real stats)
if (cfg$run_spatial_smoke) {
  coords <- as.matrix(md[, c("x_centroid","y_centroid")])
  nn <- RANN::nn2(coords, k = cfg$spatial_k + 1)$nn.idx[, -1, drop = FALSE]
  is_immune <- obj$cell_type %in% IMMUNE_TYPES
  imm_frac <- rowMeans(matrix(is_immune[nn], nrow = nrow(nn)))
  md$immune_nbr_frac <- imm_frac
  obj$immune_nbr_frac <- imm_frac
  obj$candidate_immune_aggregate <- imm_frac > cfg$immune_nbr_frac
  n_agg <- sum(obj$candidate_immune_aggregate)
  message("  spatial smoke-test: ", sum(is_immune), " immune cells; ", n_agg,
          " cells in candidate immune aggregates (>", cfg$immune_nbr_frac,
          " immune among k=", cfg$spatial_k, " neighbors)")
  p_agg <- ggplot(md[order(md$immune_nbr_frac), ],
                  aes(x_centroid, y_centroid, colour = immune_nbr_frac)) +
    geom_point(size = 0.25) + coord_fixed() +
    scale_colour_viridis_c(option = "inferno", name = "immune\nnbr frac") +
    labs(x = NULL, y = NULL,
         title = paste0(cfg$label, " — immune-neighborhood (candidate aggregates)")) +
    theme_void()
  save_png(p_agg, fig_path(cfg, "spatial_immune_aggregates.png"),
           width = 8, height = 7)
}

# ---- Persist ----------------------------------------------------------------
saveRDS(obj, cfg$obj_02)
message("  wrote ", cfg$obj_02)
write_session_info(cfg, "02")
message("== 02 done ==")
