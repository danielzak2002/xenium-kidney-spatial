# config.R — shared configuration + helpers for the Xenium kidney pipeline.
#
# Sourced by R/01_load_qc.R and R/02_cluster_annotate.R. Nothing dataset-specific
# is hard-coded in the analysis scripts: everything flows through get_config().
#
# Select the dataset with the env var XENIUM_DATASET ("preview" | "big"),
# default "preview" (the 56,510-cell PRCC section = small->big validation target).
# Run scripts from the repo root (or set XENIUM_ROOT).

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

XENIUM_ROOT <- Sys.getenv("XENIUM_ROOT", unset = getwd())

# ---- Dataset registry ------------------------------------------------------
# Cell counts / modality confirmed from metrics_summary.csv + cell_feature_matrix.h5
# during the panel-assessment step. The preview bundle on disk is the cancer
# (PRCC) section, NOT the paired non-diseased section.
.DATASETS <- list(
  preview = list(
    label       = "kidney_preview_PRCC",
    data_dir    = "kidney_10x_preview/data",
    has_protein = FALSE,
    n_cells     = 56510L,
    tissue_desc = "papillary renal cell carcinoma (PRCC) section"
  ),
  big = list(
    label       = "kidney_RCC_protein",
    data_dir    = "kidney_10x/data",
    has_protein = TRUE,   # 27-plex protein assay; see scale-up notes in 02
    n_cells     = 465534L,
    tissue_desc = "renal cell carcinoma (Stage III, T3a N1 MX) section"
  )
)

get_config <- function(dataset = Sys.getenv("XENIUM_DATASET", unset = "preview")) {
  if (!dataset %in% names(.DATASETS))
    stop("Unknown dataset '", dataset, "'. Choose one of: ",
         paste(names(.DATASETS), collapse = ", "))
  ds <- .DATASETS[[dataset]]

  data_dir <- file.path(XENIUM_ROOT, ds$data_dir)
  out      <- file.path(XENIUM_ROOT, "outputs")

  cfg <- list(
    dataset      = dataset,
    label        = ds$label,
    has_protein  = ds$has_protein,
    n_cells_meta = ds$n_cells,
    tissue_desc  = ds$tissue_desc,

    # paths
    data_dir     = data_dir,
    cells_parquet = file.path(data_dir, "cells.parquet"),
    gene_panel   = file.path(data_dir, "gene_panel.json"),
    dir_fig      = file.path(out, "figures"),
    dir_tab      = file.path(out, "tables"),
    dir_obj      = file.path(out, "objects"),     # git-ignored
    dir_log      = file.path(out, "logs"),
    obj_01       = file.path(out, "objects", paste0(ds$label, "_01_qc.rds")),
    obj_02       = file.path(out, "objects", paste0(ds$label, "_02_annotated.rds")),

    # reproducibility
    seed         = 1234L,

    # ---- QC thresholds (Xenium-specific; flag, do not droplet-style filter) --
    # Only hard removal is zero-RNA-count cells (cannot normalize/cluster).
    min_counts_keep        = 1L,
    drop_flagged           = FALSE,   # keep flagged cells; record flags only
    neg_frac_flag          = 0.02,    # per-cell neg-control (probe+codeword) fraction
    blank_frac_flag        = 0.02,    # per-cell blank/unassigned-codeword fraction
    seg_merge_count_quantile = 0.99,  # candidate segmentation merges: high counts ...
    seg_merge_area_quantile  = 0.99,  # ... AND large cell area
    lowq_mtrnr2l_frac        = 0.5,   # per-cell: counts dominated by MTRNR2L pseudogene
    lowq_complexity          = 0.2,   # per-cell: nFeature/nCount below this = low complexity

    # ---- Normalization (operator choice: LogNormalize; SCT kept as switch) ---
    norm_method        = "LogNormalize",  # or "SCT"
    norm_scale_factor  = 1e4,
    # Targeted 377-gene panel: use ALL genes as features (no HVG selection).
    use_all_features   = TRUE,
    n_variable_features = 2000L,           # only used if use_all_features = FALSE

    # ---- Dim reduction / clustering -----------------------------------------
    npcs_compute      = 50L,
    dims              = 1:30,
    cluster_method    = "leiden",   # igraph native Leiden on the SNN graph
    leiden_objective  = "modularity",
    cluster_resolution = 0.8,
    min_cluster_size  = 10L,         # merge sub-threshold clusters (drops singletons)

    # ---- Marker DE ----------------------------------------------------------
    de_only_pos   = TRUE,
    de_min_pct    = 0.10,
    de_logfc      = 0.25,
    de_top_n      = 10L,

    # ---- Spatial smoke-test (lightweight; real stats are Phase B / squidpy) -
    run_spatial_smoke = TRUE,
    spatial_k         = 20L,        # kNN on centroids
    immune_nbr_frac   = 0.60        # candidate immune-aggregate threshold
  )
  cfg$assay <- "Xenium"
  # BIG carries 15 'Deprecated Codeword' features that LoadXenium 5.5.0 cannot
  # map -> load it via the manual h5 path proactively (see load_xenium_lean).
  cfg$force_manual_load <- identical(dataset, "big")
  cfg
}

# ---- Lean Xenium loader with Deprecated-Codeword fallback ------------------
# Builds a Seurat object from cell_feature_matrix.h5 directly: keeps the gene
# panel as the "Xenium" assay, the 27-plex protein as a SEPARATE "Protein"
# assay, and DROPS the control / blank / deprecated-codeword feature classes
# (their per-cell rates are recovered from cells.parquet for QC). No FOV/polygon
# objects are built — all spatial plotting uses x/y_centroid from meta.data, so
# the lean path needs only the matrix + centroids. Used proactively for BIG
# (LoadXenium's slot.map has no 'Deprecated Codeword' entry) and as a fallback.
build_xenium_from_h5 <- function(cfg) {
  h5 <- file.path(cfg$data_dir, "cell_feature_matrix.h5")
  if (!file.exists(h5)) stop("cell_feature_matrix.h5 not found in ", cfg$data_dir)
  mat <- Read10X_h5(h5)
  if (!is.list(mat)) mat <- list(`Gene Expression` = mat)
  if (!"Gene Expression" %in% names(mat))
    stop("No 'Gene Expression' features in ", h5)
  obj <- CreateSeuratObject(counts = mat[["Gene Expression"]], assay = cfg$assay)
  if ("Protein Expression" %in% names(mat))
    obj[["Protein"]] <- CreateAssayObject(counts = mat[["Protein Expression"]])
  obj
}

load_xenium_lean <- function(cfg) {
  manual <- function(reason) {
    message("  manual h5 loader (", reason, ")"); build_xenium_from_h5(cfg)
  }
  if (isTRUE(cfg$force_manual_load)) return(manual("dataset=big / proactive"))
  tryCatch(
    LoadXenium(cfg$data_dir, assay = cfg$assay, cell.centroids = TRUE,
               molecule.coordinates = FALSE),
    error = function(e) manual(paste("LoadXenium failed:", conditionMessage(e))))
}

ensure_dirs <- function(cfg) {
  for (d in c(cfg$dir_fig, cfg$dir_tab, cfg$dir_obj, cfg$dir_log))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  invisible(cfg)
}

# Namespaced output paths so preview/big artifacts never collide.
fig_path <- function(cfg, name) file.path(cfg$dir_fig, paste0(cfg$label, "_", name))
tab_path <- function(cfg, name) file.path(cfg$dir_tab, paste0(cfg$label, "_", name))
log_path <- function(cfg, name) file.path(cfg$dir_log, paste0(cfg$label, "_", name))

save_png <- function(plot, path, width = 8, height = 6, dpi = 200) {
  ggplot2::ggsave(path, plot = plot, width = width, height = height,
                  dpi = dpi, bg = "white")
  message("  wrote ", path)
}

# ---- Curated kidney + immune marker sets -----------------------------------
# Superset of candidates; each script intersects with rownames() so the SAME
# list works for the 377-gene preview and the 405-gene custom panel.
# B-cell / plasma emphasis per project goal. NOTE: CD79B is absent from both
# panels; B-cell gate rests on MS4A1 / CD79A / CD19 / BANK1.
marker_sets <- function() list(
  Proximal_tubule   = c("LRP2","CUBN","SLC34A1","GATM","MIOX","SLC5A12","ALDOB","SLC13A3"),
  Loop_of_Henle_TAL = c("UMOD","SLC12A1","CASR","CLDN16"),
  Distal_tubule     = c("SLC12A3","TRPM6","CALB1"),
  Principal_CD      = c("AQP2","AQP3","SCNN1G","GATA3"),
  Intercalated      = c("ATP6V0D2","SLC4A1","SLC26A4","FOXI1"),
  Podocyte          = c("NPHS1","NPHS2","PTPRO","PODXL"),
  Endothelial       = c("PECAM1","FLT1","EMCN","KDR","ACKR1","VWF"),
  Stroma_mural      = c("ACTA2","PDGFRB","PDGFRA","COL1A1","DCN","NOTCH3"),
  T_cell            = c("CD3E","CD3D","CD8A","CD4","IL7R","CCL5"),
  NK_cell           = c("NKG7","GNLY","KLRD1","KLRF1"),
  B_cell            = c("MS4A1","CD79A","CD19","BANK1","CD74"),
  Plasma            = c("MZB1","DERL3","TNFRSF17","SDC1","PRDM1"),
  Myeloid           = c("CD68","CD163","LYZ","ITGAM","C1QA","C1QB","CD14"),
  DC                = c("LAMP3","CLEC9A","ITGAX","CD1C"),
  pDC               = c("LILRA4","IL3RA","CLEC4C","PLD4","GZMB","TCL1A","IRF7","JCHAIN"),
  Mast              = c("TPSAB1","TPSB2","CPA3","MS4A2","KIT","GATA2","CMA1","HPGDS"),
  Neutrophil        = c("S100A12","S100A8","S100A9","FCGR3B","CSF3R","MCEMP1","AQP9","CXCR2"),
  Proliferating     = c("MKI67","TOP2A","PCNA"),
  Tumor_RCC         = c("CA9","NDUFA4L2","VEGFA","CD70","ANGPTL4")
)

# Immune cell-type sets used by the spatial smoke-test (which clusters are immune).
IMMUNE_TYPES <- c("T_cell","NK_cell","B_cell","Plasma","Myeloid","DC",
                  "pDC","Mast","Neutrophil")

# ---- Native Leiden on a Seurat SNN graph (no python/reticulate) -------------
# Seurat's algorithm=4 routes through leiden(py) via reticulate, which is a
# fragile dependency on Apple silicon. igraph >= 1.2 ships a C implementation;
# running it directly on the SNN graph keeps clustering pure-R and deterministic.
run_leiden <- function(obj, graph_name, resolution, seed,
                       objective = "modularity", min_size = 10L, emb = NULL) {
  stopifnot(graph_name %in% names(obj@graphs))
  A <- obj[[graph_name]]                       # Seurat Graph (sparse, symmetric)
  g <- igraph::graph_from_adjacency_matrix(
    methods::as(A, "CsparseMatrix"), mode = "undirected",
    weighted = TRUE, diag = FALSE)
  set.seed(seed)
  part <- igraph::cluster_leiden(
    g, objective_function = objective, weights = igraph::E(g)$weight,
    resolution = resolution, n_iterations = -1L)
  memb <- as.integer(igraph::membership(part))

  # Merge sub-threshold clusters (incl. singleton low-count cells) into the
  # majority cluster among their SNN neighbors; iterate until stable.
  if (min_size > 1L) repeat {
    small <- as.integer(names(which(table(memb) < min_size)))
    if (!length(small)) break
    changed <- FALSE
    for (v in which(memb %in% small)) {
      nb <- as.integer(igraph::neighbors(g, v))
      nb <- nb[!(memb[nb] %in% small)]
      if (length(nb)) {
        memb[v] <- as.integer(names(which.max(table(memb[nb])))); changed <- TRUE
      }
    }
    if (!changed) break
  }

  # Residual small clusters = SNN-isolated cells; place each by nearest cell in
  # PCA space (guarantees every cell gets a real cluster, no singletons left).
  if (min_size > 1L && !is.null(emb)) {
    small <- as.integer(names(which(table(memb) < min_size)))
    if (length(small)) {
      si <- which(memb %in% small); bi <- which(!(memb %in% small))
      nn <- RANN::nn2(emb[bi, , drop = FALSE], emb[si, , drop = FALSE], k = 1)$nn.idx[, 1]
      memb[si] <- memb[bi][nn]
    }
  }

  # relabel clusters by descending size, starting at 1
  ord  <- names(sort(table(memb), decreasing = TRUE))
  relab <- setNames(seq_along(ord), ord)
  out <- factor(relab[as.character(memb)], levels = seq_along(ord))
  names(out) <- rownames(A)   # cell barcodes -> metadata assigns by name
  out
}

write_session_info <- function(cfg, tag) {
  p <- log_path(cfg, paste0("sessionInfo_", tag, ".txt"))
  writeLines(capture.output(sessionInfo()), p)
  message("  wrote ", p)
}
