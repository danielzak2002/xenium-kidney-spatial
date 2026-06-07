#!/usr/bin/env Rscript
# 11_export_anndata.R — export ONE .h5ad per dataset for Phase B (squidpy).
# Platforms are NOT merged (per-section analysis; mixing um/mm coords would be a mess).
#
#   XENIUM_DATASET=kidney_preview_PRCC Rscript R/11_export_anndata.R
#   (label keys: kidney_preview_PRCC | kidney_RCC_protein | cln_cosmx  -> pass the cfg dataset key)
#
# X = raw counts (RNA/Xenium assay), layers['lognorm'] = log-norm. obs carries the
# label as-is PLUS phase_b_label (DC -> "myeloid/DC"; CosMx T -> low-confidence),
# sample/condition/platform, IF channels (CosMx), and author label (CosMx) for ref.
# obsm['spatial'] = per-cell x/y; uns['spatial_units'] records mm (CosMx) vs um (Xenium).
# Round-trip is VERIFIED before declaring done.

Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")  # R + conda both ship libomp (macOS)
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(reticulate) })
source(file.path(Sys.getenv("XENIUM_ROOT", unset = getwd()), "R", "config.R"))
# spatial conda env python (has anndata/scanpy/squidpy). Override via SPATIAL_PYTHON.
spy <- Sys.getenv("SPATIAL_PYTHON", "")
if (!nzchar(spy) || !file.exists(spy)) {
  cand <- c("/opt/homebrew/Caskroom/miniforge/base/envs/spatial/bin/python",
            file.path(path.expand("~"), "miniforge3/envs/spatial/bin/python"))
  spy <- cand[file.exists(cand)][1]
}
use_python(spy, required = TRUE)
ad <- import("anndata"); sp <- import("scipy.sparse"); np <- import("numpy"); pd <- import("pandas")

cfg <- ensure_dirs(get_config())
message("== 11_export_anndata :: ", cfg$label, " (", cfg$platform, ") ==")
obj <- readRDS(cfg$obj_05); assay <- cfg$assay

# genes x cells dgCMatrix -> scipy CSR (cells x genes). reticulate auto-converts a
# dgCMatrix to scipy.sparse; transpose first so rows = cells, then .tocsr().
to_csr <- function(m) r_to_py(as(Matrix::t(m), "CsparseMatrix"))$tocsr()

md <- obj[[]]
celltype <- if (cfg$platform == "cosmx") as.character(obj$insitutype) else as.character(obj$final_ref_cell_type)
# phase_b_label: fold GENERIC DC -> myeloid/DC, but KEEP mregDC distinct (validated
# regulatory-DC state, a TLS player). CosMx T -> low-confidence (typing was unreliable).
pbl <- celltype
pbl[grepl("dendritic|^mDC$|^pDC$|^DC$", pbl, ignore.case = TRUE) &
    !grepl("mregDC", pbl, ignore.case = TRUE)] <- "myeloid/DC"
if (cfg$platform == "cosmx") pbl[grepl("T CD|Treg|^T_CD", pbl, ignore.case = TRUE)] <- "T cell (low-conf)"

obs <- data.frame(
  cell_type     = celltype,
  phase_b_label = pbl,
  sample        = if (cfg$multi_sample) as.character(obj[[cfg$sample_col]][, 1]) else cfg$label,
  condition     = if (!is.na(cfg$condition_col)) as.character(obj[[cfg$condition_col]][, 1]) else cfg$tissue_desc,
  platform      = cfg$platform,
  stringsAsFactors = FALSE)
if (cfg$platform == "cosmx") {
  obs$author_celltype <- as.character(obj$author_celltype)
  for (ch in c("Mean.PanCK","Max.PanCK","Mean.CD45","Max.CD45","Mean.DAPI","Max.DAPI"))
    if (ch %in% colnames(md)) obs[[ch]] <- md[[ch]]
}
rownames(obs) <- colnames(obj)

ft <- tryCatch(as.character(obj[[assay]]@meta.data$feature_type), error = function(e) NULL)
var <- data.frame(gene = rownames(obj), row.names = rownames(obj), stringsAsFactors = FALSE)
if (!is.null(ft) && length(ft) == nrow(obj)) var$feature_type <- ft

adata <- ad$AnnData(X = to_csr(GetAssayData(obj, assay = assay, layer = "counts")),
                    obs = r_to_py(obs), var = r_to_py(var))
adata$layers["lognorm"] <- to_csr(GetAssayData(obj, assay = assay, layer = "data"))
adata$obsm["spatial"]   <- np$array(as.matrix(md[, c("x_centroid", "y_centroid")]))
adata$uns["spatial_units"] <- cfg$coord_units
adata$uns["platform"]      <- cfg$platform
adata$uns["dataset"]       <- cfg$label

h5 <- file.path(cfg$dir_obj, paste0(cfg$label, ".h5ad"))
adata$write_h5ad(h5)
message("  wrote ", h5, " (", adata$n_obs, " x ", adata$n_vars, ")")

# ---- VERIFY round-trip ------------------------------------------------------
rl <- ad$read_h5ad(h5)
chk <- list(
  n_cells   = rl$n_obs == ncol(obj),
  n_genes   = rl$n_vars == nrow(obj),
  labels    = all(as.character(rl$obs[["cell_type"]]) == celltype),
  has_sample = "sample" %in% names(rl$obs), has_condition = "condition" %in% names(rl$obs),
  units     = rl$uns[["spatial_units"]] == cfg$coord_units)
set.seed(1); i <- sample(ncol(obj), 5)
coords_match <- all(abs(rl$obsm[["spatial"]][i, ] - as.matrix(md[i, c("x_centroid","y_centroid")])) < 1e-4)
counts_match <- abs(sum(rl$X) - sum(GetAssayData(obj, assay = assay, layer = "counts"))) < 1
chk$coords <- coords_match; chk$counts_total <- counts_match
message("  VERIFY: ", paste(sprintf("%s=%s", names(chk), unlist(chk)), collapse = " "))
if (!all(unlist(chk))) stop("round-trip verification FAILED")

summ <- data.frame(dataset = cfg$label, platform = cfg$platform, units = cfg$coord_units,
  n_cells = ncol(obj), n_genes = nrow(obj),
  n_samples = if (cfg$multi_sample) length(unique(obs$sample)) else 1L,
  cell_type_col = if (cfg$platform == "cosmx") "insitutype" else "final_ref_cell_type",
  roundtrip_ok = all(unlist(chk)))
write.csv(summ, tab_path(cfg, "anndata_export_summary.csv"), row.names = FALSE)
message("== 11 done ==")
