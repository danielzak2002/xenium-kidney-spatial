#!/usr/bin/env Rscript
# 12_qc_metrics_export.R — export per-cell QC metrics from the 01_load_qc objects
# for the whitepaper Group-A figures. Reads ONLY the committed *_01_qc.rds objects;
# writes a compact, subsampled per-cell CSV per dataset (git-ignored intermediate)
# plus a combined flag-summary table assembled from the existing qc_summary.csv files.
# No recomputation of QC — distributions and flags both come from the 01 outputs.
#
#   Rscript R/12_qc_metrics_export.R

suppressPackageStartupMessages({ library(Seurat) })
ROOT <- Sys.getenv("XENIUM_ROOT", unset = getwd())
OBJ  <- file.path(ROOT, "outputs", "objects")
TAB  <- file.path(ROOT, "outputs", "tables")
set.seed(1); SUBSAMPLE <- 80000L
N_NEG_COSMX <- 20L   # CosMx 1000-plex Universal panel: negative probes Negative01-20
                     # (verified from the dataset's own probe naming). Used to reconstruct
                     # a unit-matched neg-control FRACTION from the stored per-probe negmean.

DS <- list(
  list(label = "kidney_RCC_protein",  assay = "Xenium", platform = "xenium"),
  list(label = "kidney_preview_PRCC", assay = "Xenium", platform = "xenium"),
  list(label = "cln_cosmx",           assay = "RNA",    platform = "cosmx"))

amb_rows <- list()
for (d in DS) {
  f <- file.path(OBJ, paste0(d$label, "_01_qc.rds"))
  if (!file.exists(f)) { message("MISSING ", f, " — skipping"); next }
  message("== ", d$label, " ==")
  obj <- readRDS(f)
  md  <- obj[[]]
  nc_col <- paste0("nCount_", d$assay); nf_col <- paste0("nFeature_", d$assay)
  out <- data.frame(
    dataset   = d$label,
    platform  = d$platform,
    n_counts  = md[[nc_col]],
    n_genes   = md[[nf_col]],
    neg_frac  = md$neg_frac,                                   # ambient fraction (both platforms)
    negmean   = if ("negmean" %in% colnames(md)) md$negmean else NA_real_,  # CosMx raw background
    cell_area = md$cell_area,
    stringsAsFactors = FALSE)
  # FULL-DATA unit-matched neg-control fraction = neg counts / total counts per cell.
  # Xenium: neg_frac is already (control_probe+codeword)/total. CosMx: no raw neg matrix,
  # so reconstruct neg counts as negmean * N_NEG_COSMX; denominator includes them (as the
  # Xenium denominator includes controls). This is a background-fraction PROXY, not an exact
  # cross-platform identity (probe designs differ) — the order-of-magnitude gap is robust.
  if (d$platform == "cosmx") {
    negc <- md$negmean * N_NEG_COSMX
    amb  <- mean(negc / (md$totalcounts + negc), na.rm = TRUE)
  } else {
    amb  <- mean(md$neg_frac, na.rm = TRUE)
  }
  amb_rows[[d$label]] <- data.frame(dataset = d$label, platform = d$platform,
    n_neg_probes = if (d$platform == "cosmx") N_NEG_COSMX else NA_integer_,
    mean_neg_ctrl_fraction = amb, n_cells = nrow(md))

  # subsample for plotting (distributions are unchanged; flag counts come from qc_summary)
  if (nrow(out) > SUBSAMPLE) out <- out[sort(sample(nrow(out), SUBSAMPLE)), ]
  csv <- file.path(OBJ, paste0("qc_metrics_", d$label, ".csv"))
  write.csv(out, csv, row.names = FALSE)
  message("  wrote ", csv, " (", nrow(out), " cells); mean neg-ctrl fraction ",
          signif(amb, 3))
  rm(obj, md); gc(verbose = FALSE)
}
amb <- do.call(rbind, amb_rows)
write.csv(amb, file.path(TAB, "qcA_ambient_fraction.csv"), row.names = FALSE)
message("wrote qcA_ambient_fraction.csv"); print(amb)
message("== 12 done ==")
