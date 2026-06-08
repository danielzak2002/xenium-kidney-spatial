#!/usr/bin/env Rscript
# 15_insitutype_posterior.R — recover the InSituType per-cell posterior for whitepaper
# slot B3 WITHOUT desync risk. Mirrors R/10's exact pipeline (same seed, same inputs,
# same reference) and ONLY persists res$prob if the regenerated labels match the
# committed obj$insitutype EXACTLY for every cell. If they differ, writes nothing and
# the placeholder stands. Never overwrites obj_05.
#
#   XENIUM_DATASET=cln_cosmx Rscript R/15_insitutype_posterior.R

suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(InSituType) })
source(file.path(Sys.getenv("XENIUM_ROOT", unset = getwd()), "R", "config.R"))
Sys.setenv(XENIUM_DATASET = "cln_cosmx")
cfg <- ensure_dirs(get_config()); set.seed(cfg$seed)   # SAME seed as R/10
message("== 15_insitutype_posterior :: ", cfg$label, " ==")
obj <- readRDS(cfg$obj_05)
committed <- as.character(obj$insitutype)               # the labels we must reproduce

# ---- inputs IDENTICAL to R/10 ----------------------------------------------
x   <- Matrix::t(GetAssayData(obj, assay = cfg$assay, layer = "counts"))
neg <- obj$negmean
IFm <- as.matrix(obj[[c("Mean.PanCK", "Mean.CD45", "Mean.DAPI")]])
cohort <- fastCohorting(IFm)
data("ioprofiles", package = "InSituType")
hca_path <- file.path(cfg$data_dir, "Childhood onset lupus nephritis analyses", "data", "Kidney_HCA.RData")
ev <- new.env(); load(hca_path, envir = ev); kid <- as.matrix(ev$profile_matrix)
immnames <- c("CD4.T.cell","NK.cell","B.cell","CD8.T.cell","Neutrophil","NKT.cell",
              "MNP.d.Tissue.macrophage","Mast.cell","Plasmacytoid.dendritic.cell",
              "MNP.c.dendritic.cell","MNP.a.classical.monocyte.derived","MNP.b.non.classical.monocyte.derived")
sg <- intersect(rownames(ioprofiles), rownames(kid))
prof <- cbind(kid[sg, setdiff(colnames(kid), immnames)],
              ioprofiles[sg, setdiff(colnames(ioprofiles), c("fibroblast", "endothelial"))])
prof <- sweep(prof, 2, colSums(prof), "/") * 1e6
prof <- prof[intersect(rownames(prof), colnames(x)), ]
nclusts <- eval(parse(text = Sys.getenv("XENIUM_NCLUSTS", "1:5")))

res <- insitutype(x = x, neg = neg, assay_type = "rna", cohort = cohort,
                  reference_profiles = prof, n_clusts = nclusts,
                  update_reference_profiles = TRUE, n_starts = 2)
lab <- as.character(res$clust[colnames(obj)])

# ---- EXACT-match gate -------------------------------------------------------
n_match <- sum(lab == committed); n_tot <- length(committed)
message(sprintf("  label match: %d / %d (%.4f%%)", n_match, n_tot, 100 * n_match / n_tot))
if (n_match == n_tot) {
  prob <- res$prob[colnames(obj)]
  out <- data.frame(cell_id = colnames(obj), insitutype = lab, posterior = as.numeric(prob))
  write.csv(out, file.path(cfg$dir_obj, "wp_insitutype_posterior.csv"), row.names = FALSE)
  message("  EXACT MATCH -> wrote wp_insitutype_posterior.csv (zero desync). Posterior summary:")
  print(summary(out$posterior))
} else {
  message("  MISMATCH -> NOT reproducible to the bit; posterior NOT saved, placeholder stands.")
}
message("== 15 done ==")
