#!/usr/bin/env Rscript
# 10_insitutype.R — CosMx-native, contamination-aware cell typing (InSituType), to
# recover the lymphoid compartment that the engine's transfer + co-expression (09)
# could not, and to serve as the platform-matched reference capstone. The per-cell
# neg-probe background model + IF cohorting replace the gate that over-captured
# epithelium in 09. Semi-supervised so unmatched populations still fit de novo.
#
#   XENIUM_SAMPLES="<ctrl>,<sle>" XENIUM_DATASET=cln_cosmx Rscript R/10_insitutype.R  # 2-slide val
#   XENIUM_DATASET=cln_cosmx Rscript R/10_insitutype.R                                 # full
#
# Reference profiles are INDEPENDENT external references (NanoString ioprofiles for
# immune + Human Cell Atlas kidney for epithelium/stroma) — NOT derived from the
# authors' clust on this data, so the benchmark stays non-circular.

suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(InSituType); library(ggplot2) })
source(file.path(Sys.getenv("XENIUM_ROOT", unset = getwd()), "R", "config.R"))
cfg <- ensure_dirs(get_config()); set.seed(cfg$seed)
message("== 10_insitutype :: ", cfg$label, " ==")
obj <- readRDS(cfg$obj_05)

samp <- Sys.getenv("XENIUM_SAMPLES", "")
tag  <- ""
if (nzchar(samp)) { want <- trimws(strsplit(samp, ",")[[1]])
  obj <- subset(obj, cells = colnames(obj)[obj[[cfg$sample_col]][, 1] %in% want]); tag <- "_val2slide"
  message("  validation subset: ", ncol(obj), " cells from ", length(want), " slides") }

# ---- contamination-aware inputs --------------------------------------------
x   <- Matrix::t(GetAssayData(obj, assay = cfg$assay, layer = "counts"))  # cells x genes
neg <- obj$negmean                                                        # per-cell background
IFm <- as.matrix(obj[[c("Mean.PanCK", "Mean.CD45", "Mean.DAPI")]])
cohort <- fastCohorting(IFm)
message("  inputs: ", nrow(x), " cells x ", ncol(x), " genes; ",
        length(unique(cohort)), " IF cohorts; neg median ", round(median(neg), 2))

# ---- INDEPENDENT combined reference profiles (immune ioprofiles + HCA kidney)
data("ioprofiles", package = "InSituType")
hca_path <- file.path(cfg$data_dir, "Childhood onset lupus nephritis analyses", "data", "Kidney_HCA.RData")
ev <- new.env(); load(hca_path, envir = ev); kid <- as.matrix(ev$profile_matrix)
# HCA kidney's own immune columns to drop (immune comes from ioprofiles instead):
immnames <- c("CD4.T.cell","NK.cell","B.cell","CD8.T.cell","Neutrophil","NKT.cell",
              "MNP.d.Tissue.macrophage","Mast.cell","Plasmacytoid.dendritic.cell",
              "MNP.c.dendritic.cell","MNP.a.classical.monocyte.derived","MNP.b.non.classical.monocyte.derived")
sg <- intersect(rownames(ioprofiles), rownames(kid))
prof <- cbind(kid[sg, setdiff(colnames(kid), immnames)],
              ioprofiles[sg, setdiff(colnames(ioprofiles), c("fibroblast", "endothelial"))])
prof <- sweep(prof, 2, colSums(prof), "/") * 1e6
prof <- prof[intersect(rownames(prof), colnames(x)), ]                    # to panel genes
message("  reference profiles: ", nrow(prof), " panel genes x ", ncol(prof),
        " types (", paste(head(colnames(prof), 4), collapse = ","), "...)")

# ---- InSituType: semi-supervised (de-novo clusters for unmatched pops) -------
# n_clusts via env: "1:5" semi-supervised (default), "0" = SUPERVISED-ONLY — the
# on-target diagnostic for rare lymphoid (known types already in the reference;
# isolates panel/profile resolving power from any de-novo / param confound).
nclusts <- eval(parse(text = Sys.getenv("XENIUM_NCLUSTS", "1:5")))
if (identical(nclusts, 0)) tag <- paste0(tag, "_sup")
message("  insitutype n_clusts = ", paste(range(nclusts), collapse = ":"),
        if (identical(nclusts, 0)) " (SUPERVISED-only diagnostic)" else " (semi-supervised)")
res <- insitutype(x = x, neg = neg, assay_type = "rna", cohort = cohort,
                  reference_profiles = prof, n_clusts = nclusts,
                  update_reference_profiles = TRUE, n_starts = 2)
lab <- res$clust[colnames(obj)]
col <- if (grepl("_sup", tag)) "insitutype_sup" else "insitutype"
obj[[col]] <- lab
if (!nzchar(samp)) saveRDS(obj, cfg$obj_05)                 # never persist a validation SUBSET
write.csv(as.data.frame(sort(table(lab), decreasing = TRUE)),
          tab_path(cfg, paste0("insitutype_counts", tag, ".csv")), row.names = FALSE)
message("  InSituType -> ", length(unique(lab)), " types; ",
        sum(grepl("^[a-z]$", lab)), " cells in de-novo clusters")

# ---- benchmark: recall AND precision vs authors (lymphoid-focused) ----------
to_imm <- function(x) { x <- as.character(x); o <- rep(NA_character_, length(x)); m <- function(p) grepl(p, x, ignore.case = TRUE)
  o[m("B-cell|B cell")] <- "B"; o[m("plasmabl|plasma")] <- "Plasma"
  o[m("T CD8")] <- "T_CD8"; o[m("T CD4|Treg")] <- "T_CD4"; o[m("^NK$|natural killer|NK.cell")] <- "NK"
  o[m("macrophage")] <- "Macrophage"; o[m("monocyte")] <- "Monocyte"
  o[m("mDC|dendritic")] <- "DC"; o[m("pDC|plasmacytoid")] <- "pDC"; o[m("mast")] <- "Mast"; o[m("neutrophil")] <- "Neutrophil"; o }
ai <- to_imm(obj$author_celltype); oi <- to_imm(lab)
lev <- c("B","Plasma","T_CD4","T_CD8","NK","DC","pDC","Mast","Macrophage","Monocyte","Neutrophil")
# before = the 09 two-stage label, if present
o09 <- if ("immune_subtype" %in% colnames(obj[[]])) obj$immune_subtype else rep(NA, ncol(obj))
b09 <- to_imm(o09)
bench <- data.frame(immune_type = lev,
  n_author = sapply(lev, function(L) sum(ai == L, na.rm = TRUE)),
  recall_09  = sapply(lev, function(L){a<-which(ai==L); if(!length(a))NA else round(mean(b09[a]==L,na.rm=TRUE),3)}),
  recall_IST = sapply(lev, function(L){a<-which(ai==L); if(!length(a))NA else round(mean(oi[a]==L,na.rm=TRUE),3)}),
  precision_IST = sapply(lev, function(L){p<-which(oi==L); if(!length(p))NA else round(mean(ai[p]==L,na.rm=TRUE),3)}),
  n_ours_IST = sapply(lev, function(L) sum(oi == L, na.rm = TRUE)))
write.csv(bench, tab_path(cfg, paste0("insitutype_benchmark", tag, ".csv")), row.names = FALSE)
ov_r <- mean(oi[!is.na(ai)] == ai[!is.na(ai)], na.rm = TRUE)
message(sprintf("  immune recall vs authors (full denom %d cells): InSituType overall %.1f%% (09 was ~0 for lymphoid)",
                sum(!is.na(ai)), 100*ov_r))
print(bench, row.names = FALSE)
# full confusion (immune categories)
cm <- table(author = factor(ai, lev), insitutype = factor(oi, lev))
write.csv(as.data.frame.matrix(cm), tab_path(cfg, paste0("insitutype_confusion", tag, ".csv")))
write_session_info(cfg, paste0("10", tag)); message("== 10 done ==")
