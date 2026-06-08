#!/usr/bin/env Rscript
# 16_groupC_confusions.R — confusion matrices for whitepaper Group C that need labels
# living only in the _05 object (not the h5ad): (C2) reference-transfer (label c =
# final_ref_cell_type) coarse confusion vs author, where epithelial->immune errors
# concentrate. Uses R/08's coarse map. Author x reference-transfer, full data.
#
#   Rscript R/16_groupC_confusions.R

suppressPackageStartupMessages({ library(Seurat) })
ROOT <- Sys.getenv("XENIUM_ROOT", unset = getwd())
OBJ  <- file.path(ROOT, "outputs", "objects"); TAB <- file.path(ROOT, "outputs", "tables")
obj <- readRDS(file.path(OBJ, "cln_cosmx_05_refann_kidney.rds")); md <- obj[[]]

coarse <- function(x) { x <- as.character(x); o <- rep("Other", length(x)); m <- function(p) grepl(p, x, ignore.case = TRUE)
  o[m("tubule|PCT|proximal|distal|connecting|intercalat|principal|ascending|loop.of.henle|thin.limb|parietal|pelvic|urotheli|macula|papillary|collecting|epithelial")] <- "Epithelial"
  o[m("podocyte")] <- "Podocyte"
  o[m("endotheli|vasa.recta|capillary|lymphatic")] <- "Endothelial"
  o[m("fibroblast|myofibroblast|mesangial|pericyte|smooth|stroma|mural|renin")] <- "Stroma"
  o[m("macrophage|monocyte|dendritic|mregdc|myeloid|neutrophil|mast|basophil|phagocyte")] <- "Myeloid"
  o[m("B.?cell|naive B|memory B|plasmabl|plasma|MZB")] <- "B_Plasma"
  o[m("natural killer|NK cell|^NK$")] <- "NK"
  o[m("T cell|T CD|CD4|CD8|treg|regulatory|MAIT|helper|CCR7")] <- "T"
  o }
lv <- c("Epithelial","Podocyte","Endothelial","Stroma","Myeloid","B_Plasma","T","NK","Other")
ac <- coarse(md$author_celltype)
cc <- coarse(md$final_ref_cell_type)     # (c) reference transfer
cm <- table(author = factor(ac, lv), reference_transfer = factor(cc, lv))
write.csv(as.data.frame.matrix(cm), file.path(TAB, "cln_cosmx_reftransfer_coarse_confusion.csv"))
message("wrote cln_cosmx_reftransfer_coarse_confusion.csv")
# quick epithelial->immune readout
epi <- which(ac == "Epithelial"); n_epi <- length(epi)
to_imm <- sum(cc[epi] %in% c("Myeloid","B_Plasma","T","NK"))
message(sprintf("  author-Epithelial cells: %d; reference-transfer -> immune: %d (%.1f%%)",
                n_epi, to_imm, 100*to_imm/n_epi))
message("== 16 done ==")
