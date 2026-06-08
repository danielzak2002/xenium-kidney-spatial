#!/usr/bin/env Rscript
# 14_immune_benchmark_unified.R — single-source-of-truth immune benchmark for the cLN
# cohort: recall AND precision for BOTH typing passes (two-stage R/09 immune_subtype, and
# InSituType), computed with ONE consistent immune mapping from the committed _05 object.
# Supersedes the recall-only / split tables for whitepaper Group C/B4.
#
#   Rscript R/14_immune_benchmark_unified.R

suppressPackageStartupMessages({ library(Seurat) })
ROOT <- Sys.getenv("XENIUM_ROOT", unset = getwd())
OBJ  <- file.path(ROOT, "outputs", "objects"); TAB <- file.path(ROOT, "outputs", "tables")
obj <- readRDS(file.path(OBJ, "cln_cosmx_05_refann_kidney.rds"))
md  <- obj[[]]

lev <- c("B","Plasma","T_CD4","T_CD8","NK","DC","pDC","Mast","Macrophage","Monocyte","Neutrophil")
# to_imm: EXACT R/10 mapping (the one that produced the committed InSituType benchmark) —
# requires "B-cell"/"B cell" so the de-novo cluster letter "b" never leaks into B, etc.
to_imm <- function(x) { x <- as.character(x); o <- rep(NA_character_, length(x))
  m <- function(p) grepl(p, x, ignore.case = TRUE)
  o[m("B-cell|B cell")] <- "B"; o[m("plasmabl|plasma")] <- "Plasma"
  o[m("T CD8")] <- "T_CD8"; o[m("T CD4|Treg")] <- "T_CD4"; o[m("^NK$|natural killer|NK.cell")] <- "NK"
  o[m("macrophage")] <- "Macrophage"; o[m("monocyte")] <- "Monocyte"
  o[m("mDC|dendritic")] <- "DC"; o[m("pDC|plasmacytoid")] <- "pDC"
  o[m("mast")] <- "Mast"; o[m("neutrophil")] <- "Neutrophil"; o }
# two-stage (R/09) immune_subtype labels are ALREADY canonical lineage names -> identity map.
ts_map <- function(x) { x <- as.character(x); ifelse(x %in% lev, x, NA_character_) }

au  <- to_imm(md$author_celltype)
ts  <- ts_map(md$immune_subtype)   # two-stage (R/09)
ist <- to_imm(md$insitutype)       # InSituType (R/10)

rp <- function(pred) {
  data.frame(
    n_author  = sapply(lev, function(L) sum(au == L, na.rm = TRUE)),
    n_pred    = sapply(lev, function(L) sum(pred == L, na.rm = TRUE)),
    recall    = sapply(lev, function(L){a<-which(au==L); if(!length(a))NA else round(mean(!is.na(pred[a]) & pred[a]==L),3)}),
    precision = sapply(lev, function(L){p<-which(pred==L); if(!length(p))NA else round(mean(!is.na(au[p]) & au[p]==L),3)}))
}
B <- data.frame(immune_type = lev)
ts_rp <- rp(ts); ist_rp <- rp(ist)
B$n_author      <- ts_rp$n_author
B$recall_twostage    <- ts_rp$recall;    B$precision_twostage    <- ts_rp$precision; B$n_pred_twostage <- ts_rp$n_pred
B$recall_insitutype  <- ist_rp$recall;   B$precision_insitutype  <- ist_rp$precision; B$n_pred_insitutype <- ist_rp$n_pred
write.csv(B, file.path(TAB, "cln_cosmx_immune_benchmark_unified.csv"), row.names = FALSE)
cat("== UNIFIED immune benchmark (single source of truth) ==\n"); print(B, row.names = FALSE)

# ---- VERIFY: where do author Monocyte / NK cells land under InSituType? ------
rawist <- as.character(md$insitutype)
for (L in c("Monocyte","NK")) {
  a <- which(au == L)
  cat(sprintf("\nauthor %s (n=%d) -> InSituType immune-mapped:\n", L, length(a)))
  print(round(sort(prop.table(table(factor(ist[a], lev), useNA = "always")), decreasing = TRUE)[1:5], 3))
  cat("  top raw InSituType labels:\n"); print(head(sort(table(rawist[a]), decreasing = TRUE), 5))
}
cat("\n== 14 done ==\n")
