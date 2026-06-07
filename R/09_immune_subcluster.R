#!/usr/bin/env Rscript
# 09_immune_subcluster.R — STAGE 2 of a two-stage annotation for CosMx. The global
# res-0.8 pass (01-05) resolves the abundant epithelial/stromal compartments but the
# ~3% immune minority is absorbed into mixed clusters (diagnosed: 0 clusters >50%
# immune). Here we GATE the immune compartment, re-Harmonize + re-cluster it, type
# sub-clusters by CO-EXPRESSION (never count>0, which is ambient-corrupted on CosMx),
# drop epithelial carryover, fold back, and re-score immune per-type recall on the
# FULL author-immune denominator (gate misses count as misses).
#
#   XENIUM_DATASET=cln_cosmx Rscript R/09_immune_subcluster.R   (after 05)

suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(ggplot2) })
source(file.path(Sys.getenv("XENIUM_ROOT", unset = getwd()), "R", "config.R"))
options(future.globals.maxSize = 8 * 1024^3)
cfg <- ensure_dirs(get_config()); set.seed(cfg$seed)
message("== 09_immune_subcluster :: ", cfg$label, " ==")
obj <- readRDS(cfg$obj_05)
DefaultAssay(obj) <- cfg$assay
present <- function(g) intersect(g, rownames(obj))

# ---- lineage marker sets (>=2 present required; co-expression, not count>0) ---
sets <- list(
  PanImmune  = c("PTPRC","CD3D","CD3E","CD8A","CD68","CD163","MS4A1","CD79A","MZB1","NKG7","GNLY","LYZ","ITGAX","CD14","IL7R","CCL5"),
  B          = c("MS4A1","CD79A","CD79B","CD19","BANK1"),
  Plasma     = c("MZB1","TNFRSF17","DERL3","SDC1","XBP1","PRDM1"),
  T_CD4      = c("CD3D","CD3E","CD4","IL7R","FOXP3","CTLA4"),
  T_CD8      = c("CD3D","CD3E","CD8A","CD8B","GZMK"),
  NK         = c("NKG7","GNLY","KLRD1","NCAM1","KLRF1"),
  Macrophage = c("CD68","CD163","C1QA","C1QB","MRC1"),
  Monocyte   = c("LYZ","FCN1","VCAN","S100A8","S100A9","CD14"),
  DC         = c("ITGAX","CLEC9A","CD1C","FCER1A","LAMP3"),
  pDC        = c("IL3RA","PLD4","TCL1A","GZMB","IRF7"),
  Mast       = c("CPA3","MS4A2","KIT","TPSAB1","GATA2"),
  Epithelial = c("EPCAM","CDH1","KRT8","KRT18","KRT19","PAX8","SLC34A1","CUBN","UMOD","DEFB1"))
sets <- lapply(sets, present); sets <- sets[lengths(sets) >= 2]
message("  lineage sets present (>=2 genes): ", paste(names(sets), collapse = ", "))

# ---- module scores (co-expression) for gating + typing ----------------------
obj <- AddModuleScore(obj, features = sets, name = "ms_", nbin = 20, ctrl = 20, seed = cfg$seed)
msnames <- paste0("ms_", seq_along(sets)); names(msnames) <- names(sets)
imm_score <- obj[[msnames[["PanImmune"]]]][, 1]
clus <- obj$seurat_clusters

# ---- GATE (generous; recall-favoured): immune-enriched clusters INTERSECT score
cl_imm <- tapply(imm_score, clus, mean)
enriched <- names(cl_imm)[cl_imm > quantile(cl_imm, 0.55)]      # top ~45% clusters
cd45 <- obj[[cfg$if_immune_col]][, 1]; cd45_hi <- cd45 > quantile(cd45, 0.75)
thr <- quantile(imm_score, 0.50)                                # generous within-cluster floor
gate <- (as.character(clus) %in% enriched) & (imm_score > thr | cd45_hi)
message("  immune-enriched clusters: ", paste(enriched, collapse = ","),
        " | gated cells: ", sum(gate), " (", round(100*mean(gate), 1), "%)")

# ---- re-process the gated subset: re-Harmonize by slide, re-cluster ----------
sub <- subset(obj, cells = colnames(obj)[gate])
VariableFeatures(sub) <- rownames(sub)
sub <- ScaleData(sub, features = rownames(sub), verbose = FALSE)
sub <- RunPCA(sub, npcs = 30, verbose = FALSE, features = rownames(sub))
sub <- harmony::RunHarmony(sub, group.by.vars = cfg$sample_col, reduction.use = "pca",
                           dims.use = 1:30, reduction.save = "harmony", verbose = FALSE)
sub <- FindNeighbors(sub, dims = 1:30, reduction = "harmony", verbose = FALSE)
smemb <- run_leiden(sub, paste0(DefaultAssay(sub), "_snn"), 1.2, cfg$seed,
                    min_size = 30L, emb = Embeddings(sub, "harmony")[, 1:30])
sub$sub <- smemb
message("  immune subset ", ncol(sub), " cells -> ", nlevels(smemb), " sub-clusters")

# ---- type sub-clusters by CO-EXPRESSION (mean module score per sub-cluster) --
lin <- setdiff(names(sets), "PanImmune")
sc <- sapply(lin, function(L) tapply(sub[[msnames[[L]]]][, 1], smemb, mean))  # subcluster x lineage
panck <- tapply(sub[[cfg$if_epithelial_col]][, 1], smemb, mean); panck <- panck / max(panck)
assign_lin <- colnames(sc)[max.col(sc, ties.method = "first")]
# epithelial carryover: argmax is Epithelial, or PanCK-IF-high with weak immune
carry <- assign_lin == "Epithelial" | (panck > 0.6 & apply(sc[, setdiff(lin,"Epithelial"),drop=FALSE], 1, max) < 0.05)
assign_lin[carry] <- "Epithelial_carryover"
names(assign_lin) <- levels(smemb)
subtab <- data.frame(sub_cluster = levels(smemb), n = as.integer(table(smemb)),
                     lineage = assign_lin, panck_if = round(as.numeric(panck), 2),
                     top_score = round(apply(sc, 1, max), 3))
write.csv(subtab, tab_path(cfg, "immune_subclusters.csv"), row.names = FALSE)
message("  sub-cluster types: ", paste(sprintf("%s=%d", names(table(assign_lin)), table(assign_lin)), collapse = " "))

# ---- fold back: per-cell immune lineage (drop carryover) --------------------
cell_lin <- assign_lin[as.character(smemb)]; names(cell_lin) <- colnames(sub)
obj$immune_subtype <- NA_character_
keep <- cell_lin != "Epithelial_carryover"
obj$immune_subtype[names(cell_lin)[keep]] <- cell_lin[keep]
saveRDS(obj, cfg$obj_05)

# ---- immune per-type recall vs authors, FULL denominator --------------------
au <- obj$author_celltype
a2lin <- function(x) { x <- as.character(x); o <- rep(NA_character_, length(x)); m <- function(p) grepl(p, x, ignore.case = TRUE)
  o[m("B-cell|B cell")] <- "B"; o[m("plasmabl|plasma")] <- "Plasma"
  o[m("T CD8")] <- "T_CD8"; o[m("T CD4|Treg")] <- "T_CD4"; o[m("^NK$|natural killer")] <- "NK"
  o[m("macrophage")] <- "Macrophage"; o[m("monocyte")] <- "Monocyte"
  o[m("mDC")] <- "DC"; o[m("pDC")] <- "pDC"; o[m("mast")] <- "Mast"; o }
al <- a2lin(au); ol <- obj$immune_subtype
lev <- c("B","Plasma","T_CD4","T_CD8","NK","Macrophage","Monocyte","DC","pDC","Mast")
rec <- data.frame(immune_type = lev,
  n_author = sapply(lev, function(L) sum(al == L, na.rm = TRUE)),
  recall_after = sapply(lev, function(L) { a <- which(al == L); if (!length(a)) NA else round(mean(!is.na(ol[a]) & ol[a] == L), 3) }),
  n_ours = sapply(lev, function(L) sum(ol == L, na.rm = TRUE)))
write.csv(rec, tab_path(cfg, "immune_recall_twostage.csv"), row.names = FALSE)
overall <- mean(ol[which(!is.na(al))] == al[which(!is.na(al))], na.rm = TRUE)
message(sprintf("  TWO-STAGE immune recall (full denominator, %d author-immune cells): overall %.1f%%",
                sum(!is.na(al)), 100*overall))
print(rec, row.names = FALSE)
write_session_info(cfg, "09"); message("== 09 done ==")
