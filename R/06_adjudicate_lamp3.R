#!/usr/bin/env Rscript
# 06_adjudicate_lamp3.R — targeted adjudication of the LAMP3/CCR7/CD83-high, CD3+
# cluster (BIG c24 / preview c19), which the CD3 gate in 05 over-called as Treg.
# Diagnostic + transparent: runs three tests, records the decisive numbers, and
# sets the label per a fixed decision tree. Never silently overwrites — the
# evidence is written to <label>_lamp3_adjudication.csv and the chosen label
# carries flag_review where the populations are mixed.
#
#   XENIUM_DATASET=preview Rscript R/06_adjudicate_lamp3.R   (run AFTER 05)
#
# The target is found by SIGNATURE (max LAMP3), not a hard ID, so it survives
# cluster-ID shifts. B/plasma and all other clusters are untouched.

suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(ggplot2) })
source(file.path(Sys.getenv("XENIUM_ROOT", unset = getwd()), "R", "config.R"))

cfg <- ensure_dirs(get_config())
set.seed(cfg$seed)
message("== 06_adjudicate_lamp3 :: ", cfg$label, " ==")

obj  <- readRDS(cfg$obj_05)
clus <- obj$seurat_clusters; cl <- levels(clus)
ctn  <- GetAssayData(obj, assay = cfg$assay, layer = "counts")
da   <- GetAssayData(obj, assay = cfg$assay, layer = "data")
present <- function(g) intersect(g, rownames(obj))

# per-cluster percent-positive (counts) and mean (data) for a gene set
pctC  <- function(genes) { genes <- present(genes); if (!length(genes)) return(setNames(rep(NA, length(cl)), cl))
  tapply(Matrix::colSums(ctn[genes, , drop = FALSE]) > 0, clus, mean)[cl] }
meanC <- function(genes) { genes <- present(genes); if (!length(genes)) return(setNames(rep(NA, length(cl)), cl))
  tapply(Matrix::colMeans(da[genes, , drop = FALSE]), clus, mean)[cl] }

# ---- Identify target (max LAMP3) + reference Treg / CD8 clusters by signature -
lamp3_mean <- meanC("LAMP3")
tgt   <- cl[which.max(lamp3_mean)]
foxp3_pct <- pctC("FOXP3")
treg  <- cl[setdiff(order(foxp3_pct, decreasing = TRUE), which(cl == tgt))[1]]
cyto  <- pctC("NKG7") + pctC("GZMK")
cd8   <- cl[setdiff(order(cyto, decreasing = TRUE), which(cl == tgt))[1]]
gi <- function(k) if (k == "global") seq_along(clus) else which(clus == k)
message("  target=c", tgt, " (LAMP3 mean=", round(max(lamp3_mean), 2),
        ")  Treg-ref=c", treg, "  CD8-ref=c", cd8)

# ---- Assemble the evidence table (target / Treg / CD8 / global) --------------
genes <- list(LAMP3="LAMP3", CCR7="CCR7", CD83="CD83", CD3="c(CD3D,CD3E,CD3G,TRAC)",
              FOXP3="FOXP3", IL2RA="IL2RA", CTLA4="CTLA4", SELL="SELL", IL7R="IL7R",
              GZMB="GZMB", PRF1="PRF1", NKG7="NKG7", GZMK="GZMK", CD4="CD4", CD8A="CD8A")
gset <- list(LAMP3="LAMP3", CCR7="CCR7", CD83="CD83",
             CD3=c("CD3D","CD3E","CD3G","TRAC"), FOXP3="FOXP3", IL2RA="IL2RA",
             CTLA4="CTLA4", SELL="SELL", IL7R="IL7R", GZMB="GZMB", PRF1="PRF1",
             NKG7="NKG7", GZMK="GZMK", CD4="CD4", CD8A="CD8A")
ev <- do.call(rbind, lapply(names(gset), function(nm) {
  p <- pctC(gset[[nm]]); m <- meanC(gset[[nm]])
  data.frame(metric = nm,
             target_pct = round(p[tgt],3), treg_pct = round(p[treg],3),
             cd8_pct = round(p[cd8],3),    global_pct = round(mean(Matrix::colSums(ctn[present(gset[[nm]]),,drop=FALSE])>0),3),
             target_mean = round(m[tgt],2), treg_mean = round(m[treg],2),
             cd8_mean = round(m[cd8],2)) }))

# morphology / doublet
fm <- obj$flag_seg_merge; ar <- obj$cell_area; nc <- obj$nCount_Xenium
it <- gi(tgt)
seg_t <- mean(fm[it]); seg_g <- mean(fm)
# co-occurrence + within-cluster LAMP3~CD3 correlation
lamp <- ctn["LAMP3", ] > 0
cd3g <- present(c("CD3D","CD3E","CD3G","TRAC")); cd3any <- Matrix::colSums(ctn[cd3g, , drop = FALSE]) > 0
dbl  <- mean(lamp[it] & cd3any[it])
rho  <- suppressWarnings(cor(da["LAMP3", it], Matrix::colMeans(da[cd3g, , drop = FALSE])[it], method = "spearman"))
# protein (BIG)
prot <- c(CD3E.1 = NA, `HLA-DR` = NA, CD11c = NA, CD68.1 = NA)
if (cfg$has_protein && "Protein" %in% Assays(obj)) {
  pd <- GetAssayData(obj, assay = "Protein", layer = "data")
  for (pn in names(prot)) if (pn %in% rownames(pd)) prot[pn] <- mean(pd[pn, it])
}

# ---- Decision tree ----------------------------------------------------------
not_treg     <- foxp3_pct[tgt] < 0.5 * foxp3_pct[treg]
dc_high      <- lamp3_mean[tgt] > 3 * mean(meanC("LAMP3"))     # LAMP3 well above background
t_present    <- pctC("CD3E")[tgt] > 0.3
cyto_low     <- meanC(c("GZMB","PRF1","NKG7"))[tgt] < 0.5 * meanC(c("GZMB","PRF1","NKG7"))[cd8]
naive_high   <- pctC("CCR7")[tgt] > pctC("CCR7")[cd8] | pctC("IL7R")[tgt] > pctC("IL7R")[cd8]
seg_enriched <- seg_t > 2 * seg_g & median(ar[it]) > median(ar) & median(nc[it]) > median(nc)
prot_both    <- !is.na(prot["HLA-DR"]) && prot["CD3E.1"] > 1 && prot["HLA-DR"] > 1

ev_str <- sprintf(paste0("FOXP3 pct=%.2f vs Treg-ref %.2f (NOT Treg); LAMP3 mean=%.1f (%.0f%%+); ",
  "CCR7/IL7R pct=%.2f/%.2f, cytotoxic(GZMB/PRF1/NKG7) pct=%.2f/%.2f/%.2f; ",
  "seg-merge=%.4f vs %.4f, area %.0f vs %.0f, nCount %.0f vs %.0f (NOT doublet); ",
  "LAMP3~CD3 rho=%.2f (double+=%.0f%%)%s"),
  foxp3_pct[tgt], foxp3_pct[treg], lamp3_mean[tgt], 100*pctC("LAMP3")[tgt],
  pctC("CCR7")[tgt], pctC("IL7R")[tgt], pctC("GZMB")[tgt], pctC("PRF1")[tgt], pctC("NKG7")[tgt],
  seg_t, seg_g, median(ar[it]), median(ar), median(nc[it]), median(nc), rho, 100*dbl,
  if (!is.na(prot["HLA-DR"])) sprintf("; protein CD3=%.2f HLA-DR=%.2f CD11c=%.2f", prot["CD3E.1"], prot["HLA-DR"], prot["CD11c"]) else "")

if (not_treg && dc_high && t_present && rho < 0 && !seg_enriched) {
  label <- "mregDC/CCR7+ T (mixed)"; flag <- TRUE
  basis <- paste0("ADJUDICATED mixed (mregDC + CCR7+ naive/CM T; sub-cluster to resolve): ", ev_str)
} else if (not_treg && naive_high && cyto_low && !dc_high) {
  label <- "Naive/CM T"; flag <- FALSE
  basis <- paste0("ADJUDICATED naive/central-memory T: ", ev_str)
} else if (seg_enriched && dbl > 0.5 && prot_both) {
  label <- "mregDC-T doublet/conjugate"; flag <- TRUE
  basis <- paste0("ADJUDICATED doublet/conjugate: ", ev_str)
} else {
  label <- "mregDC/CCR7+ T (mixed)"; flag <- TRUE
  basis <- paste0("ADJUDICATED ambiguous -> mixed (hypotheses: mregDC vs naive/CM T): ", ev_str)
}
message("  c", tgt, " : ", label, "  (flag_review=", flag, ")")
message("    ", basis)

# ---- Write evidence + update the annotation + object ------------------------
write.csv(ev, tab_path(cfg, "lamp3_adjudication.csv"), row.names = FALSE)
ann <- read.csv(tab_path(cfg, "cluster_annotation.csv"), check.names = FALSE)
r <- which(as.character(ann$cluster) == tgt)
ann$conflict_type[r]       <- "DC_Tmix(adjudicated)"
ann$final_ref_cell_type[r] <- label
ann$resolved[r]            <- FALSE
ann$flag_review[r]         <- flag
ann$resolution_basis[r]    <- basis
if ("protein_concordant" %in% names(ann)) ann$protein_concordant[r] <- NA
write.csv(ann, tab_path(cfg, "cluster_annotation.csv"), row.names = FALSE)

lv <- levels(obj$final_ref_cell_type)
if (!label %in% lv) levels(obj$final_ref_cell_type) <- c(lv, label)
obj$final_ref_cell_type[clus == tgt] <- label
obj$flag_review[clus == tgt] <- flag
obj$final_ref_cell_type <- droplevels(obj$final_ref_cell_type)
ctc <- as.data.frame(table(final_ref_cell_type = obj$final_ref_cell_type))
write.csv(ctc, tab_path(cfg, "final_celltype_counts.csv"), row.names = FALSE)
saveRDS(obj, cfg$obj_05)

# refresh the final figures so they reflect the adjudicated label
save_png(DimPlot(obj, reduction = "umap", group.by = "final_ref_cell_type",
                 label = TRUE, repel = TRUE) +
           ggtitle(paste0(cfg$label, " — final layered cell types (immune + kidney)")),
         fig_path(cfg, "umap_final_celltypes.png"), width = 12, height = 8)
md <- obj[[]]
save_png(ggplot(md, aes(x_centroid, y_centroid, colour = final_ref_cell_type)) +
           geom_point(size = 0.2) + coord_fixed() +
           guides(colour = guide_legend(override.aes = list(size = 2), ncol = 1)) +
           labs(x = NULL, y = NULL, title = paste0(cfg$label, " — final cell types in situ")) +
           theme_void(),
         fig_path(cfg, "spatial_final_celltypes.png"), width = 12, height = 9)
write_session_info(cfg, "06")
message("== 06 done ==")
