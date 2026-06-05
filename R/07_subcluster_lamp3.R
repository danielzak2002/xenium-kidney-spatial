#!/usr/bin/env Rscript
# 07_subcluster_lamp3.R — split the under-resolved mregDC + CCR7+ T cluster
# (BIG c24 / preview c19; flagged "mixed" by 06) into its two real populations.
# LAMP3 vs CD3 are anti-correlated within it (06), so it sub-clusters cleanly.
#
#   XENIUM_DATASET=preview Rscript R/07_subcluster_lamp3.R   (run AFTER 06)
#   Run order: 01 -> 02 -> 04 -> 05 -> 06 -> 07.
#
# Target is found by signature (max mean LAMP3), so it survives cluster-ID shifts.
# Split markers: LAMP3 (mregDC) vs CD3D/CD3E (T); CCR7 is shared -> never used to
# split. Falls back to a per-cell LAMP3-vs-CD3 gate if re-clustering is degenerate.
# B/plasma and every other cluster are untouched.

suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(ggplot2) })
source(file.path(Sys.getenv("XENIUM_ROOT", unset = getwd()), "R", "config.R"))

cfg <- ensure_dirs(get_config())
set.seed(cfg$seed)
message("== 07_subcluster_lamp3 :: ", cfg$label, " ==")

obj  <- readRDS(cfg$obj_05)          # post-06 (target labelled "mregDC/CCR7+ T (mixed)")
clus <- obj$seurat_clusters
pres <- function(g) intersect(g, rownames(obj))
da   <- GetAssayData(obj, assay = cfg$assay, layer = "data")

# ---- Identify target by signature (max mean LAMP3) --------------------------
lamp3_mean <- tapply(da["LAMP3", ], clus, mean)
tgt   <- names(lamp3_mean)[which.max(lamp3_mean)]
cells <- colnames(obj)[clus == tgt]
message("  target=c", tgt, " (", length(cells), " cells, LAMP3 mean=",
        round(max(lamp3_mean), 2), ")")

# ---- Re-process the target cells in isolation -------------------------------
sub <- subset(obj, cells = cells)
DefaultAssay(sub) <- cfg$assay
VariableFeatures(sub) <- rownames(sub)
sub <- ScaleData(sub, features = rownames(sub), verbose = FALSE)
sub <- RunPCA(sub, npcs = min(30L, ncol(sub) - 1L), verbose = FALSE,
              features = rownames(sub))
dims <- 1:min(20L, ncol(Embeddings(sub, "pca")))
sub <- FindNeighbors(sub, dims = dims, verbose = FALSE)
gname <- paste0(DefaultAssay(sub), "_snn")
memb <- NULL
for (res in c(0.8, 1.0, 1.5, 2.0)) {
  m <- run_leiden(sub, gname, res, cfg$seed, min_size = 10L,
                  emb = Embeddings(sub, "pca")[, dims])
  if (nlevels(m) >= 2) { memb <- m; message("  sub-leiden res=", res, " -> ",
                          nlevels(m), " sub-clusters: ",
                          paste(as.integer(table(m)), collapse = ", ")); break }
}

# ---- Decisive split: LAMP3 (mregDC) vs CD3D/CD3E (T); CCR7 NOT used ----------
sd  <- GetAssayData(sub, assay = cfg$assay, layer = "scale.data")
cd3 <- pres(c("CD3D", "CD3E"))
zL  <- sd["LAMP3", ]; zC <- colMeans(sd[cd3, , drop = FALSE])
method <- "subcluster"; assign <- setNames(rep(NA_character_, ncol(sub)), colnames(sub))
clean <- FALSE
if (!is.null(memb)) {
  subL <- tapply(zL, memb, mean); subC <- tapply(zC, memb, mean)
  subtype <- ifelse(subL > subC, "mregDC", "CCR7+ T (naive/CM)")
  clean <- length(unique(subtype)) >= 2          # both populations present
  if (clean) assign[] <- subtype[as.character(memb)]
}
if (!clean) {                                    # fallback: per-cell LAMP3 vs CD3 gate
  method <- "per-cell gate"
  ctn <- GetAssayData(sub, assay = cfg$assay, layer = "counts")
  lp <- ctn["LAMP3", ] > 0; cp <- Matrix::colSums(ctn[cd3, , drop = FALSE]) > 0
  assign[lp & !cp] <- "mregDC"; assign[cp & !lp] <- "CCR7+ T (naive/CM)"
  assign[(lp & cp) | (!lp & !cp)] <- NA          # double-pos / double-neg -> review
  message("  re-clustering not decisive -> per-cell gate")
}
n_amb <- sum(is.na(assign))
message("  split (", method, "): mregDC=", sum(assign == "mregDC", na.rm = TRUE),
        "  CCR7+ T=", sum(assign == "CCR7+ T (naive/CM)", na.rm = TRUE),
        "  ambiguous=", n_amb)

# ---- Sanity vs markers / protein --------------------------------------------
grp <- function(lab) names(assign)[which(assign == lab)]
mm <- function(genes, cc) if (length(cc)) round(mean(colMeans(da[pres(genes), cc, drop = FALSE])), 2) else NA
message("  RNA check  mregDC: LAMP3=", mm("LAMP3", grp("mregDC")), " CD3=",
        mm(cd3, grp("mregDC")), " | T: CD3=", mm(cd3, grp("CCR7+ T (naive/CM)")),
        " IL7R=", mm("IL7R", grp("CCR7+ T (naive/CM)")), " cytotoxic=",
        mm(c("GZMB","PRF1","NKG7"), grp("CCR7+ T (naive/CM)")))
if (cfg$has_protein && "Protein" %in% Assays(obj)) {
  pd <- GetAssayData(obj, assay = "Protein", layer = "data")
  pm <- function(pn, cc) if (pn %in% rownames(pd) && length(cc)) round(mean(pd[pn, cc]), 2) else NA
  message("  protein check  mregDC: HLA-DR=", pm("HLA-DR", grp("mregDC")), " CD11c=",
          pm("CD11c", grp("mregDC")), " CD3=", pm("CD3E.1", grp("mregDC")),
          " | T: CD3=", pm("CD3E.1", grp("CCR7+ T (naive/CM)")), " HLA-DR=",
          pm("HLA-DR", grp("CCR7+ T (naive/CM)")))
}

# ---- Fold back into the object ----------------------------------------------
lab_amb <- "mregDC/CCR7+ T (mixed)"
finlab <- ifelse(is.na(assign), lab_amb, assign)
lv <- union(levels(obj$final_ref_cell_type), c("mregDC", "CCR7+ T (naive/CM)", lab_amb))
obj$final_ref_cell_type <- factor(as.character(obj$final_ref_cell_type), levels = lv)
obj$final_ref_cell_type[names(finlab)] <- finlab
obj$final_ref_cell_type <- droplevels(obj$final_ref_cell_type)
obj$flag_review[names(assign)] <- is.na(assign)          # confident split -> clear flag
obj$subpop <- NA_character_; obj$subpop[names(assign)] <- finlab
saveRDS(obj, cfg$obj_05)

# ---- Update cluster_annotation.csv: add the two sub-population rows ----------
ann <- read.csv(tab_path(cfg, "cluster_annotation.csv"), check.names = FALSE,
                stringsAsFactors = FALSE)
pr <- which(as.character(ann$cluster) == tgt)
ann$resolution_basis[pr] <- paste0("split in 07 (", method, ") -> c", tgt,
  ".mregDC / c", tgt, ".T; see sub-population rows")
ann$flag_review[pr] <- n_amb > 0
mk_row <- function(suffix, lab, n, why) { r <- ann[pr, ]; r$cluster <- paste0(tgt, suffix)
  r$n_cells <- n; r$marker_type <- "DC"; r$conflict_type <- "subpop(07)"
  r$resolved <- TRUE; r$flag_review <- FALSE; r$final_ref_cell_type <- lab
  r$resolution_basis <- why; r }
new <- rbind(
  mk_row(".mregDC", "mregDC", sum(assign == "mregDC", na.rm = TRUE),
         sprintf("LAMP3-dominant sub-population (%s split)", method)),
  mk_row(".T", "CCR7+ T (naive/CM)", sum(assign == "CCR7+ T (naive/CM)", na.rm = TRUE),
         sprintf("CD3-dominant, cytotoxic-low CCR7+ T sub-population (%s split)", method)))
if (n_amb > 0) new <- rbind(new, mk_row(".ambiguous", lab_amb, n_amb,
  "LAMP3/CD3 double-pos or double-neg -> flag_review")[, names(new)])
if (n_amb > 0) new$flag_review[nrow(new)] <- TRUE
ann <- rbind(ann, new)
write.csv(ann, tab_path(cfg, "cluster_annotation.csv"), row.names = FALSE)
ctc <- as.data.frame(table(final_ref_cell_type = obj$final_ref_cell_type))
write.csv(ctc, tab_path(cfg, "final_celltype_counts.csv"), row.names = FALSE)

# ---- Refresh the final figures ----------------------------------------------
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
write_session_info(cfg, "07")
message("== 07 done ==")
