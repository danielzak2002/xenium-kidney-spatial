#!/usr/bin/env Rscript
# 05_kidney_reference.R — Phase A step 5: KIDNEY/EPITHELIAL reference layer,
# layered on top of the immune labels from 04. Uses the KPMP/HuBMAP Azimuth
# human-kidney reference (Zenodo 10694842, CC BY 4.0; cached git-ignored) at
# annotation.l1 (16 nephron / stroma / endothelial / immune classes).
#
#   XENIUM_DATASET=preview Rscript R/05_kidney_reference.R
#
# NOTE on method: this reference's ref.Rds ships only the supervised-PCA reduction
# (the gene matrix is emptied), so it cannot be run correlation-style via SingleR.
# We therefore use Azimuth's own mapping (SCTransform + projection onto the reference
# SPCA + label transfer) and gate by its prediction score. The immune layer (04)
# stays SingleR/Monaco. The kidney label is adopted ONLY for non-immune, non-(tumour/
# LowQ) clusters; immune clusters keep their Monaco label; Tumor_RCC / LowQ_MTRNR2L /
# Mast are never overridden. Output: a single combined final_ref_cell_type.

suppressPackageStartupMessages({
  library(Seurat); library(Azimuth); library(ggplot2); library(dplyr)
})
source(file.path(Sys.getenv("XENIUM_ROOT", unset = getwd()), "R", "config.R"))

# Azimuth's label transfer routes the (already in-memory) reference + anchor weights
# through future. Pin a sequential plan (no worker copies) and lift the global-export
# guard — at 465k cells the closure exceeds any fixed cap, but nothing is duplicated.
future::plan("sequential")
options(future.globals.maxSize = 1e12)

cfg <- ensure_dirs(get_config())
set.seed(cfg$seed)
message("== 05_kidney_reference (Azimuth KPMP) :: ", cfg$label, " ==")

obj  <- readRDS(cfg$obj_04)          # carries immune ref_cell_type + marker cell_type
DefaultAssay(obj) <- cfg$assay
clus <- obj$seurat_clusters

# ---- Azimuth mapping onto the KPMP kidney reference -------------------------
if (!file.exists(file.path(cfg$ref_kidney_dir, "ref.Rds")))
  stop("Kidney reference not found in ", cfg$ref_kidney_dir,
       " — download per DATA.md (Zenodo 10694842: ref.Rds + idx.annoy).")
pred_col  <- paste0("predicted.", cfg$ref_kidney_level)
score_col <- paste0(pred_col, ".score")
azi <- RunAzimuth(obj, reference = cfg$ref_kidney_dir,
                  annotation.levels = cfg$ref_kidney_level,
                  assay = cfg$assay, verbose = FALSE)
am <- azi[[]]
obj$kidney_l1            <- am[colnames(obj), pred_col]
obj$kidney_score         <- am[colnames(obj), score_col]
obj$kidney_mapping_score <- am[colnames(obj), "mapping.score"]
obj$kidney_lowconf       <- obj$kidney_score < cfg$ref_kidney_score_min
obj$kidney_pruned        <- ifelse(obj$kidney_lowconf, NA, obj$kidney_l1)
message("  Azimuth done on ", ncol(obj), " cells; ",
        length(unique(na.omit(obj$kidney_l1))), " ", cfg$ref_kidney_level,
        " classes; low-confidence (score<", cfg$ref_kidney_score_min, "): ",
        sum(obj$kidney_lowconf))

# ---- Per-cluster consensus (kidney) + base layered label --------------------
imm_tab <- read.csv(tab_path(cfg, "cluster_reference_labels.csv"),   # 04 immune layer
                    check.names = FALSE, stringsAsFactors = FALSE)
rownames(imm_tab) <- as.character(imm_tab$cluster)
consensus <- function(x) { x <- x[!is.na(x)]; if (!length(x)) return(c(NA, 0))
  tb <- sort(table(x), decreasing = TRUE); c(names(tb)[1], tb[1] / length(x)) }
cl_levels   <- levels(clus)
marker_type <- vapply(cl_levels, function(k)
  as.character(obj$cell_type[which(clus == k)[1]]), character(1))
immune_lab  <- vapply(cl_levels, function(k)              # from 04 (immune layer)
  as.character(obj$ref_cell_type[which(clus == k)[1]]), character(1))
kcons <- t(vapply(cl_levels, function(k) consensus(obj$kidney_pruned[clus == k]), character(2)))
kid_lab <- kcons[, 1]; kid_agree <- as.numeric(kcons[, 2])
non_immune <- !(marker_type %in% IMMUNE_TYPES)
# Base label: immune subtype (04) for immune clusters; nephron/stroma/endo for
# non-immune clusters with adequate kidney consensus (coarse "Immune" not adoptable).
adopt_kidney <- non_immune & !(marker_type %in% cfg$ref_keep_marker) &
                kid_agree >= cfg$ref_consensus_min & !is.na(kid_lab) & kid_lab != "Immune"
base <- immune_lab
base[adopt_kidney] <- kid_lab[adopt_kidney]

# CosMx: the marker layer is unreliable, so identify immune clusters by Azimuth's
# "Immune" class and label them with the Monaco fine consensus (04, per-cell SingleR);
# non-immune clusters take the Azimuth nephron class. CD45 IF / PTPRC are reported as
# orthogonal validators (cols below). Azimuth+Monaco IS the resolution -> no resolver.
if (isTRUE(cfg$immune_via_reference)) {
  monaco_fine <- vapply(cl_levels, function(k) { x <- obj$singler_pruned[clus == k]
    x <- x[!is.na(x)]; if (!length(x)) NA_character_ else names(sort(table(x), decreasing = TRUE))[1] },
    character(1))
  azi_immune <- !is.na(kid_lab) & kid_lab == "Immune" & kid_agree >= cfg$ref_consensus_min
  azi_epi    <- !is.na(kid_lab) & kid_lab != "Immune" & kid_agree >= cfg$ref_consensus_min
  base <- ifelse(marker_type %in% cfg$ref_keep_marker, marker_type,
          ifelse(azi_immune & !is.na(monaco_fine), monaco_fine,
          ifelse(azi_epi, kid_lab, marker_type)))
  immune_lab <- base
  adopt_kidney <- azi_epi & !(marker_type %in% cfg$ref_keep_marker)
  message("  immune-via-reference: ", sum(azi_immune), " Azimuth-immune cluster(s) -> Monaco fine")
}

# ---- Flag the conflicts to resolve ------------------------------------------
singler_main <- imm_tab[cl_levels, "singler_main"]
concordant   <- as.logical(imm_tab[cl_levels, "concordant"])
tish <- function(x) x %in% c("CD4+ T cells", "CD8+ T cells", "T cells")
imm_conflict <- (marker_type %in% IMMUNE_TYPES) & !(marker_type %in% cfg$ref_keep_marker) &
                !is.na(concordant) & !concordant
kid_conflict <- non_immune & !(marker_type %in% cfg$ref_keep_marker) &
                kid_agree >= cfg$ref_consensus_min & !is.na(kid_lab) & kid_lab == "Immune"
ctype <- rep("none", length(cl_levels))
ctype[imm_conflict & marker_type %in% c("DC", "pDC") & tish(singler_main)] <- "DC_Treg"
ctype[imm_conflict & marker_type == "NK_cell"        & tish(singler_main)] <- "NK_CD8"
ctype[imm_conflict & ctype == "none"]                                       <- "immune_other"
ctype[kid_conflict]                                                         <- "Epi_Immune"
if (isTRUE(cfg$immune_via_reference)) ctype <- rep("none", length(cl_levels))  # gating is the resolution

# ---- Lineage-gate evidence: per-cluster percent-positive (counts) ------------
ctn <- GetAssayData(obj, assay = cfg$assay, layer = "counts")
g   <- lapply(lineage_gate_markers(), intersect, rownames(obj))
pctpos <- function(genes) { genes <- intersect(genes, rownames(ctn))
  if (!length(genes)) return(rep(NA_real_, length(cl_levels)))
  as.numeric(tapply(Matrix::colSums(ctn[genes, , drop = FALSE]) > 0, clus, mean)[cl_levels]) }
cd3 <- pctpos(g$T_lineage); ptprc <- pctpos(g$PanImmune)
lamp3 <- pctpos(intersect("LAMP3", rownames(obj)))
foxp3 <- pctpos(intersect("FOXP3", rownames(obj)))
nkpos <- pctpos(cfg$nk_gate_markers)
# Epithelial signal: RNA epithelial markers if on-panel; else fall back to the PanCK
# IF channel (per-cluster mean, scaled to [0,1] across clusters) — CosMx lacks
# FXYD2/CDH16. CD45 IF is carried as the orthogonal pan-immune (PTPRC) validator.
ifq <- function(col) { if (is.na(col) || !(col %in% colnames(obj[[]]))) return(rep(NA_real_, length(cl_levels)))
  v <- as.numeric(tapply(obj[[col]][, 1], clus, mean)[cl_levels]); v / max(v, na.rm = TRUE) }
epi <- pctpos(g$Epithelial); if (all(is.na(epi))) epi <- ifq(cfg$if_epithelial_col)
cd45_if <- ifq(cfg$if_immune_col); panck_if <- ifq(cfg$if_epithelial_col)
HI <- cfg$gate_pos_high; NEG <- cfg$gate_pos_neg
f2 <- function(x) ifelse(is.na(x), "NA", sprintf("%.2f", x))

# ---- Resolver: CD3 = primary T vs non-T split; PTPRC = immune vs epithelial.
# Never force a call — ambiguous evidence keeps the existing label, flag_review=TRUE.
final <- base; basis <- rep("", length(cl_levels))
resolved <- rep(FALSE, length(cl_levels)); flag_review <- rep(FALSE, length(cl_levels))
for (i in which(ctype != "none")) {
  if (ctype[i] == "DC_Treg") {
    if (!is.na(cd3[i]) && cd3[i] >= HI) { final[i] <- "Treg"; resolved[i] <- TRUE
      basis[i] <- sprintf("CD3-high (pct=%s), FOXP3 pct=%s -> Treg", f2(cd3[i]), f2(foxp3[i]))
    } else if (!is.na(cd3[i]) && cd3[i] < NEG && !is.na(lamp3[i]) && lamp3[i] >= HI) {
      final[i] <- "mregDC"; resolved[i] <- TRUE
      basis[i] <- sprintf("CD3-neg (pct=%s), LAMP3-high (pct=%s) -> mregDC", f2(cd3[i]), f2(lamp3[i]))
    } else { flag_review[i] <- TRUE
      basis[i] <- sprintf("CD3 ambiguous (pct=%s), LAMP3 pct=%s -> keep %s (review)", f2(cd3[i]), f2(lamp3[i]), base[i]) }
  } else if (ctype[i] == "NK_CD8") {
    if (!is.na(cd3[i]) && cd3[i] >= HI) { final[i] <- "CD8_T"; resolved[i] <- TRUE
      basis[i] <- sprintf("CD3-high (pct=%s) -> CD8_T", f2(cd3[i]))
    } else if (!is.na(cd3[i]) && cd3[i] < NEG && !is.na(nkpos[i]) && nkpos[i] >= HI) {
      final[i] <- "NK"; resolved[i] <- TRUE
      basis[i] <- sprintf("CD3-neg (pct=%s), NK-high (NCAM1/KLRD1 pct=%s) -> NK", f2(cd3[i]), f2(nkpos[i]))
    } else { flag_review[i] <- TRUE
      basis[i] <- sprintf("CD3 ambiguous (pct=%s), NK pct=%s -> keep %s (review)", f2(cd3[i]), f2(nkpos[i]), base[i]) }
  } else if (ctype[i] == "Epi_Immune") {
    if (!is.na(ptprc[i]) && ptprc[i] < NEG && !is.na(epi[i]) && epi[i] >= HI) {
      final[i] <- marker_type[i]; resolved[i] <- TRUE       # override the reference "Immune"
      basis[i] <- sprintf("PTPRC-neg (pct=%s), epithelial-high (pct=%s) -> keep %s", f2(ptprc[i]), f2(epi[i]), marker_type[i])
    } else if (!is.na(ptprc[i]) && ptprc[i] >= HI) {
      final[i] <- "Immune infiltrate/doublet"; flag_review[i] <- TRUE
      basis[i] <- sprintf("PTPRC-high (pct=%s) -> immune infiltrate/doublet (review)", f2(ptprc[i]))
    } else { flag_review[i] <- TRUE
      basis[i] <- sprintf("ambiguous (PTPRC pct=%s, epi pct=%s) -> keep %s (review)", f2(ptprc[i]), f2(epi[i]), marker_type[i]) }
  } else { flag_review[i] <- TRUE       # immune_other: flagged but no specific gate
    basis[i] <- sprintf("immune conflict (marker=%s vs %s), no gate -> keep %s (review)",
                        marker_type[i], singler_main[i], base[i]) }
}

# ---- Protein cross-check (BIG only): CD3 protein for CD3-gated clusters ------
prot_cd3 <- rep(NA_real_, length(cl_levels)); prot_conc <- rep(NA, length(cl_levels))
if (cfg$has_protein && "Protein" %in% Assays(obj)) {
  pd <- GetAssayData(obj, assay = "Protein", layer = "data")          # CLR-normalized
  # Seurat appends ".1" to protein features that collide with RNA gene names (CD3E.1).
  pn <- intersect(c("CD3E.1", "CD3E", "CD3"), rownames(pd))[1]
  if (!is.na(pn)) {
    prot_cd3 <- as.numeric(tapply(pd[pn, ], clus, mean)[cl_levels])
    med <- median(prot_cd3, na.rm = TRUE)
    prot_conc[final %in% c("Treg", "CD8_T")] <- prot_cd3[final %in% c("Treg", "CD8_T")] >  med
    prot_conc[final %in% c("mregDC", "NK")]  <- prot_cd3[final %in% c("mregDC", "NK")]  <= med
  }
}

# ---- Write the resolved cluster annotation ----------------------------------
ann <- data.frame(
  cluster = cl_levels, n_cells = as.integer(table(clus)), marker_type = marker_type,
  immune_label = immune_lab, kidney_l1 = kid_lab, kidney_agreement = round(kid_agree, 3),
  conflict_type = ctype,
  cd3_pct = round(cd3, 3), ptprc_pct = round(ptprc, 3), lamp3_pct = round(lamp3, 3),
  foxp3_pct = round(foxp3, 3), nk_pct = round(nkpos, 3), epi_pct = round(epi, 3),
  cd45_if = round(cd45_if, 3), panck_if = round(panck_if, 3),
  protein_CD3_mean = round(prot_cd3, 3), protein_concordant = prot_conc,
  resolution_basis = basis, resolved = resolved, flag_review = flag_review,
  final_ref_cell_type = final, row.names = NULL)
write.csv(ann, tab_path(cfg, "cluster_annotation.csv"), row.names = FALSE)
message("  conflicts: ", sum(ctype != "none"), " (",
        paste(names(table(ctype[ctype != "none"])), table(ctype[ctype != "none"]),
              sep = ":", collapse = " "), "); resolved: ", sum(resolved),
        "; flag_review: ", sum(flag_review))
for (i in which(ctype != "none"))
  message("    c", cl_levels[i], " [", ctype[i], "] ", base[i], " -> ", final[i], " : ", basis[i])

cl2final <- setNames(final, cl_levels)
obj$final_ref_cell_type <- factor(unname(cl2final[as.character(clus)]))
obj$flag_review <- unname(setNames(flag_review, cl_levels)[as.character(clus)])
ctc <- as.data.frame(table(final_ref_cell_type = obj$final_ref_cell_type))
write.csv(ctc, tab_path(cfg, "final_celltype_counts.csv"), row.names = FALSE)

# ---- Figures ----------------------------------------------------------------
save_png(DimPlot(obj, reduction = "umap", group.by = "final_ref_cell_type",
                 label = TRUE, repel = TRUE) +
           ggtitle(paste0(cfg$label, " — final layered cell types (immune + kidney)")),
         fig_path(cfg, "umap_final_celltypes.png"), width = 12, height = 8)

md <- obj[[]]
p_sp <- ggplot(md, aes(x_centroid, y_centroid, colour = final_ref_cell_type)) +
  geom_point(size = 0.2) + coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size = 2), ncol = 1)) +
  labs(x = NULL, y = NULL,
       title = paste0(cfg$label, " — final cell types in situ")) +
  theme_void()
save_png(p_sp, fig_path(cfg, "spatial_final_celltypes.png"), width = 12, height = 9)

# ---- Persist ----------------------------------------------------------------
saveRDS(obj, cfg$obj_05)
message("  wrote ", cfg$obj_05)
write_session_info(cfg, "05")
message("== 05 done ==")
