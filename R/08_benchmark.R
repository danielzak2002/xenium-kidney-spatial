#!/usr/bin/env Rscript
# 08_benchmark.R — score the engine against the authors' independent annotation
# (CosMx cLN: annot$clust, 35 types). Benchmarks THREE labelings vs the authors:
#   (a) cluster structure   — clusters -> majority author type (backbone validation)
#   (b) marker + lineage-gate label (obj$cell_type; 02 + 06/07) — panel-intrinsic, PRIMARY
#   (c) reference-transfer label (obj$final_ref_cell_type; 04 SingleR/Monaco + 05 Azimuth)
# The (b) vs (c) contrast is the finding. Also: immune per-type recall (Phase B), and an
# IF cross-check that corrects reference "Immune" calls which are CD45-neg + PanCK-pos
# (= real epithelium) — orthogonal evidence that clustering is sound, transfer is weak.
#
#   XENIUM_DATASET=cln_cosmx Rscript R/08_benchmark.R   (run after 07)

suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
source(file.path(Sys.getenv("XENIUM_ROOT", unset = getwd()), "R", "config.R"))
cfg <- ensure_dirs(get_config())
message("== 08_benchmark :: ", cfg$label, " ==")
obj <- readRDS(cfg$obj_05)
if (!"author_celltype" %in% colnames(obj[[]])) { message("  no author_celltype; skip"); quit(save="no") }
md <- obj[[]]

# ---- coarse lineage map (both vocabularies) ---------------------------------
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
# finer immune categories for per-type recall (Phase B)
immcat <- function(x) { x <- as.character(x); o <- rep(NA_character_, length(x)); m <- function(p) grepl(p, x, ignore.case = TRUE)
  o[m("plasmabl|plasma|TNFRSF17|MZB")] <- "Plasma"
  o[m("B.?cell|naive B|memory B|exhausted B")] <- "B"
  o[m("CD8")] <- "T_CD8"; o[m("CD4|treg|regulatory|helper")] <- "T_CD4"
  o[m("natural killer|NK cell|^NK$")] <- "NK"
  o[m("plasmacytoid|pDC")] <- "pDC"
  o[m("dendritic|mDC|mregdc|FSCN")] <- "DC"
  o[m("macrophage|monocyte|myeloid|phagocyte|neutrophil")] <- "Myeloid"
  o[m("mast|basophil")] <- "Mast"
  o }

ac <- coarse(md$author_celltype)
bc <- coarse(md$cell_type)              # (b) marker + lineage-gate
cc <- coarse(md$final_ref_cell_type)    # (c) reference transfer
kp <- ac != "Other"
ag <- function(x) mean(x[kp & x != "Other"] == ac[kp & x != "Other"])

# ---- (a) cluster-level concordance + purity ---------------------------------
cl  <- md$seurat_clusters
maj <- tapply(ac, cl, function(v) { v <- v[v != "Other"]; if (!length(v)) NA_character_ else names(sort(table(v), decreasing = TRUE))[1] })
clagree <- mean(maj[as.character(cl)][kp] == ac[kp], na.rm = TRUE)
purity  <- mean(tapply(seq_along(cl), cl, function(i) { t <- table(ac[i]); max(t)/sum(t) }))

# ---- (3) IF cross-check: reference "Immune" that is CD45-neg + PanCK-pos -----
# Per-cluster IF (scaled 0..1 across clusters). Correct clusters whose reference
# label is immune but IF says epithelium (low CD45, high PanCK).
ccc <- cc                                # corrected reference label (coarse)
n_corr <- 0
if (!is.na(cfg$if_immune_col) && cfg$if_immune_col %in% colnames(md)) {
  cd45 <- tapply(md[[cfg$if_immune_col]], cl, mean); cd45 <- cd45 / max(cd45)
  pck  <- tapply(md[[cfg$if_epithelial_col]], cl, mean); pck <- pck / max(pck)
  ref_imm <- tapply(cc == "Myeloid" | cc == "T" | cc == "B_Plasma" | cc == "NK", cl, mean) > 0.5
  bad <- names(which(ref_imm & cd45 < 0.25 & pck > 0.40))   # immune label, but epithelial IF
  if (length(bad)) { fix <- as.character(cl) %in% bad; ccc[fix] <- "Epithelial"; n_corr <- sum(fix) }
  message(sprintf("  IF cross-check: %d cell(s) in %d cluster(s) reference-called Immune but CD45-neg/PanCK-pos -> corrected to Epithelial",
                  n_corr, length(bad)))
}
agc <- mean(ccc[kp & ccc != "Other"] == ac[kp & ccc != "Other"])

message(sprintf("  (a) CLUSTER concordance: %.1f%% (purity %.1f%%)", 100*clagree, 100*purity))
message(sprintf("  (b) MARKER+gate label : %.1f%%   <-- primary CosMx label", 100*ag(bc)))
message(sprintf("  (c) REFERENCE transfer: %.1f%%  (IF-corrected %.1f%%)", 100*ag(cc), 100*agc))

write.csv(data.frame(
  labeling = c("a_cluster_structure","b_marker_lineage_gate","c_reference_transfer","c_reference_IF_corrected","within_cluster_purity"),
  agreement_vs_author = round(c(clagree, ag(bc), ag(cc), agc, purity), 4)),
  tab_path(cfg, "benchmark_summary.csv"), row.names = FALSE)

# ---- (2) immune per-type recall (author -> ours), PRIMARY label (b) ---------
ai <- immcat(md$author_celltype); bi <- immcat(md$cell_type)
lev <- c("B","Plasma","T_CD4","T_CD8","Myeloid","DC","pDC","Mast","NK")
rec <- sapply(lev, function(L) { a <- which(ai == L); if (!length(a)) NA else mean(bi[a] == L, na.rm = TRUE) })
imm <- data.frame(immune_type = lev,
                  n_author = sapply(lev, function(L) sum(ai == L, na.rm = TRUE)),
                  n_ours_b = sapply(lev, function(L) sum(bi == L, na.rm = TRUE)),
                  recall_b = round(rec, 3))
write.csv(imm, tab_path(cfg, "benchmark_immune_recall.csv"), row.names = FALSE)
message("  immune per-type recall (primary label b):"); print(imm, row.names = FALSE)

# ---- coarse confusion (primary b) + fine cross-tab --------------------------
lv <- sort(unique(c(bc, ac)))
write.csv(as.data.frame.matrix(table(author = factor(ac, lv), ours_b = factor(bc, lv))),
          tab_path(cfg, "benchmark_coarse_confusion.csv"))
fine <- as.data.frame(table(author = md$author_celltype, ours_b = md$cell_type))
fine <- fine[fine$Freq > 0, ]; write.csv(fine[order(-fine$Freq), ], tab_path(cfg, "benchmark_fine_crosstab.csv"), row.names = FALSE)

# ---- figure: three-way agreement bar ----------------------------------------
bars <- data.frame(labeling = factor(c("(a) cluster\nstructure","(b) marker +\nlineage gate","(c) reference\ntransfer","(c) ref, IF-\ncorrected"),
                     levels = c("(a) cluster\nstructure","(b) marker +\nlineage gate","(c) reference\ntransfer","(c) ref, IF-\ncorrected")),
                   agreement = c(clagree, ag(bc), ag(cc), agc))
p <- ggplot(bars, aes(labeling, agreement, fill = labeling)) +
  geom_col(width = .65, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.0f%%", 100*agreement)), vjust = -0.4) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = paste0(cfg$label, " — labelings vs authors' 35 types"),
       subtitle = "clustering is platform-robust; reference-label transfer is the weak link",
       x = NULL, y = "coarse-lineage agreement") + theme_bw()
save_png(p, fig_path(cfg, "benchmark_threeway.png"), width = 8, height = 5)
write_session_info(cfg, "08")
message("== 08 done ==")
