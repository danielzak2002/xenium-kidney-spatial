#!/usr/bin/env Rscript
# 08_benchmark.R — score OUR engine's labels against the authors' independent
# annotation (CosMx cLN: annot$clust, 35 types). Both label vocabularies are
# mapped to a COMMON COARSE LINEAGE set for the agreement metric; the fine
# cross-tab is reported separately (no raw-string scoring). Runs only when an
# author_celltype benchmark column is present.
#
#   XENIUM_DATASET=cln_cosmx Rscript R/08_benchmark.R   (run after 07)

suppressPackageStartupMessages({ library(Seurat); library(ggplot2) })
source(file.path(Sys.getenv("XENIUM_ROOT", unset = getwd()), "R", "config.R"))
cfg <- ensure_dirs(get_config())
message("== 08_benchmark :: ", cfg$label, " ==")

obj <- readRDS(cfg$obj_05)                       # carries 04/05/06/07 final labels
if (!"author_celltype" %in% colnames(obj[[]])) { message("  no author_celltype; skip"); quit(save="no") }

# ---- Common coarse lineage map (handles BOTH vocabularies via keywords) -----
coarse <- function(x) {
  x <- as.character(x); o <- rep("Other", length(x)); m <- function(p) grepl(p, x, ignore.case = TRUE)
  o[m("tubule|PCT|proximal|distal|connecting|intercalat|principal|thick.?ascending|loop.of.henle|\\bTAL\\b|thin.limb|parietal|pelvic|urotheli|macula|papillary|collecting|\\bDCT\\b|\\bCNT\\b|epithelial")] <- "Epithelial"
  o[m("podocyte")]                                                              <- "Podocyte"
  o[m("endotheli|vasa.recta|capillary|glomerular.endo|lymphatic")]              <- "Endothelial"
  o[m("fibroblast|myofibroblast|mesangial|pericyte|smooth.?muscle|stroma|mural|renin")] <- "Stroma"
  o[m("macrophage|monocyte|dendritic|\\bDC\\b|mregdc|\\bpdc\\b|\\bmdc\\b|myeloid|neutrophil|mast|basophil|granulocyte|phagocyte")] <- "Myeloid"
  o[m("\\bB\\b|B.?cell|naive B|memory B|plasmabl|plasma|\\bMZB")]                <- "B_Plasma"
  o[m("\\bNK\\b|natural killer")]                                               <- "NK"
  o[m("\\bT\\b|T.?cell|T CD|\\bCD4\\b|\\bCD8\\b|treg|regulatory|\\bMAIT\\b|gd T|Th1|Th2|Th17|helper|CCR7\\+ T")] <- "T"
  o
}
md <- obj[[]]
oc <- coarse(md$final_ref_cell_type); ac <- coarse(md$author_celltype)
keep <- oc != "Other" & ac != "Other"
agree <- mean(oc[keep] == ac[keep])

# CLUSTER-level concordance: assign each Leiden cluster its majority author coarse
# type and measure agreement. This isolates clustering quality (does the engine
# RECOVER the populations) from reference-label transfer (does labelling them work).
cl  <- md$seurat_clusters
maj <- tapply(ac, cl, function(v) { v <- v[v != "Other"]
  if (!length(v)) NA_character_ else names(sort(table(v), decreasing = TRUE))[1] })
clpred  <- maj[as.character(cl)]
clagree <- mean(clpred[keep] == ac[keep], na.rm = TRUE)
purity  <- mean(tapply(seq_along(cl), cl, function(i) { t <- table(ac[i]); max(t)/sum(t) }))
message(sprintf("  CLUSTER-level concordance: %.1f%% (purity %.1f%%) — clustering recovers populations",
                100*clagree, 100*purity))
message(sprintf("  reference-LABEL agreement:  %.1f%% over %d cells (%.1f%% mapped)",
                100*agree, sum(keep), 100*mean(keep)))
write.csv(data.frame(
  metric = c("cluster_level_concordance", "within_cluster_purity", "reference_label_agreement"),
  value  = round(c(clagree, purity, agree), 4)),
  tab_path(cfg, "benchmark_summary.csv"), row.names = FALSE)

lev <- sort(unique(c(oc, ac)))
cm  <- table(author = factor(ac, lev), ours = factor(oc, lev))
write.csv(as.data.frame.matrix(cm), tab_path(cfg, "benchmark_coarse_confusion.csv"))
# per-lineage recall (author -> ours) and our composition
recall <- diag(cm) / pmax(rowSums(cm), 1)
summ <- data.frame(lineage = lev, n_author = as.integer(rowSums(cm)),
                   n_ours = as.integer(colSums(cm)), recall = round(recall, 3))
write.csv(summ, tab_path(cfg, "benchmark_agreement.csv"), row.names = FALSE)
message("  overall coarse agreement written; per-lineage recall:")
print(summ, row.names = FALSE)

# ---- Fine cross-tab (author 35 vs our fine labels), long form ---------------
fine <- as.data.frame(table(author = md$author_celltype, ours = md$final_ref_cell_type))
fine <- fine[fine$Freq > 0, ]; fine <- fine[order(-fine$Freq), ]
write.csv(fine, tab_path(cfg, "benchmark_fine_crosstab.csv"), row.names = FALSE)

# ---- Figure: coarse confusion heatmap (row-normalized recall) ---------------
pr <- as.data.frame(prop.table(cm, 1)); colnames(pr) <- c("author","ours","frac")
p <- ggplot(pr, aes(ours, author, fill = frac)) + geom_tile() +
  geom_text(aes(label = ifelse(frac > 0.04, sprintf("%.2f", frac), "")), size = 3) +
  scale_fill_viridis_c(name = "row frac") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = sprintf("%s — coarse lineage agreement %.0f%% (vs authors' 35 types)",
                       cfg$label, 100*agree), x = "our coarse lineage", y = "author coarse lineage")
save_png(p, fig_path(cfg, "benchmark_coarse_confusion.png"), width = 8, height = 6)
write_session_info(cfg, "08")
message("== 08 done ==")
