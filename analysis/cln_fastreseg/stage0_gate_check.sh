#!/usr/bin/env bash
# stage0_gate_check.sh — reproducible STAGE 0 gate checks for the cLN FastReseg attempt.
# Read-only. Pure science. Determines whether FastReseg can run at all on the available data.
set -euo pipefail
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
CLN_RDATA="$REPO/Danaher24/data/cleaneddata.RData"
DEMOULIN="$REPO/Demoulin26/data/spatial_adata_xenium_cosmx_zenodo.h5ad"

# ---- GATE (a): does the cLN release contain a per-TRANSCRIPT table? ----
# FastReseg needs per-molecule {target, global x/y(+z), original cell assignment}.
# List every object in the cLN .RData and its class/dim; inspect column names.
Rscript -e '
e <- new.env(); load("'"$CLN_RDATA"'", envir=e)
for (n in ls(e)) { o <- get(n,e)
  cat(sprintf("%-12s %-14s %s\n", n, paste(class(o),collapse=","),
      paste(if(!is.null(dim(o))) dim(o) else length(o), collapse="x"))) }
cat("\nannot columns:\n"); print(colnames(get("annot",e)))
cat("\ncoordinate objects (per-CELL, not per-transcript):\n")
print(colnames(get("customlocs",e)))
# extract cLN gene panel (the 957-dim of raw) for GATE (b)
raw <- get("raw",e); gd <- which(dim(raw)==957); writeLines(dimnames(raw)[[gd]], "/tmp/cln_panel_genes.txt")
'
# VERDICT (a): objects are annot (cell metadata), customlocs/um/viz (per-CELL centroids),
# raw (gene x cell counts), clust (per-cell). NO per-molecule export -> GATE (a) FAILS.

# ---- GATE (b): cLN panel ∩ Demoulin CosMx panel (moot if (a) fails, but recorded) ----
conda run -n spatial python -c "
import anndata as ad
a=ad.read_h5ad('$DEMOULIN', backed='r'); dem=set(map(str,a.var_names)); a.file.close()
cln=set(l.strip() for l in open('/tmp/cln_panel_genes.txt') if l.strip())
inter=cln&dem
key=['CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2','CD4','CD8A','CD8B','PTPRC','MS4A1','CD79A','MZB1','KRT8','KRT18','EPCAM','CD68']
print('cLN',len(cln),'Demoulin',len(dem),'intersection',len(inter))
print('CD3-family + lineage in intersection:', [g for g in key if g in inter])
print('absent from both 1k panels:', [g for g in key if g not in inter])
"
