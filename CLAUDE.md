# Spatial Transcriptomics Pilot — Human Kidney (Xenium)

## Goal
Build a reproducible analysis workflow for 10x Xenium human-kidney spatial data,
validated on a tiny dataset and then scaled to a full multimodal section. Scientific
emphasis: immune / B-cell spatial biology in kidney tissue — cell typing, spatial
niches, and tertiary-lymphoid-structure / immune-aggregate detection.

## Operator
Expert computational + translational biologist. R-first (expert), Python proficient.
Assume high fluency: do NOT explain single-cell basics, normalization theory, or
clustering fundamentals. Be concise and technical. Do explicitly flag non-obvious,
Xenium-specific pitfalls.

## Environment (hard constraints — do not exceed)
- MacBook Air, Apple M4, 24 GB unified memory, arm64.
- Usable working RAM ~16–18 GB after OS/GPU. Process ONE section at a time.
- Python: conda env `spatial` (miniforge, arm64, conda-forge). Never install into base.
- R: 4.4+ (CRAN arm64 build), Seurat v5.
- Large image and per-transcript files stay on disk. Never load them wholesale.

## Datasets
Two public 10x Genomics Xenium human-kidney datasets (CC BY 4.0). Raw outputs live in
each dataset's data/ subfolder and are git-ignored (see DATA.md). Each folder also has a
readme.txt with the background pasted from the 10x dataset page — READ IT during assessment.
- kidney_10x/data/
- kidney_10x_preview/data/

The panel-assessment step establishes panel size, modality (gene vs gene+protein), and
cell count for each (from gene_panel.json + metrics_summary.csv). The dataset with the
SMALLER cell count is the validation target for the small->big gate; the larger one is
scaled to only after the loop is proven. Do NOT assume which is which — confirm from
metrics_summary.csv. Treat all data/ contents as READ-ONLY. Large image and per-transcript
files (.ome.tif, transcripts.parquet) stay on disk; load the cell-feature matrix +
coordinates only.

## Working method — SMALL → BIG is a hard gate
1. ASSESS PANELS FIRST. Parse gene_panel.json for BOTH datasets. Report panel size
   and confirm presence of immune/B-cell markers (MS4A1, CD79A, CD79B, CD19) and
   plasma markers (MZB1, DERL3) before designing any cell typing. Enumerate the
   negative-control probe and blank-codeword features.
2. Validate the full loop on TINY: load -> QC -> normalize -> cluster -> annotate ->
   neighborhood enrichment. Parametrize paths and thresholds (no hard-coding).
3. ONLY THEN point the SAME parametrized code at BIG, layering in the protein
   modality and the spatial-niche / immune-aggregate analysis.

## Xenium QC rules (this is NOT scRNA-seq QC)
- Targeted, low-plex: per-cell counts are inherently low. Do NOT hard-filter cells
  on minimum counts/genes the way you would droplet data.
- QC on negative-control-probe and blank-codeword rates; flag cells/regions with
  elevated control signal.
- Implausibly high-count cells are usually segmentation merges — flag, don't trust.
- Report segmentation quality before trusting any cell-level result.

## Multimodal (gene + protein) handling
- Keep RNA and the 27-plex protein as SEPARATE assays/layers.
- Normalize protein independently (CLR, as for CITE-seq ADTs). Do NOT co-normalize
  protein with RNA. RNA: standard log-normalize or SCT.

## Tooling order
- Phase A (R / Seurat v5): LoadXenium(..., molecule.coordinates = FALSE) to stay lean.
  Initial visualization, clustering, and annotation here.
- Phase B (Python / squidpy + spatialdata): spatialdata_io.xenium() with lazy /
  Zarr-backed images; transcripts loaded lazily or skipped. Spatial graph stats:
  neighborhood enrichment, co-occurrence, Moran's I, niche detection.
- For large in-memory operations prefer BPCells (R) or backed AnnData (Python).

## Reproducibility
- Pin everything: export the conda env; capture sessionInfo() for R.
- Deterministic seeds. Prefer scripts over notebooks for anything rerun; notebooks
  for exploration only.
- Write intermediate objects to outputs/. Treat data/ as READ-ONLY.

## Human-in-the-loop
Strategy is reviewed in a separate Claude chat. After each milestone, emit a COMPACT
summary built for pasting: metrics tables, marker lists, cluster counts, and key code
diffs — not walls of log output.

## Repository (public git repo)
This is a public, reproducible repo: CC BY 4.0 data (10x Genomics), MIT-licensed code.
Discipline:
- NEVER commit anything under */data/ or large binaries (.ome.tif, .parquet, .zarr, .h5,
  .h5ad, .rds, .mtx). The .gitignore enforces this — verify `git status` is clean of data
  before every commit.
- Commit small figures (outputs/figures/*.png) and tables (outputs/tables/*.csv) that
  document results; git-ignore large objects (outputs/objects/).
- No secrets, no .Rhistory/.RData. Use descriptive, conventional commit messages.
- Keep README.md and DATA.md accurate as the work evolves — they are the public face.

## Project layout
README.md            # public overview
DATA.md              # how to obtain the (git-ignored) raw data
CLAUDE.md            # this file
LICENSE              # MIT (code)
.gitignore
environment.yml      # Python env
setup_R.R            # R deps
kidney_10x/          # readme.txt + data/ (raw, git-ignored)
kidney_10x_preview/  # readme.txt + data/ (raw, git-ignored)
R/                   # Seurat scripts
py/                  # squidpy / spatialdata scripts
outputs/             # figures/ + tables/ (committed), objects/ (git-ignored)
notebooks/           # exploratory only
