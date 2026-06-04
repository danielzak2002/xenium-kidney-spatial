# Seurat v5 spatial workflow + RStudio MCP for Claude Code.
# Run once:  Rscript setup_R.R   (or source() inside RStudio)

options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages(c(
  "Seurat",       # v5: LoadXenium(), spatial workflows
  "SeuratObject",
  "Matrix",
  "arrow",        # read transcripts.parquet / cells.parquet
  "sf",           # cell / nucleus boundary polygons
  "data.table",
  "hdf5r",        # read cell_feature_matrix.h5
  "ggplot2",
  "patchwork",
  "remotes",
  "mcptools",     # exposes R as an MCP server to Claude Code
  "btw"           # default MCP toolset: env inspection, package docs, session info
))

# Not on CRAN:
remotes::install_github("immunogenomics/presto")  # fast Wilcoxon DE for markers
remotes::install_github("bnprks/BPCells/r")        # on-disk matrices for the big RCC run

cat("
Installed. Next steps for the RStudio MCP:

  1. Add to ~/.Rprofile so every interactive session is discoverable:
       if (interactive()) mcptools::mcp_session()

  2. Register the server with Claude Code (run in your terminal):
       claude mcp add -s user r-mcptools -- Rscript -e \"mcptools::mcp_server()\"

  3. Open RStudio (or Positron) so a live R session exists, then run  /mcp  inside
     Claude Code to confirm the connection.

  (mcptools was recently renamed from 'acquaint' and is evolving — if a command
   differs, check the mcptools README/CRAN page for the current syntax.)
")
