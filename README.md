# Spatial Transcriptomics of Human Kidney (10x Xenium)

A reproducible pilot analysis of imaging-based spatial transcriptomics in human kidney
tissue, using two public 10x Genomics Xenium datasets. The workflow covers panel
assessment, QC tailored to imaging-based ST, cell-type annotation, spatial-neighborhood
analysis, and immune / B-cell niche detection — implemented in both R (Seurat v5) and
Python (squidpy / spatialdata) to compare the two ecosystems on the same data.

## Motivation
Imaging-based spatial transcriptomics resolves single cells in situ, preserving the
tissue architecture that dissociation-based scRNA-seq discards. In kidney, that spatial
context is central to understanding the immune contexture — where infiltrating and
resident immune populations localize relative to nephron structures, and whether they
organize into aggregates or tertiary lymphoid structures. This pilot builds a clean,
reusable workflow for those questions on Xenium data.

## Datasets
Two public 10x Genomics Xenium human-kidney datasets (CC BY 4.0). Raw outputs are **not**
included in this repository — see [DATA.md](DATA.md) for download and placement.

| Folder | Dataset | Panel / modality | Cells |
|---|---|---|---|
| `kidney_10x/` | Xenium Protein FFPE Human Renal Cell Carcinoma (Stage III, T3a N1 MX) | 405 genes / 430 probes (Human Multi-Tissue & Cancer + 28-gene custom add-on, 25 boosted) **+ 27-plex protein** — gene + protein | 465,534 |
| `kidney_10x_preview/` | Human Kidney Preview — **cancer (PRCC) section only** (`hKidney_cancer_section`) | 377 genes (Human Multi-Tissue & Cancer, dev v1.5.0) — gene only | 56,510 ⟵ **validation target** |

Populated from `gene_panel.json` + `metrics_summary.csv` (panel assessment step). Two
notes: (1) the `kidney_10x_preview/data/` bundle on disk is the **cancer/PRCC section**
of the preview release (56,510 cells), not the paired non-diseased section (97,560);
(2) it has the fewer cells, so it is the small→big validation target. Both panels share
the same 377-gene Human Multi-Tissue & Cancer base and the same 20 negative-control
probes / 41 negative-control codewords; only `kidney_10x` adds the custom genes and the
27-plex protein assay.

## Approach
1. **Panel assessment** — parse each `gene_panel.json`; report panel size and confirm
   immune/B-cell and plasma markers before designing cell typing.
2. **QC for imaging-based ST** — negative-control-probe and blank-codeword rates;
   segmentation quality; flag segmentation merges (no scRNA-style count hard-filters).
3. **Clustering + annotation** (Seurat v5) — Leiden clustering, marker-based annotation.
4. **Spatial statistics** (squidpy / spatialdata) — neighborhood enrichment,
   co-occurrence, Moran's I, spatial-niche detection.
5. **Immune characterization** — localize immune populations; detect immune aggregates /
   candidate tertiary lymphoid structures.

## Reproduce
```bash
# Python (Apple-silicon arm64; miniforge recommended)
conda env create -f environment.yml
conda activate spatial
python -m ipykernel install --user --name spatial

# R (4.4+)
Rscript setup_R.R

# Data: follow DATA.md to download the two datasets into their data/ subfolders.
# Then run the R scripts in R/ and the Python scripts in py/.
```

## Repository layout
```
README.md            Public overview (this file)
DATA.md              How to obtain the (git-ignored) raw data
CLAUDE.md            Analysis plan, constraints, and conventions
environment.yml      Python environment
setup_R.R            R dependencies
R/                   Seurat scripts
py/                  squidpy / spatialdata scripts
outputs/             figures/ + tables/ (committed); objects/ (git-ignored)
notebooks/           Exploratory only
kidney_10x/          readme.txt + data/ (raw outputs, git-ignored)
kidney_10x_preview/  readme.txt + data/ (raw outputs, git-ignored)
```

## Tooling
R 4.4+ with Seurat v5; Python 3.11 with scanpy + squidpy + spatialdata. Developed on
Apple silicon (arm64).

## Results
_To be populated as analyses complete (figures in `outputs/figures/`, tables in
`outputs/tables/`)._

## License
- **Code:** MIT — see [LICENSE](LICENSE).
- **Data:** Creative Commons Attribution 4.0 (CC BY 4.0), © 10x Genomics. Attribute 10x
  Genomics as the data source when presenting results. See [DATA.md](DATA.md).

## Author
Daniel — _add a sentence and a link (LinkedIn/site) here._
