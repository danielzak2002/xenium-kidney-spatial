# Data

The raw Xenium outputs are **not** stored in this repository. They are large
(tens of GB) and released under their own license, so they are git-ignored. This
document explains how to obtain them and where to put them.

## Source & license
Both datasets are public 10x Genomics Xenium human-kidney datasets, released under
**Creative Commons Attribution 4.0 (CC BY 4.0)**, © 10x Genomics. When presenting or
publishing any results derived from them, attribute 10x Genomics as the data source.

- `kidney_10x/` — *Xenium Protein FFPE Human Renal Carcinoma* (gene + protein):
  https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-carcinoma
- `kidney_10x_preview/` — *Human Kidney Preview Data (Xenium Human Multi-Tissue and Cancer Panel)*;
  the `data/` bundle here is the **cancer (PRCC) section** (`hKidney_cancer_section`):
  https://www.10xgenomics.com/datasets/human-kidney-preview-data-xenium-human-multi-tissue-and-cancer-panel-1-standard

(Each folder's `readme.txt` holds the background pasted from the 10x dataset page.)

## Download & placement
1. Create a free account at 10xgenomics.com and accept the dataset terms.
2. Download each dataset's Xenium output bundle.
3. Extract so the 10x outputs land in the `data/` subfolder of each dataset:

```
xenium/
├── kidney_10x/
│   ├── readme.txt
│   └── data/             <- 10x output bundle extracted here
└── kidney_10x_preview/
    ├── readme.txt
    └── data/             <- 10x output bundle extracted here
```

Expected files in each `data/` folder include:
`cell_feature_matrix.h5`, `cells.parquet`, `cell_boundaries.parquet`,
`nucleus_boundaries.parquet`, `transcripts.parquet`, `morphology*.ome.tif`,
`gene_panel.json`, `experiment.xenium`, `metrics_summary.csv`.

## Size note
The bundles are dominated by morphology images and the per-transcript
`transcripts.parquet`. The analysis loads only the cell-feature matrix and cell
coordinates, so the in-memory footprint is a small fraction of the on-disk size
(see CLAUDE.md for the loader flags that keep it lean).

## Tip
`gene_panel.json` is small and useful as documentation — if you want it under version
control, copy it out of the (git-ignored) `data/` folder into `outputs/tables/` and
commit that copy.

## Reference data (cell-type annotation)
Reference-based annotation uses two references, both **outside this repository** (the
`refs/` folder and the celldex cache are git-ignored):

- **Immune layer** (`R/04_reference_annotation.R`) — the **Monaco immune** reference from
  the Bioconductor `celldex` package, fetched + cached automatically on first use
  (`celldex::fetchReference("monaco_immune", ...)`). No manual download.
- **Kidney/epithelial layer** (`R/05_kidney_reference.R`) — the **KPMP/HuBMAP Azimuth
  human-kidney** reference (Zenodo record 10694842, CC BY 4.0). Download both files into
  `refs/azimuth_kidney/` once:

  ```bash
  mkdir -p refs/azimuth_kidney
  curl -L -o refs/azimuth_kidney/ref.Rds   "https://zenodo.org/records/10694842/files/ref.Rds?download=1"
  curl -L -o refs/azimuth_kidney/idx.annoy "https://zenodo.org/records/10694842/files/idx.annoy?download=1"
  ```

## Danaher et al. 2024 — childhood-onset lupus nephritis (cLN) CosMx cohort
A second platform/cohort for cross-platform validation of the annotation engine:
**NanoString CosMx 1000-plex** spatial transcriptomics of childhood-onset lupus
nephritis kidney (8 cLN + 2 control), with the authors' processed objects and 35-type
cell annotations.

- **Source:** Zenodo record **13964258** — https://zenodo.org/records/13964258
- **License:** Creative Commons Attribution 4.0 (CC BY 4.0). Attribute Danaher et al. 2024
  and NanoString when presenting derived results.
- **Reference:** Danaher P, *et al.* "Childhood-onset lupus nephritis is characterized by
  complex interactions between kidney stroma and infiltrating immune cells", 2024.
  Code: https://github.com/Nanostring-Biostats/childhood-onset-lupus-nephritis-analyses

The raw `.RData` objects and the authors' analysis files are **not** stored here — the
entire `Danaher24/data/` folder is git-ignored (only `Danaher24/readme.txt` is tracked).
Download the Zenodo bundle and extract so the layout is:

```
Danaher24/
├── readme.txt
└── data/                                   <- git-ignored (whole folder)
    ├── cleaneddata.RData                    # processed objects + annotations
    ├── data_with_failed_slide.RData
    └── Childhood onset lupus nephritis analyses/   # authors' R code / outputs
```
