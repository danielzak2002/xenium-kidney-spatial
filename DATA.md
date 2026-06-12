# Data

Raw spatial-transcriptomics outputs are **not** stored in this repository. They are large
(GB–tens of GB) and released under their own licenses, so they are git-ignored. This document
explains how to obtain each dataset and where to put it. All committed work is reproducible from
these public sources; attribute the original providers when presenting results.

| Folder | Provider · platform | Source |
|---|---|---|
| `kidney_10x/`, `kidney_10x_preview/` | 10x Genomics · Xenium | 10x dataset pages (below) |
| `Danaher24/` | NanoString / Danaher et al. 2024 · CosMx | Zenodo `13964258` |
| `Demoulin26/` | Dumoulin et al. 2026 · Xenium 5K + CosMx 1k | Zenodo `19868428` |
| figshare ccRCC *(external, not stored)* | figshare · Xenium | figshare DOI `10.6084/m9.figshare.25685961` |

---

## 10x Genomics Xenium — RCC (`kidney_10x/`, `kidney_10x_preview/`)
Public 10x Genomics Xenium human-kidney datasets, **CC BY 4.0**, © 10x Genomics. Attribute 10x
Genomics as the data source.

- `kidney_10x/` — *Xenium Protein FFPE Human Renal Carcinoma* (clear-cell RCC; gene + 27-plex protein):
  https://www.10xgenomics.com/datasets/xenium-protein-ffpe-human-renal-carcinoma
- `kidney_10x_preview/` — *Human Kidney Preview Data (Xenium Human Multi-Tissue and Cancer Panel)*;
  the `data/` bundle here is the **cancer (PRCC) section** (`hKidney_cancer_section`):
  https://www.10xgenomics.com/datasets/human-kidney-preview-data-xenium-human-multi-tissue-and-cancer-panel-1-standard

(Each folder's `readme.txt` holds the background pasted from the 10x dataset page.)

**Download & placement**
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

Expected files in each `data/` folder include: `cell_feature_matrix.h5`, `cells.parquet`,
`cell_boundaries.parquet`, `nucleus_boundaries.parquet`, `transcripts.parquet`,
`morphology*.ome.tif`, `gene_panel.json`, `experiment.xenium`, `metrics_summary.csv`. The bundles
are dominated by morphology images and the per-transcript `transcripts.parquet`; the analysis loads
only the cell-feature matrix and cell coordinates, so the in-memory footprint is a small fraction of
the on-disk size (see CLAUDE.md for the loader flags). `gene_panel.json` is small and useful as
documentation — copy it into `outputs/tables/` if you want it under version control.

## NanoString CosMx — childhood lupus nephritis (`Danaher24/`)
A second platform/cohort: **NanoString CosMx 1000-plex** panel (957 genes in the deposited object)
spatial transcriptomics of **childhood-onset lupus nephritis** kidney — **control + SLE**, 14 tissue
sections analysed — with the authors' processed objects and cell-type annotations.

- **Source:** Zenodo record **13964258** — https://zenodo.org/records/13964258
- **License:** CC BY 4.0. Attribute Danaher et al. 2024 and NanoString when presenting results.
- **Reference:** Danaher P, *et al.* "Childhood-onset lupus nephritis is characterized by complex
  interactions between kidney stroma and infiltrating immune cells", 2024.
  Code: https://github.com/Nanostring-Biostats/childhood-onset-lupus-nephritis-analyses

The raw `.RData` objects and the authors' analysis files are **not** stored here — the entire
`Danaher24/data/` folder is git-ignored (only `Danaher24/readme.txt` is tracked). Download the Zenodo
bundle and extract so the layout is:

```
Danaher24/
├── readme.txt
└── data/                                          <- git-ignored (whole folder)
    ├── cleaneddata.RData                            # processed objects + annotations (analysed)
    ├── data_with_failed_slide.RData                 # includes the flagged failed slide
    └── Childhood onset lupus nephritis analyses/    # authors' R code / outputs
```

## Xenium 5K + CosMx 1k DKD atlas — `Demoulin26/`
The **diabetic kidney disease** spatial atlas of Dumoulin et al. 2026 (*Spatial atlas of diabetic
kidney disease reveals a B cell-rich subgroup*, Nature 2026): an **annotated** CosMx 1k + Xenium 5k
atlas as a single AnnData object (4,337,862 cells × 5,443-gene union; 16 Xenium + 48 CosMx samples),
plus a sample-diagnosis sheet. **Only the processed/annotated object is deposited — raw molecule
tables (`transcripts.parquet`) are NOT in the deposit** (so transcript-level re-segmentation of DKD
is not possible from Zenodo alone).

- **Source:** Zenodo record **19868428** — https://zenodo.org/records/19868428
- **License:** CC BY 4.0. Attribute Dumoulin et al. 2026.

Download both files into `Demoulin26/data/` (whole folder git-ignored; `readme.txt` + `ASSESSMENT.md`
are tracked):

```
Demoulin26/
├── readme.txt
├── ASSESSMENT.md
├── assess_dataset.py
└── data/                                          <- git-ignored (whole folder)
    ├── spatial_adata_xenium_cosmx_zenodo.h5ad      # annotated CosMx 1k + Xenium 5k atlas (8.7 GB)
    └── Diagnosis.xlsx                              # per-sample diagnosis / covariates
```

The object is read backed (`anndata backed='r'`, X never materialized) — see `Demoulin26/ASSESSMENT.md`
and `Demoulin26/assess_dataset.py`.

## figshare ccRCC — external, not stored
Used only for the three-cohort cross-comparison (`analysis/three_cohort_assessment/`,
`analysis/three_cohort_integration/`). It is **streamed, not committed** and not placed in any
`data/` folder.

- **Record:** *Xenium Renal Cell Carcinoma* (ProSeg demo data) — figshare DOI
  **10.6084/m9.figshare.25685961**.
- **Contents:** 10 × `transcripts.csv.gz` molecule tables (one per sample-region; 4 ccRCC patients,
  tumor + adjacent; custom 380-gene panel). The record ships **no pre-built cell × gene matrix** — it
  is derived on the fly from the transcript tables (`cell_id` + `feature_name`).
- **Access:** the HTML page blocks bots, but the figshare API and file endpoints do not — query
  `https://api.figshare.com/v2/articles/25685961` for file IDs, then stream each from
  `https://ndownloader.figshare.com/files/<id>`. License: CC BY 4.0; attribute the figshare depositor.

---

## Reference data (cell-type annotation)
Reference-based annotation uses two references, both **outside this repository** (the `refs/` folder
and the celldex cache are git-ignored):

- **Immune layer** (`R/04_reference_annotation.R`) — the **Monaco immune** reference from the
  Bioconductor `celldex` package, fetched + cached automatically on first use
  (`celldex::fetchReference("monaco_immune", ...)`). No manual download.
- **Kidney/epithelial layer** (`R/05_kidney_reference.R`) — the **KPMP/HuBMAP Azimuth human-kidney**
  reference (Zenodo record 10694842, CC BY 4.0). Download both files into `refs/azimuth_kidney/` once:

  ```bash
  mkdir -p refs/azimuth_kidney
  curl -L -o refs/azimuth_kidney/ref.Rds   "https://zenodo.org/records/10694842/files/ref.Rds?download=1"
  curl -L -o refs/azimuth_kidney/idx.annoy "https://zenodo.org/records/10694842/files/idx.annoy?download=1"
  ```
