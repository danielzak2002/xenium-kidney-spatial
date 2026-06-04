# Data

The raw Xenium outputs are **not** stored in this repository. They are large
(tens of GB) and released under their own license, so they are git-ignored. This
document explains how to obtain them and where to put them.

## Source & license
Both datasets are public 10x Genomics Xenium human-kidney datasets, released under
**Creative Commons Attribution 4.0 (CC BY 4.0)**, © 10x Genomics. When presenting or
publishing any results derived from them, attribute 10x Genomics as the data source.

- `kidney_10x/` — _dataset title + URL, fill in from the readme.txt_
- `kidney_10x_preview/` — _dataset title + URL, fill in from the readme.txt_

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
