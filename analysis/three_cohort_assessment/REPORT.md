# Three-cohort all-Xenium integration — feasibility assessment

Read-only. Assesses the figshare ccRCC dataset for a 3-cohort all-Xenium integration with the local BIG RCC (kidney_10x) and the Demoulin DKD Xenium subset; computes the shared gene space and the ProSeg (transcript-table) dependency. Pure science.

## STEP 0 — figshare ccRCC record (DOI 10.6084/m9.figshare.25685961)

**Xenium Renal Cell Carcinoma** (ProSeg demo data, *Cell Simulation as Cell Segmentation*). Pulled via the figshare API (`/v2/articles/25685961`) + `ndownloader` file endpoints.

- **Samples:** four ccRCC patients, each **tumor + adjacent** tissue (per the record description).
- **Files: 10 × `transcripts.csv.gz`** (104 MB – 5.0 GB; one molecule table per sample-region).
- **(a) 4 samples tumor+adjacent:** ✓ confirmed (description).
- **(b) cell × gene matrix:** **ABSENT** — the record ships *only* transcript tables. A matrix is derivable (the transcript table carries `cell_id`), but none is provided pre-built.
- **(c) transcript-level molecule table:** ✓ **present** — columns `transcript_id, cell_id, overlaps_nucleus, feature_name, x_location, y_location, z_location, qv` → target + x/y(/z) + cell assignment, exactly what ProSeg needs.
- **Panel:** custom **380 genes** (280 from 10x's breast base panel + 100 kidney / T-cell-subset genes) + 161 control/blank codewords (extracted from `feature_name`).

## STEP 1 — gene-panel intersection (the integration feature space)

| panel | genes |
|---|---|
| figshare ccRCC (custom) | 380 |
| BIG RCC (kidney_10x, 10x panel) | 405 |
| DKD Xenium (Demoulin var, union CosMx+Xenium) | 5443 |

- Pairwise: **figshare∩BIG = 149**, figshare∩DKD = 283, BIG∩DKD = 279.
- **THREE-WAY shared gene set = 123 genes** (`three_way_shared_genes.txt`) — the integration feature space.
- The bottleneck is **figshare∩BIG (149)**: two different ~380–405-gene custom panels share only a moderate core; DKD's 5k panel covers most of that core, so the 3-way ≈ figshare∩BIG minus a few.
- *Caveat:* DKD genes are the var **union** (CosMx+Xenium); a gene in the 3-way set that is CosMx-only on DKD would be a structural zero on DKD Xenium. The immune markers below are checked for actual Xenium measurement to avoid that trap.

## STEP 2 — immune-marker coverage across the three (load-bearing markers)

DKD column = present in var **AND measured on Xenium** (not a structural zero).

| group | gene | figshare | BIG RCC | DKD Xenium | all 3 |
|---|---|:--:|:--:|:--:|:--:|
| T-lineage | CD3D | ✓ | ✓ | — |  |
| T-lineage | CD3E | ✓ | ✓ | ✓ | **✓** |
| T-lineage | CD3G | ✓ | ✓ | ✓ | **✓** |
| T subtype | CD4 | ✓ | ✓ | ✓ | **✓** |
| T subtype | CD8A | ✓ | ✓ | ✓ | **✓** |
| T subtype | CD8B | — | — | ✓ |  |
| Treg | FOXP3 | ✓ | ✓ | ✓ | **✓** |
| Treg | IL2RA | ✓ | ✓ | ✓ | **✓** |
| Treg | CTLA4 | ✓ | ✓ | ✓ | **✓** |
| Cytotoxic | GZMB | — | ✓ | ✓ |  |
| Cytotoxic | GZMK | ✓ | ✓ | ✓ | **✓** |
| Cytotoxic | NKG7 | ✓ | ✓ | — |  |
| B | MS4A1 | ✓ | ✓ | ✓ | **✓** |
| B | CD79A | ✓ | ✓ | ✓ | **✓** |
| Plasma | MZB1 | ✓ | ✓ | ✓ | **✓** |
| Myeloid | CD68 | ✓ | ✓ | ✓ | **✓** |
| Myeloid | LYZ | ✓ | — | — |  |
| Endothelial | PECAM1 | ✓ | ✓ | ✓ | **✓** |
| Epithelial | EPCAM | ✓ | ✓ | ✓ | **✓** |
| Epithelial | KRT8 | ✓ | ✓ | — |  |

**14/20 markers survive all three:** CD3E, CD3G, CD4, CD8A, FOXP3, IL2RA, CTLA4, GZMK, MS4A1, CD79A, MZB1, CD68, PECAM1, EPCAM.

**Load-bearing Treg/cytotoxic differential is supportable in the 3-way space:** the full Treg set **FOXP3 / IL2RA / CTLA4** and cytotoxic **CD8A / GZMK** are present in all three. Drop-outs to note: **GZMB** absent from the figshare panel (CD8A+GZMK still cover cytotoxic); **CD8B** only on DKD; **CD3D and NKG7** are structural-zero on DKD Xenium (CD3D, NKG7, LYZ, KRT8 not measured) — use **CD3E/CD3G** for T-lineage.

## STEP 3 — ProSeg feasibility (transcript availability) — GO / NO-GO

| cohort | transcript / molecule table | ProSeg |
|---|---|:--:|
| figshare ccRCC | 10 × `transcripts.csv.gz` (feature_name + x/y/z + cell_id) | **GO** |
| BIG RCC (kidney_10x) | `transcripts.parquet` (1.0 GB) + `transcripts.zarr.zip` in the 10x bundle | **GO** |
| DKD (Demoulin Xenium) | **none** — the Zenodo object is a **cell-level AnnData** (X = cell×gene, obsm = centroids/X_umap/spatial_fov); no per-molecule table | **NO-GO** |

**Blocker:** uniform ProSeg across all three is **blocked by DKD** — the public Demoulin Zenodo release has no transcript/molecule table, only the segmented cell×gene AnnData. Re-segmenting DKD with ProSeg would require the **raw Demoulin Xenium output bundles** (per-sample `transcripts.parquet`), which are **not in the Zenodo deposit** and would have to be obtained elsewhere (authors / original repository). figshare and BIG RCC can be ProSeg-re-segmented today.

## Bottom line

- A 3-cohort all-Xenium integration is feasible on a **123-gene shared space** with the load-bearing immune markers intact (14/20, full Treg + CD8A/GZMK).
- **Uniform ProSeg re-segmentation is NOT currently possible** (DKD has no transcripts in the public release). Options: (i) integrate on the existing segmentations (figshare-derived + 10x + Demoulin), or (ii) ProSeg only figshare+BIG RCC and keep DKD on its native segmentation, noting the asymmetry, or (iii) obtain the raw Demoulin Xenium bundles to ProSeg all three uniformly.

*Read-only assessment; no raw files modified. figshare files were streamed for panel extraction (smallest transcript table only) and not committed.*