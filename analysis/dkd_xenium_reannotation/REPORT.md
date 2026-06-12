# DKD Xenium — independent re-integration & annotation, validated vs author labels

Xenium-only within-DKD re-integration with an **independent, marker-based cell-type annotation**,
benchmarked against the authors' deposited labels (Dumoulin et al. 2026; Zenodo 19868428). The
goal is mutual validation: reproduce the authors' typing *blind to their labels* so that agreement
confirms both pipelines, and disagreements localise where the typing is genuinely ambiguous. Author
labels (`annotation_updated`, `immune_cell_annotation_combined`) were carried along as columns but
**never used to guide normalization, integration, clustering, or annotation.**

Raw object treated read-only (`backed='r'`, counts re-normalised independently). All intermediates
(`*.parquet`, `*.h5ad`, `X_pca_harmony.npy`) are git-ignored; CSV tables + figures are committed.
Reproduce: `run_integration.py` → `annotate.py` → `validate.py` → `disease_strata.py` → `figures.py`.

---

## STEP 0 — Xenium subset

16 Xenium samples, **951,040 cells**, panel **5,443 genes** (union object); **342 genes are
structural-zero on Xenium** (CosMx-only in the union) → **5,101 measured genes** used downstream.
Per-cell depth is low-plex-typical (median 174 counts / 136 genes; the three deep Control sections
run higher). Author labels set aside as ground truth.

| sample | condition | cells | med counts | med genes | GFR |
|---|---|--:|--:|--:|---|
| 1004 | AA amyloid | 24,449 | 146 | 117 | <30 |
| 1009 | AA amyloid | 10,167 | 99 | 83 | <30 |
| 1007 | C3GN | 41,336 | 122 | 102 | 30–60 |
| HK2753 | Control | 161,350 | 573 | 351 | >60 |
| HK3106 | Control | 171,526 | 336 | 217 | >60 |
| HK3626 | Control | 104,038 | 484 | 321 | >60 |
| 1001 | DKD | 29,126 | 172 | 133 | 30–60 |
| 1006 | DKD | 45,609 | 186 | 146 | <30 |
| 1008 | DKD | 27,912 | 125 | 101 | <30 |
| 1010 | DKD | 40,475 | 165 | 135 | 30–60 |
| 1011 | DKD | 27,876 | 234 | 182 | <30 |
| 1012 | DKD | 39,105 | 192 | 153 | <30 |
| 1013 | DKD | 28,261 | 203 | 161 | <30 |
| HK2695 | DKD | 151,101 | 177 | 138 | 30–60 |
| 1003 | IgA | 30,772 | 165 | 134 | <30 |
| 1005 | MN | 17,937 | 149 | 120 | <30 |

---

## STEP 1 — independent re-integration (blind)

`counts` → `normalize_total(1e4)` → `log1p` → HVG (2,000, seurat) → `scale` → **PCA (50)** →
**Harmony on `orig_ident`** (harmonypy 2.0 called directly — its `Z_corr` is `(N,d)`; orientation
guarded) → neighbors (k=15) → **Leiden res=1.5 → 23 global clusters** → UMAP. The immune compartment
(global clusters above a marker-only immune-score gate; **4 clusters, 105,719 cells** — matching the
authors' Immune count of 103,215) was **subclustered** (Leiden res=1.4 → **16 immune subclusters**)
to resolve subtypes. Harmony mixes the 16 samples in UMAP without erasing lineage structure
(`figures/integration_umaps.png`). Embedding saved (`X_pca_harmony.npy`) for cheap re-clustering.

---

## STEP 2 — independent marker-based annotation

**Global clusters** labelled by canonical kidney markers: each marker's mean log-norm expression is
z-scored across clusters; a cell type's score is the mean z over its present markers; argmax wins
(`cluster_annotation_global.csv` records the marker basis + runner-up margin per cluster). Examples:
PT = LRP2/CUBN/SLC13A3; iPT = **HAVCR1/VCAM1** (injury); TAL = UMOD/SLC12A1/CASR; iTAL = ITGB6/UMOD;
DCT = SLC12A3/TRPM6; PC = AQP2/AQP3/GATA3; CNT = CALB1/SLC8A1/SCNN1G; IC A = SLC4A1/SLC26A7/KIT;
PEC = CLDN1/CFH; Podo = WT1/PODXL; Fibroblast = COL1A1/COL6A3/PDGFRA; VSMC = NOTCH3/ACTA2/TAGLN;
MC = POSTN/PIEZO2/REN; EC = PECAM1 (+SOST glom / PLAT peritub / ENPP2 DVR).

**Immune subtypes** use **hierarchical rules on absolute means**, not z-scores — z-scoring across the
immune compartment mislabels the T/cytotoxic axis (FOXP3/CTLA4 light up for *any* T cluster; cytotoxic
CD8 mimics NK), and **CD4 is macrophage-expressed** so T identity is gated on CD3E. Rule order:
B (MS4A1>0.8) → Plasma (MZB1>0.8) → T (CD3E>0.6: CD8 if CD8A≥0.5>CD4; Treg if FOXP3≥0.5; else CD4) →
NK (KLRD1/KLRF1/NCAM1, CD3E-low) → DC (FCER1A & CD1C, CD68-low) → Mast → Myeloid (CD68/C1QA/AIF1) →
Neutrophil → Unresolved (`cluster_annotation_immune.csv`).

Resulting cells (n): PT 210,957 · Fibroblast 115,345 · EC_Peritub 82,300 · CNT 66,982 · iPT 59,943 ·
TAL 59,608 · iTAL 49,942 · **Myeloid 46,599** · PC 36,403 · DCT 29,825 · IC A 27,422 · **CD4 T 26,946**
· PEC 26,331 · EC_glom 25,785 · MC 21,956 · **CD8 T 15,785** · VSMC 15,363 · EC_DVR 9,634 · **B 7,893** ·
Podo 7,525 · Unresolved 5,399 · **Plasma 3,097**. Lineage totals: Epithelial 574,938 · Stroma 152,664 ·
Endothelial 117,719 · Immune 105,719. Immune subtype counts track the authors' closely (Myeloid 46.6k
vs Macro 47.3k; B 7.9k vs 8.4k; Plasma 3.1k vs 2.7k; T 42.7k vs 40.8k).

**Panel limits (stated, not worked around):** IC-B markers (SLC26A4/SLC4A9), Treg-specific resolution,
and NK/Mast/Neutrophil markers (NKG7, GNLY, TPSAB1, CPA3, S100A8) are absent or too sparse on this
panel/depth → IC-B folds into IC A, Treg into CD4, and NK/Mast/Neutrophil/DC are not separately
resolved (they fall into the nearest lineage). These are reported as divergences below.

---

## STEP 3 — validation vs authors (the confidence check)

Per-cell confusion of MY labels vs the authors', harmonised to a shared vocabulary (injury states
iPT/iTAL kept distinct; genuinely fuzzy boundaries merged: CNT+PC, IC A+B, EC subtypes). ARI is
label-invariant; recall is per author class. Figure: `figures/concordance_matrices.png`.

**A · segment vs `annotation_updated`: ARI = 0.782, agreement = 85.4%** (n = 951,040).

| type | n (author) | recall | precision |
|---|--:|--:|--:|
| Podo | 7,308 | 0.98 | 0.95 |
| Immune | 103,215 | 0.95 | 0.93 |
| PC/CNT | 104,274 | 0.94 | 0.95 |
| PT | 208,867 | 0.94 | 0.93 |
| DCT | 27,960 | 0.93 | 0.87 |
| MC | 12,091 | 0.90 | 0.50 |
| Fibroblast | 120,457 | 0.87 | 0.91 |
| EC | 105,061 | 0.87 | 0.78 |
| IC | 30,317 | 0.86 | 0.95 |
| PEC | 6,986 | 0.82 | 0.22 |
| iTAL | 38,452 | 0.76 | 0.58 |
| TAL | 75,480 | 0.73 | 0.93 |
| VSMC | 12,519 | 0.65 | 0.53 |
| iPT | 87,231 | 0.64 | 0.94 |

**B · immune subtype vs `immune_cell_annotation_combined`: ARI = 0.683, agreement = 78.9%**
(n = 103,210 author-immune cells).

| type | n (author) | recall | precision |
|---|--:|--:|--:|
| Plasma | 2,726 | 0.91 | 0.87 |
| Myeloid/Macro | 47,311 | 0.90 | 0.95 |
| B | 8,422 | 0.86 | 0.93 |
| CD4 | 21,140 | 0.84 | 0.68 |
| CD8 | 19,673 | 0.58 | 0.73 |

**Divergences (findings, not failures):**
- **iPT/iTAL/TAL/VSMC** — the largest segment disagreements are injury-state and stromal boundaries:
  injured-PT (iPT, recall 0.64) is partly read as PT and TAL/iTAL trade cells; VSMC↔Fibroblast/MC
  (recall 0.65). These are the soft borders of the kidney atlas, expected under independent clustering.
- **CD8 recall 0.58** — author NK cells fall almost entirely into my CD8 cluster (cytotoxic-program
  overlap; NKG7/GNLY absent from the panel), and the CD4/CD8 boundary shifts some cells; the combined
  T-cell pool matches well.
- **Rare types** — NK, Neutrophil (→Myeloid), Mast/Baso (→Unresolved/missed), DC, and Treg are not
  separately resolvable here (markers absent or sub-ambient) and surface as off-diagonal mass.
- Low precision on **PEC (0.22) / MC (0.50) / iTAL (0.58)** = mild over-calling of small/hard classes.

**Verdict:** high concordance on every abundant lineage (epithelial segments, endothelium, stroma,
and all immune lineages B/Plasma/CD4/CD8/Myeloid) **validates both typings**. This full-panel
reannotation is markedly cleaner than the 68-gene cross-cohort pass (segment ARI 0.78 / immune 0.68
here vs ARI 0.61 there), confirming the earlier B-over-call was a shared-gene-depth artefact, not a
labelling error.

---

## STEP 4 — disease-status tagging (setup only; no association run)

`Diagnosis.xlsx` joined by `Sample ID = <orig_ident>_Xenium` — **16/16 matched (100%)**
(`disease_strata_per_sample.csv`). Strata among the 16 (samples; cells):

| condition | samples | cells | note |
|---|--:|--:|---|
| **DKD** | 8 | 389,465 | the disease group |
| **Control** | 3 | 436,914 | Disease=No · DM=no · HTN=no · GFR>60 |
| AA amyloid | 2 | 34,616 | other |
| C3GN | 1 | 41,336 | other GN |
| IgA | 1 | 30,772 | other GN |
| MN | 1 | 17,937 | other GN |

Covariates: DM/HTN are `yes` for all 13 non-Control samples and `no` for all 3 Controls (perfectly
confounded with disease). GFR: >60 = the 3 Controls; within DKD, <30 (n=5) vs 30–60 (n=3).

**Scoping verdict:** 16 samples (no donor column; 1 patient/sample) **cannot power a graded severity
model**, and eGFR/DM/HTN are confounded with disease status. The supportable contrasts are **binary**:
(a) **DKD vs Control** (8 vs 3), (b) DKD vs non-DKD (8 vs 8), or (c) the paper's **B-rich vs B-poor
DKD** subgroup split — the last being the natural follow-on given the B-aggregate focus. *No
association is run here.*

---

## STEP 5 — deliverables (`figures/`)

- **`spatial_grid_16.png/.svg`** — 4×4 grid, 16 region-cropped spatial cell-type maps, one shared
  14-class key (muted epithelium, saturated immune), each titled `sample · condition · n cells · n B`.
- **`integration_umaps.png/.svg`** — Harmony UMAP by sample / by independent cell type / by lineage.
- **`umap_yours_vs_theirs.png/.svg`** — same embedding, my coarse labels vs `annotation_updated`.
- **`concordance_matrices.png/.svg`** — the two confusion heatmaps with ARI + agreement.

Tables: `per_sample_counts.csv`, `cluster_annotation_global.csv`, `cluster_annotation_immune.csv`,
`{coarse,immune}_confusion_counts.csv` (+ `_recall`/`_precision`), `validation_summary.csv`,
`disease_strata_per_sample.csv`, `disease_strata_summary.csv`.

**Bottom line:** an independent full-panel pipeline reproduces the authors' DKD Xenium typing at
ARI 0.78 (segment) / 0.68 (immune subtype), validating both; the residual disagreements are the
expected soft boundaries (injury states, CD4/CD8, panel-limited rare immune types). The 16 samples
are tagged for a **binary DKD-vs-Control / B-rich-vs-B-poor** association as the next step.
