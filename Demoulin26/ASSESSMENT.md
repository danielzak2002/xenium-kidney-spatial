# Dumoulin et al. 2026 DKD spatial atlas — read-only assessment

## Source & object
Dumoulin et al., *Spatial atlas of diabetic kidney disease reveals a B cell-rich subgroup* (Nature 2026; Zenodo 19868428). Read-only inventory via `assess_dataset.py` (anndata backed='r', X never materialized). Raw `.h5ad`/`.xlsx` are git-ignored.

**Object:** `spatial_adata_xenium_cosmx_zenodo.h5ad` — **4,337,862 cells × 5,443 genes**, sparse `X` (float32) + `counts` layer; `obsm` = X_umap / spatial / spatial_fov; `uns` = orig_ident_colors only (no scVI/scANVI latent, no niche color map).

## Samples & design
- 60 distinct `orig_ident`; **48 CosMx + 16 Xenium = 64 platform-samples** (matches the paper's 64). 4 IDs measured on both platforms (HK2695, HK2753, HK3106, HK3626). No patient/donor column in `.obs` (paper: 58 patients) — patient mapping must come from the diagnosis sheet.
- Diagnosis sheet (64 rows): conditions DKD 29, Control 15, DM 6, DM/HTN 4, + 8 rarer GN/amyloid groups; covariates Age/Sex/Race/DM/HTN/GFR/treatments present.
- Cells: 3,386,822 CosMx + 951,040 Xenium; median 57992 cells/sample (range 1328–267652).

## Annotations present (validated reference — YES)
- **Cell types:** `annotation_updated` (20 kidney epithelial/stromal/endothelial types; immune as one coarse `Immune` bucket).
- **Immune subtypes:** `immune_cell_annotation_combined` (12): Macro, **CD8+, CD4+, B, Plasma**, Neutrophil, NK, Baso_Mast, cDC, pDC, mDC. No Treg label (T resolved only to CD4+/CD8+).
- **Niches:** `niches_annotation_based` (12 tissue niches incl. `Immune niche`) + microenvironment columns `immune_ME`/`Immune_ME_20um` (incl. **`B predom. Immune ME`**) and `iPT_iLOH_ME(_20um)`.

## Cross-platform immune resolution (the key finding)
Unlike the cLN 1k atlas (CD4/CD8 collapsed on CosMx), here **CD4+ and CD8+ are separately annotated on BOTH platforms**, as are B and Plasma. Per-platform immune cell counts (`immune_cell_annotation_combined` × `tech`):

| immune type | CosMx | Xenium |
|---|---:|---:|
| Macro | 83143 | 47311 |
| CD8+ | 27120 | 19673 |
| CD4+ | 21842 | 21140 |
| B | 7848 | 8422 |
| Neutrophil | 6497 | 1376 |
| Plasma | 2696 | 2726 |
| NK | 1279 | 1982 |
| Baso_Mast | 255 | 367 |
| cDC | 117 | 102 |
| pDC | 1 | 111 |
| mDC | 2 | 0 |

Caveat: labels exist on both platforms, but CosMx 1k immune calls inherit the platform's ambient/low-plex limits (cf. the cLN finding) — treat CosMx CD4/CD8 as author-asserted, to be validated, not assumed clean.

## B / plasma / TLS-like niche
- **`B predom. Immune ME`** = the B-cell-predominant microenvironment: ~83% B cells (B 2748, CD4+ 362, CD8+ 113, Macro 55, Plasma 22, NK 5) — a B-aggregate / TLS-like niche.
- Present in **9/60 samples**, predominantly **DKD** (then C3GN, Control, AA amyloid, DM/HTN); both platforms (Xenium-leaning). This is the spatial substrate of the paper's 'B cell-rich subgroup'.

## Data-quality flags / mismatches
- `.var` has **no panel/platform column** (only `n_cells`); the object holds 5,443 genes with no CosMx-1k vs Xenium-5k mask. CosMx cells are effectively restricted to the ~1k CosMx panel (structural zeros elsewhere), so any **cross-platform gene-level analysis must be limited to the shared/overlap panel — which is NOT pre-computed here**.
- Diagnosis `Disease` column has inconsistent casing (Yes/YES, No/NO); `Sex` missing for 229,669 cells; large `Unknown` in immune (non-immune cells) and ME columns (ME assigned only to immune-region cells).
- 60 `orig_ident` vs paper's 58 patients (replicate samples e.g. HK2844_2, HK3035_2, HK3531_2; 4 dual-platform IDs).

## Candidate next analyses (extend the existing project — NOT run here)
1. **Transfer the B-aggregate / plasma–myeloid niche workflow** (squidpy nhood-enrichment, DBSCAN aggregate delineation, in/out-of-niche composition) onto `B predom. Immune ME` across the 9 B-rich DKD samples — a direct extension of the RCC B–Treg and cLN plasma–myeloid analyses to a third disease context.
2. **Cross-platform validation of CD4/CD8 typing**: use the dual-platform IDs (HK2695/HK2753/HK3106/HK3626) to test whether CosMx-1k CD4/CD8 calls agree with paired Xenium-5k — the cleanest test of the cLN T-cell limitation.
3. **B-rich vs B-poor DKD contrast**: compare niche composition, plasma load, and myeloid state between the 9 B-niche samples and the remaining DKD samples (the paper's subgroup split), reusing the per-slide methodology.
4. **Gene-usability gate on both panels**: compute per-platform detection/enrichment for B-survival and recruitment genes (MS4A1/MZB1/TNFRSF17 vs BAFF/chemokines) to define the analyzable overlap panel before any modeling.
5. **Niche conservation across disease**: test whether the B-predominant niche composition (B–CD4 core, CD8-low) is conserved vs the RCC immunoregulatory aggregate and cLN plasma–myeloid niche.
