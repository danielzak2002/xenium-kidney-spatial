# cLN FastReseg — STAGE 0 gate: **NO-GO** (transcript table absent)

**Goal (intended):** transcript-aware re-segmentation (NanoString FastReseg) of the Danaher cLN
CosMx to recover *T-lineage* cells (not CD4/CD8 subtypes — established unrecoverable on 1k),
attacking the mechanism whereby ambient CD3-family transcripts are mis-assigned to tubular cells
by image-based segmentation, then re-typing with InSituType using a Demoulin CosMx-derived
on-organ lineage reference.

**Outcome:** STAGE 0 gate (a) **FAILS** — FastReseg cannot run on the available data. Reconnaissance
stopped here, by design. Read-only throughout; no raw data modified.

---

## Gate (a) — per-transcript table: **ABSENT (hard blocker)**

FastReseg requires a per-**molecule** table: `target` (gene) + global `x/y` (+`z`) + the original
cell assignment, so it can score each transcript against its cell's profile and trim/reassign
mis-fit molecules. The cLN deposit (Zenodo 13964258 → Figshare; `cleaneddata.RData` /
`data_with_failed_slide.RData`) contains **only cell-level objects**:

| object | class | dim | content |
|---|---|---|---|
| `annot` | data.frame | 532392 × 21 | cell metadata: `cell_ID, fov, Area, Mean/Max PanCK, Mean/Max CD45, Mean/Max DAPI, tissue, class` (cell type), `negmean, totalcounts` |
| `customlocs`, `um`, `viz` | matrix | 532392 × 2 | per-**cell** centroid coordinates (`sdimx, sdimy`) |
| `raw` | dgCMatrix | 957 genes × 532392 cells | cell × gene **count matrix** |
| `clust` | character | 532392 | per-cell cluster id |

No per-transcript / per-molecule export (`*_tx_file.csv`, AtoMx molecule table, `target + x_global +
CellId`) is present anywhere in the deposit; the analysis scripts (`1–8.R`) only ever load the
cell × gene matrix + cell centroids + IF metadata. The published repo README points to the same two
`.RData` as the complete dataset. **The transcript layer FastReseg operates on does not exist here.**

The Demoulin CosMx data (in the combined `.h5ad`) is likewise cell-level (counts layer + centroids,
no molecules), so **neither available dataset supports FastReseg.**

## Gate (b) — panel intersection: **PASS (moot)**

Recorded for completeness; it confirms the *only* blocker is the missing transcripts, not the panel.

- cLN CosMx panel: **957 genes** (1000-plex; negative-control probes already filtered from the matrix).
- cLN ∩ Demoulin (var 5443): **903 genes** (`cln_demoulin_panel_intersection.csv`).
- Present in the intersection: **CD3D, CD3E, CD3G**, CD4, CD8A, CD8B, PTPRC, MS4A1, CD79A, MZB1,
  KRT8, KRT18, EPCAM, CD68, ITGAM, LYZ, NKG7 — i.e. the T-lineage discriminators and the
  structural/epithelial + immune-lineage markers needed for both the re-segmentation profiles and
  the reference.
- Absent from **both** 1k panels: **TRAC, TRBC1/2** (TCR-constant genes; consistent with the
  Demoulin CosMx structural-zeros found in `cd4_cd8_support`). CD3D/E/G carry T-lineage instead.

So if a transcript table existed, lineage profiles could be built on the 903-gene intersection and
the mechanism would be addressable. It does not, so the analysis is blocked at the input stage.

## Decision

**NO-GO for FastReseg** on the current data. The CD3-ambient-reassignment mechanism cannot be
attacked, because trimming/reassigning ambient transcripts requires molecule-level coordinates that
the cLN release does not include. Stages 1–5 (reference profiles → FastReseg → re-type → evaluate)
were not run.

### What would unblock it
The raw CosMx per-FOV transcript files (AtoMx flat-file export / `*_tx_file.csv` with `target`,
`x_global_px`, `y_global_px`, `z`, `CellId`) for the cLN cohort — not part of the public deposit.
Would require obtaining them from the authors or a raw re-export if available.

### In-scope fallback that uses what IS available (not run; offered)
The cLN cell-level objects *do* support a **reference-free contamination assessment** without
re-segmentation: `annot` carries `Mean/Max PanCK` and `Mean/Max CD45` IF per cell plus the cell ×
gene counts, so the **CD3-family count rate in PanCK+ / epithelial cells vs CD45+ immune regions**
(the same primary metric FastReseg's effect would have been judged on) is directly measurable, and
InSituType lineage re-typing with the Demoulin on-organ reference is feasible. **Caveat:** these can
only *quantify* ambient CD3 contamination and re-type existing masks — they cannot recover
mis-segmented T-lineage cells (no transcript trimming), so they bound the problem rather than fix it.
That bound would still sharpen the *exclusion* argument for cLN T cells (we'd be measuring the actual
mechanism's footprint), consistent with the project's standing finding that cLN T-lineage signal is
weak. Flag for a follow-up decision.

---
*Reproduce:* `bash analysis/cln_fastreseg/stage0_gate_check.sh` (read-only; lists cLN objects, extracts
the panel, computes the intersection). Raw data git-ignored.
