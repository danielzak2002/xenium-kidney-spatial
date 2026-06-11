# Three-cohort all-Xenium integration — readouts A & B

Three Xenium cohorts integrated on the **123-gene shared space** (`three_cohort_assessment/`):
**RCC_big** (kidney_10x, 1 section), **RCC_figshare** (record 25685961, transcript tables →
cell×gene built natively, 10 samples), **DKD** (Demoulin Xenium subset, 16 sections).
~2.4 M cells (figshare 1.47 M · RCC_big 456 k · DKD 445 k). Read-only; pure science.

**Framing:** Harmony standardizes *typing*, not the comparison. Both readouts are computed
**per cohort** on each cohort's native cells + coordinates with the harmonized labels.
**Native segmentation on all three** (ProSeg not applied — DKD has no transcript table, so
uniform segmentation is kept; ProSeg is a separate RCC-only sensitivity, see caveats).

## Stage 1 — integration for typing (validation)

Per-cell lineage is assigned by **uniform shared-marker definitions** in every cohort (scalable
to 2.4 M cells); Harmony+Leiden on a balanced subsample (60 k/cohort) **validates** that this
typing is consistent across cohorts. Lineage counts: Epithelial 614 k, Stroma 386 k, Myeloid
372 k, Endothelial 327 k, T 266 k, B 209 k, NK 113 k, Plasma 89 k. Among T cells: Treg (FOXP3+)
59 k, cytotoxic (CD8A/GZMK+) 107 k.

Validation passes: (a) the two RCC cohorts merge in the Harmony embedding (batch baseline); (b)
the per-cohort marker dot-plot shows each lineage defined by the right markers in **every** cohort
(`INT_dotplot`); (c) the Treg:cytotoxic composition difference between RCC and DKD survives
labeling (`INT_composition`). Figures: `INT_umap`, `INT_dotplot`, `INT_composition`.

## Readout A — immunoregulatory aggregate differential (THE headline)

Per section: DBSCAN B-aggregates (ε=50, minPts=20); per aggregate, count-pooled
**Δlog₂ = log₂(Treg enrichment) − log₂(cytotoxic enrichment)** (inside vs section background),
bootstrap 95 % CI. (`readoutA_differential.csv`, `READOUT_A`.)

| cohort | aggregates / sections | Δlog₂ [95 % CI] | fold | Treg-favoring? |
|---|---|---|---|---|
| **RCC_big** | 70 / 1 | **+0.49 [+0.39, +0.59]** | ~1.4× | **yes (CI > 0)** |
| **RCC_figshare** | 95 / 5 | **+0.67 [+0.49, +0.82]** | ~1.6× | **yes (CI > 0)** |
| DKD | 62 / 13 | +0.22 [−0.16, +0.59] | ~1.2× | no (CI crosses 0) |

**Both independent RCC cohorts show a Treg-favoring bias in B-aggregates (CIs above zero); DKD is
indistinguishable from zero.** Because the two RCC cohorts are different patients, a different
custom panel, and a different lab, their agreement **de-confounds disease from batch** — the
immunoregulatory bias tracks the tumor context, not a single dataset. (Magnitudes are smaller
than the native-label ccRCC analysis [Δ≈+2.6] because the 123-gene harmonized typing is coarser
and cytotoxic here = CD8A/GZMK-gated; the **direction and the two-cohort replication** are the
contribution.)

## Readout B — endothelial / inflammatory stress near B-aggregates (cross-context axis)

**Usability gate (detection in endothelial cells), per cohort** (`readoutB_usability.csv`):

| gene | RCC_big | RCC_figshare | DKD |
|---|---|---|---|
| ANGPT2 (endo-activation) | 0.20 ✓ | 0.40 ✓ | **0.013 ✗** |
| CXCL9 | 0.13 ✓ | 0.15 ✓ | **0.005 ✗** |
| CXCL10 | 0.10 ✓ | 0.06 ✓ | **0.003 ✗** |
| HLA-DRA | 0.76 ✓ | 0.78 ✓ | 0.72 ✓ |
| PECAM1 (endo id) | 0.92 | 0.88 | 0.76 |
| VWF (endo id) | 0.62 | 0.77 | **0.00 ✗ (struct-zero)** |

**Key limitation:** ANGPT2 / CXCL9 / CXCL10 are **sub-ambient on the DKD Xenium 5k panel** — the
vascular-activation axis is only *measurable* on the two RCC cohorts; on DKD only HLA-DRA passes.

Endothelial-cell module score NEAR (≤50 µm) vs FAR (>200 µm) from B-aggregates, matched within
section (`readoutB_near_far.csv`, `readoutB_gradient.csv`, `READOUT_B`):

| cohort | endo-activation (ANGPT2) Δ(near−far) | inflammatory (CXCL9/10/HLA-DRA) Δ |
|---|---|---|
| RCC_big | −0.31 | −0.29 |
| RCC_figshare | −0.58 (2/5 sec +) | +0.03 (3/5 sec +) |
| DKD | n/a (markers gated out) | −0.08 (HLA-DRA only, 9/11 sec) |

**Result: no elevation of endothelial/inflammatory activation near B-aggregates** — if anything,
endothelial-activation is *lower* near aggregates in RCC, and the inflammatory module is null.
Combined with the DKD panel gap, **Readout B does not support a cross-context vascular-stress
niche** around B-aggregates on these panels; it is reported as an honest null with the coverage
caveat. (The within-RCC depletion may reflect that B-aggregates sit in immune-dense interstitium
displacing vasculature — associational, not tested here.)

## Caveats

- **123-gene shared space** (reduced depth vs native panels). **DKD effectively 104/123** — 19 of
  the shared genes are CosMx-only and **structural-zero on DKD Xenium** (incl. VWF, NKG7, GNLY,
  PTPRC, CCL5); the readout markers (FOXP3/IL2RA/CTLA4, CD8A/GZMK, HLA-DRA, PECAM1) survive, but
  ANGPT2/CXCL9/CXCL10 do not on DKD (Readout B limitation above).
- **Native segmentation on all three** — ProSeg was *not* applied, to keep segmentation uniform
  (DKD has no transcript table). ProSeg re-segmentation of the two RCC cohorts (which have
  transcripts) is a separate **RCC-only sensitivity check**, not part of this uniform comparison.
- **Harmony for labels/typing only.** Per-cell labels are uniform marker definitions; Harmony+Leiden
  validate cross-cohort consistency. Batch = disease is managed by the **two-RCC de-confounding**
  (both RCC cohorts independently show the Readout-A bias) + the per-cohort dot-plot.
- **No DKD patient column** (donor clustering uncontrolled). **RCC epithelium is malignant** —
  excluded from any epithelial-stress cross-cohort claim (EMT/fibrosis markers excluded by design).
- **Association, not causation**; figshare cell×gene matrices were built from the public transcript
  tables (native cell_id) and are not committed.

## Bottom line

- **Readout A replicates the immunoregulatory (Treg-favoring) aggregate bias across two independent
  RCC cohorts vs DKD** — the strongest cross-cohort result, de-confounding disease from batch.
- **Readout B is panel-limited and null**: the vascular/inflammatory-activation markers are
  sub-ambient on DKD Xenium, and where measurable (RCC) activation is not elevated near aggregates.
