# Three-cohort all-Xenium integration — reliable-gene set (updated)

Redo of the three-cohort integration on the **reliably-detected-in-all-3 gene set** (above the
per-cohort ambient floor, not just non-zero), applied to **both** the integration features **and**
the lineage definitions; with a **direct cohort-difference test** added to Readout A and Readout B
reframed. (First pass on the full 123-gene name-intersection is preserved in
`REPORT_firstpass_123gene.md`.) Cohorts: RCC_big (kidney_10x), RCC_figshare (record 25685961),
DKD (Demoulin Xenium). Read-only; pure science.

## Stage 0b — reliable gene set

For each cohort, per-gene detection is computed **per lineage**; a gene is **reliable** if its
detection in its best lineage is **≥ 2.5 %** AND **≥ 1.5× the cross-lineage median** (the median
across lineages is the non-expressing-reference "ambient" proxy — no negative-control probes
survive in the built matrices). Keeping genes reliable in **all three** cohorts:

- **Reliable-in-all-3 = 68 genes** (of 123). 55 dropped (`reliable_genes.txt`, `reliable_drops.csv`).
- **DKD drives most drops** — its Xenium 5k per-cell detection is sparse, so endothelial/NK/
  cytotoxic-adjacent genes thin out: VWF, NKG7, GNLY, IL2RA, CCL5, CD3D, ANGPT2, CXCL9/10 … all
  fail on DKD. A few fail on figshare (e.g. FOXP3 is borderline there but clears at 2.5 %).
- The lineage-aware floor correctly keeps **rare-but-specific** markers (FOXP3 in T cells: 2.6–9.9 %
  across cohorts vs 0.4–3 % in epithelium/endothelium) while dropping uniform-ambient genes.

## Stage 1b — re-type + re-integrate on the reliable set

Per-cell lineage re-derived with **only reliable markers** (uniform definitions). Availability:

| lineage | reliable markers | status |
|---|---|---|
| B | MS4A1, CD79A | robust |
| Plasma | MZB1, TNFRSF17 | robust |
| T | CD3E, CD3G | robust |
| Myeloid | CD68, CD14, CD163, AIF1, ITGAX | robust |
| Endothelial | PECAM1, EGFL7 | robust (VWF dropped on DKD) |
| Epithelial | EPCAM, CDH1 | robust |
| Stroma | PDGFRA, PDGFRB, ACTA2 | robust |
| **NK** | **KLRD1 only** | **weakened** (GNLY, NKG7 dropped on DKD) |
| **Treg gate** | **FOXP3, CTLA4** | (IL2RA dropped on DKD) |
| **cytotoxic gate** | **CD8A, GZMK** | robust |

Re-integration (Harmony on the 68-gene PCA, harmonypy called directly) **re-validates**: the two
RCC cohorts merge, lineages resolve (`INTrel_umap`), and the per-cohort dot-plot is marker-faithful
in every cohort (`INTrel_dotplot`). Lineage counts: Epithelial 561 k, Stroma 371 k, Myeloid 370 k,
Endothelial 339 k, T 291 k, B 257 k, Plasma 93 k, NK 93 k; among T: Treg 50 k, cytotoxic 120 k.

## Readout A — replication + DIRECT cohort-difference test

Per-cohort count-pooled **Δlog₂ = log₂(Treg enrich) − log₂(cytotoxic enrich)** per B-aggregate,
bootstrap 95 % CI (`readoutA_reliable.csv`); and the **difference (cohort − DKD)** bootstrapped by
resampling aggregates within cohort (`readoutA_difference.csv`, `READOUT_A_reliable`).

| cohort | aggregates / sections | Δlog₂ [95 % CI] |
|---|---|---|
| RCC_big | 77 / 1 | +0.49 [−1.01, +1.99] |
| RCC_figshare | 361 / 9 | +0.44 [−0.40, +1.27] |
| RCC_pooled | 438 / 10 | +0.55 [−0.44, +1.53] |
| DKD | 61 / 13 | +0.26 [−1.00, +1.54] |

**Direct cohort-difference test (does it exclude zero?):**

| contrast | difference in Δlog₂ [95 % CI] | excludes 0? |
|---|---|---|
| RCC_pooled − DKD | **+0.28 [−0.14, +0.72]** | **no** |
| RCC_big − DKD | +0.23 [−0.19, +0.67] | no |
| RCC_figshare − DKD | +0.18 [−0.29, +0.63] | no |

### Replication vs separation — kept distinct (do not conflate)

- **Replication (integrated, reliable-set):** all three RCC estimates are **positive and consistent
  in direction** (+0.44 to +0.55) across two independent cohorts (different patients, panels, labs),
  with DKD lower (+0.26). The Treg-favoring *direction* reproduces.
- **But the integrated difference is NOT statistically resolved:** every cohort−DKD contrast has a
  95 % CI that **includes zero**. On the reliable 68-gene space the test is **underpowered** — the
  Treg gate is reduced to FOXP3/CTLA4, DKD aggregates are few (61) with thin immune-marker detection
  (FOXP3 ~0.5 % of DKD cells), and RCC_big is a single section (very wide CI).
- **Separation evidence comes from the native-label analysis**, not this one: `bniche_dbscan`
  (full native panels, native effector-CD8 / Treg labels) gives RCC Δ ≈ **+2.6** vs DKD **+0.24**
  with **non-overlapping** bootstrap CIs. That is where the magnitude/separation lives.
- **Bottom line:** the integration **replicates the direction** across two RCC cohorts; the
  **separation** is established at native depth. The reliability filter buys cross-cohort
  comparability at the cost of power, and honestly does not by itself resolve RCC≠DKD.

## Readout B — reframed (panel limitation, not a cross-context readout)

Usability gate (endothelial-cell detection) is **decisive**: **ANGPT2, CXCL9, CXCL10 are
sub-ambient on the DKD Xenium 5k panel** (0.3–1.3 %) and were dropped from the reliable set — the
vascular/inflammatory-activation axis is **not measurable cross-cohort**. (VWF is structural-zero
on DKD.) Where measurable on the two RCC cohorts, endothelial-activation near vs far from
B-aggregates is **null / slightly depleted** (RCC_big −0.31, RCC_figshare −0.58); HLA-DRA (the only
DKD-usable inflammatory marker) shows no gradient. **Readout B is reported only as a within-RCC
null with a panel-coverage limitation; it is not a cross-context result.**

## Caveats

- **RCC_big is a single section** → very wide aggregate-bootstrap CI; lean on **figshare** for
  multi-sample RCC (9–10 sections). **DKD aggregates are sparse** (61) with thin immune detection.
- **Reliable 68-gene set** reduces depth and **weakens the Treg gate (FOXP3/CTLA4)** and **NK
  (KLRD1 only)**; this is the cost of cross-cohort reliability and the main reason the integrated
  difference is underpowered.
- **Native segmentation on all three** (uniform); ProSeg not applied (DKD has no transcripts) —
  an **RCC-only ProSeg sensitivity check is still pending**.
- **No DKD patient column** (donor clustering uncontrolled); RCC epithelium malignant (excluded
  from any epithelial-stress cross-claim). **Association, not causation.**
