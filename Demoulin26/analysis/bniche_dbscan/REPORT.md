# DKD vs RCC immunoregulatory aggregate — method-matched (DBSCAN, 1:1)

Read-only. Xenium-only (CosMx subtype imputed). Borrowed RCC pipeline verbatim from `py/phaseB_02_rcc_aggregates.py` + committed eps=50 rebuild (`py/whitepaper_recompute.py`).

## HEADLINE — burden-corrected differential

The clean, confound-free comparison is the **Treg-vs-cytotoxic bias** per aggregate: Δlog2 = log2(Treg enrichment) − log2(effector-CD8 enrichment), count-pooled, bootstrap CI. This cancels the 19× global cytotoxic-burden difference between RCC tumor and DKD kidney.

- **RCC: Δlog2 = +2.60 [95% CI +2.32, +2.84] — Treg favored over cytotoxic by ~6.1×.**
- **DKD: Δlog2 = +0.24 [95% CI -0.31, +0.72] — ~1.2× (no bias).**
- CIs do NOT overlap: the RCC immunoregulatory Treg-over-cytotoxic bias is **absent** in DKD aggregates. *This single number is immune to the burden confound and to the pseudocount issue.*

## Borrowed parameters (verbatim)

- DBSCAN(**eps=50.0**, **min_samples=20**) on B-cell coords, per section.
- Aggregate footprint = all cells within **R=50.0 um** of any member B cell.
- Enrichment = **count-pooled** log2(Σ state cells / Σ expected cells); per-aggregate k/N and Wilcoxon for direction. **Unconditional** fractions (state cells / all region cells) — matches RCC (differs from the earlier niche test's conditional FOXP3+-among-CD4+). Count-pooling (vs mean-of-per-aggregate-log2) avoids dividing empty single-aggregate rings — kills the −12.5 artifacts.
- Spatial-unit sanity: median NN cell distance 8.1 um -> microns confirmed; eps transfers.

## Aggregate delineation

- DBSCAN on B cells (independent of author niche) -> **N=32 aggregates** across 10 Xenium sections (median 52 B/agg; RCC median 346 — DKD aggregates are smaller/looser).
- Author-niche overlap: only **6/32 (19%)** overlap 'B predom. Immune ME' (>5% of footprint) — the DBSCAN delineation and the coarse author niche are **largely different units**, itself direct evidence that unit-of-analysis drives the comparison.

## Method-matched result

| axis | RCC | DKD Xenium | verdict |
|---|---|---|---|
| ** Δlog2 Treg−eff-CD8 (burden-corrected) ** | +2.60 [+2.32,+2.84] (~6.1x Treg bias) | +0.24 [-0.31,+0.72] (~1.2x) | DIVERGE |
| Treg ENRICHED in aggregate | 36/37 (log2 1.337) | 28/32 (count-pooled log2 +2.08) | REPLICATE |
| effector-CD8 EXCLUDED | 34/37 (log2 -0.937) | 3/32 (count-pooled log2 +1.84) | DIVERGE |
| Treg collar at margin | n/a | core +1.85 -> margin +1.97 | NONE |
| plasma (structural) | log2 0.106 | count-pooled log2 +0.77 | DESCRIPTIVE |

Within-aggregate Treg:effector-CD8 balance (count-pooled inside fractions): **DKD ~1.18× Treg:CD8 vs RCC ~4.8×.** RCC aggregates are Treg-dominant and cytotoxic-excluding; DKD aggregates admit cytotoxic CD8 in balance with Treg.

## Radial (core -> margin) — count-pooled, with bootstrap CI

| ring | Treg-like | eff-CD8 |
|---|---|---|
| core 0-50 | +1.85 | +2.21 |
| margin 50-100 | +1.97 | +1.81 |
| outer 100-150 | +1.82 | +1.38 |

**No Treg collar:** the margin is not enriched over the core.

## Conclusion (three-part)

1. **The earlier 'Treg absent' DIVERGE was a unit-of-analysis + metric artifact.** Under the *exact* RCC method (unconditional, count-pooled), Treg-like enrichment **replicates** (28/32, count-pooled log2 +2.08 vs RCC +1.34) and is enriched uniformly across core->margin->outer (core +1.85 / margin +1.97 / outer +1.82). The prior conditional FOXP3+-among-CD4+ niche-membership test missed it because the conditional metric and the coarse-niche unit differ from the RCC method. (Note: the earlier mean-of-per-aggregate-log2 reported a spurious margin collar; count-pooling removes it.)
2. **But the RCC discriminator does NOT replicate.** DKD aggregates **co-enrich** effector-CD8 (excluded only 3/32, count-pooled log2 +1.84). The burden-corrected **Δlog2 = +0.24 [-0.31,+0.72] (DKD) vs +2.60 [+2.32,+2.84] (RCC)** — non-overlapping CIs. RCC favors Treg over cytotoxic ~6.1×; DKD ~1.2× (no bias).
3. **Net: PARTIAL.** The DKD B-cell-rich aggregate concentrates both Treg-like and cytotoxic CD8 (an immune-dense B/plasma aggregate) but lacks the RCC immunoregulatory 'Treg-in / cytotoxic-out' architecture. The only spatial structure is a mild **cytotoxic-core gradient**: effector-CD8 is highest at the core (+2.21) and falls outward (+1.38) while Treg stays flat — the opposite of a Treg collar.

**Metric note:** raw 'exclusion vs whole-section background' is confounded by global cytotoxic burden (RCC tumor eff-CD8 bg 0.135 vs DKD kidney 0.007, ~19×). The Δlog2 differential and the radial geometry are the burden-immune readouts, and both show DKD is *less* immunoregulatory than RCC.

## Caveats

- Borrowed RCC params (eps=50.0, minPts=20, R=50.0) listed above; **N=32 aggregates**.
- **No patient column** -> donor clustering uncontrolled (aggregates within a section share a donor).
- States reconstructed by marker+ within compartment (atlas lacks native Treg/effector-CD8 labels); GNLY structural-zero on Xenium dropped from effector-CD8.
- Magnitudes are log2 vs section background; Xenium-only scope (CosMx subtype imputed).