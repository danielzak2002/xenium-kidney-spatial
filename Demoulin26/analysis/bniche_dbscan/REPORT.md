# DKD vs RCC immunoregulatory aggregate — method-matched (DBSCAN, 1:1)

Read-only. Xenium-only (CosMx subtype imputed). Borrowed RCC pipeline verbatim from `py/phaseB_02_rcc_aggregates.py` + committed eps=50 rebuild (`py/whitepaper_recompute.py`).

## Borrowed parameters (verbatim)

- DBSCAN(**eps=50.0**, **min_samples=20**) on B-cell coords, per section.
- Aggregate footprint = all cells within **R=50.0 um** of any member B cell.
- Enrichment = log2((inside_mean+1e-6)/(section background+1e-6)); per-aggregate Wilcoxon; k/N = aggregates above background. **Unconditional** fractions (state cells / all region cells) — matches RCC (differs from the earlier niche test's conditional FOXP3+-among-CD4+).
- Spatial-unit sanity: median NN cell distance 8.1 um -> microns confirmed; eps transfers.

## Aggregate delineation

- DBSCAN on B cells (independent of author niche) -> **N=32 aggregates** across 10 Xenium sections; median 52 B/agg (RCC median 346 — DKD aggregates are smaller/looser).
- Author-niche overlap: only **6/32 (19%)** overlap 'B predom. Immune ME' (>5% of footprint; mean overlap 6%). The DBSCAN delineation and the coarse author niche are **largely different units** — itself direct evidence that unit-of-analysis drives the comparison.

## Method-matched result

| axis | RCC | DKD Xenium | verdict |
|---|---|---|---|
| Treg ENRICHED in aggregate | 36/37 (log2 1.337) | 28/32 (log2 1.974) | REPLICATE |
| effector-CD8 EXCLUDED | 34/37 (log2 -0.937) | 3/32 (log2 1.992) | DIVERGE |
| Treg collar at margin | n/a | core -0.19 -> margin 0.96 | COLLAR |
| plasma (structural) | log2 0.106 | log2 0.568 | DESCRIPTIVE |

**Treg-like ENRICHED 28/32 (replicates RCC). effector-CD8 is CO-ENRICHED 29/32 (excluded only 3/32) — RCC's cytotoxic EXCLUSION does NOT replicate. Overall: PARTIAL.**

Within-aggregate Treg:effector-CD8 balance: **DKD ≈ 0.0253/0.0286 = 0.88** (co-equal) vs **RCC ≈ 0.189/0.070 = 2.7** (Treg-dominant). DKD B-aggregates admit cytotoxic CD8 in balance with Treg; RCC aggregates are Treg-dominant and cytotoxic-excluding.

## Radial (core -> margin) — the collar test

| ring | Treg-like | eff-CD8 |
|---|---|---|
| core 0-50 | -0.19 | -0.31 |
| margin 50-100 | 0.96 | 0.54 |
| outer 100-150 | 0.27 | -1.12 |

**A Treg collar IS present at the aggregate margin** despite a non-enriched core.

## Conclusion (three-part)

1. **The earlier "Treg absent" DIVERGE was a unit-of-analysis + metric artifact.** Under the *exact* RCC method (DBSCAN eps=50/minPts=20, R=50, unconditional fraction), Treg-like enrichment **replicates** (28/32, log2 +1.97 vs RCC +1.34), with a clear **margin collar** (core −0.19 → margin +0.96, still +0.27 at the outer shell). The prior conditional FOXP3+-among-CD4+ niche-membership test averaged over this collar and missed it.

2. **But the load-bearing RCC discriminator — cytotoxic exclusion — does NOT replicate.** DKD B-aggregates **co-enrich** effector-CD8 (29/32, log2 +1.99) instead of excluding it (RCC 34/37 excluded, log2 −0.94). The within-aggregate Treg:effector-CD8 ratio is **~0.9 in DKD vs ~2.7 in RCC**: RCC aggregates are Treg-dominant and cytotoxic-excluding; DKD aggregates admit cytotoxic CD8 in balance with Treg.

3. **Net: PARTIAL replication.** The DKD B-cell-rich aggregate has a Treg margin-collar (a tertiary-lymphoid-like organization) but lacks the RCC immunoregulatory "Treg-in / cytotoxic-out" architecture. The one geometric hint of segregation is the **outer shell** (Treg +0.27 vs effector-CD8 −1.12) — cytotoxic CD8 drops off beyond the collar while Treg persists.

**Metric caveat (important):** "exclusion" vs whole-section background is sensitive to global cytotoxic burden — RCC tumor effector-CD8 background is **0.135** vs DKD kidney **0.007** (19×). In a cytotoxic-flooded tumor an aggregate reads as "excluding"; in cytotoxic-sparse kidney the same aggregate reads as "concentrating." The composition-independent readouts — the **Treg:effector-CD8 ratio** and the **radial geometry** — are the fair cross-tissue comparison, and both say DKD is *less* immunoregulatory than RCC, not equivalent.

## Caveats

- Borrowed RCC params (eps=50.0, minPts=20, R=50.0) listed above; **N=32 aggregates**.
- **No patient column** -> donor clustering uncontrolled (aggregates within a section share a donor).
- States reconstructed by marker+ within compartment (atlas lacks native Treg/effector-CD8 labels); GNLY structural-zero on Xenium dropped from effector-CD8.
- Magnitudes are log2 vs section background; Xenium-only scope (CosMx subtype imputed).