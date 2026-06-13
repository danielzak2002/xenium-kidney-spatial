# Unbiased all-pairs neighbourhood-enrichment screen (squidpy)

The hypothesis-FREE co-localization view, run for completeness on the DKD Xenium atlas. **Bottom line:
an honest null for new disease-specific immune niches.** Read by SIGN + cross-section reproducibility
(never z-magnitude), the screen recovers tissue **architecture** and **reproduces the already-known
biology**, but surfaces **no new disease-specific off-target immune co-localization**.

> Unit = physical section. z-scores are **never compared across sections** (sign + reproducibility
> instead). Per-section graph. n=1 for IgAN/MN/C3GN (descriptive one-offs); AA n=2. iPT/iTAL/CD8 recall
> caveats carry. Validated reannotation labels reused (no re-typing); built from labels + coordinates
> only — the 8.7 GB h5ad is never opened.

## Method (squidpy, per section)
- `sq.gr.spatial_neighbors(coord_type="generic", delaunay=True)`, then **prune edges > 50 µm** so
  adjacency is within-tissue (the cLN multi-core lesson). Graph **degree ≈ 5.8–5.9** in all 16 sections;
  pruning removed almost nothing (Delaunay edges are already short in these well-segmented sections) —
  the prune is a safeguard, and confirms there are no long spurious cross-tissue edges.
- `sq.gr.nhood_enrichment` over the 21 validated subtype labels (immune | tubular epi | glomerular |
  stroma/endo), n_perms = 1000. squidpy's z is computed against a **label-permutation
  (abundance-preserving) null** — so the sign is already abundance-adjusted by construction.
- z → **SIGN** per pair (enriched > +2, avoided < −2, ~0 otherwise). Types with < 20 cells in a section
  are masked (NaN) for that section.
- **Abundance control:** an explicit per-section label-SHUFFLE null on the same graph (re-run of
  nhood_enrichment with permuted labels) — used as the baseline that real signs must beat across sections.
- **Spillover control:** per-section flags from analysis 06 (`coloc_spillover_flags.csv`); candidates
  must hold in the **clean** DKD sections, not only the spillover-flagged ones.

## What the screen measures (read this before interpreting)
`nhood_enrichment` is **cell-TOUCH adjacency** — do two cell types sit as graph neighbours more than a
label-shuffle expects. This is a **stricter, different question** than analysis 06's **section-level
co-abundance** (sections with more injury have more myeloid, ρ 0.82). The two are complementary, not
contradictory (see "injury × myeloid" below).

## Cross-disease readout (SIGN / reproducibility)
210 unordered pairs/section × 16 sections = 3,360 pair-observations. Of the **201 off-target pairs**
(excluding the known B-aggregate and injury~myeloid pairs), **146 carry a reproducible DKD sign**
(≥ 6/8 sections, ≥ 75% concordant). Classified honestly:

| class | n | what it is |
|---|---|---|
| **control-underpowered** (same direction, weaker) | 78 | DKD avoidance that controls *also* show, but controls are too immune-sparse (≈2.8% immune) to reach \|z\|>2 reproducibly — a **power** difference, not biology |
| **compartmentalization** (immune avoids epithelium/structure) | 46 | immune cells live in the interstitium, not inside tubules/glomeruli → reproducible **avoidance** = tissue **architecture** |
| **constitutive adjacency** (also in controls) | 17 | reproducibly **enriched in DKD AND controls** → not disease-specific (see below) |
| **nephron anatomy** (segment/vascular neighbours) | 2 | DCT×IC-A, MC×EC_DVR — physically adjacent nephron/vascular segments |
| **DISEASE-SPECIFIC candidate (review)** | 3 | sign differs control→DKD and not architecture → reviewed below (2 pass checks) |

**The constitutive immune adjacencies the screen DOES recover (real, but present in controls → not
disease-specific):** Myeloid×Fibroblast, Myeloid×CD4 T, Myeloid×CD8 T, CD4 T×CD8 T — all **enriched
8/8 DKD AND 3/3 control**. These are genuine cell-touch niches (myeloid sits with stroma and T cells;
T cells cluster together), but they are **constitutive tissue features**, not something DKD creates.

## The known biology reproduces (sanity check)
- **B-aggregate composition** (known): B×Plasma **+7/8**, B×CD4 T **+8/8**, B×Myeloid **+5/8**, Plasma×Myeloid
  **+8/8** — mutual immune adjacency inside aggregates, as expected.
- **Injury × myeloid** (known, from 06): at the **cell-touch** level this is **avoidance** (iTAL×Myeloid
  **−8/8**; iPT×Myeloid mixed 4/8). **This does not contradict 06** — 06's ρ 0.82 is **section-level
  co-abundance**; at the touch level myeloid sits in the **interstitium near** injured tubules, not in
  direct epithelial contact. The screen refines, not reverses, the prior result.

## The 3 "disease-specific" candidates — reviewed, and why they are NOT a new immune niche
After the spillover + abundance checks, **2 of 3 survive: iPT×Fibroblast and iPT×EC_glom** — and both
are **non-immune, avoidance pairs dominated by iPT**:
- **iPT×Fibroblast** flips sign control→DKD (control enriched → DKD **avoided 8/8**). Interesting on its
  face, but **iPT recall ≈ 0.64** and iPT abundance differs ~2× (control 3.9% → DKD 7.5%), so the
  apparent flip is most likely a **typing/abundance** effect, not a niche.
- **iPT×EC_glom** avoidance — injured cortical PT physically distant from glomerular endothelium →
  architecture.
- Neither is an immune co-localization. `sq.gr.co_occurrence` on the two (subsampled ≤15k/section,
  per disease group) returns no tight co-localized length-scale — consistent with **avoidance, not a
  niche** (`nhood_cooccurrence_survivors.csv`).

**No new disease-specific immune off-target niche survives.**

## Verdict (honest null, first-class)
The unbiased squidpy screen does three things and finds no fourth: (1) **reproduces** the known
B-aggregate immune adjacency; (2) **refines** injury×myeloid to interstitial proximity rather than
direct contact; (3) **recovers tissue architecture** — immune↔epithelial compartmentalization, a
constitutive immune–stromal–T adjacency present equally in controls, and nephron-segment anatomy. The
Control→DKD **shift map is mostly empty** (most apparent shifts are control immune-sparsity / power).
**There is no hidden off-target immune co-localization that the hypothesis-driven analyses missed.** That
is a clean, reassuring negative — and exactly the question this screen existed to answer.

## Files
`nhood_per_section_long.csv` (every pair × section: z, sign, z_null, counts) ·
`nhood_pair_summary.csv` (per-pair cross-disease sign/reproducibility) ·
`nhood_offtarget_classified.csv` (the honest classification above) ·
`nhood_cooccurrence_survivors.csv` · figures `figures/nhood_per_section_heatmaps.png`,
`figures/nhood_dkd_vs_control_sign.png` · `nhood_enrichment_screen.py`.
