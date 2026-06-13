# Injury ↔ immune co-localization across nephropathies — compositional addendum

Door-closing exploratory addendum to the DKD consolidation. We asked whether the B-rich DKD subgroup
tracks tubular injury. **It does not track B-lineage specifically — tubular injury co-localizes with
general, myeloid-led immune infiltration.** Builds on `reconcile_extend.py` (validated reannotation
labels; B-lineage=B+Plasma DBSCAN eps=50µm/minPts=20; near≤50µm / far>200µm; within-section
per-compartment z; permutation null positions+count-fixed; injPT genes VCAM1/HAVCR1/PROM1/SPP1).
**Compositional / descriptive; unit = sample/section; honest nulls first-class; controls marked;
no per-gene DE (gated Stage 2). Read-only.** Script `colocalization.py`; tables `coloc_*.csv`;
figures `coloc_fig_*.png`.

## Disease → sample map (all 16) and per-section spillover flags
- **DKD (8):** 1001, 1006, 1008, 1010, 1011, 1012, 1013, HK2695
- **Control (3):** HK2753, HK3106, HK3626 · **AA amyloid (2):** 1004, 1009 · **C3GN (1):** 1007 ·
  **IgAN (1):** 1003 · **MN (1):** 1005   *(IgAN=1003, MN=1005 confirmed)*
- **Spillover (non-PT-ambient Δz > 0):** SPILLOVER in **1005, 1007, 1009, HK2695**; clean elsewhere.
  The **clean, interpretable aggregate sections = 1004, 1008, 1012, HK3626** (1006 is ambient-clean
  but iPT-*depleted* near aggregates). `coloc_spillover_flags.csv`.

## PANEL A — cross-sample specificity (n = 16 samples; descriptive)
**A1 — de-circularized partial (the airtight test).** `total-immune` contains B-lineage, so the
original partial had B partly controlling for itself; re-running against **non-B immune
(myeloid+T+Unresolved)** is the clean test (`coloc_panelA_partials.csv`):

| test | ρ [95% CI] |
|---|---|
| iPT ~ **myeloid frac** | **0.82 [0.46, 0.95]** |
| iPT ~ non-B-immune frac | 0.80 [0.50, 0.94] |
| iPT ~ total-immune frac | 0.79 [0.44, 0.92] |
| iPT ~ B-lineage frac | 0.67 [0.26, 0.87] |
| iPT ~ **B-lineage \| NON-B immune** (de-circ) | **0.13 [−0.37, 0.62]** |
| iPT ~ B-lineage \| total-immune (orig) | 0.07 [−0.46, 0.60] |

B-lineage retains **no signal beyond non-B immune** (partial 0.13, CI spans 0).

**A2 — CLR closure sensitivity** (`coloc_panelA_clr.csv`; coarse 10-part composition, multiplicative
zero-replacement). On CLR coordinates: clr(iPT)~clr(**Myeloid**) **0.39** ≥ clr(iPT)~clr(B) 0.15;
non-B partial clr(iPT)~clr(B)|nonB-immune = −0.05. **Direction and ordering are preserved** (myeloid ≥
B; partial collapses). Closure biases positive part-correlations *downward* (conservative), so the
iPT~immune signals are **if anything understated**, not inflated.

**A3 — the surviving positive (stated in its own right).** **iPT ~ myeloid ρ = 0.82 [0.46, 0.95]:
tubular injury co-localizes with myeloid infiltration** across samples — honestly associational at
n=16, with B-lineage along for the ride, not driving it.

## PANEL B — disease specificity / cross-nephropathy
- **B1** (`coloc_fig_B1_disease.png`): per-sample iPT-fraction rises with both total-immune and
  myeloid fraction; **the three Controls sit together in the low-injury / low-immune corner**
  (squares), the nephropathies spread up-and-right. No disease is an outlier from the shared trend.
- **B2** (`coloc_panelB_spatial.csv`, `coloc_fig_B2_spatial.png`): near B-aggregates, **iPT cell-fraction
  is enriched in 6/11 sections** (perm p<.05: 1004, 1008, 1009, 1012, HK2695, **HK3626 [control]**);
  **myeloid in 4/11** (1001, 1008, 1009, HK3626). **Flat in MN (1005)**; **IgAN (1003) has no
  aggregates**; depleted in the B-richest 1006. The signal is present across DKD / AA-amyloid / C3GN
  **and the one aggregate-bearing control** — i.e. it is **not disease-specific and not B-specific**.
- **B3** (`coloc_fig_B3_gallery.png`): one representative section per disease — injured tubule
  (iPT/iTAL) + immune (myeloid/B/plasma) with B-aggregate hulls.

## PANEL C — which cells participate (compositional, gated)
- **C1 — epithelium near aggregates shifts to INJURED states** (`coloc_panelC1_epi_subtype.csv`,
  `coloc_fig_C1_epi_subtype.png`): **healthy PT is depleted** near aggregates (PT near−far mostly
  negative) while **iPT/iTAL are enriched**; the injured-state gate (iPT > healthy-PT **and** >
  neutral stroma/endo bystander **and** perm p<.05) **passes 6/11** (1004, 1008, 1009, 1012, HK2695,
  HK3626) — **fails in the B-richest 1006** (iPT depleted there). So injured epithelium, not
  epithelium-in-general, sits near immune aggregates — in most but not all sections.
- **C2 — the near-injury infiltrate is MYELOID-skewed** (`coloc_panelC2_immune_injury.csv`,
  `coloc_fig_C2_immune_injury.png`): immune composition within ≤30 µm of an injured-tubule (iPT/iTAL)
  cell vs >100 µm (diffuse injury → distance-to-nearest-injured-cell, not DBSCAN). **Myeloid is
  enriched near injury in 12/16 sections** (perm p<.05; pooled Δ +0.086; myeloid is 36–71% of the
  near-injury immune compartment), while **B-lineage is depleted** near injury (negative in most
  sections, e.g. 1006 −0.40). **Injury attracts myeloid, not B.**

## Bottom line (door-closing)
The injured-PT ↔ immune co-localization is **real, reproducible, and myeloid-led — not B-lineage
specific.** Cross-sample, injury tracks myeloid (ρ 0.82) and general immune infiltration; B-lineage
adds nothing beyond non-B immune (de-circularized partial 0.13, CI spans 0; CLR-robust). Spatially,
injured epithelial states concentrate near immune aggregates (6/11) and the near-injury infiltrate is
myeloid-skewed / B-depleted (12/16). **The B-rich subgroup's significance is architectural /
compositional** (the reproduced subgroup; cross-nephropathy B-dominated/mixed/absent aggregates),
**not a translational tubular-injury role.** No B-lineage→kidney-damage claim is supported.

**Caveats:** n=16 samples / 11 aggregate-bearing (9 with matched near/far); **n=1 per non-DKD
nephropathy** (IgAN/MN/C3GN) — single sections, descriptive only; iPT recall ~0.64 (label-fraction
used; reconcile showed program-score agrees cross-sample); no donor column (donor clustering
uncontrolled); spillover sections flagged; centroids + expression only (no morphology); associational,
not causal. Stage-2 per-gene DE deliberately NOT run.
