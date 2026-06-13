# Making the neighbourhood-enrichment heatmap legible (explanatory panels)

Pedagogical companion to `REPORT_nhood_enrichment.md` — **descriptive only, no new tests.** Shows a bench
biologist WHERE the z comes from and reconciles the injury×myeloid scale puzzle, worked on one
representative DKD section (**1006**; n = 45,189 cells, mean graph degree **5.86**, 132,324 ≤50 µm edges).
Feeds the organisation section of the combined writeup.

## Panels (figures committed in `figures/`)
1. **`nhood_explain_buildup.png` — graph → counts → z.** (a) The Delaunay graph pruned ≤50 µm on a region
   crop: cells coloured by type, edges drawn — "neighbours = who touches whom." (b) `sq.gr.interaction_matrix`
   raw A–B neighbour-pair **counts** (log scale) — the numerator, before any statistics. (c) the
   `nhood_enrichment` **z-matrix** — the same counts normalised against label-shuffle. Same colour key as the walkthrough.
2. **`nhood_explain_permnull.png` — the normalisation, worked twice.** Histograms of label-shuffled
   expected counts (n = 1,000) with the observed count marked:
   - **Myeloid × Fibroblast (enriched):** observed 6,127 vs null mean 4,691 → **z = +22.6**.
   - **iTAL × Myeloid (avoided):** observed 903 vs null mean 1,169 → **z = −8.3**.
   "z = how many SDs the observed touch-count sits above/below random labelling." Read SIGN/scale, not z-as-importance.
3. **`nhood_explain_cooccurrence.png` — the SCALE reconciliation.** `sq.gr.co_occurrence` ratio vs radius:
   **iTAL × Myeloid starts ~0.33 at the touch scale (z < 0) and rises through 1.0 by the interstitial
   radius (~40–50 µm), plateauing ≈1.08.** So cell-TOUCH avoidance (squidpy z) and section-level
   co-abundance (analysis 06, ρ 0.82) are **not in conflict — they are different length-scales**: myeloid
   sits in the interstitium *near* injured tubules, not in epithelial contact. **B × Plasma** is shown as a
   tight-positive contrast (ratio ~2.6–3.2, high AT the touch scale and staying high — a real adjacency niche).
4. **`nhood_explain_ripley.png` (optional warm-up).** B-cell Ripley's L vs a CSR envelope: at the aggregate
   scale (<300 µm) B's L sits far above the random expectation → B cells are clustered, not random. (Plotted
   to 400 µm; beyond ~1.5 mm L saturates as a finite-window artefact, not dispersion.)

## One-line teaching points
- The heatmap is just **(observed touches) normalised by (random-labelling expectation)** — panels 1→3.
- **z sign = touch enrichment/avoidance**; magnitude is within-section only, never compared across sections.
- **Touch ≠ co-abundance:** the same pair can be touch-avoided yet co-abundant at interstitial range —
  `co_occurrence` shows the length-scale; this is the honest reconciliation with analysis 06.

Script: `nhood_enrichment_explain.py` (reuses the per-section graph recipe + validated labels from
`nhood_enrichment_screen.py`; coords + labels only — the 8.7 GB h5ad is never opened).
