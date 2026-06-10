# cLN CD3-family ambient contamination — cell-level bound (close-out)

Read-only. Bounds ambient CD3 mis-assignment to epithelial cells; **does not recover** mis-segmented T cells (no transcript layer — see `analysis/cln_fastreseg`). Completes the cLN T-cell-unreliability account alongside the CD4/CD8-imputed finding (`cd4_cd8_support`).

## Method

- Compartments **orthogonal from IF**: PanCK+ epithelial (Mean.PanCK>P60=899, CD45<P85) vs CD45+ immune (Mean.CD45>P85=265, PanCK<P60); cross-checked vs author labels (42% of author-epithelial are IF-epithelial).
- CD3D/CD3E/CD3G (TRAC/TRBC absent from the 957 panel — confirmed). Formal ambient anchor = per-cell **negmean** (negative-probe mean; cLN retained it).
- cLN resolves only **Treg** among T cells (CD4/CD8 collapsed) — Treg is the 'true T' reference.

## The bound

- **Ambient CD3 false-positive rate in PanCK+ epithelial: ~28.3% (IF) / ~35.2% (author-epithelial)** of epithelial cells carry >=1 CD3 count.
- **T-vs-epithelial CD3 separability: ~2.3x** (Treg vs epithelial sum-mean; IF 1.4x). Modest — CD3 is detectable in epithelium at a substantial fraction of its level in T cells.
- **Epithelial CD3 sits ~7.4x above the negmean floor** (epi CD3 sum-mean 0.560 vs negmean 0.075) — i.e. epithelial CD3 is real spillover/ambient transcripts above pure background, not zero, consistent across the 14 slides (epi CD3+ rate 0.20-0.50).

## Reconciliation & interpretation

Consistent with the earlier observation (~21% epithelial CD3+, ~2x T enrichment): here 35% of epithelial cells are CD3+ with only ~2.3x T-vs-epithelial separation. A ~2.3x separability and double-digit epithelial false-positive rate mean CD3-based T-lineage calls in cLN CosMx are **swamped by ambient spillover** — the image-segmentation mis-assignment mechanism is real and quantitatively large.
**This BOUNDS, does not fix.** It strengthens the *exclusion* of reliable cLN T-cell calls: the contamination is large enough that, without transcript-level re-segmentation (unavailable), T-lineage recall cannot be trusted. Pairs with the CD4/CD8-imputation result for a complete 'cLN T cells are unreliable on this platform/segmentation' account.

## Caveats

- Cell-level only; no transcript trimming possible (FastReseg blocked). Bound, not recovery.
- IF thresholds are percentile cuts (P60/P85) on an epithelial-dominant tissue; author-label compartments agree. 'True T' limited to Treg (the only resolved T label).
- negmean is the mean of negative-probe counts retained in the cLN release.