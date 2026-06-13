# Per-sample cell-type composition by disease group

Top-level layer for the walkthrough. **Reuses the validated reannotation labels** (summary 01:
951,040 cells / 16 samples); no typing was re-run. Descriptive only — **n = 1 for IgAN, MN, C3GN**
(AA amyloid n = 2; DKD n = 8; Control n = 3). **No statistics.** Fractions are compositional
(sum-to-one); reported descriptively with a CLR sensitivity small-multiple. Raw data read-only.

## Groups (confirmed against the object)
Control {HK2753, HK3106, HK3626} · DKD {1001, 1006, 1008, 1010, 1011, 1012, 1013, HK2695} ·
IgAN {1003} · MN {1005} · AA amyloid {1004, 1009} · C3GN {1007}.

## What the labels support
- **Coarse lineage** (`my_lineage`, 4 parts): Epithelial / Immune / Stroma / Endothelial.
- **Immune subtypes** (`my_label`): Myeloid / CD4 T / CD8 T / B / Plasma. **NK and DC are not
  separately typed** in the validated labels — not reported as their own classes.
- **Epithelial subtypes** (`my_coarse`): PT / iPT / TAL / iTAL / DCT / PC-CNT / IC / Podo / PEC
  (in the CSV; the two committed dotplots are coarse + immune).
- `Unresolved` (5,399 cells, 0.6%) is excluded from fraction denominators.

## Descriptive read (group medians, % of section cells)

| group | Epithelial | Immune | Stroma | Endothelial | Myeloid | T (CD4+CD8) | B | Plasma | healthy-PT | iPT |
|---|---|---|---|---|---|---|---|---|---|---|
| Control | **74.4** | **2.8** | 10.3 | 11.0 | 2.2 | 1.0 | 0.05 | 0.02 | **37.1** | 3.9 |
| DKD | 49.9 | 14.2 | 21.0 | 15.0 | 5.6 | 6.7 | 0.85 | 0.43 | 14.2 | 7.5 |
| IgAN\* | 47.4 | 18.2 | 19.5 | 15.0 | 10.9 | 6.5 | 0.47 | 0.35 | 11.1 | **14.0** |
| MN\* | 38.3 | 24.0 | 24.5 | 13.2 | **13.3** | 7.5 | 1.33 | **1.83** | 6.5 | 14.0 |
| AA amyloid | 43.7 | 21.9 | 17.5 | 17.0 | 7.6 | 10.1 | **3.51** | 0.62 | 12.2 | 13.1 |
| C3GN\* | 37.3 | 25.8 | 24.2 | 12.7 | 8.2 | **14.2** | 2.54 | 0.90 | 7.1 | 7.8 |

\* n = 1 (single section — anecdotal, descriptive).

**Patterns (descriptive, not tested):** Control is strongly epithelial-dominant (74%) with a near-bare
immune compartment (2.8%); every disease group shows **epithelial loss with immune + stromal
expansion**, and the proximal-tubule shift is consistent throughout — **healthy PT collapses
(37%→6–14%) as injured iPT rises (3.9%→7–14%)**. Among the single-section references, **MN** is the
most immune/plasma-infiltrated glomerular disease (immune 24%, plasma 1.8% — the plasma-skew echoed in
the B-lineage analysis), **C3GN** is the most T-skewed (T 14%), and **AA amyloid** carries the highest
B-cell fraction (3.5%). DKD spans a wide spread (the dotplots show its 8 sections fanning across each
panel), and the IgAN/MN/C3GN single dots sit within or just beyond that spread — i.e. **none of the
single-section references is clearly outside the DKD range except on its own standout axis** (MN plasma,
C3GN T, AA B). The **CLR sensitivity row tracks the raw-fraction row** in every panel — the ordering is
not an artifact of closure.

## Files
`composition_by_group.csv` (tidy: group · sample · resolution · cell_type · fraction · n_cells) ·
`figures/composition_by_group_coarse.png` · `figures/composition_by_group_immune.png` ·
`composition_by_group.py`.
