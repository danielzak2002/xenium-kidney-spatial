# Does the RCC immunoregulatory B-aggregate signature replicate in the DKD B-niche?

Read-only. Niche column `immune_ME`, value `B predom. Immune ME`. Subset: immune cells in niche-containing samples (73,174 cells materialized; full X never loaded).

## Scope — why subtype claims are Xenium-only

Per `cd4_cd8_support`, CosMx 1k CD4/CD8 subtype is **reference-imputed** (CD8A AUROC 0.56, CD4 0.50 at chance), so **all Treg / effector-CD8 / mregDC state claims here use Xenium 5k only**. CosMx contributes a coarse B-predominance consistency check (Step 5) with no subtype claims.

## Power & confounds

- Xenium samples containing the niche: **n = 5** (1006, 1007, 1009, HK2695, HK3626).
- Small cohort — interpret k/n descriptively.
- **No patient column** in obs -> patient-level clustering cannot be controlled (samples may share donors).
- Measured Xenium markers (struct-zeros dropped): Treg=['FOXP3', 'IL2RA', 'CTLA4'], eff-CD8=['GZMB', 'GZMK', 'PRF1'], mregDC=['LAMP3', 'CCR7', 'FSCN1', 'CD274'].

## RCC-axis result (Xenium)

- **Treg-like (FOXP3+/IL2RA+/CTLA4+ among CD4+) ENRICHED in-niche: 1/4 samples.**
- **effector-CD8 (GZMB/GZMK/PRF1/GNLY+ among CD8+) DEPLETED in-niche: 2/4 samples.**
- mregDC-like (LAMP3/CCR7/FSCN1/CD274+ among DC): descriptive only (RCC inconclusive). DC compartments in-niche are tiny (0–1 cells/sample) — uninformative.
- Coarse composition confirms B-dominance; non-B immune pattern in `xenium_coarse_enrichment_per_sample.csv`.

### What replicates vs what diverges
- **Replicates (both platforms):** the niche is genuinely B-dominated (B log2 ~3.9–7.1 in-niche on Xenium; ~5.3–7.8 on CosMx) and **plasma is enriched in-niche** (Xenium 3/5 samples positive incl. HK3626 log2 8.1; CosMx 5/5). So the *structural* B/plasma aggregate is shared.
- **Diverges (the load-bearing RCC axis):** Treg-like is **not enriched — it trends depleted**: FOXP3+/IL2RA+/CTLA4+ fraction among CD4+ is *lower* in-niche than out-niche in 3/4 testable samples (log2 −0.85, −0.25, −0.21, +0.13). effector-CD8 is mixed (2/4 depleted, log2 spanning −0.32…+0.99). The DKD B-niche is a B/plasma aggregate **without** the RCC immunoregulatory Treg(+)/effector-CD8(−) collar. Magnitudes are modest (|log2| < 1) — read as direction, not strong effect.

## Head-to-head verdict

| axis | RCC | DKD B-niche (Xenium) | verdict |
|---|---|---|---|
| Treg ENRICHED in B-niche | 36/37 (0.97) | 1/4 (0.25) | DIVERGE |
| effector-CD8 EXCLUDED from B-niche | 34/37 (0.92) | 2/4 (0.50) | PARTIAL |
| mregDC membership | inconclusive | descriptive (see refined CSV) | DESCRIPTIVE |
| plasma | inconclusive | coarse only | DESCRIPTIVE |

**Overall, scoped to the load-bearing Treg(+)/effector-CD8(-) axis: DIVERGE.**

## Caveats

- Coarse subtype labels -> marker-refined states within compartments (FOXP3+ among CD4+, etc.); not identical to the RCC DBSCAN-aggregate spatial delineation (this is a niche-membership, not a density-aggregate, test). Direction of the Treg(+)/eff-CD8(-) axis is the comparable quantity.
- Small n and no patient control (above); CosMx excluded from all subtype claims.
- Out-niche background = all non-niche cells in the SAME sample (controls sample composition).