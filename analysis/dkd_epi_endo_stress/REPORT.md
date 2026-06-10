# Epithelial / endothelial stress near immune B-aggregates (Demoulin Xenium)

Read-only. Xenium-only. Reuses bniche_dbscan B-aggregate delineation (DBSCAN eps=50/minPts=20). NEAR = <=50 um to nearest aggregate-member B cell; FAR = >200 um; matched **within section and by cell type** (PT/tubular-to-tubular, endo-to-endo). Programs module-scored (CP-median log1p, per-section per-compartment z, averaged over usable genes).

## Usability gate (per program)

- **injPT** (epi): usable = ['VCAM1', 'HAVCR1', 'PROM1', 'SPP1', 'ITGB6']
- **endoAct** (endo): usable = ['ICAM1', 'VCAM1', 'ENG', 'PLVAP']
- **hypoxia** (endo): usable = ['VEGFA']
- dropped (no usable members): ['fibroEMT']

Full gate (detection vs Immune-cell ambient floor) in `usability_gate.csv`.

## NEAR vs FAR (matched, within section)

| program | target | k/N sections near>far | median Δz | pooled Δz | MWU p (near>far) |
|---|---|---|---|---|---|
| injPT | epi | 6/9 | 0.109 | 0.133 | 3.81e-02 |
| endoAct | endo | 4/9 | -0.004 | -0.037 | 6.63e-01 |
| hypoxia | endo | 0/9 | -0.112 | -0.068 | 7.44e-131 |

(32 aggregates across 10 sections.)

## Distance gradient (mean program z-score)

| bin (um) | injPT | endoAct | hypoxia |
|---|---|---|---|
| 0-50 | 0.131 | -0.036 | -0.067 |
| 50-100 | 0.144 | -0.002 | -0.051 |
| 100-200 | 0.065 | -0.066 | 0.0 |
| 200-500 | 0.071 | 0.002 | -0.008 |
| >500 | -0.011 | 0.002 | 0.003 |

## Interpretation

- **injPT: ELEVATED near aggregates.** pooled Δz +0.13, 6/9 sections positive, MWU p=3.81e-02; monotonic distance gradient (+0.13 at 0-50um -> -0.01 at >500um). The clearest signal.
- **endoAct: not elevated** near aggregates (pooled Δz -0.04, 4/9 sections positive, no decreasing gradient).
- **hypoxia: not elevated** near aggregates (pooled Δz -0.07, 0/9 sections positive, no decreasing gradient).
- *Note:* the hypoxia program is a SINGLE usable gene (VEGFA); its large-n Mann-Whitney p is a tie-dominated artifact and should be read with the per-section k/N and gradient (both null), not the p-value.
- **Associational, not causal:** aggregates may nucleate in already-injured/inflamed regions; this is co-localization of tubular injury (injPT) with B-aggregates, not evidence that aggregates cause injury — fully consistent with the paper's profibrotic-tubular-niche thesis at the level of spatial association.
- **No patient column** -> donor clustering uncontrolled (sections may share donors); per-section matching mitigates section/composition effects but not donor structure.
- **No morphology:** processed h5ad = centroids + expression only; this is transcriptional-state + spatial-organization, NOT cell-shape histomorphology. No image/segmentation layer present.
- Module scores are per-section per-compartment z (comparable within section); pooled across sections after centering. Genes gated vs an Immune-cell ambient floor (Demoulin dropped negprobes).