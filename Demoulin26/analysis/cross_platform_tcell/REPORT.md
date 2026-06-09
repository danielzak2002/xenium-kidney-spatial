# CD4/CD8 subtype labels: measured or imputed? (Dumoulin 2026 DKD atlas)

Read-only analysis. Object: 169,775 cells materialized (CD4/CD8 T cells + epithelial reference 'PT'). Platform column `tech`; raw counts from layer `counts`.

## Decision rule

- **T-lineage supported** on a platform if CD3-family markers sit **above the ambient floor** (epithelial reference) in *both* CD4- and CD8-labeled cells.
- **Subtype split supported** if discriminating markers (CD8A/CD8B/CD4) are above ambient **and discriminate** CD8- from CD4-labeled cells: **AUROC well off 0.5**. AUROC ~0.5 with markers at the ambient floor => the subtype labels are **reference-imputed, not measured**.

## Marker measurement (flag #1: structural zeros)

- **CosMx**: structural-zero / absent markers excluded from its tests: TRAC, TRBC2.
- **Xenium**: structural-zero / absent markers excluded from its tests: CD3D, TRAC, TRBC2, NKG7, IL7R.

## Subtype discrimination (lead metric, AUROC CD8-label vs CD4-label)

| marker | CosMx | Xenium |
|---|---|---|
| CD4 | 0.500 | 0.507 |
| CD8A | 0.564 | 0.694 |
| CD8B | 0.515 | 0.532 |
| **composite (CD8A+CD8B−CD4)** | 0.570 | 0.686 |

Lineage-marker AUROC (control — should be ~0.5):

| marker | CosMx | Xenium |
|---|---|---|
| CD3D | 0.508 | nan |
| CD3E | 0.526 | 0.627 |
| CD3G | 0.507 | 0.530 |
| TRAC | nan | nan |
| TRBC2 | nan | nan |

## Per-platform verdict

### CosMx
- T-lineage above ambient (CD3 family): **YES**.
- Subtype discriminating AUROC (CD8A/CD8B/CD4): ['0.56', '0.52', '0.50'] -> subtype split **NOT supported (≈chance -> imputed)**.

### Xenium
- T-lineage above ambient (CD3 family): **YES**.
- Subtype discriminating AUROC (CD8A/CD8B/CD4): ['0.69', '0.53', '0.51'] -> subtype split **SUPPORTED (measured)**.

## Dual-platform shared blocks (CD4 fraction of T cells)

| sample | CosMx CD4frac | Xenium CD4frac | CosMx ratio | Xenium ratio |
|---|---|---|---|---|
| HK2695 | 0.5127 | 0.503 | 1.052 | 1.012 |
| HK2753 | 0.4061 | 0.561 | 0.684 | 1.278 |
| HK3106 | 0.4559 | 0.5964 | 0.838 | 1.478 |
| HK3626 | 0.3046 | 0.566 | 0.438 | 1.304 |

(Region-level only — cells are NOT co-registered across platforms.)

## Cross-platform conclusion

**T-lineage is measured on both platforms; the CD4/CD8 SUBTYPE split is measured on Xenium but reference-imputed on CosMx.**

1. **Lineage (CD3 family) is real on both.** CD3E/CD3G sit well above the PT ambient floor in both CD4- and CD8-labeled cells on each platform (e.g. CosMx CD3E 0.16/0.13 vs ambient 0.042; Xenium CD3E 0.47/0.42 vs ambient 0.002). The "T cell" call is transcript-supported everywhere.

2. **Subtype is measured on Xenium.** CD8A discriminates CD8- from CD4-labeled cells (AUROC 0.69; composite 0.69), with CD8A detected in 41% of CD8-labeled vs 0.1% ambient, and CD4 in 23% of CD4-labeled vs 8% wrong-label / 0.3% ambient. Real per-cell signal drives the split.

3. **Subtype is NOT measured on CosMx — it is imputed.** CD8A barely separates (AUROC 0.56) and **CD4 is at chance (AUROC 0.50)**, against a uniformly elevated ambient floor (every immune marker detected in ~3–7% of PT cells — the CosMx ambient signature). CD8-labeled cells retain a faint residual CD8A bump (0.17 vs 0.046 ambient), so CD8 calls have a weak measured component; **CD4 calls have essentially none**. The CD4/CD8 partition on CosMx therefore reflects the scANVI reference prior, not 1k-plex transcripts.

4. **Dual-platform blocks corroborate.** On the same four tissue blocks, CosMx skews toward CD8 (CD4 frac 0.30–0.51) while Xenium is balanced/CD4-leaning (0.50–0.60). Same tissue, divergent CD4:CD8 — what you expect when one platform's subtype assignment is prior-driven rather than measured (region-level comparison; cells are not co-registered).

5. **Consistent with the cLN 1k finding.** As in the cLN CosMx 1000-plex cohort where CD4/CD8 collapsed, the 1k panel here cannot resolve the subtype while the 5k Xenium can. **Practical rule: trust CD4/CD8 subtype on Xenium; treat CosMx CD4/CD8 as imputed (especially CD4) and do not use it for subtype-resolved cross-platform claims** — pool to "T cell" or restrict subtype analyses to Xenium.

### Flags / caveats
- **No negative-control features in `var_names`** (CosMx Negative01–20 and Xenium NegControlProbe/BlankCodeword were dropped from the published `.h5ad`). The ambient floor is anchored on PT epithelium only; the uniformly ~3–7% CosMx immune-marker detection in PT *is* the ambient readout, but a true neg-probe rate cannot be computed here.
- **AUROC on raw zero-inflated counts is conservative** — heavy ties at zero compress it toward 0.5, so it understates CD4 (which shows a clear detection-rate gap on Xenium: 0.226 vs 0.095). Detection-rate-above-ambient is the more sensitive readout and agrees in direction: Xenium ≫ CosMx for subtype support.
- **Structural-zero exclusions (flag #1):** CD3D, NKG7, IL7R, TRAC, TRBC2 absent on Xenium; TRAC/TRBC2 absent on CosMx — excluded from the respective platform's tests, not read as zero expression.
- Ambient floor from a single epithelial reference (PT); dual-platform comparison is region-level only.
