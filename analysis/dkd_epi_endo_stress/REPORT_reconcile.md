# Injured-PT ↔ B-lineage — reconcile & extend under post-BAFF rigor

Re-examines the prior `dkd_epi_endo_stress` result (author labels + B-only aggregates: **injPT
near>far Δz +0.13, 6/9 sections, MWU p=0.038**, monotonic gradient) using the **validated
reannotation labels** (incl. `iPT`) and **B-lineage = B+Plasma** DBSCAN aggregates (eps=50/minPts=20,
inherited). injPT program genes inherited verbatim (VCAM1/HAVCR1/PROM1/SPP1/ITGB6). **Unit of
replication = sample; descriptive throughout; read-only.** Script `reconcile_extend.py`; tables
`reconcile_*.csv`; figures `reconcile_fig_primary.png`, `reconcile_fig_spatial.png`.

**Task B status:** committed **and pushed** — current commit `3129358` on `origin/main`, path
`analysis/dkd_epi_endo_stress/`. (`2955126` is its *pre-history-rewrite* SHA; content byte-identical.)

## Usability gate (injPT genes, PT-compartment detection vs Immune ambient floor — this cohort)
| gene | PT detect | immune floor | ratio | usable |
|---|--:|--:|--:|:--:|
| VCAM1 | 0.163 | 0.023 | 7.2× | ✓ |
| HAVCR1 | 0.079 | 0.003 | 24.0× | ✓ |
| PROM1 | 0.103 | 0.005 | 23.1× | ✓ |
| SPP1 | 0.684 | 0.045 | 15.1× | ✓ |
| **ITGB6** | 0.021 | 0.003 | 6.6× | **✗ (<3% detect in PT)** |

**4 of 5 pass** under reannotation PT labels; ITGB6 drops (2.1% PT detection < 3% gate). Scored on
VCAM1/HAVCR1/PROM1/SPP1. (Prior kept all 5 on the broader author-`epi` compartment.)

## PRIMARY — cross-sample correlation (n = 16 samples), per-sample table
| sample | cond | iPT frac (epi) | injPT prog (PT) | B-lin frac | n_agg | tot-imm frac | myeloid frac |
|---|---|--:|--:|--:|--:|--:|--:|
| 1001 | DKD | 0.097 | 0.156 | 0.0072 | 1 | 0.100 | 0.053 |
| 1003 | IgA | 0.295 | 0.308 | 0.0082 | 0 | 0.185 | 0.108 |
| 1004 | AA amyloid | 0.317 | 0.311 | 0.0261 | 1 | 0.224 | 0.085 |
| 1005 | MN | 0.365 | 0.338 | 0.0311 | 3 | 0.254 | 0.131 |
| 1006 | DKD | 0.289 | 0.265 | 0.0292 | 5 | 0.217 | 0.100 |
| 1007 | C3GN | 0.208 | 0.185 | 0.0342 | 4 | 0.262 | 0.082 |
| 1008 | DKD | 0.205 | 0.139 | 0.0132 | 2 | 0.119 | 0.052 |
| 1009 | AA amyloid | 0.282 | 0.316 | 0.0562 | 4 | 0.220 | 0.067 |
| 1010 | DKD | 0.122 | 0.214 | 0.0051 | 0 | 0.161 | 0.087 |
| 1011 | DKD | 0.208 | 0.125 | 0.0104 | 0 | 0.132 | 0.053 |
| 1012 | DKD | 0.081 | 0.223 | 0.0131 | 3 | 0.155 | 0.054 |
| 1013 | DKD | 0.219 | 0.299 | 0.0196 | 1 | 0.235 | 0.076 |
| HK2695 | DKD | 0.084 | 0.280 | 0.0203 | 4 | 0.143 | 0.058 |
| HK2753 | Control | 0.053 | 0.172 | 0.0006 | 0 | 0.033 | 0.017 |
| HK3106 | Control | 0.065 | 0.170 | 0.0004 | 0 | 0.031 | 0.022 |
| HK3626 | Control | 0.053 | 0.211 | 0.0081 | 3 | 0.053 | 0.026 |

**Spearman ρ (bootstrap 95% CI):**
| x | y | ρ | 95% CI |
|---|---|--:|---|
| iPT frac | B-lineage frac | **0.67** | [0.26, 0.87] |
| iPT frac | B-aggregate count | 0.18 | [−0.33, 0.64] |
| iPT frac | **total-immune frac** | **0.79** | [0.44, 0.92] |
| iPT frac | **myeloid frac** | **0.82** | [0.50, 0.96] |
| iPT frac | **B-lineage frac \| total-immune** (partial) | **0.07** | [−0.46, 0.60] |
| injPT prog (PT) | B-lineage frac | 0.58 | [0.14, 0.86] |
| injPT prog (PT) | total-immune frac | 0.69 | [0.28, 0.90] |
| injPT prog (PT) | B-lineage frac \| total-immune (partial) | 0.05 | [−0.49, 0.68] |

**SPECIFICITY VERDICT — NOT B-lineage-specific.** iPT burden correlates with B-lineage (ρ 0.67), but
it correlates **at least as strongly with total-immune (0.79) and myeloid (0.82)**, and the
**partial correlation iPT~B-lineage controlling for total-immune collapses to 0.07** (CI spans 0).
Per the pre-registered rule: report as **"tubular injury tracks general immune infiltration"**
(myeloid-led), not a B-lineage-specific association. The hard iPT-label fraction and the continuous
injPT program score **agree** (both ~0.6 with B-lineage; both partial-collapse). iPT~B-aggregate
count is weak/null (0.18).

## SECONDARY — spatial near/far injPT across the 9 aggregate-bearing sections (gated like BAFF)
| section | cond | n near/far | injPT Δz | ctrl-gene floor | non-PT ambient | nCount Δz | injPT−ctrl | iPT-frac near/far | perm p |
|---|---|---|--:|--:|--:|--:|--:|---|--:|
| 1004 | AA amyloid | 51/6483 | −0.06 | +0.59 | −0.32 | +0.74 | −0.66 | 0.57 / 0.31 | 0.000 |
| 1005 | MN | 41/3300 | −0.11 | +0.07 | +0.38 | +0.29 | −0.18 | 0.39 / 0.37 | 0.376 |
| 1006 | DKD | 37/11446 | −0.75 | −0.29 | −0.40 | −0.24 | −0.46 | 0.14 / 0.29 | 1.000 |
| 1007 | C3GN | 21/5945 | +0.10 | +0.28 | +0.11 | −0.05 | −0.19 | 0.16 / 0.21 | 0.858 |
| 1008 | DKD | 108/9332 | +1.26 | +0.90 | −0.09 | +0.65 | +0.37 | 0.81 / 0.18 | 0.000 |
| 1009 | AA amyloid | 125/1760 | +0.14 | +0.22 | +0.25 | +0.06 | −0.08 | 0.39 / 0.26 | 0.001 |
| 1012 | DKD | 66/3622 | −0.28 | +0.24 | −0.00 | −0.24 | −0.52 | 0.42 / 0.07 | 0.000 |
| HK2695 | DKD | 140/15027 | +0.10 | +0.14 | +0.42 | +0.05 | −0.04 | 0.21 / 0.08 | 0.000 |
| HK3626 | Control | 131/41229 | +0.22 | +1.01 | −0.06 | −0.17 | −0.79 | 0.25 / 0.05 | 0.000 |

**Summary:** injPT near>far **5/9**; **exceeds the control-gene floor only 1/9**; iPT-fraction
permutation p<0.05 in **6/9**; non-PT ambient also rises in 4/9. **Pooled injPT Δz = +0.069 (median
+0.097); pooled control-gene floor = +0.352; pooled non-PT ambient = +0.032.**

## RECONCILIATION vs the prior +0.13 / 6-of-9
- **Program-intensity spatial reading ATTENUATES and does NOT survive the control-gene floor → retract.**
  The injPT program near/far drops to pooled +0.07 (median +0.097, 5/9) from the prior +0.13/6-of-9,
  and — decisively — the **neutral control-gene panel rises *more* near aggregates (+0.35) than the
  injPT program (+0.07)**; the program exceeds its own inflation floor in only **1/9** sections. So
  the modest near-aggregate "injury" intensity is **generic transcript inflation in the dense
  peri-aggregate field, not an injPT-specific gradient** — the same spillover/inflation mechanism that
  retracted the stromal-BAFF claim.
- **What survives is compositional, not intensity, and is non-specific.** The **iPT cell-fraction is
  enriched near aggregates in 6/9 sections** (within-section permutation, positions+count fixed) —
  injured-PT-*labeled* cells co-localize with B-aggregates. But the PRIMARY analysis shows this
  co-localization **tracks general immune infiltration** (total-immune ρ 0.79 ≥ B-lineage 0.67;
  partial collapses to 0.07), so it is **not B-lineage-specific**.
- **Hard label vs continuous score — flagged divergence (iPT recall ~0.64).** Cross-sample, the two
  agree. **Spatially they diverge:** the iPT *fraction* is enriched near aggregates (6/9 perm) while
  the per-cell injPT *program intensity* is not (fails the control floor, 1/9). Read as **more
  injured-PT-labeled cells near aggregates (composition), not a per-cell injury-program gradient.**
- **B-rich {1006, HK2695} (descriptive color only).** The two B-rich DKD sections **disagree**: 1006
  shows injPT *depleted* near aggregates (Δz −0.75; iPT-fraction near 0.14 < far 0.29, perm p=1.0),
  HK2695 shows a weak positive (Δz +0.10). No coherent B-rich enhancement — reported, not leaned on.
- **Section-set note:** current aggregate-bearing 9 = {1004,1005,1006,1007,1008,1009,1012,HK2695,HK3626}
  vs prior {…,1013,…}; 1013 drops (single small aggregate, <10 near PT cells) and 1005 (MN) enters
  under B-lineage (B+Plasma) aggregates.

## Bottom line
Under post-BAFF rigor the injured-PT ↔ B-lineage association **is real as co-localization but is
neither B-lineage-specific nor a program-intensity gradient.** It reduces to: **tubular injury (iPT)
tracks general, myeloid-led immune infiltration across samples**, and **injured-PT-labeled cells sit
near immune aggregates compositionally**, while the **per-cell injury-program near/far rise does not
exceed neutral-gene transcript inflation** (prior +0.13 retracted). Honest attenuation, consistent
with the BAFF result that "near a dense aggregate" signals must clear the control-gene/specificity
bar before a niche reading.

**Caveats:** n=16 samples / 9 aggregate-bearing sections — underpowered, descriptive; no patient
column (donor clustering uncontrolled); iPT recall ~0.64 (both label-fraction and program-score
reported); processed h5ad = centroids + expression only (no morphology); associational, not causal.
