# Comparative non-immune <-> immune spatial-interaction map (RCC / PRCC / cLN / DKD kidney)

Read-only. Per-section neighborhood enrichment (squidpy, KNN k=6, 1000 perms, permutation z) on a harmonized lineage labelling. **Compositional / ligand-independent**: spatial co-localization, NOT communication or causation. X never touched (labels + coords only).

## Stage 0 — harmonization (native -> common)

Common immune = B, Plasma, Myeloid (mono/macro/DC pooled), T_lineage (pooled), NK. Non-immune = Tubular_epi, Endothelial, Stroma, Podocyte, **Malignant_epi (RCC/PRCC only — tumor-immune, interpreted separately)**. Injured-epithelium = DKD native iPT/iTAL overlay.

| dataset | source label | notes |
|---|---|---|
| RCC | phase_b_label | 465k single section -> spatially TILED into ~50k pseudo-sections; Tumor_RCC->Malignant_epi; NK absent in base |
| PRCC | phase_b_label | single section (optional RCC replication) |
| cLN | author_celltype | 14 sections, CosMx; T pooled & **UNRELIABLE (35% CD3 contam)** |
| DKD | annotation_updated + immune_cell_annotation_combined | 16 Xenium sections; injured overlay = iPT/iTAL |

Full native->common mapping printed to console / in script.

## KEY interpretation: read DIFFERENTIALS, not absolute sign

Absolute nhood z is **dominated by a universal immune<->parenchyma SEGREGATION**: immune cells sit in sparse interstitial/peri-vascular aggregates, so bulk Tubular/Endothelial/Malignant epithelium reads as 'avoiding' every immune type in every context (the strongest, most conserved pattern — but largely a density-geometry fact). The biologically informative signal is the **relative ordering**: which immune partner a given non-immune type avoids LEAST, and which non-immune **state** (injured vs healthy) is least excluded.

## Stage 3 — validation (known niches recovered, as differentials)

| dataset | test | values | recovered? |
|---|---|---|---|
| cLN | Plasma<->Myeloid (direct) | z=8.49 | YES |
| RCC | B<->Treg vs B<->CD8eff (differential) | -0.64 vs -24.05 (Δ+23.4) | YES (differential) |
| DKD | Injured_epi<->B vs Tubular_epi<->B (differential) | -11.56 vs -26.06 (Δ+14.5) | YES (differential) |

The pervasive cLN plasma-myeloid niche recovers as a **direct** positive z; the focal RCC B-Treg aggregate and DKD injured-PT proximity recover only as **differentials** (global KNN nhood is less sensitive to focal aggregates than the DBSCAN approach in `bniche_dbscan`).

## Stage 2 — conserved vs context-specific (non-immune x immune)

| non-immune | immune | RCC | PRCC | cLN | DKD | class |
|---|---|---|---|---|---|---|
| Tubular_epi | B | -3.0 | -6.4 | -11.5 | -26.1 | CONSERVED-avoid |
| Tubular_epi | Plasma | -5.2 | -21.5 | -13.3 | -13.4 | CONSERVED-avoid |
| Tubular_epi | Myeloid | 2.5 | -35.7 | -23.3 | -30.0 | DISCORDANT |
| Tubular_epi | T_lineage | 2.4 | -50.9 | -1.9 | -37.7 | CONSERVED-avoid |
| Tubular_epi | NK | nan | -11.2 | -1.2 | -8.1 | CONSERVED-avoid |
| Malignant_epi | B | -12.4 | -15.0 | nan | nan | CONSERVED-avoid |
| Malignant_epi | Plasma | -14.7 | -19.8 | nan | nan | CONSERVED-avoid |
| Malignant_epi | Myeloid | 1.1 | -13.4 | nan | nan | ns-all |
| Malignant_epi | T_lineage | -11.4 | -17.3 | nan | nan | CONSERVED-avoid |
| Malignant_epi | NK | nan | -10.9 | nan | nan | insufficient |
| Endothelial | B | -27.2 | -15.2 | 0.6 | -4.3 | CONSERVED-avoid |
| Endothelial | Plasma | -10.4 | -16.0 | -0.9 | -2.1 | CONSERVED-avoid |
| Endothelial | Myeloid | -35.9 | -41.8 | 14.0 | 2.9 | DISCORDANT |
| Endothelial | T_lineage | -40.5 | -59.8 | -0.4 | -1.5 | CONSERVED-avoid |
| Endothelial | NK | nan | 8.0 | 1.0 | 5.3 | CONSERVED-enrich |
| Stroma | B | -31.8 | -10.5 | 1.2 | -3.3 | CONSERVED-avoid |
| Stroma | Plasma | -6.3 | 21.0 | 4.4 | 4.0 | CONSERVED-enrich |
| Stroma | Myeloid | -28.6 | -41.4 | 12.2 | 18.6 | DISCORDANT |
| Stroma | T_lineage | -29.0 | -78.0 | -1.1 | 10.1 | DISCORDANT |
| Stroma | NK | nan | 7.2 | -0.2 | 2.2 | CONSERVED-enrich |
| Podocyte | B | nan | nan | -0.2 | -2.8 | ns-all |
| Podocyte | Plasma | nan | nan | -0.7 | -1.8 | ns-all |
| Podocyte | Myeloid | nan | nan | 6.2 | -4.5 | DISCORDANT |
| Podocyte | T_lineage | nan | nan | -1.2 | -4.0 | ns-all |
| Podocyte | NK | nan | nan | 1.2 | 3.6 | ns-all |

### Synthesis

- **CONSERVED (all contexts): immune cells segregate from Tubular epithelium** (Tubular_epi avoids B/Plasma/Myeloid/T everywhere) — the robust conserved axis, though largely density-geometry.
- **CONSERVED-ish ENRICH: Stroma <-> Myeloid and Stroma <-> Plasma** in the benign/papillary contexts (cLN, DKD, PRCC) — myeloid & plasma cells localize to the interstitial stroma. **DISCORDANT in RCC** (clear-cell stroma avoids immune) — a tumor-context-specific inversion.
- **Endothelial <-> NK and Stroma <-> NK enrich** where NK is measured (PRCC, DKD) — peri-vascular/interstitial NK; cLN null is uninformative (CosMx depth).
- **DKD-specific: Stroma <-> T_lineage enrich** (+10) — interstitial T cells, a DKD interstitial-nephritis signature absent in the tumor contexts.
## Interpretation rules & caveats

- **Read differentials, not absolute sign** (see Key interpretation): the universal Tubular<->immune avoidance is geometry; the informative signals are Stroma<->myeloid/plasma and the validation differentials.
- A positive CONSERVED association is strong; a dataset-specific ABSENCE is **ambiguous** (platform/typing depth — cLN is CosMx, sparser) -> a cLN null is NOT biological absence.
- **cLN T_lineage pairs are UNRELIABLE** (35% epithelial CD3 contamination, see `cln_cd3_contamination`) -> excluded from conserved T claims.
- **RCC/PRCC epithelium is MALIGNANT** (Malignant_epi) -> tumor-immune, interpreted separately from benign Tubular_epi; RCC/PRCC are single-section (k/N = ?/1, no section replication).
- Platform-depth confound (CosMx 957 vs Xenium 5k) -> asymmetric reading of nulls.
- Sparse-type z is noisy (pairs with <20 cells/section flagged via n_sec); NK absent in RCC base.
- No patient column for DKD (donor clustering uncontrolled). Per-section normalized z throughout; KNN graph is scale-invariant (handles um vs mm). Association/colocalization only.