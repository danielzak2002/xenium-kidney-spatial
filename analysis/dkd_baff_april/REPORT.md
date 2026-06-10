# BAFF/APRIL/TACI/BCMA/BAFF-R survival-niche assessment (Demoulin Xenium) — GATED

Read-only, Xenium-only. Standard B-cell-survival immunology. Leads with the usability gate (valuable either way): detection in expected source/target cells vs the PT-ambient floor.

## Usability gate (DECISIVE)

| gene | protein | expected cells | expect detect | PT floor detect | ratio | verdict |
|---|---|---|---|---|---|---|
| TNFSF13B | BAFF | myeloid_stromal_endo | 0.0102 | 0.0038 | 2.68 | **sub-ambient** |
| TNFSF13 | APRIL | myeloid_stromal_endo | 0.0173 | 0.0504 | 0.34 | **sub-ambient** |
| TNFRSF17 | BCMA | plasma | 0.1955 | 0.0003 | 756.27 | **usable** |
| TNFRSF13B | TACI | B | 0.0886 | 0.0002 | 528.6 | **usable** |
| TNFRSF13C | BAFF-R | B | 0.3407 | 0.001 | 343.73 | **usable** |

Usable: ligands [], receptors ['BCMA', 'TACI', 'BAFF-R'].

## Source analysis: STOPPED (both ligands fail)

**BAFF (TNFSF13B) and APRIL (TNFSF13) are sub-ambient on the Xenium 5k panel** in their expected myeloid/stromal/endothelial sources (≈PT-floor detection), **consistent with the cLN CosMx result** where both were sub-ambient. A ligand source/gradient analysis is not interpretable and was not run.

## Stage 1 (passing genes) — in- vs out-aggregate detection

| gene | protein | pop | in-aggregate | out (>200um) | log2 in/out |
|---|---|---|---|---|---|
| TNFRSF17 | BCMA | plasma | 0.2029 | 0.1929 | 0.073 |
| TNFRSF13C | BAFF-R | B | 0.3424 | 0.3402 | 0.009 |
| TNFRSF13B | TACI | B | 0.0734 | 0.1131 | -0.616 |

## Verdict & feasible follow-up

**No detectable BAFF/APRIL survival-ligand niche on this panel** (ligands sub-ambient). Receptor detection (BCMA/TACI/BAFF-R) is reported for completeness but cannot establish a local survival niche without measurable ligand.

**NOTED (not run) feasible alternative inflammatory readout** the paper emphasizes and that is more likely to be on-panel/above-floor: complement (C3, C4A, C8B), TNF-receptor-family, and granzyme programs — a candidate follow-up to characterize the aggregate microenvironment without relying on the low-abundance survival ligands.

## Caveats

- Gate vs PT-ambient floor (Demoulin dropped negprobes); sub-ambient = not interpretable, not biological absence. Xenium-only (CosMx subtype imputed). No patient column.