# DKD Xenium — B-lineage disease analysis (subgroup reproduction, references, BAFF/APRIL)

Built on the **validated reannotation** (`analysis/dkd_xenium_reannotation/`, segment ARI 0.78 /
immune ARI 0.68). **B-LINEAGE = B + Plasma.** Goal: reproduce the paper's B-rich DKD subgroup
*within* DKD, place IgAN/MN and other non-DKD as individual references, and re-test the BAFF/APRIL
ligand axis under a correctly *conditioned* design. Raw object read-only; pure science.
Reproduce: `substrate.py` → `baff_april.py` → `figures.py`.

> **Caveats up front (apply throughout).** DKD n=8 vs Control n=3, non-DKD references n=1–2 →
> **descriptive / underpowered**; any p-value is a flag, not inference. **One patient per sample**
> (no donor replication). Architecture work inherits the reannotation's typing caveats: **CD8
> recall 0.58** (author NK folds into CD8) and **iPT/iTAL** soft boundaries — immune-aggregate
> composition is robust at the B/Plasma/CD4/CD8/Myeloid level, not finer.

---

## STEP 0 — groups

Comparison groups: **DKD (8)** {1001,1006,1008,1010,1011,1012,1013,HK2695} and **Control (3)**
{HK2753,HK3106,HK3626}. Individual references (NOT pooled): **AA amyloid (2)** {1004,1009},
**C3GN (1)** {1007}, **IgAN (1)** {1003}, **MN (1)** {1005}.

## STEP 1 — per-sample B-lineage substrate

Per sample: B-lineage fraction, and B-aggregate burden from **DBSCAN on B-lineage cells
(eps=50 µm, minPts=20)** — aggregate count + aggregate cells, normalised per 10k cells and per mm²
(tissue area = occupied 50 µm grid bins). Full table: `per_sample_substrate.csv`.

| sample | cond | cells | B-lin % | B-lin/10k | n_agg | agg cells/10k | frac B-lin in agg | author B-predom ME |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| **1006** | DKD | 45,609 | 2.92 | 292 | 5 | **202** | 0.69 | 797 |
| **HK2695** | DKD | 151,101 | 2.03 | 203 | 4 | **115** | 0.57 | 1576 |
| 1013 | DKD | 28,261 | 1.96 | 196 | 1 | 7 | 0.04 | 0 |
| 1008 | DKD | 27,912 | 1.32 | 132 | 2 | 40 | 0.30 | 0 |
| 1012 | DKD | 39,105 | 1.31 | 131 | 3 | 29 | 0.22 | 0 |
| 1011 | DKD | 27,876 | 1.04 | 104 | 0 | 0 | 0.00 | 0 |
| 1001 | DKD | 29,126 | 0.72 | 72 | 1 | 9 | 0.12 | 0 |
| 1010 | DKD | 40,475 | 0.51 | 51 | 0 | 0 | 0.00 | 0 |
| HK2753 | Control | 161,350 | 0.06 | 6 | 0 | 0 | 0.00 | 0 |
| HK3106 | Control | 171,526 | 0.04 | 4 | 0 | 0 | 0.00 | 0 |
| HK3626 | Control | 104,038 | 0.81 | 81 | 3 | 34 | 0.42 | 121 |
| 1003 (IgAN) | IgA | 30,772 | 0.82 | 82 | 0 | 0 | 0.00 | 0 |
| 1005 (MN) | MN | 17,937 | 3.11 | 311 | 3 | 68 | 0.22 | 0 |
| 1004 | AA amyloid | 24,449 | 2.61 | 261 | 1 | 73 | 0.28 | 0 |
| 1009 | AA amyloid | 10,167 | 5.62 | 562 | 4 | 335 | 0.60 | 52 |
| 1007 | C3GN | 41,336 | 3.42 | 342 | 4 | 110 | 0.32 | 271 |

## STEP 2 — within-DKD B-rich/B-poor reproduction (primary), validated vs authors

The paper's "B cell-rich subgroup" is defined by **B-aggregates/TLS**, not raw B fraction (which
mis-ranks the large HK2695 section). On **B-aggregate burden** the DKD distribution has a clean ~3×
gap between the top-2 (1006 = 202, HK2695 = 115 per 10k) and #3 (1008 = 40); cutting in that gap
(≥75) gives **B-rich DKD = {1006, HK2695}**, B-poor = the other six. (The data alone is between a
top-2 and a top-4 break; the authors' niche pins it to the top-2 — see validation.)

**Validation vs authors:** the authors' `B predom. Immune ME` niche is present in DKD **only** in
1006 (797 cells) and HK2695 (1,576) — zero in the other six DKD. Our B-rich membership is therefore
**100% concordant (8/8)** with the authors' B-predominant niche. Reproduction confirmed.
(`dkd_subgroup_split.csv`, `figures/dkd_subgroup.png`).

## STEP 3 — non-DKD references (descriptive) + DKD vs Control

**DKD vs Control** (`dkd_vs_control_test.csv`, `figures/dkd_vs_control.png`): B-lineage fraction is
higher in DKD (median 1.32% vs 0.06%, Mann-Whitney **p=0.049**, *descriptive*); B-aggregate burden
trends higher but n.s. (p=0.35) — driven by the two B-rich samples. Two of three Controls have
**zero** B-aggregates (HK3626 carries a few). Individual non-DKD references overlaid as labelled
points.

**IgAN (1003) — FOCUS.** B-lineage 0.82%, **zero B-aggregates** (none reach DBSCAN density). On the
B-rich/B-poor spectrum IgAN sits firmly at the **B-poor / no-TLS** end here — its B-lineage
inflammatory pattern looks **distinct from B-rich DKD**: scattered B/plasma, no organised B
follicular structure in this section. (n=1, exploratory.)

**MN (1005) — FOCUS.** B-lineage 3.11% (high) and it **does form B-aggregates** (3; 68 agg cells/10k),
but their **composition is distinct from B-rich DKD**: in/around the aggregate (≤50 µm) MN is
B 11% · **Plasma 9% · Myeloid 12%** · CD4 T 25% · CD8 T 15% — a *plasma/myeloid/T-mixed* aggregate,
versus B-rich DKD which is strongly **B-dominated** (HK2695 B 36% / Plasma 2% / Myeloid 4%; 1006
B 38% / Plasma 2% / Myeloid 5%). So MN resembles B-rich DKD in *forming aggregates* but differs in
*type* — less a B-follicular niche, more a plasma/myeloid-leaning infiltrate. (n=1, exploratory.)

**Lighter references.** AA amyloid: 1009 is strikingly B-lineage-high (5.6%) with the heaviest
aggregate burden of any sample (335/10k) — an amyloid-associated B/plasma infiltrate; 1004 is milder
(2.6%, 73/10k). C3GN (1007): B-lineage 3.4%, 110 agg cells/10k, and the authors flag a modest
B-predom niche (271) — the one non-DKD that most resembles the B-rich pattern.
(`aggregate_composition.csv`.)

## STEP 4 — per-participant B-lineage figures (all 16)

`figures/b_lineage_gallery_16.png` — one region-cropped panel per sample, B (blue) + Plasma (orange)
highlighted on grey tissue, **B-lineage aggregate hulls** drawn (DBSCAN eps=50/minPts=20), titled
`sample · group · B-lineage% · agg cells`. Ordered **B-rich DKD → IgAN, MN (placed adjacent for
direct contrast) → B-poor DKD → Control → other one-offs**. Visual read: B-rich DKD show compact B
follicular aggregates; IgAN has none; MN has looser mixed aggregates; Controls are essentially bare.

## STEP 5 — BAFF/APRIL conditioned re-assessment (verdict)

Both ligands and all three receptors are **on-panel and measured** on Xenium (not CosMx-only
structural-zeros): BAFF/TNFSF13B, APRIL/TNFSF13, BCMA/TNFRSF17, TACI/TNFRSF13B, BAFF-R/TNFRSF13C.
Detection conditioned on **producers** (not a tissue-wide average) and on **space**
(`baff_detection_by_celltype.csv`, `baff_near_far.csv`, `baff_receptors.csv`, `figures/baff_april_panel.png`).

**BAFF (TNFSF13B) → GO (overturns the prior global null).**
- **Producer-enriched:** myeloid detect **3.1% vs 0.39% epithelial floor = 7.9×** (global rate was
  0.62% — the average had washed it out). Fibroblast 0.88% (2.2×).
- **Receptor axis real and aggregate-concentrated:** BAFF-R on B cells **35.7%**, and **elevated in
  B-aggregates (35.4% near vs 22.1% far, 1.6×)**.
- ~~**Peri-aggregate ligand source = stroma:** stromal BAFF enriched near B-aggregates (2.1–2.5×).~~
  **⚠ RETRACTED** — this localized-stromal claim did **not** survive an ambient/spillover stress-test
  (`REPORT_baff_ambient.md`): the near/far BAFF rise is **non-specific** (non-producer epithelium 2.15×
  and endothelium 2.86× rise just like stromal 2.44×, density flat) → a **local spillover field**, not
  stromal production. Myeloid BAFF is tissue-wide (0.66–1.14× near).
- Verdict (revised): a measurable BAFF **producer** signal exists (**myeloid tissue-wide**, 7.9× the
  epithelial floor; **anchored** — see `REPORT_baff_receptor.md`: cell-intrinsic 24× per-transcriptome,
  reproducible 16/16) and B/plasma cells **constitutively** express the receptors; but there is **no
  localized peri-aggregate stromal production** AND **no aggregate-specific receptor concentration** —
  the BAFF-R "1.6×" did **not** survive nCount/control-gene rigor (`REPORT_baff_receptor.md`). Net: a
  real tissue-wide myeloid BAFF producer with constitutive B-lineage receptors, **no localized niche**.

**APRIL (TNFSF13) → NO-GO (global null confirmed under conditioning).**
- **Not producer-specific:** myeloid 5.4% is **below** the epithelial floor 6.8% → APRIL is broad/
  ambient, not a myeloid-niche signal. Its receptor BCMA on plasma is real (17%) but **not**
  aggregate-concentrated (0.37× near). No spatial niche signal.

---

## Bottom line

The B-rich DKD subgroup **reproduces** from our validated labels (B-rich = {1006, HK2695}, 100%
concordant with the authors' B-predom niche). IgAN **lacks** B-aggregates (distinct from B-rich DKD);
MN **forms** aggregates but of a **plasma/myeloid-mixed**, non-B-follicular type. The conditioned
ligand re-test **splits the prior global null**: **BAFF has a real tissue-wide myeloid producer
signal** (7.9× the epithelial floor) and B cells express BAFF-R, whereas **APRIL remains an ambient
null**. NOTE: a follow-up stress-test (`REPORT_baff_ambient.md`) **retracted the localized
peri-aggregate stromal-BAFF claim** as ambient/spillover — keep myeloid tissue-wide BAFF, drop the
spatial-niche-production reading. All comparisons are descriptive (n=8 vs 3; references n=1–2; one patient
per sample) and inherit the CD8/iPT typing caveats.
