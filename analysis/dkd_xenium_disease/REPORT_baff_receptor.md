# BAFF axis — receptor nCount-control + myeloid producer anchor

Completes the BAFF assessment with the same nCount rigor that retracted the stromal-ligand claim
(`REPORT_baff_ambient.md`): test the **receptors** (BAFF-R/TNFRSF13C, BCMA/TNFRSF17, TACI/TNFRSF13B)
for genuine aggregate upregulation, and **anchor** the myeloid BAFF (TNFSF13B) producer. APRIL is
absent on this tissue, so BAFF is the relevant ligand for BCMA/TACI. Reuses the `dkd_xenium_disease`
object (validated reannotation labels + per-cell receptor/ligand counts + nCount + a neutral
control-gene panel TPT1/PPIA/YWHAZ/TMSB10/UBB + myeloid markers). Read-only. Scripts:
`baff_receptor_anchor.py`; figures `figures/baff_receptor_control.png`, `figures/baff_myeloid_anchor.png`.

Aggregate "inside" = the target cell is a DBSCAN B-lineage aggregate member (eps=50/minPts=20);
"outside" = same cell type, not an aggregate member. Tested in the four aggregate-bearing sections
(B-rich DKD 1006, HK2695; C3GN 1007; AA amyloid 1009).

---

## PART A — receptors, within-target-cell-type, nCount + control-gene controlled

For each receptor we compare detection in its target cell type inside vs outside aggregates, and ask
whether any rise **exceeds** (i) nCount inflation and (ii) the same inside/outside ratio for neutral
control genes in those same cells. **Pooled** (`receptor_aggregate_control.csv`):

| receptor (target) | n in/out | det in % | det out % | in/out ratio | nCount ratio | norm-expr ratio | control-gene ratio | verdict |
|---|---|--:|--:|--:|--:|--:|--:|---|
| BAFF-R (B) | 3295 / 1424 | 36.6 | 26.7 | 1.37 | 1.12 | 1.29 | **1.28** | nCount/ambient inflation |
| TACI (B+Plasma) | 3458 / 2930 | 7.4 | 6.7 | 1.10 | 1.12 | 0.99 | 1.52 | null |
| BCMA (Plasma) | 163 / 1506 | 14.7 | 12.1 | 1.22 | 0.92 | 1.30 | **1.34** | nCount/ambient inflation |

**No receptor is upregulated beyond the generic inflation baseline.** BAFF-R's in-aggregate rise
(norm-expr 1.29×) **equals** the neutral control-gene inflation (1.28×) — i.e., aggregate-resident B
cells detect *modestly more of everything*, and BAFF-R rides that, not exceeds it. BCMA's rise
(1.30×) is **below** its control baseline (1.34×). TACI is flat. **Per-section** the only receptor with
any real signal is BAFF-R, and it is **inconsistent**: it exceeds control in 1006 and 1007 but sits
at/below control in HK2695 and 1009 (and the strongest single result, BCMA 2.8× in 1007/C3GN, does
not reproduce in the B-rich DKD sections). **The earlier "BAFF-R aggregate-elevated 1.6×" claim does
not survive** — it was transcript-count inflation.

What is real: B/plasma cells **constitutively** express the receptors at substantial baseline
(BAFF-R ~27–37% of B cells, BCMA ~12–17% of plasma) — just **not aggregate-concentrated**.

## PART B — myeloid BAFF producer anchor

| check | result | reading |
|---|---|---|
| **nCount control** | myeloid BAFF per-transcriptome **23.7×** epithelium; myeloid have *fewer* counts (median 113 vs 419) | not a count artifact — opposite direction |
| **nCount-matched fold** | **23.9×** (raw detection fold was 7.9×) | matching removes epithelium's count advantage → fold rises |
| **Identity** | BAFF⁺ myeloid are *richer* in CD68/CD14/CD163/ITGAX/AIF1/C1QA than BAFF⁻ (e.g. CD68 24% vs 16%, C1QA 53% vs 40%) | bona fide, cell-intrinsic — activated-macrophage phenotype, not doublets/ambient |
| **Reproducibility** | detected in **all 16/16 samples**, median **2.94%**, IQR 2.2–3.3%, range 1.3–5.7% | stable anchor, not driven by 1–2 sections |

**Verdict: myeloid BAFF is a RELIABLE, anchored producer** — cell-intrinsic (per-transcriptome 24×
epithelium despite fewer counts), bona fide myeloid identity, reproducible across every sample.

---

## Synthesis

**Receptors wash out under nCount control; the myeloid producer is rock-solid.** The defensible
statement the data supports:

> **Myeloid cells produce BAFF tissue-wide (robustly anchored: cell-intrinsic, ~24× epithelium
> per-transcriptome, reproducible across all 16 samples), and B/plasma cells constitutively express
> the BAFF receptors — but there is NO aggregate-specific concentration of either ligand or receptor.
> The axis reduces to tissue-wide myeloid-derived BAFF with constitutive B-lineage receptors; it does
> not form a spatially localized aggregate/TLS niche.**

This closes the BAFF story consistently across all three tests: the global null was wrong about the
*producer* (myeloid BAFF is real and strong), but right about *localization* — neither the
peri-aggregate stromal ligand (`REPORT_baff_ambient.md`) nor the in-aggregate receptor concentration
(here) survives transcript-density/control-gene rigor. Honest result: a real tissue-wide BAFF
producer, no localized B-aggregate niche.

**Caveats:** receptors and ligand are sparsely detected (BCMA inside n=163 plasma → that cell is the
weakest); pooled across 4 aggregate-bearing sections; aggregate membership inherits the B-lineage
DBSCAN definition and the reannotation typing caveats.
