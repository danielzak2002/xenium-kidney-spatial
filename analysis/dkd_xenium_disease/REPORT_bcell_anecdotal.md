# B-lineage mechanistic signatures across nephropathies — anecdotal

**Extends summary 02; does not re-derive it.** This layers B-lineage *mechanism* onto the B-lineage
*burden* established in 02 (per-sample B-lin%, aggregate counts, hull composition, per-participant
gallery — all pulled, none recomputed).

> **READ THIS FIRST.** Every non-DKD condition here is **n = 1 participant** (IgAN 1003, MN 1005,
> C3GN 1007; AA amyloid n = 2). **No statistics, no p-values, nothing is hypothesis-*tested*.** These
> are **descriptive, single-section, hypothesis-*generating*** observations, explicitly framed as
> *suggestive of a different disease process* and nothing more. DKD is the only condition with a
> usable n (16 sections / B-rich n = 2). Treat each cross-disease contrast as an anecdote to be tested
> in a powered cohort, not a result.

Validated reannotation labels (ARI 0.78 segment / 0.68 immune vs authors). Raw data read-only.

---

## Preflight

**(a) Structural anchors — usable.** Glomerulus = Podo (7,525) + MC (21,956) + EC_glom (25,785) +
PEC (26,331); tubule = PT/iPT/TAL/iTAL/DCT/CNT/PC/IC-A; vessel = EC_Peritub/EC_DVR/VSMC;
interstitial = >30 µm from any of these. Glomerular types are reliably typed (Podo recall ≈ 0.98),
so the peri-glomerular axis is trustworthy. **iPT recall ≈ 0.64 carries** — the injured-tubule arm of
the damage-coupling layer is conservative (injury under-called, distances biased long).

**(b) Gene usability gate (detection in the correct producer cell vs ambient floor; PASS = ≥3% in
producer AND ≥2× floor).** Conditioning matters — the loose preflight over-counted:

| readout | gene | producer % | ambient % | fold | call |
|---|---|---|---|---|---|
| Ig isotype | **IGHG1 (IgG)** | 62.5 (plasma) | 0.38 | 163× | **PASS** |
| Ig isotype | IGHA1 (IgA) | 0.97 (plasma) | 0.01 | 88× | **FAIL — sub-floor** |
| Ig isotype | IGHM / IGHG2 / IGKC / IGLC / JCHAIN | — | — | — | **NOT ON XENIUM PANEL** (0/951k) |
| plasma-state | MZB1 / TNFRSF17 / PRDM1 | 57 / 17 / 16 | <0.8 | 19–276× | PASS |
| B-state | TCL1A / BANK1 / CD27 | 3.3 / 30 / 7.0 | <0.9 | 13–99× | PASS |
| B-state | CD83 | 7.3 | 3.4 | 2.1× | weak-PASS |
| TLS | **CCL19** (T-zone/FRC) | 4.5 (immune) | 0.08 | 52× | **PASS** |
| TLS | CXCL13 / CXCR5 / CCR7 | 0.4 / 0.7 / 0.6 | <0.03 | 12–61× | specific but <1% → **presence-only** |
| TLS | CR2 | 0.9 | 3.4 | 0.3× | **FAIL — ambient-saturating** |

> **Two consequences worth stating plainly.** (1) **The IgA-IgAN axis cannot be tested on this panel.**
> IGHA1 is plasma-specific (88×) but detected in <1% of plasma — and IGHM/IGKC/IGLC are simply absent
> from the Xenium gene set. Only the **IgG** arm of the isotype hypothesis is measurable. (2) **The
> canonical follicular-TLS marker CXCL13 is below the detection floor** (0.42%); only the T-zone
> chemokine **CCL19** is quantitative. CXCL13/CXCR5/CCR7 are reported as presence counts only.

---

## Per-condition table — new mechanistic layers

| condition (sample) | B:Plasma | plasma IgG% (IgA sub-gate) | aggregate? CCL19% / CXCL13+ cells | B localization (periglom / perivasc) | B→injury vs B→glom (µm) | plasma maturity MZB1/BCMA |
|---|---|---|---|---|---|---|
| **DKD B-rich** 1006 | **6.4 (B-skew)** | 53 (1.7) | yes · 6.7% / 2 | 0.28 / 0.38 | **94 vs 34** (glom) | 71 / 25 |
| **DKD B-rich** HK2695 | 2.0 | 72 (0.3) | yes · **10.5% / 303** | 0.43 / 0.24 | **56 vs 26** (glom) | 41 / 9 |
| DKD B-poor (n=6) | 0.8–5.3 | 23–82 | mostly small/none | periglom 0.12–0.31 | ~ balanced–glom | 60–79 / 16–28 |
| Control (n=3) | high (tiny n) | 35–54 | sparse | **periglom 0.04–0.16** (low) | far from both | 15–81 / 18–27 |
| **IgAN** 1003 | 1.4 | 59 (2.8 — best, still sub-floor) | **NO aggregate** · — / 0 | **0.13 / 0.49 (peri-vascular)** | 37 vs 45 (~balanced) | **88 / 37 (most mature)** |
| **MN** 1005 | **0.73 (plasma-skew)** | **69 (0.9 — IgG-dominant)** | yes · **27.4%** / 13 (mixed) | 0.29 / 0.34 | **46 vs 29** (glom) | 63 / 27 |
| AA amyloid 1004/1009 | 4.8 / 6.1 (B-skew) | 24 / 45 | small · ~1 CXCL13 | periglom 0.18–0.21 | ~ balanced | **low: 25–38 / 8–15** |
| C3GN 1007 | 2.8 | 59 | yes · 25% / 40 | periglom 0.15 | 60 vs 53 | 72 / 15 |

(Full per-sample CSVs: `bcell_content.csv`, `bcell_isotype.csv`, `bcell_tls.csv`,
`bcell_localization.csv`, `bcell_damage_coupling.csv`, `bcell_state.csv`. Hull composition pulled from
02's `aggregate_composition.csv`: MN B11/Pl9/My12/CD4T25 vs B-rich DKD B36–38/Pl2/My4–5.)

---

## Anecdotal read — how IgAN vs MN vs DKD *appear* to differ in B-lineage mechanism

**These three sections look like three different B-lineage programs — a hypothesis to test, not a result.**

- **DKD (B-rich)** is the only *organized-lymphoid* picture: B-predominant (B:Plasma 2–6), the **only
  sections carrying follicular CXCL13** (HK2695, 303+ cells) on top of T-zone CCL19, with B sitting
  **near glomeruli/aggregates and far from injured tubules** (B→injury 56–94 µm ≫ B→glom 26–34 µm)
  while myeloid hugs injury (~30 µm) — i.e. a **TLS-like** structure, orthogonal to the injury→myeloid
  axis of summary 06.
- **MN** is the most *antibody-mediated-looking*: the **only plasma-skewed** glomerular disease
  (B:Plasma 0.73), **IgG-dominant / IgA-absent plasma** (69% IGHG1+), B/plasma oriented toward
  glomeruli (B→glom 29 µm), with aggregates present but **mixed/T-zone (CCL19 27%, low CXCL13)** rather
  than follicular — consistent with an **IgG glomerular process** plus disorganized infiltrate.
- **IgAN** is *diffuse and unorganized*: B-lineage present but with **no aggregate at all** (no TLS,
  no CCL19 region), B sitting **peri-vascular not peri-glomerular** (0.49 vs 0.13), yet the **most
  mature plasma** of any section (MZB1 88%, BCMA 37%) and the numerically highest IgA — pointing to a
  **scattered plasma/antibody** process rather than an organized lymphoid one.

> **Hard caveat on the IgAN headline:** the defining IgA axis **is not measurable** on this panel
> (IGHA1 sub-floor; IGHM/IGKC/IGLC absent), so "IgAN looks plasma-diffuse, not TLS-organized" rests on
> architecture + plasma-maturity, **not** on isotype. AA amyloid separately stands out with the
> **least-mature plasma** (BCMA 8–15%) and unusual B activation (CD83 20% in 1004) — also n=1–2, also
> untested.

**Every line above is n=1 per non-DKD condition, descriptive, hypothesis-generating, and NOT
statistically tested.**

---

## Method (one paragraph)

Validated reannotation labels + coordinates (`dkd_xenium_reannotation/cells.parquet`); gated gene
counts re-extracted per-cell from the Demoulin et al. 2026 DKD Xenium h5ad (backed read, Xenium rows
only, aligned on (orig_ident, cell_id)). B-lineage burden, aggregate counts (DBSCAN eps=50 µm,
minPts=20) and hull composition are **pulled from 02** (`per_sample_substrate.csv`,
`aggregate_composition.csv`), not recomputed. New layers: B:Plasma from labels; Ig isotype as detection
fraction **within plasma** (gated); TLS markers among cells in the ≤50 µm aggregate region (same
definition as 02); structural localization by positional nearest-anchor (≤30 µm else interstitial);
damage×B coupling as per-cell median distance to nearest injured-tubule vs nearest glomerular cell;
state as identity-conditioned detection fraction. No permutation, no test — descriptive distributions
only. Figures: `figures/bcell_fig_A_split_isotype.png`, `…_B_tls_state.png`,
`…_C_localization_coupling.png`, `…_D_glom_crops.png`.
