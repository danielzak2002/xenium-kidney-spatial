# Spatial Transcriptomics of B-lineage Immune Organization in Kidney Disease

Imaging-based spatial transcriptomics (10x **Xenium**, NanoString **CosMx**) across **tumor,
autoimmune, and metabolic** kidney disease. The work spans panel assessment, QC tailored to
imaging-based ST, cell-type annotation, spatial-neighborhood and B-cell-aggregate analysis, and
cross-platform / cross-disease comparison — implemented in both R (Seurat v5) and Python
(scanpy / squidpy). Public data only (CC BY 4.0); raw outputs are git-ignored (see [DATA.md](DATA.md)).

## Motivation
Imaging-based spatial transcriptomics resolves single cells in situ, preserving the tissue
architecture that dissociation-based scRNA-seq discards. In kidney, that spatial context is central
to the immune contexture — where infiltrating and resident immune populations localize relative to
nephron structures, and whether they organize into aggregates or tertiary-lymphoid-structure (TLS)
–like niches. The recurring question across these datasets is **B-lineage organization**: where B
and plasma cells aggregate, what they co-organize with, and whether that differs by disease — tested
where possible across biological replicates and across the two imaging platforms.

## Datasets
Public datasets (CC BY 4.0). Raw outputs are **not** committed — see [DATA.md](DATA.md) for download
and placement. Panel sizes and counts below were confirmed from each folder's `gene_panel.json` /
`metrics_summary.csv` / deposited object, not assumed.

| Folder | Disease context | Platform · panel | Scope (confirmed) |
|---|---|---|---|
| `kidney_10x/` | RCC, **clear-cell** ("neoplastic clear cells") — *discovery* | Xenium · 405 genes / 430 probes (377-gene Human Multi-Tissue & Cancer base + custom) **+ 27-plex protein** | 1 section · **465,534** cells |
| `kidney_10x_preview/` | RCC, **papillary (PRCC)** — *independent validation* | Xenium · 377-gene Human Multi-Tissue & Cancer (gene only) | 1 section · **56,510** cells |
| `Danaher24/` | **Childhood lupus nephritis** (autoimmune); Danaher et al. 2024 | CosMx · **957-plex** | **14 sections** (control + SLE) |
| `Demoulin26/` | **Diabetic kidney disease** (metabolic); Dumoulin et al. 2026 | Xenium ~5K-plex **+** CosMx ~1k-plex (one deposited object, **5,443**-gene union; CosMx restricted to its ~1k panel) | **16 Xenium** (951,040 cells) **+ 48 CosMx** samples; processed/annotated AnnData only — **raw molecule tables not deposited** |
| figshare `25685961` *(external, not committed)* | ccRCC — *cross-cohort comparison* | Xenium · custom **380-gene** (280 breast base + 100 kidney / T-cell) | 4 patients, tumor + adjacent; **streamed**, not stored in-repo |

Notes: `kidney_10x_preview/data/` on disk is the **cancer/PRCC section** of the preview release
(56,510 cells), not the paired non-diseased section — and it has the fewer cells, so it served as the
small→big validation target. The `Demoulin26` folder name follows the repo; the paper is Dumoulin et
al. 2026 (Zenodo). `Danaher24` ships the authors' CosMx object + R code; the cell × gene matrix is in
the deposited `.RData`.

## Approach
1. **Panel assessment** — parse each `gene_panel.json`; report panel size and confirm immune / B-cell
   and plasma markers before designing cell typing.
2. **QC for imaging-based ST** — negative-control-probe / blank-codeword rates; segmentation quality;
   flag segmentation merges and ambient/spillover (no scRNA-style count hard-filters).
3. **Clustering + annotation** — Leiden clustering, marker-based annotation (Seurat v5 / scanpy);
   where a deposit provides labels, an independent re-typing is validated against them.
4. **Spatial statistics** (squidpy) — neighborhood enrichment, co-occurrence, and DBSCAN
   B-cell-aggregate delineation with in/around-aggregate composition.
5. **Cross-platform / cross-disease comparison** — harmonized lineage labels, shared-gene-space
   integration for label consistency, with effect sizes computed per cohort on native panels.

## Results
Pure science; figures and full numbers live in each analysis's `REPORT.md` and `figures/`. Headline
findings (each appropriately scoped to what the data support):

- **Tumor (RCC).** In clear-cell RCC, B-cell aggregates are **regulatory-T-enriched and
  cytotoxic-CD8-excluded** (native-panel DBSCAN: Δlog₂(Treg − effector-CD8) ≈ **+2.6**, ~6× Treg
  bias, n = 35 aggregates). The direction **replicates** in an independent papillary-RCC section
  (B–Treg +1.68, CD8-exclusion −1.42; p = 0.031 each) — structure conserved across patient, subtype
  and panel. *(`Demoulin26/analysis/bniche_dbscan/`, `analysis/interaction_map/`.)*
- **Autoimmune (cLN, CosMx).** ~**one-third of epithelial cells carry ambient CD3** (≈28–35% CD3⁺)
  with only ~**2.3× T-vs-epithelial** separability — CD3-based T-lineage calls are swamped by
  segmentation spillover, so **T-cell calls are unreliable on this deposit** (transcript-level
  re-segmentation was not possible — no molecule table). A platform/segmentation QC caution.
  *(`analysis/cln_cd3_contamination/`, `analysis/cln_fastreseg/`.)*
- **Metabolic (DKD, Xenium ~5K).** An independent, label-blind **full-panel reannotation reproduces
  the authors' typing** (segment ARI **0.78** / 85% agreement; immune-subtype ARI **0.68**), and the
  authors' **B-cell-rich DKD subgroup is independently reproduced** (2 of 8 DKD sections, **100%
  concordant** with the authors' B-predominant niche). *(`analysis/dkd_xenium_reannotation/`,
  `analysis/dkd_xenium_disease/`.)*
- **Cross-disease.** B-cell / TLS-like aggregates **recur** across tumor, metabolic and other
  nephropathies, but their **composition is context-specific**: B-dominated and Treg-favoring in RCC;
  B-dominated but Treg/CD8-balanced in DKD (Δlog₂ ≈ +0.24, no Treg bias); plasma–myeloid in cLN;
  mixed in membranous nephropathy and absent in IgA nephropathy. The disease-specific signal is
  aggregate **composition, not the mere presence** of aggregates. *(`analysis/interaction_map/`,
  `Demoulin26/analysis/bniche_replication/`, `analysis/dkd_xenium_disease/`.)*
- **BAFF survival axis (DKD).** Myeloid cells are a **robust, cell-intrinsic, tissue-wide BAFF
  source** (≈8× the epithelial floor; ~24× per-transcriptome; reproducible across all 16 sections)
  and B/plasma cells **constitutively express the receptors** — but under ambient / transcript-density
  controls **neither ligand nor receptor concentrates at B-aggregates**. A tissue-wide producer, not
  a localized niche (an earlier peri-aggregate stromal-BAFF reading was retracted under the controls).
  *(`analysis/dkd_xenium_disease/REPORT_baff_ambient.md`, `…/REPORT_baff_receptor.md`.)*
- **Methods finding.** Cross-panel imaging-platform integration is **gated by whether enough
  reliably-detected shared genes exist to type the cell type of interest**. On the three-cohort
  **68-gene** shared space, B cells could not be typed reliably (over-called), so the integration is
  used for **label consistency only**, with differentials computed **per cohort on native panels**.
  *(`analysis/three_cohort_assessment/`, `analysis/three_cohort_integration/`.)*

## Repository layout
```
README.md                          Public overview (this file)
DATA.md                            How to obtain the (git-ignored) raw data
CLAUDE.md                          Analysis plan, constraints, and conventions
spatial_kidney_report_narrative.md Long-form narrative behind the figure report
environment.yml / setup_R.R        Python env / R dependencies
R/  py/                            Seurat scripts / scanpy·squidpy scripts
outputs/                           figures/ + tables/ + report HTMLs (committed); objects/ (ignored)
notebooks/                         Exploratory only

analysis/                          Cross-dataset analyses (each: scripts, CSV tables, REPORT.md, figures/)
  cln_cd3_contamination/             cLN CosMx — epithelial CD3 ambient-contamination bound (T unreliability)
  cln_fastreseg/                     cLN — transcript-level re-segmentation gate (NO-GO: no molecule table)
  dkd_baff_april/                    DKD — initial BAFF/APRIL/receptor survival-niche assessment (gated)
  dkd_cosmx_aggregates/              DKD — B-aggregate composition, CosMx vs Xenium (cross-platform)
  dkd_epi_endo_stress/               DKD — epithelial/endothelial stress near immune B-aggregates
  dkd_xenium_reannotation/           DKD — independent full-panel reannotation, validated vs author labels
  dkd_xenium_disease/                DKD — B-lineage disease analysis (B-rich subgroup, references, BAFF stress-tests)
  interaction_map/                   RCC/PRCC/cLN/DKD — comparative neighborhood-enrichment map
  three_cohort_assessment/           figshare ccRCC + RCC + DKD — integration feasibility / shared-gene space
  three_cohort_integration/          68-gene three-cohort integration (label consistency; readouts A & B)

presentation/figures/              Slide-grade figure builders (figstyle.py + build_*.py)
  integration_walkthrough/           Integration walkthrough panels

kidney_10x/  kidney_10x_preview/   RCC Xenium — readme.txt + data/ (git-ignored)
Danaher24/                         cLN CosMx — readme.txt + data/ (git-ignored)
Demoulin26/                        DKD atlas — ASSESSMENT.md, assess_dataset.py, data/ (git-ignored),
                                     analysis/ {bniche_dbscan, bniche_replication, cross_platform_tcell}
```
Committed deliverables are PNG figures, CSV tables, and `REPORT.md` files; large objects, raw data,
embedded-base64/vector artifacts, and run logs are git-ignored (see `.gitignore` and CLAUDE.md).

## Reproduce
```bash
# Python (Apple-silicon arm64; miniforge recommended)
conda env create -f environment.yml
conda activate spatial
python -m ipykernel install --user --name spatial

# R (4.4+)
Rscript setup_R.R

# Data: follow DATA.md to download each dataset into its data/ subfolder, then run the
# scripts in R/ , py/ , and the per-analysis scripts under analysis/ .
```

## Tooling
R 4.4+ with Seurat v5; Python 3.11 with scanpy + squidpy + spatialdata (plus harmonypy,
scikit-learn, scipy). Developed on Apple silicon (arm64); one section processed at a time.

## License
- **Code:** MIT — see [LICENSE](LICENSE).
- **Data:** Creative Commons Attribution 4.0 (CC BY 4.0); attribute the original providers
  (10x Genomics; NanoString / Danaher et al.; Dumoulin et al.) when presenting results. See
  [DATA.md](DATA.md).

## Author
Daniel ([@danielzak2002](https://github.com/danielzak2002)).
