# Cross-Platform Imaging Spatial Transcriptomics of B-Lineage Immune Organization in Human Kidney

*Technical whitepaper / reference document. This is intentionally richer and more detailed than a manuscript — it documents every stage, the metrics that justified each decision, and the limits encountered. Figure/table slots are marked `[FIGURE: …]` / `[TABLE: …]`; many require generation from the data (see the companion generation pass). Where a metric is summarized, an underlying representative plot is requested alongside it so the derivation is visible, not just asserted.*

---

## 1. Overview — what we asked and what we found

This study asks how B-lineage immune cells — B cells and antibody-secreting plasma cells — and the immune populations around them are spatially organized in human kidney, and whether that organization differs between a tumor context and an autoimmune-inflammatory context. It uses three publicly available imaging-based spatial transcriptomics datasets across two platforms:

- **RCC (Xenium)** — a renal cell carcinoma section, ~465,000 cells, 405-gene panel plus a 27-plex cell-surface protein co-assay. *Discovery substrate.*
- **PRCC preview (Xenium)** — an independent papillary renal cell carcinoma section, ~56,000 cells, 377-gene panel. *Independent replication.*
- **cLN (CosMx)** — a published childhood lupus-nephritis atlas, ~532,000 cells across 14 tissue regions (controls plus ISN class III/IV/IV+V disease, 8 patients), 957-gene panel. *Contrasting disease context.*

Three findings emerged, each tested across biological replicates and, where the metric allowed, across platforms:

1. **Tumor (RCC):** B cells form an *immunoregulatory* lymphoid aggregate — co-organized with regulatory T cells and mature-regulatory dendritic cells, while excluding cytotoxic CD8 T cells — replicated in an independent papillary-RCC section.
2. **Lupus nephritis (cLN):** plasma cells form aggregates into which myeloid cells are recruited — present across disease slides, absent in controls.
3. **Shared:** plasma cells self-aggregate in both diseases; what *differs* is myeloid recruitment to those aggregates (lupus only).

**Honest framing up front.** Imaging spatial transcriptomics resolved tissue architecture and cell–cell co-organization robustly and reproducibly here. It could *not* reliably type low-RNA-content immune cells (T/NK) on the larger panel, nor read the low-expression secreted/inducible signaling genes (chemokines, interferon-response) a mechanistic account would need. The findings that survived validation are structural-spatial; mechanism remains out of reach on these data. This document foregrounds those limits — they define where the platforms are trustworthy.

---

## 2. Platforms: powers and limitations, observed directly

The two platforms behaved differently in ways that shaped what could be concluded.

**Panel and breadth.** Xenium sections used 377–405-gene panels; CosMx used 957 genes. The larger CosMx panel is the reason the lupus cohort could be considered for B-survival-axis genes at all — but breadth came with costs.

**Sensitivity and background.** CosMx carried substantial ambient (off-target/diffusion) background; on the Xenium sections this was far less pronounced. Ambient background is the single most consequential platform difference in this study, and it is quantified directly in Section 3.

**Orthogonal modalities.** The RCC Xenium section carried a 27-plex protein co-assay; both platforms carried immunofluorescence channels (PanCK for epithelium, CD45 for immune), used throughout for validation.

**Central empirical lesson — recoverability scales with abundance, at two levels.**

- *Cell-type level:* high-RNA cells (plasma, myeloid, B) typed reliably on both platforms; low-RNA cells (T, NK) typed cleanly on Xenium but collapsed on CosMx, where CD4/CD8 could not be separated at all.
- *Gene-class level:* abundant structural markers were readable; low-expression secreted/inducible genes (chemokines, interferon-response, the BAFF ligand) sat at or below ambient and were unreadable on CosMx.

**Best use, on this evidence.** Xenium: discovery-grade immune typing (including rare states) plus an orthogonal protein assay, at a smaller panel. CosMx: breadth and a large multi-sample disease cohort, at the cost of ambient contamination and loss of low-signal cells/genes.

---

## 3. Quality control across datasets

Imaging spatial transcriptomics QC differs from droplet scRNA-seq: cells are segmented from images, transcripts are counted from a targeted panel, and dedicated negative/control probes report ambient background directly. The QC philosophy here was **flag, don't filter** — only zero-count cells were dropped; cells with high negative-probe signal, blank/control-codeword signal, or segmentation-merge signatures were flagged and tracked, so QC effects on downstream results stay auditable.

**Metrics tracked, and why each matters** (real distributions requested per dataset):

- **Transcripts (counts) per cell** — basic sensitivity; determines whether a cell can be typed at all. Low-count cells are where typing fails first, and are the root of the T-cell problem (small, low-RNA cells). [FIGURE: counts-per-cell distribution, per dataset — violin/histogram]
- **Genes detected per cell** — transcriptional complexity; complements counts. [FIGURE: genes-per-cell distribution, per dataset]
- **Negative-probe signal** (negative-probe fraction; CosMx per-cell `negmean`) — the *direct* ambient readout and the single most important metric here. It corrupts low-expression genes and low-signal cell types, and is why the lupus cohort needed contamination-aware typing. [FIGURE: negative-probe / negmean distribution, per dataset]
- **Cell area / segmentation morphology** — flags segmentation merges/doublets (abnormally large, high-count cells); basis of the seg-merge flag. [FIGURE: cell-area distribution, per dataset]
- **Flag rates** — fraction dropped (zero-count) vs flagged (negative/blank/seg-merge), per dataset. [TABLE: per-dataset QC flag summary]

**Per-dataset interpretation** (summary text to accompany the plots):

- *RCC / PRCC preview (Xenium):* sensitivity and complexity per cell; low ambient background; clean QC profile supporting discovery-grade typing.
- *cLN (CosMx):* larger panel and higher per-cell breadth, but notably higher ambient `negmean` — the empirical basis for the contamination-aware typing in Section 5. The negative-probe distribution is the plot to study here.

**The comparative panel is the most informative single view** — the three datasets side-by-side on counts, genes, and especially ambient — because it shows the ambient difference *empirically* rather than by assertion, and that difference, more than any other factor, determined which analyses were feasible on which platform. [FIGURE: comparative QC panel across the three datasets — counts, genes, ambient]

This section grounds the platform-limitations argument (Section 2) in actual data distributions.

---

## 4. Pipeline and tools — what each was for, and what it bought us

Three tool layers were used in sequence, each addressing a problem the previous could not.

### 4.1 Seurat (R) — the typing backbone

**What it did.** Quality control, normalization, multi-sample integration (Harmony, for the 14-slide cLN cohort), graph-based clustering (Leiden), and the marker-and-reference annotation engine (marker-set scoring; SingleR/Monaco for immune reference; Azimuth for the kidney reference). It was the backbone for both Xenium sections and the first-pass typing of the CosMx cohort.

**What it bought us.** A fast, mature on-ramp with native imaging-ST loaders and a broad annotation ecosystem; unsupervised clustering recovered the major tissue populations well on every dataset (the platform-robust part of the pipeline). [FIGURE: clustering UMAP per dataset, colored by major lineage] [FIGURE: Harmony integration before/after for the cLN cohort — slides mixing without erasing biology]

**Where it reached its limit.** Reference transfer degraded across platforms, and at global resolution it could not resolve the minority immune compartment on CosMx — motivating the next layer.

### 4.2 InSituType — contamination-aware Bayesian typing

**What it did.** A cell typer built for imaging spatial data that explicitly models per-cell background from the negative-probe signal and incorporates the immunofluorescence channels as cohort covariates, assigning cells against reference profiles by posterior probability.

**What it bought us.** It is what recovered the high-RNA immune compartment (B, plasma, macrophage) on the ambient-heavy CosMx data, at usable recall *and* precision — which generic clustering, marker co-expression, and a two-stage sub-clustering could not achieve. The per-cell background model is the mechanism: it discounts ambient signal cell-by-cell rather than thresholding a contaminated count. [FIGURE: InSituType assignment-confidence/posterior distribution] [FIGURE: immune recall before (generic/two-stage) vs after (InSituType)]

**Where it reached its limit.** Even with the CD4 profile explicitly present, it could not assign CD4 T cells on the 957-gene panel — an intrinsic signal limit, not a tuning gap (Section 5).

### 4.3 squidpy — the spatial analysis layer

**What it did.** Built per-section spatial neighbor graphs (Delaunay), computed permutation-based neighborhood-enrichment, distance-resolved co-occurrence, and density-based (DBSCAN) aggregate delineation — operating on cell coordinates.

**What it bought us.** This is where the platforms' actual strength — tissue architecture — was exploited, and where every biological finding was generated. It is scverse-native and scaled to the ~465k–532k-cell sections. [FIGURE: example spatial neighbor graph on a tissue crop]

**Parameterization that mattered.** Spatial graphs were built *per physical section* (never across merged objects); distance scales were set per platform (micrometres for Xenium, millimetres for CosMx); metrics were normalized (z-scores, log2 fold, aggregation rate) so results compare across platforms despite differing units and segmentation.

---

## 5. Annotation — the hardest part, told straight (with the evidence for each pivot)

Cell typing was the most failure-prone stage. Each failure and pivot is shown, not just stated.

**Cross-platform reference transfer degrades.** Benchmarked against the published author annotations on the lupus cohort: unsupervised clustering structure recovered ~82% of assignments and marker-based labeling ~74%, but droplet-snRNA reference transfer reached only ~50%, with many low-confidence calls and epithelial cells mislabeled as immune. [FIGURE: three-way benchmark bar chart — clustering 82% / marker 74% / reference-transfer 50%] [FIGURE: confusion of reference-transfer labels vs author labels, highlighting epithelial→immune errors]

**The minority immune compartment was under-resolved and ambient-corrupted.** On CosMx the immune compartment was ~3% of cells in an epithelium-dominated section; at global resolution it never separated, and count-based marker calling was corrupted by ambient (a rare pDC marker, ambient in tubule, drove a ~34,000-cell over-call). [FIGURE: ambient contamination illustration — a lineage marker (e.g. CD3E) detected across tubular/epithelial cells] [FIGURE: where author-immune cells land across global clusters — the under-resolution]

**Contamination-aware typing recovered the high-RNA immune compartment.** On the full 14-slide cohort (single source of truth: `cln_cosmx_immune_benchmark_unified.csv`), InSituType recovered the high-RNA immune types at usable recall *and* precision — plasma 0.95 recall / 0.70 precision, macrophage 0.71 / 0.69, B 0.59 / 0.47, DC 0.78 / 0.31 — and recovered B and DC, which the two-stage co-expression pass (R/09) missed entirely (0 recall). The two-stage pass reached comparable *recall* on a few types (monocyte 0.76, NK 0.59, macrophage 0.58) only by gross over-calling — predicting 22k–81k cells for compartments of under 900–9,700 true cells, i.e. precision of 0.03 (monocyte), 0.005 (NK), 0.001 (CD8). Precision is the dimension InSituType won on, and it is the dimension that matters for a downstream spatial analysis. (A 2-slide validation *pilot* scored higher — e.g. macrophage 0.99/0.94 — but the full-run numbers above are the reference.) [TABLE: InSituType per-type recall/precision — `cln_cosmx_immune_benchmark_unified.csv`] [FIGURE: per-type recall AND precision, two-stage vs InSituType]

**Low-RNA T cells hit an intrinsic wall.** Four approaches — generic clustering, marker co-expression, semi-supervised, and pure supervised InSituType with the CD4 profile present — all failed; CD4 recall was zero even supervised, T-lineage recall ~16–22%, most T cells lost to epithelial/background. A genuine panel-and-ambient limit for low-signal cells, confirmed rather than assumed. [FIGURE: T-cell confusion matrix — author T subtypes vs assigned labels, showing loss to epithelial/background]

**Validation discipline.** Orthogonal checks throughout: the CD45/PanCK immunofluorescence channels; the published author annotations as an external benchmark; and a gene-usability gate (detection and cell-type enrichment over ambient) run *before* building any analysis on a gene — which, e.g., flagged the BAFF ligand as ambient-level and prevented a spurious signaling analysis. [FIGURE: gene-usability gate — expected-cell-type enrichment vs ambient, with positive controls (MS4A1/MZB1) and failures (BAFF/chemokines) marked]

**Practical resolution.** Analyses were scoped to the compartments that typed reliably on each platform — the full immune repertoire on Xenium; B/plasma/myeloid on CosMx — and explicitly not extended to unreliable compartments.

---

## 6. Spatial findings

All spatial statistics use normalized metrics so results compare across platforms. At these cell counts, significance is near-universal; effect size, direction, and reproducibility across replicates carry the findings. **Each summary metric is shown alongside a representative underlying spatial plot, so the derivation is visible.**

### 6.1 An immunoregulatory B-cell aggregate in renal carcinoma

B cells form discrete aggregates with a consistent composition: regulatory T cells enriched within them, cytotoxic effector CD8 T cells excluded. Neighborhood enrichment: B×Treg strongly positive, B×effector-CD8 strongly negative. Formal aggregate delineation confirmed this per-aggregate — Tregs enriched in 32/35 aggregates, effector-CD8 excluded in 34/35; mature-regulatory DCs were members, and the mregDC–CCR7⁺T pairing recurred across ~38 separated foci (de-risked against shared-annotation-origin artifact). Germinal-center markers are off-panel, so this is *not* a classic effector TLS; the composition reads as an immunoregulatory/tolerogenic aggregate.

- [FIGURE: RCC neighborhood-enrichment z-score heatmap — *summary metric*]
- [FIGURE: RCC delineated aggregates colored by cell type — *summary*]
- [FIGURE: representative aggregate crop, cell-type coloring beside MS4A1 / FOXP3 / LAMP3 / CD8A transcript overlays — *underlying derivation, marker-backed*]
- [FIGURE: representative crop of the mregDC–CCR7⁺T foci — *underlying derivation*]

**Independent replication (PRCC).** The load-bearing effects replicated in an independent papillary-RCC section — different patient, subtype, B-cell subtype, and panel: B–Treg enrichment held (present in every delineated aggregate, p = 0.031), effector-CD8 exclusion held (absent from every aggregate, p = 0.031). Directional replication (z scales with n); structure conserved while the dominant B-cell phenotype shifts. [TABLE: BIG vs preview replication summary] [FIGURE: side-by-side BIG vs preview aggregate composition]

### 6.2 A plasma–myeloid niche in lupus nephritis — shown per slide

Plasma cells form aggregates into which myeloid cells are recruited: plasma×myeloid neighborhood enrichment was positive in every plasma-bearing slide (z +4 to +17), and myeloid were enriched inside plasma aggregates (log2 +0.45 to +2.57). The niche recurred across five independent slides/patients; CD45 immunofluorescence confirmed both compartments sat in genuine immune regions, not ambient-mislabeled tubule. Plasma infiltration distinguished disease from controls (controls essentially plasma-free). No clean ISN-class gradient was demonstrable — the niche appeared across III/IV/IV+V and magnitude tracked per-slide plasma load, not class — and with one IV+V region and multiple regions per patient, class comparisons are descriptive only.

**Per-slide breakout (the most informative view — variability and consistency across specimens).** For *each* plasma-bearing slide (SP21_213 [III], SP19_1139 [IV], SP20_642 [IV+V], SP20_10838 [III], SP18_8471 [III]), and one or more representative controls for contrast, generate an individual figure set rather than a collective plot:

- [FIGURE per slide: plasma + myeloid spatial scatter on the tissue]
- [FIGURE per slide: delineated plasma aggregates with myeloid overlaid]
- [FIGURE per slide: representative high-magnification crop of one plasma–myeloid focus — *underlying derivation of the enrichment z*]
- [FIGURE: across-slide summary — per-slide plasma×myeloid z and myeloid-in-aggregate log2, with the per-slide points shown (not just a mean), so consistency and variability are both visible — *summary metric, points-level*]
- [TABLE: per-slide plasma count, n aggregates, aggregation rate, plasma×myeloid z, myeloid log2-enrich]

### 6.3 The cross-context contrast, tested symmetrically

A naive comparison suggested plasma "disperse" in tumor and "aggregate" in lupus — but that compared two different measurements. Re-analyzing the RCC sections with the *identical* plasma metrics used in lupus resolved it:

- **Plasma self-aggregation is shared** — stronger in RCC (aggregation rate 0.83, 0.44) than lupus (0.10–0.40). The earlier "dispersal" was specifically plasma position relative to *B* aggregates, not plasma aggregation in general.
- **Myeloid recruitment is the disease-specific feature.** RCC plasma×myeloid was strongly negative (z −57, −13) with no myeloid enrichment in plasma aggregates (log2 ≈ 0); lupus was positive on both.

[TABLE: RCC vs cLN plasma–myeloid comparison] [FIGURE: side-by-side representative plasma foci — RCC (myeloid-segregated) vs cLN (myeloid-infiltrated) — *underlying derivation of the contrast*]

### 6.4 Myeloid state in the plasma niche — a bounded hypothesis

A gene-usability gate found the secreted/inducible "recruitment" genes (chemokines, interferon-response) unreadable above ambient on CosMx — only structural macrophage markers usable — so a mechanistic account of recruitment is unsupported. Among usable markers, plasma-niche myeloid leaned weakly toward a complement-producing, mannose-receptor-positive, antigen-presenting state (C1QB up in 5/5 slides; MRC1, C1QA, HLA-DRA also positive; effect sizes log2 0.1–0.3). Because complement is central to lupus-nephritis pathology, a complement-producing macrophage state associating with the antibody-producing plasma niche is a coherent hypothesis — offered strictly as correlative and hypothesis-generating. A static snapshot cannot establish recruitment or its direction. [FIGURE: myeloid-state log2 enrichment in-vs-out of niche, per gene, with per-slide points] [FIGURE: usability gate result — which candidate genes passed vs failed]

---

## 7. Limitations, and where imaging spatial transcriptomics earns its place

**Where it was strong.** These platforms resolved tissue architecture, cell–cell co-organization, and niche composition robustly — and those findings replicated across biological replicates and, for shared metrics, across platforms. The combination of an unbiased neighborhood statistic, formal aggregate delineation, and orthogonal protein/IF validation produced conclusions that held up. This is the regime where imaging spatial transcriptomics delivers: *where* cells are relative to one another, and what co-organizes.

**Where it was weak or out of reach, on these data.**

- *Low-RNA cell types* (T, NK) were not reliably typeable on the larger, ambient-prone panel — confirmed across four typing strategies.
- *Low-expression gene classes* (chemokines, interferon-response, BAFF) sat at ambient and were unreadable on CosMx, foreclosing molecular signaling analyses.
- *Off-panel biology* — germinal-center markers were absent, so some cell-state questions could not be posed.
- *Causation and directionality* cannot be inferred from a static snapshot, regardless of panel quality.
- *Cross-platform quantitative comparison* requires normalized metrics and per-section graphs; absolute distances are not comparable, and multi-patient cohorts demand patient-aware, descriptive treatment of class effects.

**Practical guidance.** For discovery-grade immune typing and orthogonal validation, the smaller Xenium panel with its protein co-assay was the stronger substrate. For breadth and large disease cohorts, CosMx was enabling — provided analyses are scoped to high-abundance cells and genes, paired with contamination-aware typing and orthogonal IF validation. Within those bounds the methods supported reproducible, defensible biology; pushed beyond them they fail quietly, and the discipline that matters most is checking — with orthogonal data and usability gates — before building a conclusion on a label or a gene.
