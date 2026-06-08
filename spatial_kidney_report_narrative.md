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
- *cLN (CosMx):* larger panel and higher per-cell breadth, but markedly higher ambient — the empirical basis for the contamination-aware typing in Section 5. On a single unit-matched metric (negative-probe counts ÷ total counts per cell), the mean neg-control fraction is ≈ 0.6% on cLN versus ≈ 9e-5 on RCC-Xenium — about **65× higher**. (The cLN fraction is reconstructed as `negmean × 20` — the CosMx 1000-plex panel's 20 negative probes — since the raw negative-probe matrix is not stored; it is a background-fraction proxy, and the order-of-magnitude gap is robust.)

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

**Contamination-aware typing recovered the high-RNA immune compartment.** On the full 14-slide cohort (single source of truth: `cln_cosmx_immune_benchmark_unified.csv`), InSituType recovered the high-RNA immune types at usable recall *and* precision (recall / precision): plasma 0.952 / 0.698, macrophage 0.712 / 0.689, B 0.594 / 0.472, DC 0.784 / 0.313 — and it recovered B and DC, which the two-stage co-expression pass (R/09) missed entirely (0 recall). The two-stage pass reached comparable *recall* on a few types (monocyte 0.764, NK 0.590, macrophage 0.575) only by gross over-calling — predicting 22k–81k cells for compartments of under 900–9,700 true cells, i.e. precision of 0.030 (monocyte), 0.005 (NK), 0.001 (CD8). Precision is the dimension InSituType won on, and the dimension that matters for a downstream spatial analysis. (All numbers are the full 14-slide run; the earlier 2-slide validation subset is a pilot only and its inflated precisions are not cited.) [TABLE: InSituType per-type recall/precision — `cln_cosmx_immune_benchmark_unified.csv`] [FIGURE: per-type recall AND precision, two-stage vs InSituType]

**Low-RNA T cells hit an intrinsic wall.** Four approaches — generic clustering, marker co-expression, semi-supervised, and pure supervised InSituType with the CD4 profile present — all failed. Two distinct facts, often conflated, both hold: (i) **zero cells were assigned to CD4 specifically** (T_CD4 recall = 0 in the per-type table, even supervised — CD4 vs CD8 is not separable on this panel); and (ii) at the coarse T-lineage level, only a minority of author T cells were retained *as T at all* — author CD8 28%, author CD4/Treg just 11% (and that 11% lands as CD8/other-T, never CD4), with the remainder lost to de-novo and epithelial background. Weighted T-lineage retention is ~16%. A genuine panel-and-ambient limit for low-signal cells, confirmed rather than assumed. [FIGURE: T-cell confusion matrix — author T subtypes vs assigned labels, showing loss to de-novo / epithelial background]

**Validation discipline.** Orthogonal checks throughout: the CD45/PanCK immunofluorescence channels; the published author annotations as an external benchmark; and a gene-usability gate (detection and cell-type enrichment over ambient) run *before* building any analysis on a gene. On CosMx the plasma/B markers MZB1, CD79A, BCMA and the myeloid marker CD68 passed comfortably (3–16× ambient), while the BAFF ligand (TNFSF13B), CCL2 and CCL19 failed (≤2× ambient) — flagging the BAFF ligand as ambient-level prevented a spurious signaling analysis. MS4A1 is a *weak* pass (~2× ambient, ~9% detection), kept but not leaned on. [FIGURE: gene-usability gate — expected-cell-type enrichment vs ambient, with positive controls (MS4A1/MZB1) and failures (BAFF/chemokines) marked]

**Practical resolution.** Analyses were scoped to the compartments that typed reliably on each platform — the full immune repertoire on Xenium; B/plasma/myeloid on CosMx — and explicitly not extended to unreliable compartments.

---

## 6. Spatial findings

All spatial statistics use normalized metrics so results compare across platforms. At these cell counts, significance is near-universal; effect size, direction, and reproducibility across replicates carry the findings. **Each summary metric is shown alongside a representative underlying spatial plot, so the derivation is visible.**

### 6.1 An immunoregulatory B-cell aggregate in renal carcinoma

B cells form discrete aggregates with a consistent composition: regulatory T cells enriched within them, cytotoxic effector-memory CD8 T cells excluded. Neighborhood enrichment: B×Treg strongly positive (z +28), B×effector-CD8 strongly negative (z −117). The excluded population is specifically the **effector-memory CD8 T-cell** cluster (the boxed pair in the heatmap); a second, smaller generic CD8 cluster — labeled `CD8_T`, less clearly differentiated — is grouped with it as "effector CD8" for the composition analysis. Formal DBSCAN delineation (eps = 50 µm) recovered **37 B-cell aggregates** and confirmed the composition per-aggregate: **Tregs enriched in 36 / 37 aggregates** (log2 +1.34, p < 1e-4), **effector-CD8 excluded in 34 / 37** (log2 −0.94, p < 1e-4). Mature-regulatory DCs were members (mregDC log2 +1.40), and the mregDC–CCR7⁺T pairing recurred across ~38 separated foci (de-risked against shared-annotation-origin artifact). Germinal-center markers are off-panel, so this is *not* a classic effector TLS; the composition reads as an immunoregulatory/tolerogenic aggregate.

- [FIGURE: RCC neighborhood-enrichment z-score heatmap — *summary metric*]
- [FIGURE: RCC delineated aggregates colored by cell type — *summary*]
- [FIGURE: representative aggregate crop, cell-type coloring beside MS4A1 / FOXP3 / LAMP3 / CD8A transcript overlays — *underlying derivation, marker-backed*]
- [FIGURE: representative crop of the mregDC–CCR7⁺T foci — *underlying derivation*]

**Independent replication (PRCC).** *Only the two load-bearing effects* replicated in an independent papillary-RCC section — different patient, subtype, B-cell subtype, and panel: **B–Treg enrichment held** (preview log2 +1.68, present in every delineated aggregate, p = 0.031) and **effector-CD8 exclusion held** (preview log2 −1.42, absent from every aggregate, p = 0.031). Directional replication (z scales with n); structure conserved while the dominant B-cell phenotype shifts. The other two compartments did **not** replicate and are *outside* the replicated claim: **mregDC sign-flipped** (RCC +1.40 → preview −1.49) and **plasma sign-flipped** (RCC +0.11 → preview +1.12) — both underpowered in the preview (mregDC n = 94, and the preview plasma trend is non-significant, p = 0.156), so they are reported as inconclusive, not as replicated members of the niche. [TABLE: BIG vs preview replication summary] [FIGURE: side-by-side BIG vs preview aggregate composition]

### 6.2 A plasma–myeloid niche in lupus nephritis — shown per slide

In disease, plasma cells form spatial aggregates with myeloid cells enriched inside them. **What separates disease from control is aggregate formation, not neighborhood enrichment.** Two facts must be kept distinct:

- **Disease-specific (the load-bearing signal): plasma aggregate formation + myeloid enrichment inside aggregates.** Plasma infiltrate and form DBSCAN aggregates only in disease (six plasma-bearing slides across III/IV/IV+V form 1–18 aggregates each); the four control regions form **zero** plasma aggregates (controls are essentially plasma-free). Inside those disease aggregates, myeloid are enriched over the slide background (log2 +0.45 to +2.57). Because controls form no aggregates, the myeloid-in-aggregate metric is intrinsically disease-only, and aggregation is gap-safe (DBSCAN eps = 50 µm cannot bridge the mm gaps between tissue cores).
- **NOT disease-specific (supporting/exploratory only): the plasma×myeloid neighborhood-enrichment z.** It is positive on every plasma-bearing slide (z +4 to +17) — but it is *also* positive on a control region (SMI0016c_SP2219, 18 plasma, z = 5.5, inside the disease range). Sparse plasma sitting near an abundant myeloid compartment yields a positive z even without a niche, so the z alone does **not** demonstrate disease specificity; it is reported as supporting context, not as the discriminating statistic.

CD45 immunofluorescence confirmed both compartments sat in genuine immune regions, not ambient-mislabeled tubule. No clean ISN-class gradient was demonstrable — aggregation appeared across III/IV/IV+V and magnitude tracked per-slide plasma load, not class — and with one IV+V region and multiple regions per patient, class comparisons are descriptive only.

*Methods note (per-section graph).* The cLN sections are multi-core (8–43 cores/slide separated by mm gaps), so a naive per-slide Delaunay graph would span cores. The **primary** plasma×myeloid z values reported here are computed on the **per-section graph** — Delaunay with all cross-core edges (>50 µm) pruned, so adjacency is within-core only. As a sensitivity check, the unpruned graph gives essentially identical values (max |Δz| = 0.65 across all slides; e.g. SP21_213 17.15 vs 16.97, control 5.49 vs 5.48), because the permutation null uses the same graph — the conclusion is robust to cross-core edge handling. The aggregation analysis is intrinsically gap-safe (50 µm DBSCAN cannot bridge mm gaps).

**Per-slide breakout (the most informative view — variability and consistency across specimens).** For *each* plasma-bearing slide (SP21_213 [III], SP19_1139 [IV], SP20_642 [IV+V], SP20_10838 [III], SP19_4061 [IV], SP18_8471 [III]), and a representative control for contrast (SMI0016C_SP17SP19: 1,375 myeloid, ~0 plasma, 0 aggregates), an individual figure set is generated — faceted per tissue core to avoid the inter-core whitespace — rather than a collective plot:

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

A gene-usability gate found the secreted/inducible "recruitment" genes (chemokines, interferon-response) unreadable above ambient on CosMx — only structural macrophage markers usable — so a mechanistic account of recruitment is unsupported. Among the usable markers there is a *weak* signal, and the honest reading is that **the per-slide spread exceeds the mean lean**: the effect sizes are small (log2 ≈ 0.1–0.3) and most genes are inconsistent slide-to-slide. The one marker that is reasonably consistent is **C1QB** (complement; inside > outside in ~5/5 slides, but still small). MRC1 (mannose receptor) has a positive *mean* but is positive in only a minority of slides — it should not be read as an established niche feature. C1QA and HLA-DRA are weakly positive on average. Because complement is central to lupus-nephritis pathology, a complement-producing (C1q⁺) macrophage state associating with the antibody-producing plasma niche is a coherent hypothesis — but it is offered **strictly as correlative and hypothesis-generating**, resting mainly on C1QB; a static snapshot cannot establish recruitment or its direction, and the spread across slides means even the lean is provisional. [FIGURE: myeloid-state log2 enrichment in-vs-out of niche, per gene, with per-slide points] [FIGURE: usability gate result — which candidate genes passed vs failed]

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
