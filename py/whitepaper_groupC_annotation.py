#!/usr/bin/env python
"""
whitepaper_groupC_annotation.py — Group C (Section 5, annotation challenges) figures.

Inputs (committed/derived artifacts only):
  outputs/tables/cln_cosmx_benchmark_summary.csv             three-way agreement (R/08)
  outputs/tables/cln_cosmx_reftransfer_coarse_confusion.csv  author x reference-transfer (R/16)
  outputs/tables/cln_cosmx_immune_benchmark_unified.csv      recall/precision table (R/14)
  outputs/objects/cln_cosmx.h5ad                             counts + author_celltype + InSituType

Outputs -> outputs/figures/whitepaper/:
  qcC_threeway_benchmark.png   C1: clustering 82 / marker 74 / reference-transfer 50
  qcC_reftransfer_confusion.png C2: reference-transfer vs author (epithelial->immune)
  qcC_ambient_cd3e.png         C3: CD3E detected across tubular/epithelial cells (ambient)
  qcC_insitutype_table.png     C4: per-type recall/precision table (plot = B4)
  qcC_tcell_confusion.png      C5: author T subtypes vs InSituType (loss to epithelial/bg)
  qcC_gene_usability_gate.png  C6: expected-cell-type enrichment vs ambient (MS4A1/MZB1 pass, BAFF/chemokines fail)

  conda run -n spatial python py/whitepaper_groupC_annotation.py
"""
import os, re, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.environ.get("XENIUM_ROOT", os.getcwd())
OBJ = os.path.join(ROOT, "outputs/objects"); TAB = os.path.join(ROOT, "outputs/tables")
FIG = os.path.join(ROOT, "outputs/figures/whitepaper"); os.makedirs(FIG, exist_ok=True)

def coarse(x):
    s = str(x); l = s.lower()
    if len(s) == 1 and s.isalpha(): return "de-novo"
    if re.search(r"tubule|epitheli|podocyte|intercalated|principal|urotheli|henle|nephron|pelvic|connecting|progenitor|proximal", l): return "Epithelial"
    if re.search(r"endotheli|vasa.recta|capillary|glomerul", l): return "Endothelial"
    if re.search(r"fibroblast|myofibroblast|stroma|mural|mesangial", l): return "Stroma"
    if re.search(r"plasma|\bb[- ]?cell|memory b|naive b", l): return "B/Plasma"
    if re.search(r"natural killer|\bnk\b", l): return "NK"
    if re.search(r"cd8|cd4|treg|regulatory t|\bt cell|t cd|ccr7\+ t|effector memory", l): return "T"
    if re.search(r"monocyte|macrophage|myeloid|\bdc\b|mdc|mregdc|dendritic|mast|neutrophil|pdc|plasmacytoid", l): return "Myeloid"
    return "Other"

# ============================================================================
# C1 — three-way benchmark
# ============================================================================
bs = pd.read_csv(os.path.join(TAB, "cln_cosmx_benchmark_summary.csv")).set_index("labeling")["agreement_vs_author"]
# neutral grayscale (NOT the dataset palette — blue/green/red stay RCC/PRCC/cLN doc-wide);
# darker = higher agreement
rows = [("Clustering structure", bs["a_cluster_structure"], "#3a3a3a"),
        ("Marker + lineage gate", bs["b_marker_lineage_gate"], "#8a8a8a"),
        ("Reference transfer\n(snRNA-seq)", bs["c_reference_transfer"], "#cccccc")]
fig, ax = plt.subplots(figsize=(7, 5))
for i, (n, v, c) in enumerate(rows):
    ax.bar(i, v, color=c); ax.text(i, v + 0.01, f"{v*100:.0f}%", ha="center", fontsize=12, weight="bold")
ax.set_xticks(range(3)); ax.set_xticklabels([r[0] for r in rows]); ax.set_ylim(0, 1)
ax.set_ylabel("coarse-lineage agreement vs author labels")
ax.set_title("cLN: three labelings benchmarked vs the authors' annotation\n"
             "unsupervised structure transfers; reference transfer degrades", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcC_threeway_benchmark.png"), dpi=150); plt.close(fig)
print("wrote qcC_threeway_benchmark.png")

# ============================================================================
# C2 — reference-transfer vs author confusion (row-normalized)
# ============================================================================
cm = pd.read_csv(os.path.join(TAB, "cln_cosmx_reftransfer_coarse_confusion.csv"), index_col=0)
cmn = cm.div(cm.sum(1).replace(0, np.nan), axis=0)
fig, ax = plt.subplots(figsize=(8.5, 7))
im = ax.imshow(cmn.values, cmap="Reds", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(cm.shape[1])); ax.set_xticklabels(cm.columns, rotation=40, ha="right")
ax.set_yticks(range(cm.shape[0])); ax.set_yticklabels(cm.index)
ax.set_xlabel("reference-transfer label"); ax.set_ylabel("author label")
for i in range(cm.shape[0]):
    for j in range(cm.shape[1]):
        v = cmn.values[i, j]
        if v == v and v > 0.01:
            ax.text(j, i, f"{v*100:.0f}", ha="center", va="center", fontsize=7,
                    color="white" if v > 0.5 else "#333")
# highlight epithelial->myeloid cell
ei = list(cm.index).index("Epithelial"); mj = list(cm.columns).index("Myeloid")
ax.add_patch(plt.Rectangle((mj-0.5, ei-0.5), 1, 1, fill=False, ec="#1f77b4", lw=2.5))
fig.colorbar(im, ax=ax, fraction=0.046, label="row fraction")
ax.set_title("Reference transfer vs author (row-normalized %): author-epithelial cells\n"
             "mislabeled immune (35.8% -> Myeloid, boxed); no B/Plasma/T/NK calls at all", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcC_reftransfer_confusion.png"), dpi=150); plt.close(fig)
print("wrote qcC_reftransfer_confusion.png")

# ============================================================================
# load cLN h5ad for C3 / C5 / C6
# ============================================================================
a = ad.read_h5ad(os.path.join(OBJ, "cln_cosmx.h5ad"))
au = a.obs["author_celltype"].astype(str).values
ist = a.obs["cell_type"].astype(str).values     # InSituType label
au_c = np.array([coarse(v) for v in au]); ist_c = np.array([coarse(v) for v in ist])
def counts(gene):
    j = a.var_names.get_loc(gene)
    return np.asarray(a.X[:, j].todense()).ravel() if hasattr(a.X, "todense") else np.asarray(a.X[:, j]).ravel()

# ---- C3 ambient CD3E: detected across non-T (tubular/epithelial) cells ------
cd3e = counts("CD3E")
order = ["T", "Myeloid", "B/Plasma", "NK", "Endothelial", "Stroma", "Epithelial", "de-novo"]
detect = [float((cd3e[au_c == g] > 0).mean()) if (au_c == g).any() else np.nan for g in order]
meanex = [float(cd3e[au_c == g].mean()) if (au_c == g).any() else np.nan for g in order]
fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
cols = ["#d62728" if g == "T" else "#999" for g in order]
ax[0].bar(range(len(order)), detect, color=cols)
ax[0].set_xticks(range(len(order))); ax[0].set_xticklabels(order, rotation=30, ha="right")
ax[0].set_ylabel("fraction CD3E+ cells"); ax[0].set_title("CD3E detection by author lineage", fontsize=10)
ax[1].bar(range(len(order)), meanex, color=cols)
ax[1].set_xticks(range(len(order))); ax[1].set_xticklabels(order, rotation=30, ha="right")
ax[1].set_ylabel("mean CD3E counts/cell"); ax[1].set_title("CD3E mean expression by author lineage", fontsize=10)
epi_det = detect[order.index("Epithelial")]
fig.suptitle(f"Ambient contamination — the T-cell marker CD3E is detected in {epi_det*100:.0f}% of "
             f"epithelial/tubular cells (ambient), not only in T cells (red)", fontsize=11)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcC_ambient_cd3e.png"), dpi=150); plt.close(fig)
print("wrote qcC_ambient_cd3e.png")

# ---- C5 T-cell confusion: author T subtypes vs InSituType coarse -------------
def tsub(x):
    l = str(x).lower()
    if "cd8" in l: return "author CD8 T"
    if "cd4" in l or "treg" in l or "regulatory" in l or "helper" in l: return "author CD4/Treg"
    if re.search(r"\bt cell|t-cell|mait|gamma", l): return "author T (other)"
    return None
tlab = np.array([tsub(v) for v in au], dtype=object)
tm = tlab != None
dest_order = ["T", "de-novo", "Epithelial", "Myeloid", "NK", "B/Plasma", "Endothelial", "Stroma", "Other"]
rows_t = [r for r in ["author CD8 T", "author CD4/Treg", "author T (other)"] if (tlab == r).sum() > 0]
M = np.zeros((len(rows_t), len(dest_order)))
for ri, rt in enumerate(rows_t):
    sel = tlab == rt
    if sel.sum() == 0: continue
    dist = pd.Series(ist_c[sel]).value_counts(normalize=True)
    for cj, d in enumerate(dest_order):
        M[ri, cj] = dist.get(d, 0.0)
fig, ax = plt.subplots(figsize=(10, 3.6))
im = ax.imshow(M, cmap="Purples", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(len(dest_order))); ax.set_xticklabels(dest_order, rotation=30, ha="right")
ax.set_yticks(range(len(rows_t))); ax.set_yticklabels([f"{r}\n(n={int((tlab==r).sum())})" for r in rows_t])
for i in range(len(rows_t)):
    for j in range(len(dest_order)):
        if M[i, j] > 0.01:
            ax.text(j, i, f"{M[i,j]*100:.0f}", ha="center", va="center", fontsize=8,
                    color="white" if M[i, j] > 0.5 else "#333")
ax.set_xlabel("InSituType assignment (coarse)")
fig.colorbar(im, ax=ax, fraction=0.025, label="row fraction")
ax.set_title("T-cell loss on CosMx (author T subtypes → InSituType, coarse)\n"
             "CD8 28% kept as T; CD4/Treg only 11% (as CD8/other-T, never CD4); "
             "rest lost to de-novo / epithelial background", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcC_tcell_confusion.png"), dpi=150); plt.close(fig)
print("wrote qcC_tcell_confusion.png")

# ---- C6 gene-usability gate: expected-cell-type enrichment vs ambient --------
# gene -> expected high lineage (coarse, by author labels). Pass = detected AND >=2x ambient.
GENES = [("MS4A1", "B/Plasma", "B-cell control"), ("CD79A", "B/Plasma", "B-cell control"),
         ("MZB1", "B/Plasma", "plasma control"), ("TNFRSF17", "B/Plasma", "plasma (BCMA) control"),
         ("CD68", "Myeloid", "myeloid control"),
         ("TNFSF13B", "Myeloid", "BAFF ligand"), ("CCL2", "Myeloid", "chemokine"),
         ("CCL19", "Stroma", "chemokine")]
amb_mask = au_c == "Epithelial"
rowsg = []
for g, exp, note in GENES:
    if g not in a.var_names:
        rowsg.append(dict(gene=g, note=note, det=np.nan, log2=np.nan, usable=False, absent=True)); continue
    c = counts(g); tgt = au_c == exp
    det = float((c[tgt] > 0).mean()); mt = float(c[tgt].mean()); ma = float(c[amb_mask].mean())
    log2 = float(np.log2((mt + 1e-6) / (ma + 1e-6)))
    rowsg.append(dict(gene=g, note=note, det=round(det, 3), log2=round(log2, 2),
                      usable=(det >= 0.05 and log2 >= 1.0), absent=False))
gdf = pd.DataFrame(rowsg); gdf.to_csv(os.path.join(TAB, "qcC_gene_usability_gate.csv"), index=False)
fig, ax = plt.subplots(figsize=(10, 5))
xs = np.arange(len(gdf))
cols = ["#2ca02c" if u else "#c44e52" for u in gdf.usable]
ax.bar(xs, gdf.log2.fillna(0), color=cols)
ax.axhline(1.0, color="k", ls="--", lw=1, label="usability threshold (2× ambient)")
ax.set_xticks(xs); ax.set_xticklabels([f"{r.gene}\n{r.note}" for r in gdf.itertuples()], fontsize=8)
ax.set_ylabel("log2(mean in expected lineage / ambient epithelial)")
for xi, r in zip(xs, gdf.itertuples()):
    tag = "PASS" if r.usable else "FAIL"
    ax.text(xi, (r.log2 if r.log2 == r.log2 else 0) + 0.1, f"{tag}\ndet {r.det:.2f}" if r.det==r.det else "absent",
            ha="center", fontsize=7, color=("#2ca02c" if r.usable else "#c44e52"))
ax.set_title("Gene-usability gate (cLN CosMx): B/plasma + CD68 markers pass; "
             "BAFF ligand & chemokines fail (≤2× ambient)", fontsize=11)
ax.legend(loc="upper right", fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcC_gene_usability_gate.png"), dpi=150); plt.close(fig)
print("wrote qcC_gene_usability_gate.png")

# ---- C4 InSituType recall/precision TABLE (plot is B4) -----------------------
u = pd.read_csv(os.path.join(TAB, "cln_cosmx_immune_benchmark_unified.csv"))
show = u[["immune_type", "n_author", "recall_twostage", "precision_twostage",
          "recall_insitutype", "precision_insitutype"]].copy()
show.columns = ["type", "n", "rec(2-stage)", "prec(2-stage)", "rec(InSituType)", "prec(InSituType)"]
fig, ax = plt.subplots(figsize=(9, 4.2)); ax.axis("off")
tb = ax.table(cellText=show.round(3).fillna("—").values, colLabels=show.columns,
              cellLoc="center", loc="center")
tb.auto_set_font_size(False); tb.set_fontsize(8.5); tb.scale(1, 1.35)
for j in range(len(show.columns)):
    tb[0, j].set_facecolor("#4c72b0"); tb[0, j].set_text_props(color="white", weight="bold")
ax.set_title("InSituType per-type recall & precision vs two-stage (single source of truth: R/14)\n"
             "T_CD4 recall = 0: no cell assigned to CD4 specifically (CD4/CD8 not separable on this panel)",
             fontsize=9.5, pad=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcC_insitutype_table.png"), dpi=150); plt.close(fig)
print("wrote qcC_insitutype_table.png")
print("== Group C done ==")
