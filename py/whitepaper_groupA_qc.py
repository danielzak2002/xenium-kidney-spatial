#!/usr/bin/env python
"""
whitepaper_groupA_qc.py — Group A (Section 3, QC) whitepaper figures.

Inputs (committed/derived artifacts only — nothing fabricated):
  outputs/objects/qc_metrics_<label>.csv   per-cell QC metrics (from R/12, 01_qc objects)
  outputs/tables/<label>_qc_summary.csv    flag counts (from R/01_load_qc)

Outputs -> outputs/figures/whitepaper/:
  qcA_<label>_dist.png        per-dataset: counts / genes / ambient / cell-area distributions
  qcA_comparative.png         three datasets side-by-side: counts, genes, ambient
  qcA_flag_summary.png        per-dataset dropped(zero) vs flagged bars
  outputs/tables/qcA_flag_summary.csv   per-dataset QC flag summary table

  conda run -n spatial python py/whitepaper_groupA_qc.py
"""
import os, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.environ.get("XENIUM_ROOT", os.getcwd())
OBJ  = os.path.join(ROOT, "outputs/objects"); TAB = os.path.join(ROOT, "outputs/tables")
FIG  = os.path.join(ROOT, "outputs/figures/whitepaper"); os.makedirs(FIG, exist_ok=True)

# label -> (display name, platform, color)
DS = [("kidney_RCC_protein", "RCC (Xenium)",         "xenium", "#4c72b0"),
      ("kidney_preview_PRCC", "PRCC preview (Xenium)", "xenium", "#55a868"),
      ("cln_cosmx",           "cLN (CosMx)",          "cosmx",  "#c44e52")]
NEG_FLAG_XEN = 0.02      # config: neg_frac_flag (Xenium per-cell neg-control fraction)
NEG_FLAG_COS = 1.0       # config: cosmx_neg_ratio_flag (negmean*nGenes/nCount)

data = {lab: pd.read_csv(os.path.join(OBJ, f"qc_metrics_{lab}.csv")) for lab, *_ in DS}

def clipped(x, q=0.995):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    hi = np.quantile(x, q); return x[x <= hi]

# ---------------------------------------------------------------------------
# per-dataset distribution figures
# ---------------------------------------------------------------------------
for lab, name, plat, col in DS:
    d = data[lab]
    fig, ax = plt.subplots(1, 4, figsize=(18, 4.2))
    ax[0].hist(clipped(d.n_counts), bins=60, color=col, alpha=0.85)
    ax[0].set_xlabel("transcripts / cell"); ax[0].set_ylabel("cells")
    ax[0].set_title(f"counts (median {d.n_counts.median():.0f})", fontsize=10)
    ax[1].hist(clipped(d.n_genes), bins=60, color=col, alpha=0.85)
    ax[1].set_xlabel("genes detected / cell")
    ax[1].set_title(f"genes (median {d.n_genes.median():.0f})", fontsize=10)
    # ambient: Xenium neg-probe fraction; CosMx per-cell negmean (raw background)
    if plat == "cosmx":
        amb = clipped(d.negmean); ax[2].hist(amb, bins=60, color=col, alpha=0.85)
        ax[2].set_xlabel("per-cell negmean (background)")
        ax[2].set_title(f"ambient — negmean (median {np.median(d.negmean.dropna()):.2f})", fontsize=10)
    else:
        amb = d.neg_frac.values
        ax[2].hist(amb[amb <= 0.05], bins=60, color=col, alpha=0.85)
        ax[2].axvline(NEG_FLAG_XEN, color="k", ls="--", lw=1, label="flag thr 0.02")
        ax[2].set_xlabel("neg-control fraction / cell"); ax[2].legend(fontsize=8)
        ax[2].set_title(f"ambient — neg fraction (>0 in {100*(amb>0).mean():.1f}% cells)", fontsize=10)
    ax[3].hist(clipped(d.cell_area), bins=60, color=col, alpha=0.85)
    ax[3].set_xlabel("cell area" + (" (px/µm²)" if plat == "xenium" else " (CosMx units)"))
    ax[3].set_title(f"cell area (median {d.cell_area.median():.0f})", fontsize=10)
    fig.suptitle(f"{name} — per-cell QC distributions  (n shown = {len(d):,})", fontsize=13)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, f"qcA_{lab}_dist.png"), dpi=150); plt.close(fig)
    print(f"wrote qcA_{lab}_dist.png")

# ---------------------------------------------------------------------------
# comparative panel: counts, genes, ambient (the key view)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(16, 5))
names = [name for _, name, *_ in DS]; cols = [c for *_, c in DS]
# counts & genes: violins (log y) of real distributions
for k, (col, ttl) in enumerate([("n_counts", "transcripts / cell"), ("n_genes", "genes / cell")]):
    vals = [np.log10(data[lab][col].clip(lower=1)) for lab, *_ in DS]
    parts = ax[k].violinplot(vals, showmedians=True, showextrema=False)
    for b, c in zip(parts["bodies"], cols):
        b.set_facecolor(c); b.set_alpha(0.7)
    ax[k].set_xticks([1, 2, 3]); ax[k].set_xticklabels(names, rotation=15, ha="right")
    ax[k].set_ylabel(f"log10({ttl})"); ax[k].set_title(ttl, fontsize=11)
# ambient: mean neg-fraction (log scale) — the empirical ambient difference
qc = {lab: pd.read_csv(os.path.join(TAB, f"{lab}_qc_summary.csv")) for lab, *_ in DS}
amb_mean = [float(qc[lab]["mean_neg_frac"].iloc[0]) for lab, *_ in DS]
bars = ax[2].bar(range(3), amb_mean, color=cols, alpha=0.85)
ax[2].set_yscale("log"); ax[2].set_xticks(range(3)); ax[2].set_xticklabels(names, rotation=15, ha="right")
ax[2].set_ylabel("mean neg-control fraction / cell (log)")
ax[2].set_title("ambient background — the decisive platform difference", fontsize=11)
for b, v in zip(bars, amb_mean):
    ax[2].text(b.get_x() + b.get_width()/2, v, f"{v:.1e}", ha="center", va="bottom", fontsize=9)
fold = amb_mean[2] / amb_mean[0]
fig.suptitle(f"Comparative QC across datasets — CosMx ambient ≈ {fold:.0f}× the RCC-Xenium level", fontsize=13)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcA_comparative.png"), dpi=150); plt.close(fig)
print("wrote qcA_comparative.png")

# ---------------------------------------------------------------------------
# flag-summary table + figure (from qc_summary.csv — full-data counts)
# ---------------------------------------------------------------------------
rows = []
for lab, name, plat, col in DS:
    q = qc[lab].iloc[0]
    loaded = int(q["n_cells_loaded"]); kept = int(q["n_cells_kept"])
    rows.append(dict(dataset=name, platform=plat, n_loaded=loaded, n_kept=kept,
                     dropped_zero_count=loaded - kept,
                     flag_negative=int(q["n_flag_neg"]),
                     flag_blank=(int(q["n_flag_blank"]) if pd.notna(q["n_flag_blank"]) else 0),
                     flag_seg_merge=int(q["n_flag_seg_merge"]),
                     flag_lowq=int(q["n_flag_lowq"])))
flag = pd.DataFrame(rows)
flag.to_csv(os.path.join(TAB, "qcA_flag_summary.csv"), index=False)
print("\nQC flag summary:\n", flag.to_string(index=False))

cats = ["dropped_zero_count", "flag_negative", "flag_blank", "flag_seg_merge", "flag_lowq"]
catcol = {"dropped_zero_count": "#444", "flag_negative": "#c44e52", "flag_blank": "#dd8452",
          "flag_seg_merge": "#8172b3", "flag_lowq": "#937860"}
fig, ax = plt.subplots(figsize=(9, 5))
x = np.arange(len(DS)); w = 0.16
for i, cat in enumerate(cats):
    ax.bar(x + (i - 2) * w, flag[cat] + 0.5, w, color=catcol[cat], label=cat.replace("_", " "))
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels([r["dataset"] for r in rows], rotation=15, ha="right")
ax.set_ylabel("cells (log; flagged, not dropped)")
ax.set_title("QC flags per dataset — flag-don't-filter (only zero-count cells dropped)", fontsize=11)
ax.legend(fontsize=8, ncol=2)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "qcA_flag_summary.png"), dpi=150); plt.close(fig)
print("wrote qcA_flag_summary.png + qcA_flag_summary.csv")
print("== Group A done ==")
