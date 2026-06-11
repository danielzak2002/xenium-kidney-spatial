#!/usr/bin/env python
"""make_index.py — INDEX.md (talk order + captions) + contact-sheet montage of all figures."""
import os, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.image as mpimg
import figstyle as fs
OUT=fs.OUT

# talk order (SPEC §7) with one-line scientific captions
TALK=[
 ("A1","Spatial cell-type maps (ccRCC, cLN, DKD): orient the audience in each tissue."),
 ("A3","CD4/CD8 subtype is measured on Xenium but imputed on CosMx (CD4 AUROC≈0.5)."),
 ("B2","cLN: ambient CD3 mis-assigned to epithelium (~35% epithelial CD3+, ~2.3× T-vs-epi)."),
 ("B3","Usability gate: which markers clear vs fail the ambient floor (BAFF/APRIL, fibroEMT fail)."),
 ("B1","★ Aggregate marker overlays: Treg-around/cytotoxic-excluded (ccRCC) vs cytotoxic-mixed (DKD)."),
 ("C1","DBSCAN B-cell density aggregates delineated over tissue (ccRCC, DKD)."),
 ("C2","Per-aggregate Treg vs effector-CD8 enrichment: ccRCC vs DKD."),
 ("C3","★ Burden-corrected differential (HEADLINE): Treg-over-cytotoxic bias is tumor-specific."),
 ("C4","Count-pooled radial profile: Treg flat across rings, mild cytotoxic-core gradient (no collar)."),
 ("C5","Comparative non-immune×immune neighborhood enrichment; absolute z=geometry, biology=differential."),
 ("D1","Conserved B/plasma scaffold with context-specific immune wiring (schematic)."),
 ("D2","Platform capability matrix: what each panel can and cannot establish."),
 ("A2","cLN CosMx immune typing benchmark (InSituType recall/precision vs author labels)."),
 ("B4","DKD injured-PT program elevated near B-aggregates (Δz +0.13, 6/9 sections, p=0.038)."),
]

# INDEX.md
lines=["# Presentation figures — talk order & captions\n",
 "Pure-science figures for the spatial-kidney work. Hero (★) = B1, C3. PNG@300 + SVG per ID.\n",
 "**Take-home:** rigor/platform first, then biology as proof — a conserved B/plasma scaffold acquires "
 "context-specific wiring: immunoregulatory (Treg-in/cytotoxic-out) in tumor, ABSENT in non-malignant disease.\n",
 "## Talk order"]
groups=[("1",["A1"]),("2",["A3","B2","B3"]),("3",["B1","C1","C2"]),("4",["C3","C4"]),
        ("5",["C5"]),("6",["D1"]),("7 (+caveats)",["D2"])]
cap={k:v for k,v in TALK}
for g,ids in groups:
    lines.append(f"\n**Block {g}**")
    for i in ids: lines.append(f"- **{i}** — {cap.get(i,'')}")
lines.append("\n## Supporting (referenced as needed)")
for i in ["A2","B4"]: lines.append(f"- **{i}** — {cap.get(i,'')}")
lines.append("\n## Caveats baked into the figures")
lines+=["- cLN T-lineage UNRELIABLE (35% CD3 contamination) — excluded from conserved-T claims.",
 "- ccRCC stroma–immune inversion in C5 is PROVISIONAL (not tile-verified).",
 "- B-aggregate niches are associational/colocalization, not communication or causation.",
 "- No patient column for DKD (donor clustering uncontrolled); RCC epithelium is malignant (separate).",
 "- D1 is a programmatic schematic — may want manual vector polish."]
open(os.path.join(OUT,"INDEX.md"),"w").write("\n".join(lines))
print("wrote INDEX.md")

# contact sheet montage in talk order
order=[i for _,ids in groups for i in ids]+["A2","B4"]
n=len(order); ncol=4; nrow=(n+ncol-1)//ncol
fig,axes=plt.subplots(nrow,ncol,figsize=(4*ncol,3*nrow))
for ax,fid in zip(axes.ravel(),order):
    p=os.path.join(OUT,f"{fid}.png")
    if os.path.exists(p):
        cp=cap.get(fid,'').replace("★","*")[:46]
        ax.imshow(mpimg.imread(p)); ax.set_title(f"{fid}  {cp}",fontsize=8)
    ax.axis("off")
for ax in axes.ravel()[n:]: ax.axis("off")
fig.suptitle("Spatial-kidney presentation — contact sheet (talk order)",fontsize=16,fontweight="bold")
fig.tight_layout(); fig.savefig(os.path.join(OUT,"CONTACT_SHEET.png"),dpi=150,bbox_inches="tight"); plt.close(fig)
print("wrote CONTACT_SHEET.png")
print("== index done ==")
