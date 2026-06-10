#!/usr/bin/env python
"""
dkd_baff_april.py — is a local B/plasma survival-ligand niche (BAFF/APRIL -> TACI/BAFF-R/BCMA)
detectable around the Demoulin Xenium B-aggregates? GATED on usability (decisive — reported
regardless of outcome).

Genes: TNFSF13B (BAFF), TNFSF13 (APRIL), TNFRSF17 (BCMA), TNFRSF13B (TACI), TNFRSF13C (BAFF-R).
Standard B-cell-survival immunology; pure science.

STAGE 0 gate (Xenium): detection + mean in the EXPECTED source/target cells vs PT-ambient floor.
  BCMA->plasma; TACI/BAFF-R->B; BAFF/APRIL->myeloid/stromal/endothelial. VERDICT per gene.
STAGE 1 (only for genes that pass; if BOTH ligands fail, STOP the source analysis & say so):
  ligand source cell types; ligand+ enrichment near/in B-aggregates; receptor in- vs out-aggregate.

Read-only; backed h5ad; only the 5 gene columns materialized.
"""
import os, warnings
warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
H5=os.path.join(REPO,"Demoulin26/data/spatial_adata_xenium_cosmx_zenodo.h5ad")
OUT=os.path.join(REPO,"analysis/dkd_baff_april"); os.makedirs(OUT,exist_ok=True)
EPS=50.0; MINPTS=20; R=50.0; DKD_COLOR="#6A3D9A"
def hdr(s): print("\n"+"="*78+"\n"+s+"\n"+"="*78)

GENES={"TNFSF13B":"BAFF","TNFSF13":"APRIL","TNFRSF17":"BCMA","TNFRSF13B":"TACI","TNFRSF13C":"BAFF-R"}
LIGANDS=["TNFSF13B","TNFSF13"]; RECEPTORS=["TNFRSF17","TNFRSF13B","TNFRSF13C"]

# ============================================================================
hdr("STAGE 0 — load Xenium, materialize the 5 ligand/receptor columns")
a=ad.read_h5ad(H5, backed="r")
present=[g for g in GENES if g in set(map(str,a.var_names))]
print("present in var:", [(g,GENES[g]) for g in present], "| absent:", [g for g in GENES if g not in present])
xen=(a.obs["tech"].astype(str).values=="Xenium"); xidx=np.where(xen)[0]
sub=a[xidx, present].to_memory()
ann=a.obs["annotation_updated"].astype(str).values[xidx]
imm=a.obs["immune_cell_annotation_combined"].astype(str).values[xidx]
samp=a.obs["orig_ident"].astype(str).values[xidx]
xy=np.asarray(a.obsm["spatial"],float)[xidx]
a.file.close()
C=sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
cnt=pd.DataFrame(C,columns=present)

# expected source/target populations
DCs=[l for l in np.unique(imm) if l.upper().endswith("DC")]
pop={"plasma":(imm=="Plasma"), "B":(imm=="B"),
     "myeloid_stromal_endo":(np.isin(imm,["Macro"]+DCs) | (ann=="Fibroblast") | np.char.startswith(ann.astype(str),"EC_"))}
floor=(ann=="PT")  # PT-ambient floor
EXPECT={"TNFSF13B":"myeloid_stromal_endo","TNFSF13":"myeloid_stromal_endo",
        "TNFRSF17":"plasma","TNFRSF13B":"B","TNFRSF13C":"B"}
print("expected-pop sizes:", {k:int(v.sum()) for k,v in pop.items()}, "| PT floor:", int(floor.sum()))

# ============================================================================
hdr("STAGE 0 (gate) — detection/mean in expected cells vs PT-ambient floor")
def det(mask,g): v=cnt[g].values[mask]; return float((v>0).mean()), float(v.mean())
rows=[]
for g in GENES:
    if g not in present:
        rows.append(dict(gene=g,protein=GENES[g],expected_pop=EXPECT[g],verdict="absent",
                         expect_detect=np.nan,floor_detect=np.nan,ratio=np.nan)); continue
    ep=pop[EXPECT[g]]; ed,em=det(ep,g); fd,fm=det(floor,g)
    ratio=ed/max(fd,1e-9)
    verdict="usable" if (ed>=0.03 and ratio>=2.0) else ("sub-ambient" if ed>0 else "absent")
    rows.append(dict(gene=g,protein=GENES[g],expected_pop=EXPECT[g],n_expect=int(ep.sum()),
        expect_detect=round(ed,4),expect_mean=round(em,4),floor_detect=round(fd,4),
        floor_mean=round(fm,4),ratio=round(ratio,2),verdict=verdict))
gate=pd.DataFrame(rows); gate.to_csv(os.path.join(OUT,"usability_gate.csv"),index=False)
print(gate.to_string(index=False))
passed=set(gate.loc[gate.verdict=="usable","gene"])
ligands_pass=[g for g in LIGANDS if g in passed]; receptors_pass=[g for g in RECEPTORS if g in passed]
print(f"\nUSABLE: ligands={[GENES[g] for g in ligands_pass]}  receptors={[GENES[g] for g in receptors_pass]}")

# ============================================================================
# STAGE 1 — only for genes that pass
# ============================================================================
ran_stage1=False
if passed:
    hdr("STAGE 1 — aggregate distance + in/out enrichment (for genes that pass)")
    # B-aggregate distance per section (reuse machinery)
    dist=np.full(len(xidx),np.inf)
    for s in np.unique(samp):
        sm=samp==s; Bm=sm&(imm=="B")
        if Bm.sum()<MINPTS: continue
        Bc=xy[Bm]; cl=DBSCAN(eps=EPS,min_samples=MINPTS).fit(Bc).labels_; members=Bc[cl!=-1]
        if len(members)==0: continue
        tree=cKDTree(members); d,_=tree.query(xy[sm]); dist[np.where(sm)[0]]=d
    in_agg=dist<=R; out_agg=dist>200
    rows1=[]
    for g in passed:
        ep=pop[EXPECT[g]]
        di=float((cnt[g].values[ep&in_agg]>0).mean()) if (ep&in_agg).sum()>=10 else np.nan
        do=float((cnt[g].values[ep&out_agg]>0).mean()) if (ep&out_agg).sum()>=10 else np.nan
        rows1.append(dict(gene=g,protein=GENES[g],pop=EXPECT[g],
            detect_in_agg=round(di,4) if di==di else np.nan,
            detect_out_agg=round(do,4) if do==do else np.nan,
            log2_in_vs_out=round(float(np.log2((di+1e-3)/(do+1e-3))),3) if (di==di and do==do) else np.nan))
    s1=pd.DataFrame(rows1); print(s1.to_string(index=False)); s1.to_csv(os.path.join(OUT,"stage1_in_vs_out.csv"),index=False)
    ran_stage1=True
    # ligand source (only if a ligand passed)
    if ligands_pass:
        hdr("STAGE 1 — ligand source cell types")
        for g in ligands_pass:
            by=pd.Series(cnt[g].values).groupby(ann).apply(lambda v:(v>0).mean()).sort_values(ascending=False)
            print(f"{GENES[g]} detection by annotation_updated (top): {dict(by.head(6).round(3))}")
    else:
        print("\nBOTH LIGANDS FAIL the gate -> source (BAFF/APRIL-producing) analysis STOPPED by design.")

# ============================================================================
hdr("FIGURE — usability (always) + in/out if ran")
fig,axes=plt.subplots(1,2,figsize=(12,4.6))
ax=axes[0]; x=np.arange(len(GENES)); gl=list(GENES)
ed=[gate.loc[gate.gene==g,"expect_detect"].iloc[0] for g in gl]
fd=[gate.loc[gate.gene==g,"floor_detect"].iloc[0] for g in gl]
ax.bar(x-0.2,ed,0.4,label="expected cells",color=DKD_COLOR)
ax.bar(x+0.2,fd,0.4,label="PT ambient floor",color="#bbbbbb")
ax.axhline(0.03,color="k",ls=":",lw=1,label="3% gate")
ax.set_xticks(x); ax.set_xticklabels([f"{GENES[g]}\n({gate.loc[gate.gene==g,'verdict'].iloc[0]})" for g in gl],fontsize=8,rotation=20)
ax.set_ylabel("detection rate"); ax.set_title("Usability: expected-cell detection vs PT floor"); ax.legend(fontsize=8)
ax2=axes[1]
if ran_stage1 and len(s1):
    xx=np.arange(len(s1))
    ax2.bar(xx-0.2,s1.detect_in_agg.values,0.4,label="in-aggregate",color=DKD_COLOR)
    ax2.bar(xx+0.2,s1.detect_out_agg.values,0.4,label="out (>200um)",color="#bbbbbb")
    ax2.set_xticks(xx); ax2.set_xticklabels([GENES[g] for g in s1.gene],fontsize=9)
    ax2.set_ylabel("detection in expected cells"); ax2.set_title("In- vs out-aggregate (passing genes)"); ax2.legend(fontsize=8)
else:
    ax2.axis("off"); ax2.text(0.5,0.5,"No genes passed the gate\n(no Stage 1)",ha="center",va="center")
fig.suptitle("BAFF/APRIL/BCMA/TACI/BAFF-R survival-niche assessment (Demoulin Xenium)")
fig.tight_layout(); fig.savefig(os.path.join(OUT,"fig_baff_april_usability.png"),dpi=150); plt.close(fig)
print("saved fig_baff_april_usability.png")

# ============================================================================
hdr("WRITING REPORT")
L=[];W=L.append
W("# BAFF/APRIL/TACI/BCMA/BAFF-R survival-niche assessment (Demoulin Xenium) — GATED\n")
W("Read-only, Xenium-only. Standard B-cell-survival immunology. Leads with the usability gate "
  "(valuable either way): detection in expected source/target cells vs the PT-ambient floor.\n")
W("## Usability gate (DECISIVE)\n")
W("| gene | protein | expected cells | expect detect | PT floor detect | ratio | verdict |")
W("|---|---|---|---|---|---|---|")
for _,r in gate.iterrows():
    W(f"| {r.gene} | {r.protein} | {r.expected_pop} | {r.get('expect_detect')} | {r.get('floor_detect')} | {r.get('ratio')} | **{r.verdict}** |")
W(f"\nUsable: ligands {[GENES[g] for g in ligands_pass]}, receptors {[GENES[g] for g in receptors_pass]}.\n")
if not ligands_pass:
    W("## Source analysis: STOPPED (both ligands fail)\n")
    W("**BAFF (TNFSF13B) and APRIL (TNFSF13) are sub-ambient on the Xenium 5k panel** in their "
      "expected myeloid/stromal/endothelial sources (≈PT-floor detection), **consistent with the "
      "cLN CosMx result** where both were sub-ambient. A ligand source/gradient analysis is not "
      "interpretable and was not run.")
if ran_stage1 and len(s1):
    W("\n## Stage 1 (passing genes) — in- vs out-aggregate detection\n")
    W("| gene | protein | pop | in-aggregate | out (>200um) | log2 in/out |")
    W("|---|---|---|---|---|---|")
    for _,r in s1.iterrows():
        W(f"| {r.gene} | {GENES[r.gene]} | {r['pop']} | {r.detect_in_agg} | {r.detect_out_agg} | {r.log2_in_vs_out} |")
W("\n## Verdict & feasible follow-up\n")
if ligands_pass:
    W("At least one ligand is measurable; partial survival-niche evidence — see Stage 1 table.")
else:
    W("**No detectable BAFF/APRIL survival-ligand niche on this panel** (ligands sub-ambient). "
      "Receptor detection (BCMA/TACI/BAFF-R) is reported for completeness but cannot establish a "
      "local survival niche without measurable ligand.")
W("\n**NOTED (not run) feasible alternative inflammatory readout** the paper emphasizes and that "
  "is more likely to be on-panel/above-floor: complement (C3, C4A, C8B), TNF-receptor-family, and "
  "granzyme programs — a candidate follow-up to characterize the aggregate microenvironment without "
  "relying on the low-abundance survival ligands.\n")
W("## Caveats\n")
W("- Gate vs PT-ambient floor (Demoulin dropped negprobes); sub-ambient = not interpretable, not "
  "biological absence. Xenium-only (CosMx subtype imputed). No patient column.")
open(os.path.join(OUT,"REPORT.md"),"w").write("\n".join(L))
print("wrote REPORT.md")
print("\n== TASK C done ==")
