#!/usr/bin/env python
"""Per-sample cell-type composition by disease GROUP -- top-level layer for the walkthrough.

REUSES the VALIDATED reannotation labels (summary 01: 951,040-cell / 16-sample object). Does NOT
re-run typing. Descriptive only; n=1 per non-DKD condition; NO statistics. Compositional closure
(fractions sum-to-one) noted, with a CLR sensitivity small-multiple. Raw data read-only.
"""
import os, warnings
warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
RE=f"{REPO}/analysis/dkd_xenium_reannotation"; OUT=f"{REPO}/analysis/dkd_xenium_disease"; FIG=f"{OUT}/figures"
os.makedirs(FIG,exist_ok=True)

# ---- groups (confirmed against the object) ----
GROUPS={"Control":["HK2753","HK3106","HK3626"],
        "DKD":["1001","1006","1008","1010","1011","1012","1013","HK2695"],
        "IgAN":["1003"],"MN":["1005"],"AA amyloid":["1004","1009"],"C3GN":["1007"]}
GORDER=["Control","DKD","IgAN","MN","AA amyloid","C3GN"]
NSING={g for g,s in GROUPS.items() if len(s)==1}  # n=1 groups (no median line)
samp2grp={s:g for g,ss in GROUPS.items() for s in ss}
GC={"Control":"#7f7f7f","DKD":"#6A3D9A","IgAN":"#1f77b4","MN":"#E31A1C","AA amyloid":"#FF7F00","C3GN":"#33A02C"}

c=pd.read_parquet(f"{RE}/cells.parquet").reset_index(drop=True)
c["sample"]=c.orig_ident.astype(str)
assert set(c["sample"].unique())==set(samp2grp), "sample set mismatch vs declared groups"
c["group"]=c["sample"].map(samp2grp)

# ---- resolution definitions (reuse 01 label columns; recompute no typing) ----
LINEAGE=["Epithelial","Immune","Stroma","Endothelial"]               # my_lineage
IMMUNE=["Myeloid","CD4 T","CD8 T","B","Plasma"]                      # my_label (NK/DC not separately typed)
EPI=["PT","iPT","TAL","iTAL","DCT","PC/CNT","IC","Podo","PEC"]        # my_coarse epithelial subtypes

def frac_table(catcol,cats,exclude_unresolved=True):
    """per-sample fraction of TOTAL section cells in each category."""
    d=c if not exclude_unresolved else c[c.my_label!="Unresolved"]
    tot=d.groupby("sample",observed=True).size()
    ct=d.groupby(["sample",catcol],observed=True).size().unstack(fill_value=0)
    fr=ct.div(tot,axis=0)
    for k in cats:
        if k not in fr.columns: fr[k]=0.0
    fr=fr[cats]; fr.index.name="sample"
    return fr

lin=frac_table("my_lineage",LINEAGE); imm=frac_table("my_label",IMMUNE); epi=frac_table("my_coarse",EPI)

# ---- assemble tidy per-sample CSV (all three resolutions) ----
def tidy(fr,res):
    t=fr.reset_index().melt(id_vars="sample",var_name="cell_type",value_name="fraction")
    t["resolution"]=res; t["group"]=t["sample"].map(samp2grp)
    t["n_cells_sample"]=t["sample"].map(c.groupby("sample",observed=True).size())
    return t
TIDY=pd.concat([tidy(lin,"coarse_lineage"),tidy(imm,"immune_subtype"),tidy(epi,"epithelial_subtype")],ignore_index=True)
TIDY=TIDY[["group","sample","resolution","cell_type","fraction","n_cells_sample"]].sort_values(
    ["resolution","cell_type","group","sample"])
TIDY.to_csv(f"{OUT}/composition_by_group.csv",index=False)
print("wrote composition_by_group.csv", TIDY.shape)

# ---- CLR (sensitivity; closure-aware). pseudocount on counts, clr over the resolution's parts ----
def clr_table(catcol,cats):
    d=c[c.my_label!="Unresolved"]
    ct=d.groupby(["sample",catcol],observed=True).size().unstack(fill_value=0)
    for k in cats:
        if k not in ct.columns: ct[k]=0
    ct=ct[cats].astype(float)+0.5                      # multiplicative-ish pseudocount
    logp=np.log(ct); clr=logp.sub(logp.mean(axis=1),axis=0)
    return clr
lin_clr=clr_table("my_lineage",LINEAGE); imm_clr=clr_table("my_label",IMMUNE)

# ---- plotting helpers ----
def jitter(n,w=0.16): return (np.random.RandomState(0).rand(n)-0.5)*2*w if n>1 else np.zeros(n)
def dot_panel(ax,fr,ct,ylab,as_pct=True):
    for gi,grp in enumerate(GORDER):
        ss=GROUPS[grp]; vals=np.array([fr.loc[s,ct] for s in ss if s in fr.index])
        if as_pct: vals=vals*100
        xj=gi+jitter(len(vals))
        ax.scatter(xj,vals,s=46,c=GC[grp],edgecolor="white",linewidth=0.6,zorder=3,alpha=0.95)
        if len(vals)>1:  # group median bar only where n>1
            ax.plot([gi-0.26,gi+0.26],[np.median(vals)]*2,color=GC[grp],lw=2.4,zorder=2)
    ax.set_xticks(range(len(GORDER)))
    ax.set_xticklabels([g+("*" if g in NSING else "") for g in GORDER],fontsize=7.5,rotation=35,ha="right")
    ax.set_title(ct,fontsize=11); ax.set_ylabel(ylab,fontsize=8); ax.margins(x=0.08)
    ax.grid(axis="y",ls=":",alpha=0.4)

NOTE=("each dot = one sample · bar = group median (n>1 only) · *n=1 (single section, descriptive)"
      "  —  ANECDOTAL for non-DKD · NO statistics")

# ===================== FIG 1: coarse lineage (+ CLR sensitivity row) =====================
fig,axes=plt.subplots(2,len(LINEAGE),figsize=(4.0*len(LINEAGE),8.4),sharex=True)
for j,ct in enumerate(LINEAGE):
    dot_panel(axes[0,j],lin,ct,"% of section cells",as_pct=True)
    # CLR row
    ax=axes[1,j]
    for gi,grp in enumerate(GORDER):
        ss=GROUPS[grp]; vals=np.array([lin_clr.loc[s,ct] for s in ss if s in lin_clr.index])
        xj=gi+jitter(len(vals)); ax.scatter(xj,vals,s=40,c=GC[grp],edgecolor="white",linewidth=0.6,zorder=3,alpha=0.95)
        if len(vals)>1: ax.plot([gi-0.26,gi+0.26],[np.median(vals)]*2,color=GC[grp],lw=2.2,zorder=2)
    ax.axhline(0,color="#bbb",lw=0.8,ls="--")
    ax.set_xticks(range(len(GORDER))); ax.set_xticklabels([g+("*" if g in NSING else "") for g in GORDER],fontsize=7.5,rotation=35,ha="right")
    ax.set_title(f"{ct} — CLR",fontsize=9.5,color="#555"); ax.set_ylabel("CLR (sensitivity)",fontsize=8); ax.grid(axis="y",ls=":",alpha=0.4)
fig.suptitle("Coarse-lineage composition per sample by disease group  (top = raw fraction · bottom = CLR sensitivity, closure-aware)",fontsize=13)
fig.text(0.5,0.005,NOTE,ha="center",fontsize=9,color="#a33")
fig.tight_layout(rect=[0,0.02,1,0.97]); fig.savefig(f"{FIG}/composition_by_group_coarse.png",dpi=170,bbox_inches="tight"); plt.close(fig)
print("  [fig] composition_by_group_coarse.png")

# ===================== FIG 2: immune-subtype drill-down (+ CLR row) =====================
fig,axes=plt.subplots(2,len(IMMUNE),figsize=(3.4*len(IMMUNE),8.4),sharex=True)
for j,ct in enumerate(IMMUNE):
    dot_panel(axes[0,j],imm,ct,"% of section cells",as_pct=True)
    ax=axes[1,j]
    for gi,grp in enumerate(GORDER):
        ss=GROUPS[grp]; vals=np.array([imm_clr.loc[s,ct] for s in ss if s in imm_clr.index])
        xj=gi+jitter(len(vals)); ax.scatter(xj,vals,s=40,c=GC[grp],edgecolor="white",linewidth=0.6,zorder=3,alpha=0.95)
        if len(vals)>1: ax.plot([gi-0.26,gi+0.26],[np.median(vals)]*2,color=GC[grp],lw=2.2,zorder=2)
    ax.axhline(0,color="#bbb",lw=0.8,ls="--")
    ax.set_xticks(range(len(GORDER))); ax.set_xticklabels([g+("*" if g in NSING else "") for g in GORDER],fontsize=7.5,rotation=35,ha="right")
    ax.set_title(f"{ct} — CLR",fontsize=9.5,color="#555"); ax.set_ylabel("CLR (sensitivity)",fontsize=8); ax.grid(axis="y",ls=":",alpha=0.4)
fig.suptitle("Immune-subtype composition per sample by disease group  (fraction of all section cells; CLR over immune parts)   ·   NK/DC not separately typed in the validated labels",fontsize=12)
fig.text(0.5,0.005,NOTE,ha="center",fontsize=9,color="#a33")
fig.tight_layout(rect=[0,0.02,1,0.97]); fig.savefig(f"{FIG}/composition_by_group_immune.png",dpi=170,bbox_inches="tight"); plt.close(fig)
print("  [fig] composition_by_group_immune.png")

# ---- console: group-median fraction tables ----
def med_by_group(fr,cats,label):
    print(f"\n=== {label}: group-median % of section cells (n=1 groups* = the single value) ===")
    rows=[]
    for grp in GORDER:
        ss=[s for s in GROUPS[grp] if s in fr.index]
        rows.append([grp+("*" if grp in NSING else "")]+[round(float(np.median([fr.loc[s,ct] for s in ss]))*100,2) for ct in cats])
    print(pd.DataFrame(rows,columns=["group"]+cats).to_string(index=False))
med_by_group(lin,LINEAGE,"coarse lineage")
med_by_group(imm,IMMUNE,"immune subtype")
med_by_group(epi,EPI,"epithelial subtype")
print("\n== composition_by_group done ==")
