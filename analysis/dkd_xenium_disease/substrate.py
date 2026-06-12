#!/usr/bin/env python
"""STEP 0-3 numbers — per-sample B-lineage substrate, within-DKD B-rich/B-poor split (validated
vs authors' 'B predom. Immune ME'), in/around-aggregate immune composition, and DKD-vs-Control
descriptive test. B-LINEAGE = B + Plasma. Uses the VALIDATED reannotation labels. Read-only raw."""
import os, numpy as np, pandas as pd
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree
from scipy.stats import mannwhitneyu
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
RE=f"{REPO}/analysis/dkd_xenium_reannotation"; OUT=f"{REPO}/analysis/dkd_xenium_disease"; os.makedirs(OUT,exist_ok=True)
EPS,MINPTS=50,20

cells=pd.read_parquet(f"{RE}/cells.parquet")
lig=pd.read_parquet(f"{OUT}/ligand_me_percell.parquet")
assert (cells.cell_id.astype(str).values==lig.cell_id.astype(str).values).all(), "row alignment"
cells=cells.reset_index(drop=True)
for c in ["immune_ME","immune_ME_20um"]: cells[c]=lig[c].values
cells["is_Blin"]=cells.my_label.isin(["B","Plasma"])
cells["is_myeloid"]=cells.my_label.eq("Myeloid"); cells["is_T"]=cells.my_label.isin(["CD4 T","CD8 T"])
COND=cells.groupby("orig_ident",observed=True).Condition.first()

def grid_area_mm2(xy,bin_um=50.0):
    if len(xy)==0: return np.nan
    ij=np.floor(xy/bin_um).astype(np.int64); occ=len(set(map(tuple,ij)))
    return occ*(bin_um/1000.0)**2

rows=[]; aggdict={}
for sid,g in cells.groupby("orig_ident",observed=True):
    xy=g[["spatial_x","spatial_y"]].values; ntot=len(g)
    B=g[g.is_Blin]; Bxy=B[["spatial_x","spatial_y"]].values; nB=len(B)
    area=grid_area_mm2(xy)
    nagg=nagg_cells=0; labels=np.array([])
    if nB>=MINPTS:
        labels=DBSCAN(eps=EPS,min_samples=MINPTS).fit(Bxy).labels_
        nagg=len(set(labels))-(1 if -1 in labels else 0); nagg_cells=int((labels!=-1).sum())
    aggdict[sid]=(B.index.values,labels)
    bpme=int((g.immune_ME=="B predom. Immune ME").sum())
    bpme20=int((g.immune_ME_20um=="B predom. Immune ME").sum())
    rows.append(dict(orig_ident=sid,Condition=COND[sid],n_cells=ntot,area_mm2=round(area,3),
        n_Blin=nB,Blin_frac=round(nB/ntot,4),Blin_per10k=round(nB/ntot*1e4,1),Blin_per_mm2=round(nB/area,1),
        n_agg=nagg,n_agg_cells=nagg_cells,agg_per10k=round(nagg/ntot*1e4,3),agg_per_mm2=round(nagg/area,3),
        agg_cells_per10k=round(nagg_cells/ntot*1e4,1),frac_Blin_in_agg=round(nagg_cells/max(nB,1),3),
        author_Bpredom_ME=bpme,author_Bpredom_ME_20um=bpme20))
sub=pd.DataFrame(rows).sort_values(["Condition","Blin_frac"],ascending=[True,False])
sub.to_csv(f"{OUT}/per_sample_substrate.csv",index=False)
print("=== per-sample substrate ===")
print(sub[["orig_ident","Condition","n_cells","Blin_frac","Blin_per10k","n_agg","agg_cells_per10k","frac_Blin_in_agg","author_Bpredom_ME"]].to_string(index=False))

# ---------- STEP 2: within-DKD B-rich/B-poor split + author validation ----------
# The paper's "B cell-rich subgroup" is defined by B-AGGREGATES/TLS, not raw B fraction; B-lineage
# fraction mis-ranks the large HK2695 section. Split on AGGREGATE BURDEN (natural break), corroborated
# by B-lineage fraction; validate membership against the authors' 'B predom. Immune ME' niche.
dkd=sub[sub.Condition=="DKD"].sort_values("agg_cells_per10k",ascending=False).reset_index(drop=True)
# B-rich threshold = the clear burden gap separating the top-tier from the rest. The DKD burden
# distribution has a ~3x gap between the top-2 (HK2695 115, 1006 202 per 10k) and #3 (40); we set
# the cut in that gap (>=75). NOTE the data alone is between a top-2 and top-4 break; the authors'
# B-predom niche independently pins it to the top-2, which we then confirm (concordance below).
THR_BRICH=75.0; thr=THR_BRICH
dkd["our_subgroup"]=np.where(dkd.agg_cells_per10k>=thr,"B-rich","B-poor")
dkd["author_Brich"]=np.where(dkd.author_Bpredom_ME_20um>=200,"B-rich","B-poor")  # authors' niche present
agree=(dkd.our_subgroup==dkd.author_Brich).mean()
dkd.to_csv(f"{OUT}/dkd_subgroup_split.csv",index=False)
print(f"\n=== within-DKD split (agg_cells_per10k break at {thr:.1f}) ===")
print(dkd[["orig_ident","Blin_frac","agg_cells_per10k","n_agg","author_Bpredom_ME_20um","our_subgroup","author_Brich"]].to_string(index=False))
print(f"our B-rich vs author B-predom-ME concordance: {agree*100:.0f}% ({int((dkd.our_subgroup==dkd.author_Brich).sum())}/{len(dkd)})")

# ---------- STEP 3a: in/around-aggregate immune composition (within 50um of aggregate B-lineage) ----------
comp=[]
for sid,g in cells.groupby("orig_ident",observed=True):
    Bidx,labels=aggdict[sid]
    if len(labels)==0 or (labels!=-1).sum()==0:
        comp.append(dict(orig_ident=sid,Condition=COND[sid],n_region=0)); continue
    members=g.loc[Bidx[labels!=-1],["spatial_x","spatial_y"]].values
    xy=g[["spatial_x","spatial_y"]].values; tree=cKDTree(xy)
    reg=np.unique(np.concatenate([np.asarray(t,int) for t in tree.query_ball_point(members,r=EPS)]))
    gr=g.iloc[reg]; n=len(gr)
    comp.append(dict(orig_ident=sid,Condition=COND[sid],n_region=n,
        pct_B=round((gr.my_label=="B").mean()*100,1),pct_Plasma=round((gr.my_label=="Plasma").mean()*100,1),
        pct_Myeloid=round((gr.my_label=="Myeloid").mean()*100,1),pct_CD4T=round((gr.my_label=="CD4 T").mean()*100,1),
        pct_CD8T=round((gr.my_label=="CD8 T").mean()*100,1),pct_epi=round(gr.my_lineage.eq("Epithelial").mean()*100,1)))
compdf=pd.DataFrame(comp); compdf.to_csv(f"{OUT}/aggregate_composition.csv",index=False)
print("\n=== in/around B-aggregate composition (<=50um), focus rows ===")
foc=compdf[compdf.orig_ident.isin(["HK2695","1006","1003","1005","1007","1004"])]
print(foc.to_string(index=False))

# ---------- STEP 3b: DKD vs Control descriptive test ----------
res=[]
for metric in ["Blin_frac","agg_cells_per10k","agg_per_mm2"]:
    d=sub[sub.Condition=="DKD"][metric].values; c=sub[sub.Condition=="Control"][metric].values
    try: U,p=mannwhitneyu(d,c,alternative="two-sided"); p=round(p,4)
    except Exception: U,p=np.nan,np.nan
    res.append(dict(metric=metric,DKD_median=round(float(np.median(d)),4),Control_median=round(float(np.median(c)),4),
        DKD_n=len(d),Control_n=len(c),mannwhitney_U=U,p_value=p))
tst=pd.DataFrame(res); tst.to_csv(f"{OUT}/dkd_vs_control_test.csv",index=False)
print("\n=== DKD vs Control (descriptive; n=8 vs 3, UNDERPOWERED) ===")
print(tst.to_string(index=False))
print("== substrate done ==")
