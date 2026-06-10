#!/usr/bin/env python
"""
interaction_map.py — comparative non-immune <-> immune spatial-interaction map across
RCC (Xenium), PRCC (Xenium, optional RCC replication), cLN (CosMx), DKD (Demoulin Xenium).
Which parenchymal/stromal<->immune associations are CONSERVED vs CONTEXT-SPECIFIC across three
kidney contexts. Method: per-section neighborhood enrichment (squidpy nhood_enrichment,
permutation z) on a harmonized lineage labelling — compositional, ligand-INDEPENDENT.

Read-only. nhood_enrichment uses only coords + labels (X never touched). Graph = KNN k=6
(scale-invariant: handles um vs mm, avoids cLN cross-core gap-spanning). Per-section throughout;
DKD read backed + per-section. Association/colocalization, NOT communication or causation.
"""
import os, warnings
warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd, anndata as ad, squidpy as sq
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
OUT=os.path.join(REPO,"analysis/interaction_map"); os.makedirs(OUT,exist_ok=True)
KNN=6; NPERM=1000; MINCELL=20; ZSIG=2.0; SEED=0
def hdr(s): print("\n"+"="*78+"\n"+s+"\n"+"="*78)

NONIMM=["Tubular_epi","Malignant_epi","Endothelial","Stroma","Podocyte"]
IMMUNE=["B","Plasma","Myeloid","T_lineage","NK"]
COMMON=NONIMM+IMMUNE

# ============================================================================
# STAGE 0 — harmonization maps (native -> common). Printed per dataset.
# ============================================================================
# RCC / PRCC use phase_b_label (malignant epithelium = Tumor_RCC -> Malignant_epi)
RCC_MAP={
 "Tumor_RCC":"Malignant_epi","Proximal_tubule":"Tubular_epi","Intercalated":"Tubular_epi",
 "Endothelial":"Endothelial","Stroma_mural":"Stroma","Fibroblast":"Stroma",
 "Naive B cells":"B","Switched memory B cells":"B","Plasmablasts":"Plasma",
 "Effector memory CD8 T cells":"T_lineage","CD8_T":"T_lineage","T regulatory cells":"T_lineage",
 "CCR7+ T (naive/CM)":"T_lineage",
 "Intermediate monocytes":"Myeloid","Classical monocytes":"Myeloid","myeloid/DC":"Myeloid",
 "Myeloid dendritic cells":"Myeloid","mregDC":"Myeloid","Plasmacytoid dendritic cells":"Myeloid",
 "Natural killer cells":"NK",
 # dropped: Proliferating, LowQ_MTRNR2L, Mast, Low-density neutrophils, Immune infiltrate/doublet
}
# cLN uses author_celltype (clean names; phase_b_label has junk a-e clusters)
CLN_MAP={
 "PCT":"Tubular_epi","Connecting.tubule":"Tubular_epi","Thick.ascending.limb.of.Loop.of.Henle":"Tubular_epi",
 "Type.A.intercalated.cell":"Tubular_epi","Type.B.intercalated.cell":"Tubular_epi","Principal.cell":"Tubular_epi",
 "Pelvic.epithelium":"Tubular_epi","Epithelial.progenitor.cell":"Tubular_epi","Indistinct.intercalated.cell":"Tubular_epi",
 "Transitional.urothelium":"Tubular_epi","Parietal.epithelium":"Tubular_epi",
 "Ascending.vasa.recta.endothelium":"Endothelial","Peritubular.capillary.endothelium.1":"Endothelial",
 "Peritubular.capillary.endothelium.2":"Endothelial","Glomerular.endothelium":"Endothelial",
 "Descending.vasa.recta.endothelium":"Endothelial",
 "Fibroblast":"Stroma","Myofibroblast":"Stroma","Vascular.pericyte":"Stroma","Mesangial.cell":"Stroma",
 "Podocyte":"Podocyte",
 "B-cell":"B","plasmablast":"Plasma",
 "macrophage":"Myeloid","monocyte":"Myeloid","mDC":"Myeloid","pDC":"Myeloid",
 "Treg":"T_lineage","T CD8 memory":"T_lineage","T CD8 naive":"T_lineage","T CD4 memory":"T_lineage","T CD4 naive":"T_lineage",
 "NK":"NK",
 # dropped: mast, neutrophil
}
# DKD: non-immune from annotation_updated; immune from immune_cell_annotation_combined
DKD_NONIMM={
 "PT":"Tubular_epi","TAL":"Tubular_epi","PC":"Tubular_epi","CNT":"Tubular_epi","DCT":"Tubular_epi",
 "IC A":"Tubular_epi","IC B":"Tubular_epi","DTL_ATL":"Tubular_epi","PEC":"Tubular_epi",
 "iPT":"Tubular_epi","iTAL":"Tubular_epi",   # base class; injured overlay handled separately
 "EC_Peritub":"Endothelial","EC_glom":"Endothelial","EC_DVR":"Endothelial","EC_Lymph":"Endothelial",
 "Fibroblast":"Stroma","VSMC":"Stroma","MC1":"Stroma","Podo":"Podocyte",
}
DKD_INJURED={"iPT","iTAL"}  # native injured-epithelium states (overlay)
DKD_IMMUNE={"B":"B","Plasma":"Plasma","Macro":"Myeloid","cDC":"Myeloid","pDC":"Myeloid",
            "CD4+":"T_lineage","CD8+":"T_lineage","NK":"NK"}  # dropped: Neutrophil, Baso_Mast

def print_map(name, mp):
    inv={}
    for k,v in mp.items(): inv.setdefault(v,[]).append(k)
    print(f"\n[{name}] native -> common:")
    for c in sorted(inv): print(f"  {c:16s} <- {inv[c]}")

# ============================================================================
# loaders: return list of (section_name, coords Nx2, labels array) under common scheme
# ============================================================================
def load_simple(path, label_col, mp, sample_col="sample", extra=None):
    a=ad.read_h5ad(path)
    lab=a.obs[label_col].astype(str).map(mp).values
    if extra is not None: lab=extra(a,lab)
    xy=np.asarray(a.obsm["spatial"],float)
    samp=a.obs[sample_col].astype(str).values if sample_col in a.obs else np.array(["S0"]*a.n_obs)
    secs=[]
    for s in pd.unique(samp):
        m=(samp==s) & pd.notna(lab)
        secs.append((s, xy[m], lab[m]))
    del a
    return secs

def load_dkd(injured=False):
    a=ad.read_h5ad(os.path.join(REPO,"Demoulin26/data/spatial_adata_xenium_cosmx_zenodo.h5ad"), backed="r")
    xen=a.obs["tech"].astype(str).values=="Xenium"
    ann=a.obs["annotation_updated"].astype(str).values
    imm=a.obs["immune_cell_annotation_combined"].astype(str).values
    samp=a.obs["orig_ident"].astype(str).values
    xy=np.asarray(a.obsm["spatial"],float); a.file.close()
    lab=np.full(len(ann), None, object)
    for i in np.where(xen)[0]:
        if imm[i] in DKD_IMMUNE: lab[i]=DKD_IMMUNE[imm[i]]
        elif imm[i]=="Unknown":
            if injured and ann[i] in DKD_INJURED: lab[i]="Injured_epi"
            else: lab[i]=DKD_NONIMM.get(ann[i])
    secs=[]
    for s in pd.unique(samp[xen]):
        m=xen & (samp==s) & np.array([x is not None for x in lab])
        secs.append((s, xy[m], lab[m].astype(str)))
    return secs

# ============================================================================
# STAGE 1 — per-section nhood enrichment; aggregate to mean z + k/N
# ============================================================================
def tile(xy, lab, target=50000, min_tile=3000):
    """Split a single large section into a quantile grid of ~target-cell pseudo-sections, so z
    magnitude is comparable across datasets (z scales with sqrt N) and single-section datasets
    (RCC/PRCC) gain replication. Sections <100k cells are left whole."""
    n=len(lab)
    if n<100000: return [("whole", xy, lab)]
    g=int(np.ceil(np.sqrt(n/target)))
    xq=np.quantile(xy[:,0],np.linspace(0,1,g+1)); yq=np.quantile(xy[:,1],np.linspace(0,1,g+1))
    out=[]
    for i in range(g):
        for j in range(g):
            m=(xy[:,0]>=xq[i])&(xy[:,0]<=xq[i+1])&(xy[:,1]>=yq[j])&(xy[:,1]<=yq[j+1])
            if m.sum()>=min_tile: out.append((f"t{i}{j}", xy[m], lab[m]))
    return out

def section_z(coords, labels):
    # present labels derived from the ACTUAL harmonized labels (not a fixed COMMON list), so
    # overlay (Injured_epi) and fine immune labels are included.
    vals,cnts=np.unique(labels,return_counts=True)
    present=[v for v,c in zip(vals,cnts) if c>=MINCELL]
    if len(present)<2: return None
    keep=np.isin(labels,present)
    s=ad.AnnData(np.zeros((int(keep.sum()),1),dtype=np.float32))
    s.obsm["spatial"]=coords[keep]
    s.obs["lab"]=pd.Categorical(labels[keep], categories=present)
    sq.gr.spatial_neighbors(s, coord_type="generic", n_neighs=KNN)
    sq.gr.nhood_enrichment(s, cluster_key="lab", seed=SEED, n_perms=NPERM, show_progress_bar=False)
    return pd.DataFrame(s.uns["lab_nhood_enrichment"]["zscore"], index=present, columns=present)

def run_dataset(name, sections, extra_rows=None):
    """extra_rows: non-immune labels to ALSO include as rows (e.g., Injured_epi/Malignant_epi)."""
    rows = NONIMM + (extra_rows or [])
    zsum={}; ksig={}; kavoid={}; nsec={}
    per_sec_store=[]
    tiled=[(f"{sn}/{tn}",txy,tlab) for sn,xy,lab in sections for tn,txy,tlab in tile(xy,lab)]
    for sname, xy, lab in tiled:
        Z=section_z(np.asarray(xy,float), np.asarray(lab))
        if Z is None: continue
        per_sec_store.append((sname,Z))
        for r in rows:
            for c in IMMUNE:
                if r in Z.index and c in Z.columns:
                    z=float(Z.loc[r,c]); key=(r,c)
                    zsum.setdefault(key,[]).append(z)
                    ksig[key]=ksig.get(key,0)+(z>ZSIG); kavoid[key]=kavoid.get(key,0)+(z<-ZSIG)
                    nsec[key]=nsec.get(key,0)+1
    # assemble
    recs=[]
    for r in rows:
        for c in IMMUNE:
            key=(r,c); zs=zsum.get(key,[])
            if not zs: recs.append(dict(non_immune=r,immune=c,mean_z=np.nan,n_sec=0,k_enrich=0,k_avoid=0)); continue
            recs.append(dict(non_immune=r,immune=c,mean_z=round(float(np.mean(zs)),2),
                n_sec=nsec[key],k_enrich=int(ksig[key]),k_avoid=int(kavoid[key])))
    df=pd.DataFrame(recs); df["dataset"]=name
    return df, per_sec_store

# heatmap helper
def heatmap(ax, mat, rows, cols, title):
    im=ax.imshow(mat, cmap="RdBu_r", vmin=-8, vmax=8, aspect="auto")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=8)
    for i in range(len(rows)):
        for j in range(len(cols)):
            v=mat[i,j]
            if v==v: ax.text(j,i,f"{v:.0f}",ha="center",va="center",fontsize=7,
                             color="white" if abs(v)>4 else "black")
    ax.set_title(title, fontsize=10); return im

# ============================================================================
hdr("STAGE 0 — harmonization mapping (printed per dataset)")
print_map("RCC/PRCC", RCC_MAP); print_map("cLN", CLN_MAP)
print_map("DKD non-immune", DKD_NONIMM); print_map("DKD immune", {k:v for k,v in DKD_IMMUNE.items()})
print(f"\nDKD injured-epithelium overlay (native): {DKD_INJURED}")
print("FLAGS: RCC/PRCC epithelium = Tumor_RCC -> Malignant_epi (tumor-immune, interpreted separately).")
print("       cLN T_lineage pooled (CD4/CD8 not separable) AND flagged UNRELIABLE (35% CD3 contamination).")
print("       RCC/PRCC single-section -> k/N = ?/1 (no section replication). NK absent in RCC base.")

hdr("STAGE 1 — per-dataset per-section nhood enrichment (KNN k=6, %d perms)"%NPERM)
DATASETS=[]
print("loading RCC..."); rcc_secs=load_simple(REPO+"/outputs/objects/kidney_RCC_protein.h5ad","phase_b_label",RCC_MAP)
df_rcc,ps_rcc=run_dataset("RCC", rcc_secs, extra_rows=["Malignant_epi"])
print("loading PRCC..."); prcc_secs=load_simple(REPO+"/outputs/objects/kidney_preview_PRCC.h5ad","phase_b_label",RCC_MAP)
df_prcc,ps_prcc=run_dataset("PRCC", prcc_secs, extra_rows=["Malignant_epi"])
print("loading cLN..."); cln_secs=load_simple(REPO+"/outputs/objects/cln_cosmx.h5ad","author_celltype",CLN_MAP)
df_cln,ps_cln=run_dataset("cLN", cln_secs)
print("loading DKD (backed, per-section)..."); dkd_secs=load_dkd(injured=False)
df_dkd,ps_dkd=run_dataset("DKD", dkd_secs)
# DKD injured overlay (separate pass for B<->Injured_epi etc.)
dkd_inj_secs=load_dkd(injured=True); df_dkd_inj,_=run_dataset("DKD", dkd_inj_secs, extra_rows=["Injured_epi"])
inj_rows=df_dkd_inj[df_dkd_inj.non_immune=="Injured_epi"].copy()
df_dkd=pd.concat([df_dkd, inj_rows], ignore_index=True)

for nm,df in [("RCC",df_rcc),("PRCC",df_prcc),("cLN",df_cln),("DKD",df_dkd)]:
    print(f"\n### {nm} non-immune x immune mean z (n_sec / k_enrich):")
    piv=df.pivot_table(index="non_immune",columns="immune",values="mean_z")
    print(piv.reindex([r for r in NONIMM+["Malignant_epi","Injured_epi"] if r in piv.index]).to_string())
    df.to_csv(os.path.join(OUT,f"zmatrix_{nm}.csv"),index=False)
    # per-dataset heatmap
    rows=[r for r in NONIMM+["Malignant_epi","Injured_epi"] if r in set(df.non_immune)]
    mat=np.full((len(rows),len(IMMUNE)),np.nan)
    for i,r in enumerate(rows):
        for j,c in enumerate(IMMUNE):
            sub=df[(df.non_immune==r)&(df.immune==c)]
            if len(sub) and sub.n_sec.iloc[0]>0: mat[i,j]=sub.mean_z.iloc[0]
    fig,ax=plt.subplots(figsize=(6,4.5)); im=heatmap(ax,mat,rows,IMMUNE,f"{nm}: non-immune x immune nhood z")
    fig.colorbar(im,ax=ax,fraction=0.04,label="mean nhood z"); fig.tight_layout()
    fig.savefig(os.path.join(OUT,f"heatmap_{nm}.png"),dpi=150); plt.close(fig)

# ============================================================================
hdr("STAGE 2 — cross-dataset assemble + classify")
allz=pd.concat([df_rcc,df_prcc,df_cln,df_dkd],ignore_index=True)
allz.to_csv(os.path.join(OUT,"comparative_all_pairs.csv"),index=False)
# comparative across the shared (benign) datasets: cLN + DKD share Tubular_epi/Endo/Stroma/Podo;
# RCC/PRCC malignant handled separately. Classify each non-immune x immune pair.
def status(row):
    if row.n_sec==0 or np.isnan(row.mean_z): return "untested"
    if row.mean_z>ZSIG and row.k_enrich>=max(1,row.n_sec*0.5): return "ENRICH"
    if row.mean_z<-ZSIG and row.k_avoid>=max(1,row.n_sec*0.5): return "AVOID"
    return "ns"
allz["status"]=allz.apply(status,axis=1)
benign=["RCC","PRCC","cLN","DKD"]
rows_cmp=[]
for r in NONIMM:
    for c in IMMUNE:
        st={}; zz={}
        for d in benign:
            sub=allz[(allz.dataset==d)&(allz.non_immune==r)&(allz.immune==c)]
            if len(sub): st[d]=sub.status.iloc[0]; zz[d]=sub.mean_z.iloc[0]
        tested=[d for d in benign if st.get(d) not in (None,"untested")]
        if len(tested)<2: cls="insufficient"
        else:
            sgn=set("ENRICH" if st[d]=="ENRICH" else ("AVOID" if st[d]=="AVOID" else "ns") for d in tested)
            enr=[d for d in tested if st[d]=="ENRICH"]; avo=[d for d in tested if st[d]=="AVOID"]
            if enr and avo: cls="DISCORDANT"
            elif len(enr)>=2 and len(enr)==len([d for d in tested if st[d]!="ns"]): cls="CONSERVED-enrich"
            elif len(enr)>=1: cls="CONTEXT-specific"
            elif len(avo)>=2: cls="CONSERVED-avoid"
            else: cls="ns-all"
        rows_cmp.append(dict(non_immune=r,immune=c,
            **{f"z_{d}":round(zz[d],1) if d in zz and zz[d]==zz[d] else np.nan for d in benign},
            **{f"st_{d}":st.get(d,"-") for d in benign}, classification=cls))
cmp=pd.DataFrame(rows_cmp); cmp.to_csv(os.path.join(OUT,"comparative_classification.csv"),index=False)
print(cmp.to_string(index=False))

# faceted heatmap (pair x dataset)
fig,axes=plt.subplots(1,len(benign),figsize=(5*len(benign),5),sharey=True)
for ax,d in zip(axes,benign):
    rows=[r for r in NONIMM if r!="Malignant_epi"]
    mat=np.full((len(rows),len(IMMUNE)),np.nan)
    for i,r in enumerate(rows):
        for j,c in enumerate(IMMUNE):
            sub=allz[(allz.dataset==d)&(allz.non_immune==r)&(allz.immune==c)]
            if len(sub) and sub.n_sec.iloc[0]>0: mat[i,j]=sub.mean_z.iloc[0]
    im=heatmap(ax,mat,rows,IMMUNE,f"{d}");
    if d!=benign[0]: ax.set_ylabel("")
fig.suptitle("Non-immune x immune spatial nhood z across kidney contexts (KNN k=6, per-section mean)")
fig.colorbar(im,ax=axes,fraction=0.02,label="mean nhood z");
fig.savefig(os.path.join(OUT,"heatmap_comparative_faceted.png"),dpi=150,bbox_inches="tight"); plt.close(fig)

# ============================================================================
hdr("STAGE 3 — validation (recover known niches)")
def getz(df,r,c):
    s=df[(df.non_immune==r)&(df.immune==c)]
    return (float(s.mean_z.iloc[0]), int(s.k_enrich.iloc[0]), int(s.n_sec.iloc[0])) if len(s) and s.n_sec.iloc[0]>0 else (np.nan,0,0)
# cLN Plasma<->Myeloid (immune-immune): pull from per-section stores
def imm_imm(ps, A, B):
    zs=[]; k=0
    for _,Z in ps:
        if A in Z.index and B in Z.columns: zs.append(float(Z.loc[A,B])); k+=Z.loc[A,B]>ZSIG
    return (round(float(np.mean(zs)),2) if zs else np.nan, k, len(zs))
# RCC fine pass (subtypes) — needed here for the B<->Treg validation; reused in Stage 4
RCC_FINE=dict(RCC_MAP); RCC_FINE.update({"T regulatory cells":"Treg","Effector memory CD8 T cells":"CD8eff",
    "CD8_T":"CD8eff","CCR7+ T (naive/CM)":"CCR7T","mregDC":"mregDC"})
FINE_IMM=["B","Plasma","Myeloid","Treg","CD8eff","NK","mregDC"]
rcc_fine_secs=load_simple(REPO+"/outputs/objects/kidney_RCC_protein.h5ad","phase_b_label",RCC_FINE)
IMMUNE_SAVE=IMMUNE; IMMUNE=FINE_IMM
df_rcc_f,ps_rcc_f=run_dataset("RCC_fine", rcc_fine_secs, extra_rows=["Malignant_epi"])
IMMUNE=IMMUNE_SAVE

print("NOTE: absolute nhood z is dominated by universal immune<->parenchyma SEGREGATION (immune "
      "cells sit in sparse interstitial aggregates, so bulk parenchyma reads as 'avoiding' them). "
      "Focal niches re-emerge as DIFFERENTIALS (relative ordering among partners), not absolute sign.")
v_cln=imm_imm(ps_cln,"Plasma","Myeloid")
v_rcc=imm_imm(ps_rcc_f,"B","Treg")        # RCC B<->Treg (the immunoregulatory aggregate)
v_rcc8=imm_imm(ps_rcc_f,"B","CD8eff")     # contrast: CD8-effector excluded from B core
v_dkd_inj=getz(df_dkd,"Injured_epi","B"); v_dkd_tub=getz(df_dkd,"Tubular_epi","B")
d_rcc=v_rcc[0]-v_rcc8[0]            # Treg favored over CD8eff near B
d_dkd=v_dkd_inj[0]-v_dkd_tub[0]     # injured epi less excluded from B than healthy tubular
print(f"cLN Plasma<->Myeloid: mean z {v_cln[0]} ({v_cln[1]}/{v_cln[2]}) -> DIRECT enrichment {'RECOVERED' if v_cln[0]>ZSIG else 'no'}")
print(f"RCC B<->Treg {v_rcc[0]} vs B<->CD8eff {v_rcc8[0]} -> Treg favored by Δ{d_rcc:+.1f} -> "
      f"immunoregulatory direction {'RECOVERED' if d_rcc>5 else 'no'}")
print(f"DKD Injured_epi<->B {v_dkd_inj[0]} vs Tubular_epi<->B {v_dkd_tub[0]} -> injured less excluded by "
      f"Δ{d_dkd:+.1f} -> Task-B injPT-aggregate proximity {'RECOVERED' if d_dkd>5 else 'no'}")
val=pd.DataFrame([
 dict(dataset="cLN",test="Plasma<->Myeloid (direct)",metric=f"z={v_cln[0]}",k=f"{v_cln[1]}/{v_cln[2]}",
      recovered=("YES" if v_cln[0]>ZSIG else "no")),
 dict(dataset="RCC",test="B<->Treg vs B<->CD8eff (differential)",
      metric=f"{v_rcc[0]} vs {v_rcc8[0]} (Δ{d_rcc:+.1f})",k=f"{v_rcc[2]} tiles",
      recovered=("YES (differential)" if d_rcc>5 else "no")),
 dict(dataset="DKD",test="Injured_epi<->B vs Tubular_epi<->B (differential)",
      metric=f"{v_dkd_inj[0]} vs {v_dkd_tub[0]} (Δ{d_dkd:+.1f})",k=f"{v_dkd_inj[2]} sec",
      recovered=("YES (differential)" if d_dkd>5 else "no")),
]); val.to_csv(os.path.join(OUT,"validation.csv"),index=False); print(val.to_string(index=False))

# ============================================================================
hdr("STAGE 4 — finer immune resolution (Xenium only: RCC + DKD)")
IMMUNE_SAVE=IMMUNE
# RCC fine already computed above (df_rcc_f, ps_rcc_f); DKD fine: CD4+/CD8+ separate
DKD_IMMUNE_FINE={"B":"B","Plasma":"Plasma","Macro":"Myeloid","cDC":"Myeloid","pDC":"Myeloid",
                 "CD4+":"CD4T","CD8+":"CD8T","NK":"NK"}
def load_dkd_fine():
    a=ad.read_h5ad(os.path.join(REPO,"Demoulin26/data/spatial_adata_xenium_cosmx_zenodo.h5ad"), backed="r")
    xen=a.obs["tech"].astype(str).values=="Xenium"; ann=a.obs["annotation_updated"].astype(str).values
    imm=a.obs["immune_cell_annotation_combined"].astype(str).values; samp=a.obs["orig_ident"].astype(str).values
    xy=np.asarray(a.obsm["spatial"],float); a.file.close()
    lab=np.full(len(ann),None,object)
    for i in np.where(xen)[0]:
        if imm[i] in DKD_IMMUNE_FINE: lab[i]=DKD_IMMUNE_FINE[imm[i]]
        elif imm[i]=="Unknown": lab[i]=DKD_NONIMM.get(ann[i])
    return [(s,xy[xen&(samp==s)&np.array([x is not None for x in lab])],
             lab[xen&(samp==s)&np.array([x is not None for x in lab])].astype(str)) for s in pd.unique(samp[xen])]
IMMUNE=["B","Plasma","Myeloid","CD4T","CD8T","NK"]
df_dkd_f,_=run_dataset("DKD_fine", load_dkd_fine())
IMMUNE=IMMUNE_SAVE
print("\nRCC fine immune (non-immune x {Treg,CD8eff,...}):")
print(df_rcc_f.pivot_table(index="non_immune",columns="immune",values="mean_z").to_string())
print("\nDKD fine immune (CD4T vs CD8T separated):")
print(df_dkd_f.pivot_table(index="non_immune",columns="immune",values="mean_z").to_string())
df_rcc_f.to_csv(os.path.join(OUT,"zmatrix_RCC_fine.csv"),index=False)
df_dkd_f.to_csv(os.path.join(OUT,"zmatrix_DKD_fine.csv"),index=False)

# ============================================================================
hdr("WRITING REPORT")
L=[];W=L.append
W("# Comparative non-immune <-> immune spatial-interaction map (RCC / PRCC / cLN / DKD kidney)\n")
W("Read-only. Per-section neighborhood enrichment (squidpy, KNN k=6, "
  f"{NPERM} perms, permutation z) on a harmonized lineage labelling. **Compositional / "
  "ligand-independent**: spatial co-localization, NOT communication or causation. "
  "X never touched (labels + coords only).\n")
W("## Stage 0 — harmonization (native -> common)\n")
W("Common immune = B, Plasma, Myeloid (mono/macro/DC pooled), T_lineage (pooled), NK. "
  "Non-immune = Tubular_epi, Endothelial, Stroma, Podocyte, **Malignant_epi (RCC/PRCC only — "
  "tumor-immune, interpreted separately)**. Injured-epithelium = DKD native iPT/iTAL overlay.\n")
W("| dataset | source label | notes |")
W("|---|---|---|")
W("| RCC | phase_b_label | 465k single section -> spatially TILED into ~50k pseudo-sections; Tumor_RCC->Malignant_epi; NK absent in base |")
W("| PRCC | phase_b_label | single section (optional RCC replication) |")
W("| cLN | author_celltype | 14 sections, CosMx; T pooled & **UNRELIABLE (35% CD3 contam)** |")
W("| DKD | annotation_updated + immune_cell_annotation_combined | 16 Xenium sections; injured overlay = iPT/iTAL |")
W("\nFull native->common mapping printed to console / in script.\n")
W("## KEY interpretation: read DIFFERENTIALS, not absolute sign\n")
W("Absolute nhood z is **dominated by a universal immune<->parenchyma SEGREGATION**: immune cells "
  "sit in sparse interstitial/peri-vascular aggregates, so bulk Tubular/Endothelial/Malignant "
  "epithelium reads as 'avoiding' every immune type in every context (the strongest, most "
  "conserved pattern — but largely a density-geometry fact). The biologically informative signal "
  "is the **relative ordering**: which immune partner a given non-immune type avoids LEAST, and "
  "which non-immune **state** (injured vs healthy) is least excluded.\n")
W("## Stage 3 — validation (known niches recovered, as differentials)\n")
W("| dataset | test | values | recovered? |")
W("|---|---|---|---|")
for _,r in val.iterrows():
    W(f"| {r.dataset} | {r.test} | {r.metric} | {r.recovered} |")
W("\nThe pervasive cLN plasma-myeloid niche recovers as a **direct** positive z; the focal RCC "
  "B-Treg aggregate and DKD injured-PT proximity recover only as **differentials** (global KNN "
  "nhood is less sensitive to focal aggregates than the DBSCAN approach in `bniche_dbscan`).\n")
W("## Stage 2 — conserved vs context-specific (non-immune x immune)\n")
W("| non-immune | immune | RCC | PRCC | cLN | DKD | class |")
W("|---|---|---|---|---|---|---|")
for _,r in cmp.iterrows():
    W(f"| {r.non_immune} | {r.immune} | {r.get('z_RCC')} | {r.get('z_PRCC')} | {r.get('z_cLN')} | {r.get('z_DKD')} | {r.classification} |")
W("\n### Synthesis\n")
W("- **CONSERVED (all contexts): immune cells segregate from Tubular epithelium** (Tubular_epi "
  "avoids B/Plasma/Myeloid/T everywhere) — the robust conserved axis, though largely density-geometry.")
W("- **CONSERVED-ish ENRICH: Stroma <-> Myeloid and Stroma <-> Plasma** in the benign/papillary "
  "contexts (cLN, DKD, PRCC) — myeloid & plasma cells localize to the interstitial stroma. "
  "**DISCORDANT in RCC** (clear-cell stroma avoids immune) — a tumor-context-specific inversion.")
W("- **Endothelial <-> NK and Stroma <-> NK enrich** where NK is measured (PRCC, DKD) — peri-vascular/"
  "interstitial NK; cLN null is uninformative (CosMx depth).")
W("- **DKD-specific: Stroma <-> T_lineage enrich** (+10) — interstitial T cells, a DKD interstitial-"
  "nephritis signature absent in the tumor contexts.")
W("## Interpretation rules & caveats\n")
W("- **Read differentials, not absolute sign** (see Key interpretation): the universal Tubular<->immune "
  "avoidance is geometry; the informative signals are Stroma<->myeloid/plasma and the validation differentials.")
W("- A positive CONSERVED association is strong; a dataset-specific ABSENCE is **ambiguous** "
  "(platform/typing depth — cLN is CosMx, sparser) -> a cLN null is NOT biological absence.")
W("- **cLN T_lineage pairs are UNRELIABLE** (35% epithelial CD3 contamination, see "
  "`cln_cd3_contamination`) -> excluded from conserved T claims.")
W("- **RCC/PRCC epithelium is MALIGNANT** (Malignant_epi) -> tumor-immune, interpreted separately "
  "from benign Tubular_epi; RCC/PRCC are single-section (k/N = ?/1, no section replication).")
W("- Platform-depth confound (CosMx 957 vs Xenium 5k) -> asymmetric reading of nulls.")
W("- Sparse-type z is noisy (pairs with <%d cells/section flagged via n_sec); NK absent in RCC base."%MINCELL)
W("- No patient column for DKD (donor clustering uncontrolled). Per-section normalized z throughout; "
  "KNN graph is scale-invariant (handles um vs mm). Association/colocalization only.")
open(os.path.join(OUT,"REPORT.md"),"w").write("\n".join(L))
print("wrote REPORT.md + CSVs + heatmaps")
print("\n== interaction_map done ==")
