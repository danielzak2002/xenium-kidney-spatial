#!/usr/bin/env python
"""Unbiased all-pairs neighbourhood-enrichment screen (squidpy) on the DKD Xenium atlas.

The hypothesis-FREE co-localization view we had not run. Per PHYSICAL SECTION: Delaunay graph pruned
to <=50 µm edges (within-tissue adjacency; the cLN multi-core lesson), then sq.gr.nhood_enrichment
over the validated reannotation subtype labels. We read SIGN + cross-section REPRODUCIBILITY, never
z-magnitude across sections. squidpy's z is computed against a label-permutation (abundance-preserving)
null; we add an explicit per-section label-SHUFFLE null and a spillover cross-check (06's flags).

RIGOR: unit = section; z never compared across sections; per-section graph; spillover + abundance
cross-checks mandatory; n=1 non-DKD = descriptive one-offs; iPT/iTAL/CD8 recall caveats carry. An honest
null (off-diagonal structure is abundance/adjacency-driven) is a first-class outcome. Reuse labels only;
memory-safe (labels+coords from cells.parquet — the 8.7 GB h5ad is never opened). Read-only raw.
"""
import os, warnings, time
warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd, anndata as ad, squidpy as sq, scipy.sparse as sp
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
RE=f"{REPO}/analysis/dkd_xenium_reannotation"; DIS=f"{REPO}/analysis/dkd_xenium_disease"; FIG=f"{DIS}/figures"
os.makedirs(FIG,exist_ok=True)
EDGE=50.0; MINC=20; ZT=2.0; NPERM=1000; SEED=0
rng=np.random.default_rng(SEED)

# fixed label order (immune | tubular epi | glomerular | stroma/endo) + walkthrough palette for ticks
CATS=["B","Plasma","Myeloid","CD4 T","CD8 T","PT","iPT","TAL","iTAL","DCT","CNT","PC","IC A",
      "Podo","MC","PEC","EC_glom","Fibroblast","VSMC","EC_Peritub","EC_DVR"]
PAL={"B":"#1f77b4","Plasma":"#ff7f0e","Myeloid":"#2ca02c","CD4 T":"#7d3ac1","CD8 T":"#e7298a",
     "PT":"#8a8a8a","iPT":"#c79a5b","TAL":"#8a8a8a","iTAL":"#c79a5b","DCT":"#8a8a8a","CNT":"#8a8a8a",
     "PC":"#8a8a8a","IC A":"#8a8a8a","Podo":"#e7559a","MC":"#e7559a","PEC":"#e7559a","EC_glom":"#e7559a",
     "Fibroblast":"#3fb6a8","VSMC":"#c8a91f","EC_Peritub":"#8c6d5c","EC_DVR":"#8c6d5c"}
GROUPS={"Control":["HK2753","HK3106","HK3626"],
        "DKD":["1001","1006","1008","1010","1011","1012","1013","HK2695"],
        "IgAN":["1003"],"MN":["1005"],"AA amyloid":["1004","1009"],"C3GN":["1007"]}
samp2grp={s:g for g,ss in GROUPS.items() for s in ss}
ORDER=["HK2753","HK3106","HK3626","1001","1006","1008","1010","1011","1012","1013","HK2695",
       "1003","1005","1004","1009","1007"]
# 06 spillover flags (per-section; only a subset assessed)
SF=pd.read_csv(f"{REPO}/analysis/dkd_epi_endo_stress/coloc_spillover_flags.csv")
SPILL={str(r.section):(r.spillover_nonPT_ambient=="SPILLOVER") for _,r in SF.iterrows()}

# known pairs to EXCLUDE from "new" reporting (already established elsewhere)
def key(a,b): return tuple(sorted((a,b)))
KNOWN=set()
for a in ["B","Plasma"]:
    for b in ["B","Plasma","Myeloid","CD4 T","CD8 T"]: KNOWN.add(key(a,b))   # B-aggregate composition
for a in ["iPT","iTAL"]: KNOWN.add(key(a,"Myeloid"))                          # injury~myeloid (06)

c=pd.read_parquet(f"{RE}/cells.parquet"); c["sample"]=c.orig_ident.astype(str)
c=c[c.my_label!="Unresolved"]

# ---------------- per-section: graph + nhood_enrichment + label-shuffle null ----------------
Z={}; ZN={}; deg={}; ncell={}; counts={}
for s in ORDER:
    d=c[c["sample"]==s]; XY=d[["spatial_x","spatial_y"]].values; n=len(d)
    a=ad.AnnData(np.zeros((n,1),np.float32),
                 obs=pd.DataFrame({"ct":pd.Categorical(d.my_label.values,categories=CATS)}))
    a.obsm["spatial"]=XY
    sq.gr.spatial_neighbors(a,coord_type="generic",delaunay=True)
    Dco=a.obsp["spatial_distances"].tocoo(); k=Dco.data<=EDGE
    a.obsp["spatial_connectivities"]=sp.coo_matrix((np.ones(k.sum()),(Dco.row[k],Dco.col[k])),
                                                   shape=(n,n)).tocsr()
    deg[s]=a.obsp["spatial_connectivities"].getnnz()/n; ncell[s]=n
    cc=d.my_label.value_counts(); counts[s]={t:int(cc.get(t,0)) for t in CATS}
    sq.gr.nhood_enrichment(a,cluster_key="ct",n_perms=NPERM,seed=SEED,show_progress_bar=False)
    z=a.uns["ct_nhood_enrichment"]["zscore"].astype(float)
    # label-shuffle null on the SAME graph (abundance-preserving): expect ~0
    a.obs["ctn"]=pd.Categorical(rng.permutation(d.my_label.values),categories=CATS)
    sq.gr.nhood_enrichment(a,cluster_key="ctn",n_perms=NPERM,seed=SEED,show_progress_bar=False)
    zn=a.uns["ctn_nhood_enrichment"]["zscore"].astype(float)
    # mask types with < MINC cells this section -> NaN
    bad=[i for i,t in enumerate(CATS) if counts[s][t]<MINC]
    for M in (z,zn):
        M[bad,:]=np.nan; M[:,bad]=np.nan
    Z[s]=z; ZN[s]=zn
    print(f"  [{samp2grp[s]:10s}] {s:7s} n={n:>6} deg={deg[s]:.1f} "
          f"types>={MINC}={len(CATS)-len(bad)}/{len(CATS)}")

def sign(M): return np.where(np.isnan(M),np.nan,np.where(M>ZT,1.0,np.where(M<-ZT,-1.0,0.0)))
S={s:sign(Z[s]) for s in ORDER}; SN={s:sign(ZN[s]) for s in ORDER}

# ---------------- per-section z matrices -> long CSV ----------------
rows=[]
iu=[(i,j) for i in range(len(CATS)) for j in range(i+1,len(CATS))]
for s in ORDER:
    for i,j in iu:
        rows.append(dict(section=s,group=samp2grp[s],spillover=SPILL.get(s,"unassessed"),
            type_a=CATS[i],type_b=CATS[j],z=round(Z[s][i,j],2),sign=S[s][i,j],
            z_null=round(ZN[s][i,j],2),n_a=counts[s][CATS[i]],n_b=counts[s][CATS[j]]))
LONG=pd.DataFrame(rows); LONG.to_csv(f"{DIS}/nhood_per_section_long.csv",index=False)
print(f"\nwrote nhood_per_section_long.csv ({len(LONG)} pair-rows)")

# ---------------- cross-disease SIGN/reproducibility per pair ----------------
def grp_sign(pair_ij,grp,Smaps):
    i,j=pair_ij; ss=[s for s in GROUPS[grp] if not np.isnan(Smaps[s][i,j])]
    if not ss: return (0,0,0,np.nan)   # n_usable, n_enr, n_avo, dominant
    v=np.array([Smaps[s][i,j] for s in ss])
    nE=int((v==1).sum()); nA=int((v==-1).sum())
    dom=1.0 if nE>nA else (-1.0 if nA>nE else 0.0)
    return (len(ss),nE,nA,dom)
agg=[]
for i,j in iu:
    a,b=CATS[i],CATS[j]; rec=dict(type_a=a,type_b=b,known=key(a,b) in KNOWN)
    nu,nE,nA,dom=grp_sign((i,j),"DKD",S); rec.update(DKD_n=nu,DKD_enr=nE,DKD_avo=nA,DKD_dom=dom,
        DKD_repro=round(max(nE,nA)/nu,2) if nu else np.nan)
    cu,cE,cA,cdom=grp_sign((i,j),"Control",S); rec.update(Ctrl_n=cu,Ctrl_enr=cE,Ctrl_avo=cA,Ctrl_dom=cdom)
    for one in ["IgAN","MN","C3GN"]:
        sname=GROUPS[one][0]; rec[one]=S[sname][i,j]
    nu2,e2,a2,aad=grp_sign((i,j),"AA amyloid",S); rec["AA"]=aad
    # null reproducibility (abundance baseline): max same-sign fraction under label-shuffle, DKD
    nun,nEn,nAn,_=grp_sign((i,j),"DKD",SN); rec["DKD_null_repro"]=round(max(nEn,nAn)/nun,2) if nun else np.nan
    agg.append(rec)
AGG=pd.DataFrame(agg)

# candidate = DKD reproducible (>=6/8 usable same NONZERO sign) AND Control->DKD shift
def shift(r):
    if r.DKD_dom==0 or r.DKD_n<6 or r.DKD_repro<0.75: return False
    # control shifts: control dominant sign differs from DKD, or control not reproducible (<2/3)
    ctrl_repro = max(r.Ctrl_enr,r.Ctrl_avo)/r.Ctrl_n if r.Ctrl_n else 0
    return (r.Ctrl_dom!=r.DKD_dom) or (ctrl_repro<0.67)
AGG["candidate"]=AGG.apply(shift,axis=1)
AGG["offtarget"]=~AGG.known
AGG.to_csv(f"{DIS}/nhood_pair_summary.csv",index=False)

# ---------------- HONEST classification of every OFF-TARGET reproducible pair ----------------
# The unbiased screen mostly recovers tissue ARCHITECTURE, not new niches. Classify each off-target
# pair so the report cannot mistake architecture/constitutive adjacency for a disease-specific finding.
IMMUNE={"B","Plasma","Myeloid","CD4 T","CD8 T"}
EPI_STRUCT={"PT","iPT","TAL","iTAL","DCT","CNT","PC","IC A","Podo","MC","PEC","EC_glom",
            "Fibroblast","VSMC","EC_Peritub","EC_DVR"}
def spill_check(i,j):
    clean=[S[s][i,j] for s in GROUPS["DKD"] if SPILL.get(s)==False and not np.isnan(S[s][i,j])]
    flag =[S[s][i,j] for s in GROUPS["DKD"] if SPILL.get(s)==True  and not np.isnan(S[s][i,j])]
    return clean,flag
def classify(r):
    a,b,dom=r.type_a,r.type_b,r.DKD_dom; nimm=(a in IMMUNE)+(b in IMMUNE)
    if dom==0 or r.DKD_n<6 or r.DKD_repro<0.75: return "no reproducible DKD sign"
    if dom>0 and r.Ctrl_dom==dom:               return "constitutive adjacency (also in controls)"
    if dom<0 and nimm==1 and ((a in EPI_STRUCT) or (b in EPI_STRUCT)):
                                                 return "compartmentalization (immune avoids epi/structure)"
    if dom<0 and nimm==2:                        return "immune–immune avoidance (architecture)"
    if dom>0 and nimm==0:                        return "nephron anatomy (segment/vascular neighbours)"
    if r.Ctrl_dom==dom:                          return "control-underpowered (same direction, weaker)"
    return "DISEASE-SPECIFIC candidate (review)"
OT=AGG[~AGG.known].copy(); OT["class"]=OT.apply(classify,axis=1)
# disease-specific survivors must also be sign-consistent in CLEAN DKD sections + null not reproducible
def survives(r):
    if r["class"]!="DISEASE-SPECIFIC candidate (review)": return False
    i,j=CATS.index(r.type_a),CATS.index(r.type_b); cs,fs=spill_check(i,j)
    clean_ok=(len(cs)>=2 and any(x==r.DKD_dom for x in cs) and all(x in (0,r.DKD_dom) for x in cs))
    return clean_ok and (r.DKD_null_repro<0.6)
OT["disease_specific_survivor"]=OT.apply(survives,axis=1)
OT["dom"]=np.where(OT.DKD_dom>0,"enriched",np.where(OT.DKD_dom<0,"avoided","~0"))
OT=OT.sort_values(["disease_specific_survivor","class","dom"],ascending=[False,True,True])
OT[["type_a","type_b","dom","DKD_n","DKD_repro","Ctrl_dom","DKD_null_repro","class","disease_specific_survivor"]
   ].to_csv(f"{DIS}/nhood_offtarget_classified.csv",index=False)
print("\n=== OFF-TARGET reproducible structure — honest classification (counts) ===")
print(OT[OT["class"]!="no reproducible DKD sign"]["class"].value_counts().to_string())
survivors=[(r.type_a,r.type_b) for _,r in OT.iterrows() if r.disease_specific_survivor]
print(f"\nDISEASE-SPECIFIC off-target survivors (beyond known B-aggregate / injury~myeloid): "
      f"{survivors if survivors else 'NONE — honest null'}")
print("\nillustrative constitutive (also-in-control) immune adjacencies the screen recovers:")
for a,b in [("Myeloid","Fibroblast"),("Myeloid","CD4 T"),("Myeloid","CD8 T"),("CD4 T","CD8 T")]:
    r=OT[((OT.type_a==a)&(OT.type_b==b))|((OT.type_a==b)&(OT.type_b==a))]
    if len(r): rr=r.iloc[0]; print(f"  {a}×{b}: DKD enriched {int(max(rr.DKD_enr,rr.DKD_avo))}/{int(rr.DKD_n)}, control enriched {int(rr.Ctrl_enr)}/{int(rr.Ctrl_n)} -> constitutive")

print("\n=== graph degree per section (post-prune) ===")
print(pd.DataFrame({"section":ORDER,"group":[samp2grp[s] for s in ORDER],
    "n_cells":[ncell[s] for s in ORDER],"degree":[round(deg[s],1) for s in ORDER]}).to_string(index=False))

np.save(f"{DIS}/nhood_Z.npy",{s:Z[s] for s in ORDER},allow_pickle=True)
print(f"\nsurvivors: {survivors}")
print("== screen core done; see figures step ==")

# ============================ FIGURES ============================
def heat(ax,M,title,vlim=8):
    im=ax.imshow(M,cmap="RdBu_r",vmin=-vlim,vmax=vlim)
    ax.set_xticks(range(len(CATS))); ax.set_yticks(range(len(CATS)))
    ax.set_xticklabels(CATS,rotation=90,fontsize=6); ax.set_yticklabels(CATS,fontsize=6)
    for tl,t in zip(ax.get_xticklabels(),CATS): tl.set_color(PAL[t])
    for tl,t in zip(ax.get_yticklabels(),CATS): tl.set_color(PAL[t])
    ax.set_title(title,fontsize=10); return im

# per-section heatmaps (representative: one per group + both DKD subgroups)
REP=[("HK3626","Control"),("1006","DKD B-rich"),("1008","DKD B-poor"),
     ("1003","IgAN n=1"),("1005","MN n=1"),("1007","C3GN n=1"),("1004","AA n=1"),("HK2695","DKD B-rich*")]
fig,axes=plt.subplots(2,4,figsize=(21,11))
for ax,(s,ttl) in zip(axes.ravel(),REP):
    im=heat(ax,Z[s],f"{s} · {ttl}  (deg {deg[s]:.1f})")
fig.suptitle("Per-section neighbourhood-enrichment z (squidpy; Delaunay pruned <=50µm). SIGN is read, NOT magnitude; z is NOT compared across sections. *spillover-flagged",fontsize=12)
fig.colorbar(im,ax=axes,fraction=0.013,pad=0.01,label="z (within-section only)")
fig.savefig(f"{FIG}/nhood_per_section_heatmaps.png",dpi=150,bbox_inches="tight"); plt.close(fig)
print("  [fig] nhood_per_section_heatmaps.png")

# DKD vs Control SIGN-reproducibility differential (signed fraction), + shift
def signed_frac(grp,Smaps):
    M=np.full((len(CATS),len(CATS)),np.nan)
    for i,j in iu:
        ss=[s for s in GROUPS[grp] if not np.isnan(Smaps[s][i,j])]
        if not ss: continue
        v=np.array([Smaps[s][i,j] for s in ss]); f=(np.sum(v==1)-np.sum(v==-1))/len(ss)
        M[i,j]=f; M[j,i]=f
    return M
DKDf=signed_frac("DKD",S); CTf=signed_frac("Control",S)
fig,axes=plt.subplots(1,3,figsize=(22,7))
for ax,(M,ttl) in zip(axes,[(CTf,"Control sign-reproducibility (n=3)"),(DKDf,"DKD sign-reproducibility (n=8)"),
                            (DKDf-CTf,"Shift: DKD − Control")]):
    im=ax.imshow(M,cmap="RdBu_r",vmin=-1,vmax=1)
    ax.set_xticks(range(len(CATS))); ax.set_yticks(range(len(CATS)))
    ax.set_xticklabels(CATS,rotation=90,fontsize=6); ax.set_yticklabels(CATS,fontsize=6)
    for tl,t in zip(ax.get_xticklabels(),CATS): tl.set_color(PAL[t])
    for tl,t in zip(ax.get_yticklabels(),CATS): tl.set_color(PAL[t])
    ax.set_title(ttl,fontsize=11)
fig.suptitle("SIGN reproducibility (signed fraction of sections enriched − avoided), NOT z-magnitude. Right = the honest Control→DKD shift. Diagonal omitted.",fontsize=12)
fig.colorbar(im,ax=axes,fraction=0.012,pad=0.01,label="enriched(+1) ↔ avoided(−1) fraction")
fig.savefig(f"{FIG}/nhood_dkd_vs_control_sign.png",dpi=160,bbox_inches="tight"); plt.close(fig)
print("  [fig] nhood_dkd_vs_control_sign.png")

# ============================ co_occurrence on survivors ============================
COOC=[]
if survivors:
    SUBN=15000; iv=np.linspace(0,200,21)
    for s in ORDER:
        d=c[c["sample"]==s]
        if len(d)>SUBN: d=d.sample(SUBN,random_state=SEED)
        a=ad.AnnData(np.zeros((len(d),1),np.float32),
                     obs=pd.DataFrame({"ct":pd.Categorical(d.my_label.values,categories=CATS)}))
        a.obsm["spatial"]=d[["spatial_x","spatial_y"]].values
        try:
            sq.gr.co_occurrence(a,cluster_key="ct",interval=iv,show_progress_bar=False)
            occ=a.uns["ct_co_occurrence"]["occ"]; ivc=(iv[:-1]+iv[1:])/2
            for ta,tb in survivors:
                ia,ib=CATS.index(ta),CATS.index(tb)
                if counts[s][ta]<MINC or counts[s][tb]<MINC: continue
                curve=occ[ia,ib,:]
                pk=int(np.nanargmax(curve));
                # length-scale: last distance where ratio still >1 (co-localized)
                above=np.where(curve>1.0)[0]; ls=ivc[above.max()] if len(above) else np.nan
                COOC.append(dict(section=s,group=samp2grp[s],pair=f"{ta}×{tb}",
                    peak_dist=round(float(ivc[pk]),0),peak_ratio=round(float(curve[pk]),2),
                    colocalized_to_um=round(float(ls),0) if not np.isnan(ls) else np.nan))
        except Exception as e:
            print(f"   co_occurrence skip {s}: {e}")
    CO=pd.DataFrame(COOC); CO.to_csv(f"{DIS}/nhood_cooccurrence_survivors.csv",index=False)
    print(f"\n=== co_occurrence length-scales (survivors), median by group ===")
    if len(CO):
        print(CO.groupby(["pair","group"]).agg(peak_um=("peak_dist","median"),
            coloc_um=("colocalized_to_um","median"),n=("section","size")).to_string())
else:
    print("\nno surviving off-target candidates -> co_occurrence skipped (honest null is the result)")
print("\n== nhood_enrichment_screen done ==")
