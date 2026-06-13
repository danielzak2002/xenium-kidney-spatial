#!/usr/bin/env python
"""Injury<->immune co-localization across diseases — compositional/descriptive addendum.
Builds on reconcile_extend.py (validated reannotation labels; B-lineage=B+Plasma DBSCAN eps=50/
minPts=20; near<=50/far>200; within-section per-compartment z; permutation null positions+count
fixed; injPT genes VCAM1/HAVCR1/PROM1/SPP1). Door-closing: does the B-rich subgroup track tubular
injury? -> tracks general, myeloid-led infiltration, not B-lineage specifically. Read-only.
NO per-gene DE (gated Stage 2). Unit = sample/section; honest nulls first-class."""
import os, warnings
warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree, ConvexHull
from scipy.stats import spearmanr
from matplotlib.patches import Polygon
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
RE=f"{REPO}/analysis/dkd_xenium_reannotation"; OUT=f"{REPO}/analysis/dkd_epi_endo_stress"
EPS,MINPTS,NEAR,FAR=50.0,20,50.0,200.0
INJ_NEAR,INJ_FAR=30.0,100.0      # immune-to-injured-tubule proximity (diffuse injury; reported)
rng=np.random.default_rng(0)
DISCOL={"DKD":"#6A3D9A","Control":"#888888","AA amyloid":"#FF7F00","C3GN":"#33A02C","IgA":"#1f77b4","MN":"#E31A1C"}
def hdr(s): print("\n"+"="*78+"\n"+s+"\n"+"="*78)
def save(fig,n): fig.savefig(f"{OUT}/{n}.png",dpi=170,bbox_inches="tight"); plt.close(fig); print("  [ok]",n)

cells=pd.read_parquet(f"{RE}/cells.parquet").reset_index(drop=True)
samp=cells.orig_ident.astype(str).values; cond=cells.groupby("orig_ident",observed=True).Condition.first()
ml=cells.my_label.values; lin=cells.my_lineage.values; XY=cells[["spatial_x","spatial_y"]].values
SAMPLES=sorted(pd.unique(samp))
epi=(lin=="Epithelial"); imm=(lin=="Immune")
isB=np.isin(ml,["B","Plasma"]); nonBimm=imm&~isB; ismye=(ml=="Myeloid")
isiPT=(ml=="iPT"); isiTAL=(ml=="iTAL"); isPT=(ml=="PT"); isTAL=(ml=="TAL")
injured=isiPT|isiTAL
OTHEREPI=epi&~np.isin(ml,["PT","iPT","TAL","iTAL"])

# ---------------- disease -> sample map + spillover flags (from reconcile spatial) ----------------
hdr("disease -> sample map (all 16) + per-section spillover (non-PT ambient) flags")
dmap=pd.DataFrame({"sample":SAMPLES,"Condition":[cond[s] for s in SAMPLES]})
print(dmap.groupby("Condition")["sample"].apply(list).to_string())
sp=pd.read_csv(f"{OUT}/reconcile_spatial_per_section.csv")
sp["spillover_nonPT_ambient"]=np.where(sp.nonPT_amb_dz>0,"SPILLOVER","clean")
print("\nper-section spillover flag:")
print(sp[["section","Condition","nonPT_amb_dz","spillover_nonPT_ambient"]].to_string(index=False))
sp.to_csv(f"{OUT}/coloc_spillover_flags.csv",index=False)

# ============================ PANEL A — cross-sample specificity ============================
hdr("PANEL A — de-circularized partial + CLR closure sensitivity + foreground myeloid")
P=pd.read_csv(f"{OUT}/reconcile_per_sample.csv")
# add non-B immune fraction
nbf={s:(nonBimm&(samp==s)).sum()/(samp==s).sum() for s in SAMPLES}
P["nonBimm_frac"]=P["sample"].astype(str).map(nbf)
def boot_spear(x,y,n=2000):
    x=np.asarray(x,float); y=np.asarray(y,float); base=spearmanr(x,y)[0]; bs=[]
    for _ in range(n):
        i=rng.integers(0,len(x),len(x))
        if len(np.unique(x[i]))>2 and len(np.unique(y[i]))>2: bs.append(spearmanr(x[i],y[i])[0])
    return base,*np.nanpercentile(bs,[2.5,97.5])
def partial(x,y,z,n=2000):  # Spearman partial x~y|z via rank residuals
    def pr(a,b,c):
        ra,rb,rc=[pd.Series(v).rank().values for v in (a,b,c)]
        return np.corrcoef(ra-np.polyval(np.polyfit(rc,ra,1),rc), rb-np.polyval(np.polyfit(rc,rb,1),rc))[0,1]
    base=pr(x,y,z); bs=[]
    for _ in range(n):
        i=rng.integers(0,len(x),len(x))
        try: bs.append(pr(np.asarray(x)[i],np.asarray(y)[i],np.asarray(z)[i]))
        except Exception: pass
    return base,*np.nanpercentile(bs,[2.5,97.5])
rowsA=[]
for yn in ["Blin_frac","myeloid_frac","total_imm_frac","nonBimm_frac"]:
    r,lo,hi=boot_spear(P.iPT_frac_epi,P[yn]); rowsA.append(dict(test=f"iPT ~ {yn}",rho=round(r,3),lo=round(lo,3),hi=round(hi,3)))
r,lo,hi=partial(P.iPT_frac_epi,P.Blin_frac,P.nonBimm_frac); rowsA.append(dict(test="iPT ~ Blin | NON-B immune (de-circ)",rho=round(r,3),lo=round(lo,3),hi=round(hi,3)))
r,lo,hi=partial(P.iPT_frac_epi,P.Blin_frac,P.total_imm_frac); rowsA.append(dict(test="iPT ~ Blin | total-immune (orig)",rho=round(r,3),lo=round(lo,3),hi=round(hi,3)))
A=pd.DataFrame(rowsA); A.to_csv(f"{OUT}/coloc_panelA_partials.csv",index=False)
print(A.to_string(index=False))

# A2 — CLR closure sensitivity (coarse 10-part composition; multiplicative zero-replacement)
PARTS={"iPT":isiPT,"PT_healthy":isPT,"TAL":isTAL,"iTAL":isiTAL,"otherEpi":OTHEREPI,
       "Endothelial":(lin=="Endothelial"),"Stroma":(lin=="Stroma"),
       "B_lineage":isB,"Myeloid":ismye,"T_NK_other":(nonBimm&~ismye)}
comp=pd.DataFrame({p:[ (m&(samp==s)).sum() for s in SAMPLES] for p,m in PARTS.items()},index=SAMPLES).astype(float)
prop=comp.div(comp.sum(1),axis=0)
# multiplicative zero replacement
for s in prop.index:
    z=prop.loc[s]==0
    if z.any():
        delta=0.5/comp.loc[s].sum(); prop.loc[s,z]=delta; prop.loc[s,~z]=prop.loc[s,~z]*(1-z.sum()*delta)
clr=np.log(prop).sub(np.log(prop).mean(1),axis=0)
nonBimm_clrpart=np.log(prop[["Myeloid","T_NK_other"]].sum(1))-np.log(prop).mean(1)  # amalgamated log-ratio
rowsA2=[]
for yn in ["B_lineage","Myeloid","T_NK_other"]:
    r,lo,hi=boot_spear(clr["iPT"],clr[yn]); rowsA2.append(dict(test=f"clr(iPT) ~ clr({yn})",rho=round(r,3),lo=round(lo,3),hi=round(hi,3)))
r,lo,hi=partial(clr["iPT"].values,clr["B_lineage"].values,nonBimm_clrpart.values)
rowsA2.append(dict(test="clr(iPT) ~ clr(B) | nonB-immune (amalg)",rho=round(r,3),lo=round(lo,3),hi=round(hi,3)))
A2=pd.DataFrame(rowsA2); A2.to_csv(f"{OUT}/coloc_panelA_clr.csv",index=False)
print("\nCLR closure sensitivity:"); print(A2.to_string(index=False))
print("DIRECTION/ORDERING check: clr myeloid >= clr B? ", float(A2.iloc[1].rho)>=float(A2.iloc[0].rho),
      "| non-B partial collapses (|rho|<0.3)? ", abs(float(A2.iloc[3].rho))<0.3)

# ============================ distances (aggregates + injury) ============================
agg_dist=np.full(len(cells),np.inf); inj_dist=np.full(len(cells),np.inf); n_agg={}; agg_hull={}
for s in SAMPLES:
    sm=samp==s; Bc=XY[sm&isB]
    if (sm&isB).sum()>=MINPTS:
        cl=DBSCAN(eps=EPS,min_samples=MINPTS).fit(Bc).labels_; mem=Bc[cl!=-1]; n_agg[s]=len(set(cl[cl!=-1]))
        if len(mem): agg_dist[np.where(sm)[0]]=cKDTree(mem).query(XY[sm])[0]; agg_hull[s]=[Bc[cl==k] for k in set(cl) if k!=-1]
    else: n_agg[s]=0
    inj=XY[sm&injured]
    if len(inj)>=10: inj_dist[np.where(sm)[0]]=cKDTree(inj).query(XY[sm])[0]
AGGSEC=[s for s in SAMPLES if n_agg[s]>0]

# helper: near/far fraction of label among ref cells, with permutation null (positions+count fixed)
def nf_frac(labelmask,refmask,distarr,s,nt=NEAR,ft=FAR,nperm=1000):
    cm=refmask&(samp==s); d=distarr[cm]; near=d<=nt; far=d>ft; lab=labelmask[cm]
    if near.sum()<10 or far.sum()<10: return None
    on=lab[near].mean(); of_=lab[far].mean(); n=len(lab); k=int(lab.sum()); pn=np.empty(nperm)
    for i in range(nperm):
        p=np.zeros(n,bool); p[rng.choice(n,k,replace=False)]=True; pn[i]=p[near].mean()
    return dict(near=round(float(on),4),far=round(float(of_),4),delta=round(float(on-of_),4),
                perm_p=round(float((pn>=on).mean()),3),n_near=int(near.sum()),n_far=int(far.sum()))

# ============================ PANEL B2 — spatial co-localization by disease (iPT + myeloid) ============================
hdr("PANEL B2 — iPT-frac AND myeloid-frac near/far aggregates, by disease")
rowsB=[]
for s in AGGSEC:
    ip=nf_frac(isiPT,epi,agg_dist,s); my=nf_frac(ismye,np.ones(len(cells),bool),agg_dist,s)
    if ip is None: continue
    amb=float(sp.loc[sp.section==s,"nonPT_amb_dz"].iloc[0]) if (sp.section==s).any() else np.nan
    rowsB.append(dict(section=s,Condition=cond[s],
        iPT_near=ip["near"],iPT_far=ip["far"],iPT_perm_p=ip["perm_p"],
        mye_near=my["near"] if my else np.nan,mye_far=my["far"] if my else np.nan,mye_perm_p=my["perm_p"] if my else np.nan,
        nonPT_amb_dz=round(amb,3),spillover=("SPILLOVER" if amb>0 else "clean")))
B=pd.DataFrame(rowsB); B.to_csv(f"{OUT}/coloc_panelB_spatial.csv",index=False); print(B.to_string(index=False))

# ============================ PANEL C1 — epithelial-subtype near/far + gate ============================
hdr("PANEL C1 — epithelial-subtype near/far aggregates + injured-state gate")
SUB={"PT":isPT,"iPT":isiPT,"TAL":isTAL,"iTAL":isiTAL,"otherEpi":OTHEREPI}
rowsC1=[]
for s in AGGSEC:
    rec={"section":s,"Condition":cond[s]}
    for nm,m in SUB.items():
        r=nf_frac(m,epi,agg_dist,s); rec[f"{nm}_d"]=r["delta"] if r else np.nan; rec[f"{nm}_p"]=r["perm_p"] if r else np.nan
    neutral=(lin=="Stroma")|(lin=="Endothelial")   # neutral non-epi bystander (NOT immune; immune IS the aggregate)
    nonepi=nf_frac(neutral,np.ones(len(cells),bool),agg_dist,s)
    rec["nonepi_amb_d"]=nonepi["delta"] if nonepi else np.nan
    # GATE: injured > healthy (iPT>PT) AND injured rise exceeds neutral-bystander (stroma/endo) pile-up
    rec["gate_iPT"]=bool((rec["iPT_d"]>rec["PT_d"]) and (rec["iPT_d"]>rec.get("nonepi_amb_d",0)) and rec["iPT_p"]<0.05)
    rowsC1.append(rec)
C1=pd.DataFrame(rowsC1); C1.to_csv(f"{OUT}/coloc_panelC1_epi_subtype.csv",index=False)
print(C1[["section","Condition","PT_d","iPT_d","TAL_d","iTAL_d","otherEpi_d","nonepi_amb_d","gate_iPT"]].to_string(index=False))
print(f"\niPT injured-state gate PASS: {int(C1.gate_iPT.sum())}/{len(C1)} sections")

# ============================ PANEL C2 — immune-subtype near/far INJURY ============================
hdr(f"PANEL C2 — immune composition near/far INJURED TUBULE (<= {INJ_NEAR} / > {INJ_FAR} um)")
INJSEC=[s for s in SAMPLES if np.isfinite(inj_dist[samp==s]).any() and (imm&(samp==s)).sum()>50]
rowsC2=[]
for s in INJSEC:
    rec={"section":s,"Condition":cond[s]}
    ok=False
    for nm,m in [("Myeloid",ismye),("B_lineage",isB),("CD4 T",ml=="CD4 T"),("CD8 T",ml=="CD8 T")]:
        r=nf_frac(m,imm,inj_dist,s,nt=INJ_NEAR,ft=INJ_FAR)
        if r is None: rec[f"{nm}_near"]=np.nan; rec[f"{nm}_d"]=np.nan; rec[f"{nm}_p"]=np.nan
        else: rec[f"{nm}_near"]=r["near"]; rec[f"{nm}_d"]=r["delta"]; rec[f"{nm}_p"]=r["perm_p"]; ok=True
    if ok: rowsC2.append(rec)
C2=pd.DataFrame(rowsC2); C2.to_csv(f"{OUT}/coloc_panelC2_immune_injury.csv",index=False)
print(C2[["section","Condition","Myeloid_near","Myeloid_d","Myeloid_p","B_lineage_d","CD4 T_d","CD8 T_d"]].to_string(index=False))
k_my=int((C2["Myeloid_p"]<0.05).sum()); print(f"\nmyeloid enriched near injured tubule (perm p<.05): {k_my}/{len(C2)} sections; "
      f"pooled myeloid Δ(near-far) = {C2['Myeloid_d'].mean():+.3f}")

# ============================ FIGURES ============================
hdr("FIGURES")
# PANEL A (partials + CLR)
fig,axes=plt.subplots(1,2,figsize=(13,4.8))
ax=axes[0]; y=np.arange(len(A)); ax.barh(y,A.rho,color=["#6A3D9A","#2ca02c","#999","#888","#2ca02c","#d6a000"][:len(A)],edgecolor="black")
ax.errorbar(A.rho,y,xerr=[A.rho-A.lo,A.hi-A.rho],fmt="none",ecolor="black",capsize=3)
ax.axvline(0,color="gray",lw=1); ax.set_yticks(y); ax.set_yticklabels(A.test,fontsize=8.5); ax.invert_yaxis()
ax.set_xlabel("Spearman rho (95% CI)"); ax.set_title("A1 — de-circularized specificity\n(raw fractions, n=16)",fontsize=10)
for sp_ in ["top","right"]: ax.spines[sp_].set_visible(False)
ax=axes[1]; y=np.arange(len(A2)); ax.barh(y,A2.rho,color=["#2ca02c","#6A3D9A","#888","#888"][:len(A2)],edgecolor="black")
ax.errorbar(A2.rho,y,xerr=[A2.rho-A2.lo,A2.hi-A2.rho],fmt="none",ecolor="black",capsize=3)
ax.axvline(0,color="gray",lw=1); ax.set_yticks(y); ax.set_yticklabels(A2.test,fontsize=8.5); ax.invert_yaxis()
ax.set_xlabel("Spearman rho (95% CI)"); ax.set_title("A2 — CLR closure sensitivity\n(direction/ordering robust)",fontsize=10)
for sp_ in ["top","right"]: ax.spines[sp_].set_visible(False)
fig.suptitle("PANEL A — iPT~myeloid leads; B-lineage not specific beyond non-B immune (descriptive, n=16)",fontsize=12)
fig.tight_layout(); save(fig,"coloc_fig_A_specificity")
# PANEL B1
fig,axes=plt.subplots(1,2,figsize=(13,5.2))
for ax,(xn,xl) in zip(axes,[("total_imm_frac","total-immune fraction"),("myeloid_frac","myeloid fraction")]):
    for _,r in P.iterrows():
        c=DISCOL.get(r.Condition,"#444"); mk="s" if r.Condition=="Control" else "o"
        ax.scatter(r[xn],r.iPT_frac_epi,s=85,c=c,marker=mk,edgecolor="black",zorder=3)
        ax.annotate(r["sample"],(r[xn],r.iPT_frac_epi),fontsize=6.5,xytext=(4,2),textcoords="offset points")
    rr=spearmanr(P[xn],P.iPT_frac_epi)[0]; ax.set_xlabel(xl); ax.set_ylabel("iPT fraction of epithelium")
    ax.set_title(f"iPT ~ {xl}  (rho={rr:.2f})",fontsize=11)
    for sp_ in ["top","right"]: ax.spines[sp_].set_visible(False)
axes[1].legend(handles=[plt.Line2D([],[],marker=("s" if k=="Control" else "o"),ls="",mfc=v,mec="black",ms=8,label=k) for k,v in DISCOL.items()],fontsize=8,loc="upper left")
fig.suptitle("PANEL B1 — tubular injury vs immune infiltration by disease (controls = squares, low-low corner)",fontsize=12)
fig.tight_layout(); save(fig,"coloc_fig_B1_disease")
# PANEL B2
fig,ax=plt.subplots(figsize=(12,5)); x=np.arange(len(B))
ax.bar(x-0.2,B.iPT_near-B.iPT_far,0.4,label="iPT Δ(near−far)",color="#6A3D9A",edgecolor="black")
ax.bar(x+0.2,B.mye_near-B.mye_far,0.4,label="myeloid Δ(near−far)",color="#2ca02c",edgecolor="black")
ax.axhline(0,color="gray",lw=1)
for xi,r in zip(x,B.itertuples()):
    if r.iPT_perm_p<0.05: ax.text(xi-0.2,(r.iPT_near-r.iPT_far)+0.01,"*",ha="center",color="#6A3D9A",fontsize=12)
    if r.mye_perm_p<0.05: ax.text(xi+0.2,(r.mye_near-r.mye_far)+0.01,"*",ha="center",color="#2ca02c",fontsize=12)
    if r.spillover=="SPILLOVER": ax.text(xi,-0.02,"⚠",ha="center",color="#d6a000",fontsize=11)
ax.set_xticks(x); ax.set_xticklabels([f"{r.section}\n{r.Condition}" for r in B.itertuples()],fontsize=8)
ax.set_ylabel("near − far fraction"); ax.legend(fontsize=9)
ax.set_title("PANEL B2 — co-localization near B-aggregates by disease (* perm p<.05; ⚠ spillover section)",fontsize=11)
for sp_ in ["top","right"]: ax.spines[sp_].set_visible(False)
fig.tight_layout(); save(fig,"coloc_fig_B2_spatial")
# PANEL B3 gallery — one rep section per disease
reps={"DKD":"1006","AA amyloid":"1009","C3GN":"1007","MN":"1005","IgA":"1003","Control":"HK3626"}
LC={"iPT":"#6A3D9A","iTAL":"#9467bd","Myeloid":"#2ca02c","B":"#1f77b4","Plasma":"#ff7f0e"}
fig,axes=plt.subplots(2,3,figsize=(16,10))
for ax,(dis,s) in zip(axes.ravel(),reps.items()):
    sm=samp==s; ax.scatter(XY[sm,0],XY[sm,1],s=1,c="#ececec",linewidths=0,rasterized=True)
    for nm,c in LC.items():
        m=sm&(ml==nm) if nm in("iPT","iTAL","Myeloid","B","Plasma") else sm&(ml==nm)
        ax.scatter(XY[m,0],XY[m,1],s=(2 if nm in("iPT","iTAL") else 6),c=c,linewidths=0,rasterized=True)
    for P_ in agg_hull.get(s,[]):
        if len(P_)>=4:
            try: h=ConvexHull(P_); ax.add_patch(Polygon(P_[h.vertices],closed=True,fill=False,edgecolor="#08306b",lw=1.3))
            except Exception: pass
    ax.set_aspect("equal"); ax.axis("off"); ax.set_title(f"{dis} — {s}",color=DISCOL.get(dis,"#333"),fontsize=11)
fig.legend(handles=[plt.Line2D([],[],marker="o",ls="",mfc=c,mec="none",ms=8,label=k) for k,c in LC.items()]+
           [plt.Line2D([],[],marker="s",ls="",mfc="none",mec="#08306b",ms=10,label="B-aggregate")],
           loc="lower center",ncol=6,frameon=False,fontsize=10,bbox_to_anchor=(0.5,-0.01))
fig.suptitle("PANEL B3 — injured tubule (iPT/iTAL) + immune (myeloid/B/plasma) with B-aggregate hulls, one section per disease",fontsize=13)
fig.subplots_adjust(bottom=0.06,top=0.95,wspace=0.04,hspace=0.12); save(fig,"coloc_fig_B3_gallery")
# PANEL C1 epithelial subtype near/far
fig,ax=plt.subplots(figsize=(12,5)); x=np.arange(len(C1)); w=0.2
for off,col,nm in [(-1.5*w,"#c7c7c7","PT (healthy)"),(-0.5*w,"#6A3D9A","iPT"),(0.5*w,"#9467bd","iTAL"),(1.5*w,"#8c564b","stroma/endo ambient")]:
    key={"PT (healthy)":"PT_d","iPT":"iPT_d","iTAL":"iTAL_d","stroma/endo ambient":"nonepi_amb_d"}[nm]
    ax.bar(x+off,C1[key],w,label=nm,color=col,edgecolor="black",linewidth=0.3)
ax.axhline(0,color="gray",lw=1)
for xi,r in zip(x,C1.itertuples()):
    if getattr(r,"gate_iPT"): ax.text(xi,max(r.iPT_d,r.iTAL_d,0)+0.03,"✓",ha="center",color="#2ca02c",fontsize=12)
ax.set_xticks(x); ax.set_xticklabels([f"{r.section}\n{r.Condition}" for r in C1.itertuples()],fontsize=8)
ax.set_ylabel("near − far fraction (among epithelium)"); ax.legend(fontsize=8.5,ncol=4)
ax.set_title("PANEL C1 — near B-aggregates the epithelium shifts to INJURED states (healthy PT down, iPT/iTAL up; ✓=injured gate)",fontsize=10.5)
for sp_ in ["top","right"]: ax.spines[sp_].set_visible(False)
fig.tight_layout(); save(fig,"coloc_fig_C1_epi_subtype")
# PANEL C2 myeloid near-injury
fig,ax=plt.subplots(figsize=(12,5)); x=np.arange(len(C2))
ax.bar(x,C2["Myeloid_d"],color=[DISCOL.get(c,"#444") for c in C2.Condition],edgecolor="black")
ax.axhline(0,color="gray",lw=1)
for xi,r in zip(x,C2.itertuples()):
    if getattr(r,"Myeloid_p")<0.05: ax.text(xi,getattr(r,"Myeloid_d")+0.005,"*",ha="center",fontsize=12,color="#d62728")
ax.set_xticks(x); ax.set_xticklabels([f"{r.section}\n{r.Condition}" for r in C2.itertuples()],fontsize=7.5)
ax.set_ylabel("myeloid Δ fraction (near−far injured tubule)")
ax.set_title(f"PANEL C2 — near-injury infiltrate is myeloid-skewed (* perm p<.05); near<={INJ_NEAR}/far>{INJ_FAR}um to iPT/iTAL",fontsize=11)
for sp_ in ["top","right"]: ax.spines[sp_].set_visible(False)
fig.tight_layout(); save(fig,"coloc_fig_C2_immune_injury")
print("\n== colocalization done ==")
