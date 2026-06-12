#!/usr/bin/env python
"""Finish the BAFF axis with nCount rigor. PART A: receptors (BAFF-R/TNFRSF13C, BCMA/TNFRSF17,
TACI/TNFRSF13B) in their target cell type INSIDE vs OUTSIDE B-aggregates, controlled for nCount
inflation AND a neutral control-gene panel. PART B: anchor the myeloid BAFF (TNFSF13B) producer
(nCount control, cell-intrinsic identity, per-sample reproducibility, nCount-matched fold). APRIL is
absent so BAFF is the relevant ligand for BCMA/TACI. Read-only raw; reuse dkd_xenium_disease object."""
import os, numpy as np, pandas as pd
from sklearn.cluster import DBSCAN
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
RE=f"{REPO}/analysis/dkd_xenium_reannotation"; OUT=f"{REPO}/analysis/dkd_xenium_disease"
EPS,MINPTS=50,20
CTRL=["TPT1","PPIA","YWHAZ","TMSB10","UBB"]; MK=["CD68","CD14","CD163","ITGAX","AIF1","C1QA"]
AGG=["1006","HK2695","1007","1009"]

cells=pd.read_parquet(f"{RE}/cells.parquet").reset_index(drop=True)
lig=pd.read_parquet(f"{OUT}/ligand_me_percell.parquet").reset_index(drop=True)
nc=pd.read_parquet(f"{OUT}/ncount_percell.parquet").reset_index(drop=True)
ex=pd.read_parquet(f"{OUT}/genes_extra_percell.parquet").reset_index(drop=True)
for d in (lig,nc,ex): assert (cells.cell_id.astype(str).values==d.cell_id.astype(str).values).all()
for g in ["TNFSF13B","TNFRSF17","TNFRSF13B","TNFRSF13C"]: cells[g]=lig[g].values
for g in CTRL+MK: cells[g]=ex[g].values
cells["nCount"]=nc.nCount.values
cells["is_Blin"]=cells.my_label.isin(["B","Plasma"])

# aggregate membership (DBSCAN on B-lineage per section), only in aggregate-bearing sections
cells["agg_member"]=False
for sid in AGG:
    idx=cells.index[(cells.orig_ident==sid)&cells.is_Blin]
    xy=cells.loc[idx,["spatial_x","spatial_y"]].values
    if len(xy)>=MINPTS:
        lab=DBSCAN(eps=EPS,min_samples=MINPTS).fit(xy).labels_
        cells.loc[idx[lab!=-1],"agg_member"]=True

def rate(df,g): return (df[g]>0).mean()
def norm_expr(df,g): return (df[g]/df.nCount.clip(lower=1)).mean()

# ================= PART A — receptors, within-target, nCount + control-gene controlled =================
RECs=[("BAFF-R","TNFRSF13C",["B"]),("TACI","TNFRSF13B",["B","Plasma"]),("BCMA","TNFRSF17",["Plasma"])]
def targmask(df,types): return df.my_label.isin(types)
rowsA=[]
AGc=cells[cells.orig_ident.isin(AGG)]
for setname,scope in [("POOLED",AGc)]+[(s,AGc[AGc.orig_ident==s]) for s in AGG]:
    for rname,gene,types in RECs:
        tt=scope[targmask(scope,types)]
        ins=tt[tt.agg_member]; out=tt[~tt.agg_member]
        if len(ins)<10 or len(out)<10:
            rowsA.append(dict(sample_set=setname,receptor=rname,n_in=len(ins),n_out=len(out),note="too few")); continue
        det_in,det_out=rate(ins,gene),rate(out,gene); det_ratio=det_in/max(det_out,1e-9)
        ncr=ins.nCount.mean()/max(out.nCount.mean(),1e-9)
        ne_in,ne_out=norm_expr(ins,gene),norm_expr(out,gene); ne_ratio=ne_in/max(ne_out,1e-9)
        ctrl_ratios=[rate(ins,c)/max(rate(out,c),1e-9) for c in CTRL]
        ctrl_med=float(np.median(ctrl_ratios))
        verdict=("aggregate-upregulated" if (det_ratio>1.2 and ne_ratio>1.15 and det_ratio>1.2*ctrl_med)
                 else ("nCount/ambient-inflation" if det_ratio>1.2 else "null"))
        rowsA.append(dict(sample_set=setname,receptor=rname,n_in=len(ins),n_out=len(out),
            det_in_pct=round(det_in*100,2),det_out_pct=round(det_out*100,2),det_ratio=round(det_ratio,2),
            nCount_ratio=round(ncr,2),normExpr_ratio=round(ne_ratio,2),ctrl_gene_ratio=round(ctrl_med,2),
            verdict=verdict))
A=pd.DataFrame(rowsA); A.to_csv(f"{OUT}/receptor_aggregate_control.csv",index=False)
pd.set_option("display.width",260)
print("=== PART A — receptors INSIDE vs OUTSIDE aggregates (nCount + control-gene controlled) ===")
print(A[A["sample_set"]=="POOLED"].to_string(index=False))
print("\n--- per-sample ---")
print(A[A["sample_set"]!="POOLED"][["sample_set","receptor","det_in_pct","det_out_pct","det_ratio","normExpr_ratio","ctrl_gene_ratio","verdict"]].to_string(index=False))

# ================= PART B — anchor the myeloid BAFF producer =================
BAFF="TNFSF13B"; mye=cells[cells.my_label=="Myeloid"]; epi=cells[cells.my_lineage=="Epithelial"]
# 1. nCount control: BAFF normalized expr myeloid vs epithelium; BAFF vs control genes within myeloid
b1=dict(
  myeloid_BAFF_det=round(rate(mye,BAFF)*100,2), epi_BAFF_det=round(rate(epi,BAFF)*100,2),
  myeloid_nCount_med=int(np.median(mye.nCount)), epi_nCount_med=int(np.median(epi.nCount)),
  myeloid_BAFF_normExpr=float(f"{norm_expr(mye,BAFF):.2e}"), epi_BAFF_normExpr=float(f"{norm_expr(epi,BAFF):.2e}"),
  normExpr_fold=round(norm_expr(mye,BAFF)/max(norm_expr(epi,BAFF),1e-12),1))
# BAFF normExpr vs control-gene normExpr within myeloid (BAFF is specific, not generic high-count)
b1["BAFF_vs_ctrl_in_myeloid_detfold"]=round(rate(mye,BAFF)/max(np.median([rate(mye,c) for c in CTRL]),1e-9),3)
pd.DataFrame([b1]).to_csv(f"{OUT}/myeloid_anchor_ncount.csv",index=False)
print("\n=== PART B.1 — myeloid BAFF nCount control ===")
for k,v in b1.items(): print(f"  {k}: {v}")

# 4. nCount-matched fold (restrict myeloid & epithelium to overlapping nCount bins, recompute detection fold)
lo,hi=np.percentile(np.r_[mye.nCount.values,epi.nCount.values],[5,95])
mm=mye[(mye.nCount>=lo)&(mye.nCount<=hi)]; ee=epi[(epi.nCount>=lo)&(epi.nCount<=hi)]
# bin-match: within deciles of nCount, average per-bin detection
bins=np.quantile(np.r_[mm.nCount,ee.nCount],np.linspace(0,1,11))
def binned_det(df):
    d=[];
    for i in range(10):
        m=df[(df.nCount>=bins[i])&(df.nCount<bins[i+1])]
        if len(m)>30: d.append(rate(m,BAFF))
    return np.mean(d) if d else np.nan
matchfold=binned_det(mm)/max(binned_det(ee),1e-9)
print(f"\n=== PART B.4 — nCount-matched myeloid/epithelial BAFF detection fold = {matchfold:.1f}x (raw 7.9x) ===")

# 2. identity: BAFF+ vs BAFF- myeloid marker profile
mp=mye[mye[BAFF]>0]; mn=mye[mye[BAFF]==0]
idrows=[]
for g in MK:
    idrows.append(dict(marker=g,BAFFpos_det=round(rate(mp,g)*100,1),BAFFneg_det=round(rate(mn,g)*100,1)))
idt=pd.DataFrame(idrows); idt["BAFFpos_n"]=len(mp); idt["BAFFneg_n"]=len(mn)
idt.to_csv(f"{OUT}/myeloid_baff_identity.csv",index=False)
print("\n=== PART B.2 — BAFF+ vs BAFF- myeloid identity (marker detection %) ===")
print(idt.to_string(index=False))

# 3. reproducibility: per-sample myeloid BAFF detection across all 16
rep=cells[cells.my_label=="Myeloid"].groupby("orig_ident",observed=True).apply(
    lambda g: pd.Series(dict(Condition=g.Condition.iloc[0],n_myeloid=len(g),BAFF_det_pct=round(rate(g,BAFF)*100,2)))).reset_index()
rep=rep.sort_values("BAFF_det_pct",ascending=False)
rep.to_csv(f"{OUT}/myeloid_baff_per_sample.csv",index=False)
print("\n=== PART B.3 — per-sample myeloid BAFF detection (all 16) ===")
print(rep.to_string(index=False))
print(f"\nmyeloid BAFF detection: median {rep.BAFF_det_pct.median():.2f}%, IQR {rep.BAFF_det_pct.quantile(.25):.2f}-{rep.BAFF_det_pct.quantile(.75):.2f}%, range {rep.BAFF_det_pct.min():.2f}-{rep.BAFF_det_pct.max():.2f}% across {len(rep)} samples")

# ================= figures =================
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys; sys.path.insert(0,f"{REPO}/presentation/figures"); import figstyle as fs
FIG=f"{OUT}/figures"; os.makedirs(FIG,exist_ok=True)
def save(fig,name):
    fig.savefig(f"{FIG}/{name}.png",dpi=300,bbox_inches="tight"); fig.savefig(f"{FIG}/{name}.svg",bbox_inches="tight")
    plt.close(fig); print(f"  [ok] {name}")
VCOL={"aggregate-upregulated":"#2ca02c","nCount/ambient-inflation":"#d6a000","null":"#999999"}

# --- FIG 1: receptors, in/out ratio vs control-gene & nCount inflation ---
P=A[A["sample_set"]=="POOLED"].set_index("receptor")
fig,axes=plt.subplots(1,2,figsize=(13.5,5.4))
ax=axes[0]; recs=["BAFF-R","TACI","BCMA"]; x=np.arange(len(recs))
ax.bar(x,[P.loc[r,"det_ratio"] for r in recs],0.5,color=[VCOL.get(P.loc[r,"verdict"],"#999") for r in recs],edgecolor="#333",label="receptor in/out ratio")
ax.plot(x,[P.loc[r,"ctrl_gene_ratio"] for r in recs],"D",ms=11,color="black",label="neutral control-gene ratio (generic inflation)")
ax.plot(x,[P.loc[r,"nCount_ratio"] for r in recs],"s",ms=9,mfc="none",mec="#555",label="nCount in/out ratio")
ax.axhline(1.0,color="#888",lw=1)
for xi,r in zip(x,recs): ax.text(xi,P.loc[r,"det_ratio"]+0.04,P.loc[r,"verdict"].split("/")[0],ha="center",fontsize=8,fontweight="bold")
ax.set_xticks(x); ax.set_xticklabels([f"{r}\n({'B' if r=='BAFF-R' else ('B+Pl' if r=='TACI' else 'Plasma')})" for r in recs])
ax.set_ylabel("in-aggregate / outside ratio"); ax.legend(frameon=False,fontsize=8.5,loc="upper right")
ax.set_title("(a) receptor in/out ratio sits AT the neutral control-gene inflation\n-> generic transcript inflation, not specific aggregate upregulation",fontsize=10)
for s in ["top","right"]: ax.spines[s].set_visible(False)
ax=axes[1]; ps=A[(A["sample_set"]!="POOLED")&(A.receptor=="BAFF-R")]
xx=np.arange(len(ps)); w=0.36
ax.bar(xx-w/2,ps.det_ratio,w,color="#1f77b4",label="BAFF-R in/out ratio")
ax.bar(xx+w/2,ps.ctrl_gene_ratio,w,color="#cccccc",label="control-gene ratio")
ax.axhline(1.0,color="#888",lw=1)
ax.set_xticks(xx); ax.set_xticklabels(ps.sample_set,fontsize=9); ax.set_ylabel("in/out ratio")
ax.legend(frameon=False,fontsize=8.5)
ax.set_title("(b) BAFF-R per section: exceeds control only in 1006/1007,\nat/below control in HK2695/1009 -> inconsistent, not reproducible",fontsize=10)
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.suptitle("PART A -- BAFF receptors are NOT aggregate-upregulated beyond nCount/control-gene inflation",fontsize=13,y=1.01)
fig.subplots_adjust(top=0.84,wspace=0.22); save(fig,"baff_receptor_control")

# --- FIG 2: myeloid BAFF anchor ---
fig,axes=plt.subplots(1,3,figsize=(17,5.3))
ax=axes[0]
ax.bar([0,1],[b1["myeloid_BAFF_normExpr"],b1["epi_BAFF_normExpr"]],color=["#9467bd","#c7c7c7"],edgecolor="#333")
ax.set_xticks([0,1]); ax.set_xticklabels([f"Myeloid\n(nCount med {b1['myeloid_nCount_med']})",f"Epithelial\n(nCount med {b1['epi_nCount_med']})"])
ax.set_ylabel("BAFF per-transcriptome (count/nCount)")
ax.set_title(f"(1) cell-intrinsic, NOT a count artifact\nmyeloid {b1['normExpr_fold']}x epithelium per-transcriptome\n(myeloid have FEWER counts; nCount-matched {matchfold:.0f}x)",fontsize=9.5)
for s in ["top","right"]: ax.spines[s].set_visible(False)
ax=axes[1]; xx=np.arange(len(MK)); w=0.4
ax.bar(xx-w/2,idt.BAFFpos_det,w,label=f"BAFF+ myeloid (n={len(mp)})",color="#6a3d9a")
ax.bar(xx+w/2,idt.BAFFneg_det,w,label=f"BAFF- myeloid (n={len(mn)})",color="#c9b8de")
ax.set_xticks(xx); ax.set_xticklabels(MK,rotation=40,ha="right",fontsize=9); ax.set_ylabel("marker detection %")
ax.legend(frameon=False,fontsize=8.5)
ax.set_title("(2) BAFF+ myeloid are bona fide myeloid\n(richer CD68/CD14/C1QA/AIF1 than BAFF- -> activated macro)",fontsize=9.5)
for s in ["top","right"]: ax.spines[s].set_visible(False)
ax=axes[2]; rp=rep.sort_values("BAFF_det_pct",ascending=False)
cc={"DKD":fs.DATASET["DKD"],"Control":"#888"}
cols=[cc.get(c,"#b8860b") for c in rp.Condition]
ax.bar(np.arange(len(rp)),rp.BAFF_det_pct,color=cols,edgecolor="#333",linewidth=0.3)
ax.axhline(rp.BAFF_det_pct.median(),ls="--",color="#444",lw=1); ax.text(len(rp)-0.5,rp.BAFF_det_pct.median()+0.1,f"median {rp.BAFF_det_pct.median():.1f}%",ha="right",fontsize=8.5)
ax.set_xticks(np.arange(len(rp))); ax.set_xticklabels(rp.orig_ident,rotation=90,fontsize=7); ax.set_ylabel("myeloid BAFF detection %")
ax.set_title(f"(3) reproducible across ALL 16 samples\n(range {rp.BAFF_det_pct.min():.1f}-{rp.BAFF_det_pct.max():.1f}%, not 1-2 sections)",fontsize=9.5)
ax.legend(handles=[Line2D([],[],marker="s",ls="",mfc=cc["DKD"],mec="none",ms=8,label="DKD"),Line2D([],[],marker="s",ls="",mfc="#888",mec="none",ms=8,label="Control"),Line2D([],[],marker="s",ls="",mfc="#b8860b",mec="none",ms=8,label="other")],frameon=False,fontsize=8,loc="upper right")
for s in ["top","right"]: ax.spines[s].set_visible(False)
fig.suptitle("PART B -- myeloid BAFF is a RELIABLE anchored producer (cell-intrinsic, bona fide myeloid, reproducible)",fontsize=13,y=1.01)
fig.subplots_adjust(top=0.82,wspace=0.28); save(fig,"baff_myeloid_anchor")
print("== baff_receptor_anchor done ==")
