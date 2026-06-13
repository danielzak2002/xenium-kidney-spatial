#!/usr/bin/env python
"""Re-examine the injured-PT <-> B-lineage association under post-BAFF rigor.
RECONCILE-AND-EXTEND of dkd_epi_endo_stress.py (prior commit, author-labels + B-only aggregates:
injPT near>far Δz +0.13, 6/9, p=0.038). HERE: VALIDATED REANNOTATION labels (incl. iPT) +
B-lineage (B+Plasma) DBSCAN aggregates (eps=50/minPts=20). injPT program genes inherited verbatim
(VCAM1/HAVCR1/PROM1/SPP1/ITGB6). Unit of replication = SAMPLE; descriptive throughout. Read-only.

PRIMARY  = cross-sample Spearman (iPT burden ~ B-lineage burden) + specificity vs total-immune.
SECONDARY= spatial near/far iPT across aggregate-bearing sections, gated exactly like the BAFF
           controls (control-gene floor, nCount, non-PT ambient, within-section permutation null)."""
import os, warnings
warnings.filterwarnings("ignore"); os.environ.setdefault("KMP_DUPLICATE_LIB_OK","TRUE")
import numpy as np, pandas as pd
from sklearn.cluster import DBSCAN
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import sys; sys.path.insert(0,"/Users/danie/ClaudeCode/pilot_analyses/xenium/presentation/figures")
try: import figstyle as fs; DKDC=fs.DATASET["DKD"]
except Exception: DKDC="#6A3D9A"
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
RE=f"{REPO}/analysis/dkd_xenium_reannotation"; DIS=f"{REPO}/analysis/dkd_xenium_disease"
OUT=f"{REPO}/analysis/dkd_epi_endo_stress"
EPS,MINPTS,NEAR,FAR=50.0,20,50.0,200.0
INJPT=["VCAM1","HAVCR1","PROM1","SPP1","ITGB6"]
CTRL=["TPT1","PPIA","YWHAZ","TMSB10","UBB"]
rng=np.random.default_rng(0)
def hdr(s): print("\n"+"="*78+"\n"+s+"\n"+"="*78)

# ---------------- load (validated reannotation labels + coords) ----------------
cells=pd.read_parquet(f"{RE}/cells.parquet").reset_index(drop=True)
inj=pd.read_parquet(f"{DIS}/injpt_genes_percell.parquet").reset_index(drop=True)
nc=pd.read_parquet(f"{DIS}/ncount_percell.parquet").reset_index(drop=True)
ex=pd.read_parquet(f"{DIS}/genes_extra_percell.parquet").reset_index(drop=True)
for d in (inj,nc,ex): assert (cells.cell_id.astype(str).values==d.cell_id.astype(str).values).all()
for g in INJPT: cells[g]=inj[g].values
for g in CTRL: cells[g]=ex[g].values
cells["nCount"]=nc.nCount.values
PTc=cells.my_label.isin(["PT","iPT"]).values        # PT compartment
epi=(cells.my_lineage=="Epithelial").values
nonPTepi=epi&~PTc
imm=(cells.my_lineage=="Immune").values
isB=cells.my_label.isin(["B","Plasma"]).values       # B-lineage
isiPT=(cells.my_label=="iPT").values
ismye=(cells.my_label=="Myeloid").values
samp=cells.orig_ident.astype(str).values
cond=cells.groupby("orig_ident",observed=True).Condition.first()
SAMPLES=sorted(pd.unique(samp))

# ---------------- usability gate (PT vs Immune ambient floor; inherit prior gate) ----------------
hdr("USABILITY GATE — injPT genes, detection in PT compartment vs Immune floor (this cohort)")
def det(mask,g): v=cells[g].values[mask]; return float((v>0).mean())
gate=[]; usable=[]
for g in INJPT:
    td=det(PTc,g); fd=det(imm,g); ok=(td>=0.03) and (td>=2*fd)
    gate.append(dict(gene=g,PT_detect=round(td,4),immune_floor=round(fd,4),ratio=round(td/max(fd,1e-9),2),usable=ok))
    if ok: usable.append(g)
gate_df=pd.DataFrame(gate); gate_df.to_csv(f"{OUT}/reconcile_usability_gate.csv",index=False)
print(gate_df.to_string(index=False)); print("usable:",usable)

# ---------------- CP-median log1p; injPT + control module raw scores (per cell) ----------------
med=np.median(cells.nCount.values[cells.nCount.values>0])
def lognorm(g): return np.log1p(cells[g].values/np.where(cells.nCount.values>0,cells.nCount.values,1)*med)
inj_raw=np.mean([lognorm(g) for g in usable],axis=0)
ctrl_raw=np.mean([lognorm(g) for g in CTRL],axis=0)

# ---------------- B-lineage DBSCAN aggregates per section; distance for epi cells ----------------
dist=np.full(len(cells),np.inf); n_agg={}; agg_members={}
for s in SAMPLES:
    sm=samp==s; Bm=sm&isB; Bc=cells.loc[Bm,["spatial_x","spatial_y"]].values
    if Bm.sum()<MINPTS: n_agg[s]=0; continue
    cl=DBSCAN(eps=EPS,min_samples=MINPTS).fit(Bc).labels_; mem=Bc[cl!=-1]
    n_agg[s]=len(set(cl[cl!=-1])); agg_members[s]=mem
    if len(mem)==0: continue
    cellmask=sm&(epi)   # distance computed for epithelial cells (PT + nonPT)
    d,_=cKDTree(mem).query(cells.loc[cellmask,["spatial_x","spatial_y"]].values)
    dist[np.where(cellmask)[0]]=d
near=dist<=NEAR; far=dist>FAR
def grid_area(s,b=50.0):
    xy=cells.loc[samp==s,["spatial_x","spatial_y"]].values
    return len(set(map(tuple,np.floor(xy/b).astype(int))))*(b/1000)**2

# ============================ PRIMARY — cross-sample correlation ============================
hdr("PRIMARY — per-sample burdens + cross-sample Spearman")
rows=[]
for s in SAMPLES:
    sm=samp==s; n=sm.sum()
    e=sm&epi; ptm=sm&PTc
    rows.append(dict(sample=s,Condition=cond[s],n_cells=int(n),
        iPT_frac_epi=round((isiPT[e].sum())/max(e.sum(),1),4),
        prog_score_PT=round(float(inj_raw[ptm].mean()),4),
        Blin_frac=round(isB[sm].sum()/n,4), n_agg=n_agg[s], agg_per_mm2=round(n_agg[s]/grid_area(s),3),
        total_imm_frac=round(imm[sm].sum()/n,4), myeloid_frac=round(ismye[sm].sum()/n,4)))
P=pd.DataFrame(rows); P.to_csv(f"{OUT}/reconcile_per_sample.csv",index=False)
print(P.to_string(index=False))

def boot_spear(x,y,n=2000):
    x=np.asarray(x,float); y=np.asarray(y,float); base=spearmanr(x,y)[0]
    bs=[]
    for _ in range(n):
        i=rng.integers(0,len(x),len(x))
        if len(np.unique(x[i]))>2 and len(np.unique(y[i]))>2: bs.append(spearmanr(x[i],y[i])[0])
    lo,hi=np.nanpercentile(bs,[2.5,97.5]); return base,lo,hi
def partial_spear(x,y,z,n=2000):  # Spearman partial corr x~y | z (rank residuals)
    def pr(a,b,c):
        ra,rb,rc=[pd.Series(v).rank().values for v in (a,b,c)]
        rx=ra-np.polyval(np.polyfit(rc,ra,1),rc); ry=rb-np.polyval(np.polyfit(rc,rb,1),rc)
        return np.corrcoef(rx,ry)[0,1]
    base=pr(x,y,z); bs=[]
    for _ in range(n):
        i=rng.integers(0,len(x),len(x))
        try: bs.append(pr(np.asarray(x)[i],np.asarray(y)[i],np.asarray(z)[i]))
        except Exception: pass
    lo,hi=np.nanpercentile(bs,[2.5,97.5]); return base,lo,hi

cors=[]
for xn,x in [("iPT_frac_epi",P.iPT_frac_epi),("prog_score_PT",P.prog_score_PT)]:
    for yn,y in [("Blin_frac",P.Blin_frac),("n_agg",P.n_agg),("total_imm_frac",P.total_imm_frac),("myeloid_frac",P.myeloid_frac)]:
        r,lo,hi=boot_spear(x,y); cors.append(dict(x=xn,y=yn,rho=round(r,3),ci_lo=round(lo,3),ci_hi=round(hi,3)))
# partials controlling for total immune
for xn,x in [("iPT_frac_epi",P.iPT_frac_epi),("prog_score_PT",P.prog_score_PT)]:
    r,lo,hi=partial_spear(x,P.Blin_frac,P.total_imm_frac)
    cors.append(dict(x=xn,y="Blin_frac | total_imm",rho=round(r,3),ci_lo=round(lo,3),ci_hi=round(hi,3)))
C=pd.DataFrame(cors); C.to_csv(f"{OUT}/reconcile_correlations.csv",index=False)
print("\n=== cross-sample Spearman (n=%d samples; bootstrap 95%% CI) ===" % len(P)); print(C.to_string(index=False))

# ============================ SECONDARY — spatial near/far (BAFF-style controls) ============================
hdr("SECONDARY — near/far injPT program across aggregate-bearing sections (gated)")
AGGSEC=[s for s in SAMPLES if n_agg[s]>0]
def nf_z(raw,compmask,s):    # within-(section,compartment) z, then near-mean - far-mean
    cm=compmask&(samp==s); v=raw[cm]; sd=v.std()
    nn=near[cm]; ff=far[cm]
    if sd==0 or nn.sum()<10 or ff.sum()<10: return None
    z=(v-v.mean())/sd; return float(z[nn].mean()-z[ff].mean()), int(nn.sum()), int(ff.sum())
def perm_frac(s,nperm=2000):  # iPT fraction near/far among epithelial cells; positions+count fixed
    cm=epi&(samp==s); lab=isiPT[cm]; nn=near[cm]; ff=far[cm]
    if nn.sum()<10 or ff.sum()<10: return None
    obs_near=lab[nn].mean(); obs_far=lab[ff].mean(); n=len(lab); k=int(lab.sum())
    pn=np.empty(nperm)
    for i in range(nperm):
        p=np.zeros(n,bool); p[rng.choice(n,k,replace=False)]=True; pn[i]=p[nn].mean()
    pval=float((pn>=obs_near).mean())
    return round(obs_near,4),round(obs_far,4),round(float(obs_near-obs_far),4),round(pval,3)
rows=[]
for s in AGGSEC:
    prog=nf_z(inj_raw,PTc,s)              # PT injPT program
    ctl=nf_z(ctrl_raw,PTc,s)              # control-gene floor (PT)
    amb=nf_z(inj_raw,nonPTepi,s)          # non-PT epithelial ambient
    ncz=nf_z(cells.nCount.values.astype(float),PTc,s)  # density
    pf=perm_frac(s)
    if prog is None: continue
    rows.append(dict(section=s,Condition=cond[s],n_near=prog[1],n_far=prog[2],
        prog_dz=round(prog[0],3), ctrl_floor_dz=round(ctl[0],3) if ctl else np.nan,
        nonPT_amb_dz=round(amb[0],3) if amb else np.nan, nCount_dz=round(ncz[0],3) if ncz else np.nan,
        prog_minus_ctrl=round(prog[0]-(ctl[0] if ctl else 0),3),
        iPTfrac_near=pf[0] if pf else np.nan,iPTfrac_far=pf[1] if pf else np.nan,
        iPTfrac_perm_p=pf[3] if pf else np.nan))
S=pd.DataFrame(rows); S.to_csv(f"{OUT}/reconcile_spatial_per_section.csv",index=False)
print(S.to_string(index=False))
k_prog=int((S.prog_dz>0).sum()); n_prog=len(S)
k_exceed=int((S.prog_minus_ctrl>0).sum())
k_perm=int((S.iPTfrac_perm_p<0.05).sum())
amb_pos=int((S.nonPT_amb_dz>0).sum())
print(f"\nprog near>far: {k_prog}/{n_prog} | exceeds control floor: {k_exceed}/{n_prog} | "
      f"iPT-frac perm p<.05: {k_perm}/{n_prog} | non-PT ambient also up: {amb_pos}/{n_prog}")
print(f"pooled prog Δz = {S.prog_dz.mean():+.3f} (median {S.prog_dz.median():+.3f}); "
      f"pooled control floor = {S.ctrl_floor_dz.mean():+.3f}; pooled non-PT ambient = {S.nonPT_amb_dz.mean():+.3f}")

# ============================ figures ============================
hdr("FIGURES")
def save(fig,name): fig.savefig(f"{OUT}/{name}.png",dpi=200,bbox_inches="tight"); plt.close(fig); print("  [ok]",name)
CC={"DKD":DKDC,"Control":"#888","AA amyloid":"#FF7F00","C3GN":"#33A02C","IgA":"#1f77b4","MN":"#E31A1C"}
# PRIMARY
fig,axes=plt.subplots(1,3,figsize=(17,5.2))
for ax,(yn,ylab) in zip(axes[:2],[("Blin_frac","B-lineage fraction"),("total_imm_frac","total-immune fraction")]):
    for _,r in P.iterrows(): ax.scatter(r[yn],r.iPT_frac_epi,s=70,c=CC.get(r.Condition,"#444"),edgecolor="black",zorder=3)
    rr=C[(C.x=="iPT_frac_epi")&(C.y==yn)].iloc[0]
    ax.set_xlabel(ylab); ax.set_ylabel("iPT fraction of epithelium")
    ax.set_title(f"iPT ~ {ylab}\nSpearman rho={rr.rho} [{rr.ci_lo},{rr.ci_hi}]",fontsize=11)
    for s_ in ["top","right"]: ax.spines[s_].set_visible(False)
ax=axes[2]; lab=["iPT~Blin","iPT~Bagg","iPT~totImm","iPT~myel","iPT~Blin|totImm"]
keys=[("iPT_frac_epi","Blin_frac"),("iPT_frac_epi","n_agg"),("iPT_frac_epi","total_imm_frac"),("iPT_frac_epi","myeloid_frac"),("iPT_frac_epi","Blin_frac | total_imm")]
vals=[C[(C.x==k[0])&(C.y==k[1])].iloc[0] for k in keys]
x=np.arange(len(lab)); ax.bar(x,[v.rho for v in vals],color=["#6A3D9A","#6A3D9A","#999","#999","#2ca02c"],edgecolor="black")
ax.errorbar(x,[v.rho for v in vals],yerr=[[v.rho-v.ci_lo for v in vals],[v.ci_hi-v.rho for v in vals]],fmt="none",ecolor="black",capsize=4)
ax.axhline(0,color="gray",lw=1); ax.set_xticks(x); ax.set_xticklabels(lab,rotation=30,ha="right",fontsize=9)
ax.set_ylabel("Spearman rho (95% CI)"); ax.set_title("specificity: iPT vs B-lineage vs general infiltration",fontsize=11)
for s_ in ["top","right"]: ax.spines[s_].set_visible(False)
fig.suptitle("PRIMARY — cross-sample iPT~B-lineage correlation + specificity control (n=%d samples, descriptive)"%len(P),fontsize=13)
fig.tight_layout(); save(fig,"reconcile_fig_primary")
# SECONDARY
fig,ax=plt.subplots(figsize=(12,5.4)); x=np.arange(len(S))
ax.bar(x-0.27,S.prog_dz,0.27,label="injPT program (PT)",color="#6A3D9A",edgecolor="black")
ax.bar(x,S.ctrl_floor_dz,0.27,label="control-gene floor (PT)",color="#cccccc",edgecolor="black")
ax.bar(x+0.27,S.nonPT_amb_dz,0.27,label="non-PT ambient",color="#E69F00",edgecolor="black")
ax.axhline(0,color="gray",lw=1); ax.axhline(0.133,ls="--",color="#6A3D9A",lw=1.2,label="prior pooled +0.13")
for xi,r in zip(x,S.itertuples()):
    if r.iPTfrac_perm_p==r.iPTfrac_perm_p: ax.text(xi,max(r.prog_dz,0)+0.02,("*" if r.iPTfrac_perm_p<0.05 else ""),ha="center",fontsize=12,color="#d62728")
ax.set_xticks(x); ax.set_xticklabels([f"{r.section}\n{r.Condition}" for r in S.itertuples()],fontsize=8)
ax.set_ylabel("near − far Δz (within section, PT)"); ax.legend(fontsize=8.5,ncol=2)
ax.set_title("SECONDARY — per-section near/far injPT vs control-gene floor vs non-PT ambient\n(* = iPT-fraction permutation p<0.05; reconcile vs prior +0.13/6-of-9)",fontsize=11)
for s_ in ["top","right"]: ax.spines[s_].set_visible(False)
fig.tight_layout(); save(fig,"reconcile_fig_spatial")
print("\n== reconcile_extend done ==")
