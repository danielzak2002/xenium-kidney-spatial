#!/usr/bin/env python
"""STEP 2 — independent marker-based annotation of the Leiden clusters (global + immune
subclusters), BLIND to author labels. For each cluster we z-score each marker's mean log-norm
expression ACROSS clusters, then score each candidate cell type as the mean z over its present
markers and take the argmax. Marker basis (top contributing markers + z) is recorded per cluster.
Writes cluster_annotation_{global,immune}.csv and adds my_label / my_lineage / my_immune_label
to cells.parquet."""
import os, numpy as np, pandas as pd
REPO="/Users/danie/ClaudeCode/pilot_analyses/xenium"
OUT=f"{REPO}/analysis/dkd_xenium_reannotation"

# ---- canonical marker panels (kidney atlas + immune) ----
GLOBAL={
 "PT":["LRP2","CUBN","SLC5A12","SLC34A1","SLC13A3","MIOX"],
 "iPT":["VCAM1","HAVCR1","DCDC2","SPP1","ITGB6"],                 # injured PT
 "TAL":["UMOD","SLC12A1","CASR","CLDN16","ESRRB"],
 "iTAL":["UMOD","PROM1","ITGB6","DCDC2"],                         # injured TAL (UMOD+ injury)
 "DTL_ATL":["AQP1","CRYAB","AKR1B1","SLC44A5","CLDN4"],
 "DCT":["SLC12A3","TRPM6","KLHL3"],
 "CNT":["SLC8A1","CALB1","SCNN1G"],                              # connecting tubule (CALB1/SLC8A1-high)
 "PC":["AQP2","AQP3","FXYD4","GATA3"],
 "IC A":["SLC4A1","SLC26A7","ATP6V0D2","KIT","AQP6"],            # IC-B markers (SLC26A4/SLC4A9) absent from panel -> IC-B not separable

 "Podo":["NPHS2","PODXL","PTPRO","WT1","NPHS1","CLIC5"],
 "PEC":["CLDN1","CFH","RBP7","ALDH1A2"],
 "EC_glom":["EHD3","PLVAP","SOST","TBX3"],
 "EC_Peritub":["PECAM1","EMCN","IGFBP5","PLAT"],
 "EC_DVR":["PECAM1","SLC14A1","AQP1","ENPP2"],
 "EC_Lymph":["PROX1","CCL21","MMRN1","PDPN","FLT4"],
 "Fibroblast":["PDGFRA","DCN","COL1A1","COL1A2","COL6A3","C7"],
 "VSMC":["ACTA2","MYH11","TAGLN","NOTCH3","CNN1"],
 "MC":["PDGFRB","REN","POSTN","GATA3","PIEZO2"],                  # mesangial
 "Immune":["PTPRC","CD3E","CD68","LYZ","MS4A1","C1QA","CD2"],
}
# Immune is labeled by HIERARCHICAL RULES on absolute mean log-norm expression (NOT z-score):
# z-scoring across the immune compartment mislabels the T/cytotoxic axis (FOXP3/CTLA4 light up
# for any T cluster vs B/Myeloid; cytotoxic CD8 looks NK-like). CD4 is also macrophage-expressed,
# so T identity is gated on CD3E. Thresholds are on log1p(1e4-norm) means, tuned to this panel.
def annotate_immune_rules(means_csv):
    M=pd.read_csv(means_csv,index_col=0)
    def gv(cl,g): return float(M.loc[cl,g]) if g in M.columns else 0.0
    rows=[]
    for cl in M.index:
        ms4a1=gv(cl,"MS4A1"); cd79a=gv(cl,"CD79A"); mzb1=gv(cl,"MZB1"); derl3=gv(cl,"DERL3")
        cd3e=gv(cl,"CD3E"); cd8a=gv(cl,"CD8A"); cd4=gv(cl,"CD4"); foxp3=gv(cl,"FOXP3"); ctla4=gv(cl,"CTLA4")
        cd68=gv(cl,"CD68"); c1qa=gv(cl,"C1QA"); aif1=gv(cl,"AIF1"); cd163=gv(cl,"CD163")
        klrd1=gv(cl,"KLRD1"); klrf1=gv(cl,"KLRF1"); ncam1=gv(cl,"NCAM1")
        fcer1a=gv(cl,"FCER1A"); cd1c=gv(cl,"CD1C"); lilra4=gv(cl,"LILRA4")
        mast=max(gv(cl,"TPSAB1"),gv(cl,"CPA3"),gv(cl,"MS4A2")); kit=gv(cl,"KIT")
        neut=max(gv(cl,"S100A8"),gv(cl,"FCGR3B"),gv(cl,"CSF3R"))
        if ms4a1>0.8:                                   lab,basis="B",f"MS4A1={ms4a1:.2f},CD79A={cd79a:.2f}"
        elif mzb1>0.8:                                  lab,basis="Plasma",f"MZB1={mzb1:.2f},DERL3={derl3:.2f}"
        elif cd3e>0.6:                                  # T lineage gated on CD3E
            if cd8a>=0.5 and cd8a>cd4:                  lab,basis="CD8 T",f"CD3E={cd3e:.2f},CD8A={cd8a:.2f}>CD4={cd4:.2f}"
            elif foxp3>=0.5 and foxp3>0.4*cd4:          lab,basis="Treg",f"CD3E={cd3e:.2f},FOXP3={foxp3:.2f},CTLA4={ctla4:.2f}"
            else:                                       lab,basis="CD4 T",f"CD3E={cd3e:.2f},CD4={cd4:.2f},CD8A={cd8a:.2f},FOXP3={foxp3:.2f}"
        elif klrd1>0.4 or klrf1>0.3 or ncam1>0.4:       lab,basis="NK",f"KLRD1={klrd1:.2f},KLRF1={klrf1:.2f},CD3E={cd3e:.2f}"
        elif fcer1a>0.5 and cd1c>0.3 and cd68<0.4:      lab,basis="DC",f"FCER1A={fcer1a:.2f},CD1C={cd1c:.2f},LILRA4={lilra4:.2f}"
        elif mast>0.5:                                  lab,basis="Mast_Baso",f"mast(TPSAB1/CPA3/MS4A2)={mast:.2f},KIT={kit:.2f}"
        elif cd68>0.2 or c1qa>0.6 or aif1>0.5:          lab,basis="Myeloid",f"CD68={cd68:.2f},C1QA={c1qa:.2f},CD163={cd163:.2f},AIF1={aif1:.2f}"
        elif neut>0.5:                                  lab,basis="Neutrophil",f"neut(S100A8/FCGR3B/CSF3R)={neut:.2f}"
        else:                                           lab,basis="Unresolved",f"CD3E={cd3e:.2f},CD68={cd68:.2f} (low signal)"
        rows.append(dict(cluster=str(cl),label=lab,marker_basis=basis))
    return pd.DataFrame(rows)
LINEAGE={**{k:"Epithelial" for k in ["PT","iPT","TAL","iTAL","DTL_ATL","DCT","CNT","PC","IC A","IC B","Podo","PEC"]},
         **{k:"Endothelial" for k in ["EC_glom","EC_Peritub","EC_DVR","EC_Lymph"]},
         **{k:"Stroma" for k in ["Fibroblast","VSMC","MC"]}, "Immune":"Immune"}

def annotate(means_csv,panels):
    M=pd.read_csv(means_csv,index_col=0)                         # clusters x genes
    present={t:[g for g in gs if g in M.columns] for t,gs in panels.items()}
    Z=(M-M.mean(0))/(M.std(0)+1e-9)                              # z across clusters
    rows=[]
    for cl in M.index:
        best=None
        for t,gs in present.items():
            if not gs: continue
            s=float(Z.loc[cl,gs].mean())
            top=Z.loc[cl,gs].sort_values(ascending=False)
            topm=", ".join(f"{g}({Z.loc[cl,g]:+.1f})" for g in top.index[:4])
            if best is None or s>best[1]: best=(t,s,topm,gs)
        # runner-up margin (confidence)
        scores=sorted(((float(Z.loc[cl,gs].mean()),t) for t,gs in present.items() if gs),reverse=True)
        margin=scores[0][0]-scores[1][0] if len(scores)>1 else np.nan
        rows.append(dict(cluster=str(cl),label=best[0],score=round(best[1],3),margin=round(margin,3),
                         n_markers=len(best[3]),marker_basis=best[2]))
    return pd.DataFrame(rows)

ga=annotate(f"{OUT}/global_marker_means.csv",GLOBAL)
ga["lineage"]=ga.label.map(LINEAGE)
ga.to_csv(f"{OUT}/cluster_annotation_global.csv",index=False)
print("=== GLOBAL cluster annotation ==="); print(ga.to_string(index=False))

ia=annotate_immune_rules(f"{OUT}/immune_marker_means.csv")
ia.to_csv(f"{OUT}/cluster_annotation_immune.csv",index=False)
print("\n=== IMMUNE subcluster annotation (rule-based) ==="); print(ia.to_string(index=False))

# ---- map labels onto cells ----
cells=pd.read_parquet(f"{OUT}/cells.parquet")
gmap=dict(zip(ga.cluster,ga.label)); lmap=dict(zip(ga.cluster,ga.lineage))
imap=dict(zip(ia.cluster,ia.label))
cells["global_label"]=cells.leiden.astype(str).map(gmap)
cells["my_lineage"]=cells.leiden.astype(str).map(lmap)
il=cells.immune_leiden.dropna().astype(int).astype(str).map(imap)
cells["my_immune_label"]=pd.Series(np.nan,index=cells.index,dtype=object)
cells.loc[il.index,"my_immune_label"]=il.values
# final per-cell label: immune cells -> subtype; else global epithelial/EC/stroma label
cells["my_label"]=np.where(cells.is_immune_cluster & cells.my_immune_label.notna(),
                           cells.my_immune_label, cells.global_label)
# refine lineage for immune subtypes
IMM_LIN={"B":"Immune","Plasma":"Immune","CD4 T":"Immune","CD8 T":"Immune","Treg":"Immune",
         "Myeloid":"Immune","NK":"Immune","DC":"Immune","Mast_Baso":"Immune","Neutrophil":"Immune"}
cells.loc[cells.my_immune_label.notna(),"my_lineage"]="Immune"
cells.to_parquet(f"{OUT}/cells.parquet")
print("\n=== final my_label counts ===")
print(cells.my_label.value_counts(dropna=False).to_string())
print("\n=== my_lineage counts ==="); print(cells.my_lineage.value_counts(dropna=False).to_string())
print("== annotate done ==")
