import json, os, subprocess
import pandas as pd, numpy as np
d=json.load(open('/tmp/fig.json')); fs=sorted(d['files'],key=lambda x:x['size'])
MARK={'CA9','NDUFA4L2','VEGFA'}   # ccRCC tumor markers (highly specific to clear-cell)
rows=[]
for k,f in enumerate(fs,1):
    sid=f"figS{k:02d}"; gz=f"/tmp/cls_{sid}.csv.gz"
    print(f"[{sid}] download {f['size']/1e6:.0f}MB ...",flush=True)
    subprocess.run(["curl","-sL","--max-time","1800",f['download_url'],"-o",gz],check=True)
    tot=0; mk=0; cells=set()
    for ch in pd.read_csv(gz,usecols=['cell_id','feature_name'],chunksize=4_000_000,low_memory=False):
        ch=ch[ch.cell_id>0]; tot+=len(ch); mk+=int(ch.feature_name.isin(MARK).sum())
    os.remove(gz)
    rate=mk/max(tot,1)
    rows.append(dict(sample=sid,size_mb=round(f['size']/1e6),assigned_tx=tot,ccRCC_marker_tx=mk,ccRCC_rate=round(rate,5)))
    print(f"  [{sid}] assigned {tot:,} tx | CA9/NDUFA4L2/VEGFA rate {rate:.5f}",flush=True)
df=pd.DataFrame(rows)
# classify tumor vs adjacent: kmeans-ish split on rate (natural gap) — tumor = high
r=df.ccRCC_rate.values; thr=(r.max()+r.min())/2  # provisional; refine to the largest gap
sr=np.sort(r); gaps=np.diff(sr); gi=int(np.argmax(gaps)); thr=(sr[gi]+sr[gi+1])/2
df['region']=np.where(df.ccRCC_rate>=thr,'tumor','adjacent')
df.to_csv('/Users/danie/ClaudeCode/pilot_analyses/xenium/analysis/three_cohort_integration/figshare_tumor_classification.csv',index=False)
print("\n=== classification (threshold %.5f at largest gap) ===" % thr)
print(df.sort_values('ccRCC_rate',ascending=False).to_string(index=False))
print("\ntumor samples:",list(df[df.region=='tumor']['sample']))
print("== classify done ==")
