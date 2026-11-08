# -*- coding: utf-8 -*-
"""조성 최적화 시연: PEO:LiTFSI 비중(EO:Li) 탐색 → ECM 채점 (막 + 양극 이온경로 전송선 + 계면[가정])
실행: PYTHONPATH=. python3 optimize_demo.py"""
import numpy as np, pandas as pd, json
from sklearn.ensemble import HistGradientBoostingRegressor
from load_db import load
from features import build
d = load(); d = d[~d.has_add & d["Salt1 SMILES"].astype(str).str.contains(r"\[Li\+\]")].copy()
X = build(d); y = d.logS.values
m = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.05, random_state=0).fit(X, y)
ref = d[(d.poly=="PEO")&(d.salt=="LiTFSI")].iloc[0]
EO = 44.05e-3  # kg/mol (EO 반복단위)
grid = np.round(np.arange(4, 41, 1.0), 1)          # EO:Li 4~40
cell = dict(L_sep=20e-4, L_cath=60e-4, eps=0.30, R_int=40.0, i=1.5e-3)  # cm, -, Ω·cm²[가정], A/cm² (0.5C, 3 mAh/cm²)
def tlm(sig):   # 전송선 모델: 양극 이온저항 Ri, 전하이동(계면) Rct → 유효저항 = sqrt(Ri*Rct)*coth(sqrt(Ri/Rct))
    Ri = cell["L_cath"]/(sig*cell["eps"]**1.5); Rct = cell["R_int"]
    nu = np.sqrt(Ri/Rct); return np.sqrt(Ri*Rct)/np.tanh(nu)
rows=[]
for T in [25, 60]:
    for r in grid:
        mol = 1.0/(r*EO)
        dd = d.iloc[[0]].copy(); dd["S1 SMILES"]=ref["S1 SMILES"]; dd["Salt1 SMILES"]=ref["Salt1 SMILES"]
        dd["m"]=mol; dd["Mw"]=6e5; dd["T_K"]=T+273.15
        sig = 10**m.predict(build(dd))[0]
        Rb = cell["L_sep"]/sig; Rc = tlm(sig); tot = Rb+Rc
        rows.append({"T_C":T,"EO_Li":r,"sigma_pred":sig,"R_sep":Rb,"R_cath_TLM":Rc,"ASR_total":tot,"dV_0p5C":cell["i"]*tot})
r = pd.DataFrame(rows)
import os
out_dir = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(out_dir, exist_ok=True)
r.to_csv(os.path.join(out_dir, "optimize_demo.csv"), index=False)
out={}
for T in [25,60]:
    s=r[r.T_C==T]; b=s.loc[s.dV_0p5C.idxmin()]
    out[f"{T}C_best_EO_Li"]=float(b.EO_Li); out[f"{T}C_best_sigma"]=float(f"{b.sigma_pred:.2e}"); out[f"{T}C_best_dV"]=round(float(b.dV_0p5C),3)
    out[f"{T}C_sep_share_at_best"]=round(float(b.R_sep/b.ASR_total),2)
# 실측 비교: DB에서 PEO-LiTFSI 고분자량(>=1e5) 60°C 근처(55~65) 최고 σ의 EO:Li
x=d[(d.poly=="PEO")&(d.salt=="LiTFSI")&(d.Mw>=1e5)&(d.T_K.between(328,338))].copy()
x["EO_Li"]=1/(x.m*EO)
best=x.loc[x.logS.idxmax()]
out["measured_60C_best_EO_Li"]=round(float(best.EO_Li),1); out["measured_60C_best_sigma"]=float(f"{10**best.logS:.2e}"); out["measured_rows_55_65C"]=int(len(x))
out["measured_60C_EO_Li_of_top5"]=sorted([round(float(v),1) for v in x.nlargest(5,"logS").EO_Li])
print(json.dumps(out,indent=1))
