# -*- coding: utf-8 -*-
"""원인 분석: (1) PEC/PTMC 실측 비중 경향 (2) Tg 특징 추가 시 오차 변화 (3) PEO 결정 영역(<60°C) 잔차"""
import numpy as np, pandas as pd, json
from sklearn.ensemble import HistGradientBoostingRegressor
from scipy.stats import spearmanr
from load_db import load
from features import build
d = load(); d = d[~d.has_add & d["Salt1 SMILES"].astype(str).str.contains(r"\[Li\+\]")].copy()
X = build(d); y = d.logS.values
X2 = X.copy(); X2["Tg"] = pd.to_numeric(d["Tg (oC)"], errors="coerce")
mk = lambda: HistGradientBoostingRegressor(max_iter=400, learning_rate=0.05, random_state=0)
res = {}
for pp in ["PEC","PTMC","PEO"]:
    s = d[(d.poly==pp)&(d.salt=="LiTFSI")].assign(Tb=lambda z:(z.T_K/5).round())
    rs = [spearmanr(g.m, g.logS).correlation for _, g in s.groupby(["DOI","Tb"]) if g.m.nunique()>=3]
    rs = [r for r in rs if np.isfinite(r)]
    res[f"{pp}_measured_trend_sigma_vs_salt_median_rank"] = round(float(np.median(rs)),2) if rs else None
    res[f"{pp}_m_range"] = [round(float(s.m.min()),2), round(float(s.m.max()),2)]
    res[f"{pp}_Tg_median_C"] = None if s["Tg (oC)"].isna().all() else float(pd.to_numeric(s["Tg (oC)"],errors="coerce").median())
for pp in ["PEC"]:
    te = ((d.poly==pp)&(d.salt=="LiTFSI")).values
    for name, XX in [("no_Tg",X),("with_Tg",X2)]:
        m = mk().fit(XX[~te], y[~te]); yp = m.predict(XX[te])
        res[f"{pp}_MAE_{name}"] = round(float(np.mean(np.abs(yp-y[te]))),2)
# PEO-LiTFSI 온도 구간별 잔차 (논문 단위 제외 학습)
peo = ((d.poly=="PEO")&(d.salt=="LiTFSI")).values
dois = d.loc[peo,"DOI"].unique(); errs=[]
for doi in dois:
    te = (peo & (d.DOI==doi).values)
    m = mk().fit(X[~te], y[~te]); errs.append(pd.DataFrame({"T":d.T_K[te]-273.15,"err":m.predict(X[te])-y[te]}))
e = pd.concat(errs)
res["PEO_LiTFSI_paper_out_MAE_below60C"] = round(float(e[e["T"]<60].err.abs().mean()),2)
res["PEO_LiTFSI_paper_out_MAE_above60C"] = round(float(e[e["T"]>=60].err.abs().mean()),2)
res["PEO_LiTFSI_n_below_above"] = [int((e["T"]<60).sum()), int((e["T"]>=60).sum())]
print(json.dumps(res, ensure_ascii=False, indent=1))
