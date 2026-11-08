# -*- coding: utf-8 -*-
"""교차 검증: A:B로 학습 → (D:B) 또는 (A:C)를 한 번도 안 보여 주고 예측 → 실측과 비교
실행: PYTHONPATH=. python3 transfer_test.py"""
import numpy as np, pandas as pd, json
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import GroupKFold
from scipy.stats import spearmanr
from load_db import load
from features import build

d = load()
d = d[~d.has_add & d["Salt1 SMILES"].astype(str).str.contains(r"\[Li\+\]")].copy()
X = build(d); y = d["logS"].values
feat = list(X.columns)
def model(): return HistGradientBoostingRegressor(max_iter=400, learning_rate=0.05, max_leaf_nodes=31, random_state=0)
out = {}
# 1) 기준 성능: 논문 단위 5-fold (같은 논문이 학습·시험에 동시에 들어가지 않게)
gkf = GroupKFold(n_splits=5); pred = np.zeros(len(y))
for tr, te in gkf.split(X, y, groups=d["DOI"].fillna("na")):
    m = model().fit(X.iloc[tr], y[tr]); pred[te] = m.predict(X.iloc[te])
out["baseline_paper_split_MAE"] = float(np.mean(np.abs(pred - y)))
# 2) 측정 잡음 바닥: 같은 PEO-LiTFSI 조성(m 0.1 단위 반올림)·온도(5 K 반올림)에서 논문 간 차이
p = d[(d.poly=="PEO")&(d.salt=="LiTFSI")].copy()
p["mb"] = p["m"].round(1); p["Tb"] = (p["T_K"]/5).round()*5
g = p.groupby(["mb","Tb"]).filter(lambda x: x.DOI.nunique()>=2).groupby(["mb","Tb"])
noise = g.apply(lambda x: (x.logS - x.logS.median()).abs().median())
out["noise_floor_between_papers_median_abs_dev"] = float(noise.median()); out["noise_floor_groups"] = int(len(noise))
# 3) 교차 사례
cases = {
 "A:C 염 유사체  PEO:LiFSI (B=LiTFSI와 같은 이미드 계열)": ("PEO","LiFSI"),
 "A:C 염 유사체  PEO:LiBETI (이미드, 더 긴 불소 사슬)": ("PEO","LiBETI"),
 "A:C 염 다른계열 PEO:LiTfO (설폰산염)": ("PEO","LiTfO"),
 "D:B 고분자 유사체 P(EO-MO):LiTFSI (PEO 유사도 0.47)": ("[Cu]OCOCCOCC[Au]","LiTFSI"),
 "D:B 고분자 유사체 PEC:LiTFSI (유사도 0.38)": ("PEC","LiTFSI"),
 "D:B 고분자 다른계열 PTMC:LiTFSI (유사도 0.29)": ("PTMC","LiTFSI"),
}
rows = []
for name, (pp, ss) in cases.items():
    te = (d.poly==pp)&(d.salt==ss)
    m = model().fit(X[~te], y[~te]); yp = m.predict(X[te]); yt = y[te.values]
    sub = d[te].assign(pred=yp)
    # 비중(몰랄농도) 경향: 같은 논문·같은 온도 묶음 안에서 순위 상관
    rhos = []
    for _, s in sub.assign(Tb=(sub.T_K/5).round()).groupby(["DOI","Tb"]):
        if s.m.nunique()>=3:
            r = spearmanr(s.m, s.logS).correlation; rp = spearmanr(s.m, s.pred).correlation
            if np.isfinite(r) and np.isfinite(rp): rhos.append(spearmanr(s.logS, s.pred).correlation)
    # 최적 비중 위치 일치: 각 묶음의 σ 최대 몰랄농도 실측 vs 예측
    opt = []
    for _, s in sub.assign(Tb=(sub.T_K/5).round()).groupby(["DOI","Tb"]):
        if s.m.nunique()>=4:
            opt.append((s.loc[s.logS.idxmax(),"m"], s.loc[s.pred.idxmax(),"m"]))
    rows.append({"case":name, "rows":int(te.sum()), "papers":int(d[te].DOI.nunique()),
                 "MAE_log":round(float(np.mean(np.abs(yp-yt))),2), "bias_log":round(float(np.mean(yp-yt)),2),
                 "fold_err":round(float(10**np.mean(np.abs(yp-yt))),1),
                 "trend_rank_corr_median":round(float(np.median(rhos)),2) if rhos else None, "n_series":len(rhos),
                 "opt_m_same":f"{sum(abs(a-b)<0.35 for a,b in opt)}/{len(opt)}" if opt else None})
out["cases"] = rows
print(json.dumps(out, ensure_ascii=False, indent=1))
import os
out_dir = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(out_dir, exist_ok=True)
pd.DataFrame(rows).to_csv(os.path.join(out_dir, "transfer_results.csv"), index=False)
