"""External test: ML trained on CheMixHub (Li, no additive) predicts the KIT BPE-LiTFSI sigma(T)
that was never in the training set (DOI absent). BPE approximated as PEO repeat unit + LiTFSI.
Run: PYTHONPATH=. python3 external_test.py"""
import os, numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from scipy.stats import spearmanr
from load_db import load
from features import build

OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "results")
os.makedirs(OUT_DIR, exist_ok=True)

d = load()
d = d[~d.has_add & d["Salt1 SMILES"].astype(str).str.contains(r"\[Li\+\]")].copy()
X = build(d); y = d["logS"].values
m = HistGradientBoostingRegressor(max_iter=400, learning_rate=0.05, max_leaf_nodes=31, random_state=0).fit(X, y)

tmpl = d[(d.poly == "PEO") & (d.salt == "LiTFSI")].iloc[0]
G = pd.read_csv(os.path.join(OUT_DIR, "eis_anchor_sigma_by_T.csv")); G = G[G.cation == "Li"]
out = []
for Mw in (1e5, 1e6):              # BPE Mw not given in README -> two bounds
    rows = []
    for _, g in G.iterrows():
        r = tmpl.copy(); r["m"] = 1/(g.EO*0.04405); r["T_K"] = g.T_C + 273.15; r["Mw"] = Mw
        rows.append(r)
    Q = pd.DataFrame(rows); pred = m.predict(build(Q))
    out.append(G.assign(Mw=Mw, logS_meas=np.log10(G.sigma), logS_pred=pred))
E = pd.concat(out); E["err"] = E.logS_pred - E.logS_meas
E.to_csv(os.path.join(OUT_DIR, "external_test_results.csv"), index=False)
pd.set_option("display.width", 200)
for (Mw, dirn), g in E.groupby(["Mw", "direction"]):
    best_meas = g.loc[g.groupby("T_C").logS_meas.idxmax()][["T_C", "EO"]].set_index("T_C").EO
    best_pred = g.loc[g.groupby("T_C").logS_pred.idxmax()][["T_C", "EO"]].set_index("T_C").EO
    rho = [spearmanr(s.logS_meas, s.logS_pred).correlation for _, s in g.groupby("T_C")]
    print(f"Mw={Mw:.0e} {dirn}: MAE={g.err.abs().mean():.2f} log (x{10**g.err.abs().mean():.1f}), bias={g.err.mean():+.2f}, "
          f"MAE T>=55C={g[g.T_C>=55].err.abs().mean():.2f}, MAE T<55C={g[g.T_C<55].err.abs().mean():.2f}, "
          f"best-EO match {int((best_meas==best_pred).sum())}/{len(best_meas)}, EO-rank rho median={np.nanmedian(rho):.2f}")
    print("   best EO meas:", dict(best_meas), " pred:", dict(best_pred))
print(E[(E.Mw == 1e6) & (E.direction == "heating")].pivot(index="T_C", columns="EO", values=["logS_meas", "logS_pred"]).round(2))
