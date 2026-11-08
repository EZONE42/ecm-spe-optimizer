# -*- coding: utf-8 -*-
"""CheMixHub 고분자 전해질 DB 로드·정규화 (Bradford et al. 2023, CC BY 4.0)"""
import os
import pandas as pd, numpy as np
RAW = os.environ.get("CHEMIXHUB_CSV", os.path.join(os.path.dirname(__file__), "..", "data", "PolymerElectrolyteData.csv"))
POLY = {"[Cu]OCC[Au]":"PEO", "[Cu]OC(=O)OCCC[Au]":"PTMC", "[Cu]OC(=O)OCC[Au]":"PEC"}
SALT = {"O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F.[Li+]":"LiTFSI",
        "O=S(=O)([O-])C(F)(F)F.[Li+]":"LiTfO", "[Li+].[O-][Cl+3]([O-])([O-])[O-]":"LiClO4"}
def load():
    d = pd.read_csv(RAW, low_memory=False)
    d["poly"] = d["S1 SMILES"].map(POLY).fillna(d["S1 SMILES"])
    s = d["Salt1 SMILES"].astype(str)
    d["salt"] = s.map(SALT)
    d.loc[d["salt"].isna() & s.str.contains(r"\[N-\]\(S\(=O\)\(=O\)F\)|FS\(=O\)\(=O\)\[N-\]S\(=O\)\(=O\)F|O=S\(=O\)\(F\)\[N-\]S\(=O\)\(=O\)F", regex=True), "salt"] = "LiFSI"
    d["salt"] = d["salt"].fillna(d["Salt 1"])
    d["T_K"] = pd.to_numeric(d["Temperature (oC)"], errors="coerce") + 273.15
    d["logS"] = pd.to_numeric(d["log Conductivity (S/cm)"], errors="coerce")
    d["m"] = pd.to_numeric(d["Salt1 Molality (mol salt/kg polymer)"], errors="coerce")
    d["Mw"] = pd.to_numeric(d["S1 Mn or Mw"], errors="coerce")
    d["has_add"] = d[["Solvent 2 (S2)","Solvent 3 (S3)","Solvent 4 (S4)","Salt 2","Inorganic Material 1 (IM1)"]].notna().any(axis=1)
    d["add_name"] = d["Solvent 2 (S2)"].fillna(d["Inorganic Material 1 (IM1)"])
    d["add_wt"] = pd.to_numeric(d["S2 Weight %"], errors="coerce").fillna(pd.to_numeric(d["IM1 Weight %"], errors="coerce"))
    d = d.dropna(subset=["T_K","logS"])
    return d
if __name__ == "__main__":
    d = load()
    print("rows", len(d), "papers", d.DOI.nunique())
    b = d[~d.has_add]
    g = b.groupby(["poly","salt"]).agg(rows=("logS","size"), papers=("DOI","nunique"),
        m_min=("m","min"), m_max=("m","max"), nm=("m","nunique"), Tmin=("T_K","min"), Tmax=("T_K","max"))
    print(g.sort_values("rows", ascending=False).head(15).to_string())
    a = d[d.has_add & (d.poly=="PEO") & (d.salt=="LiTFSI")]
    print("\nPEO-LiTFSI+additive rows", len(a), "papers", a.DOI.nunique())
    print(a.groupby("add_name").agg(rows=("logS","size"), papers=("DOI","nunique"), wt_n=("add_wt","nunique"), wt_min=("add_wt","min"), wt_max=("add_wt","max")).sort_values("rows",ascending=False).head(12).to_string())
    print("\nunmapped salt SMILES sample containing F and N-:")
    print(d.loc[d["Salt 1"].astype(str).str.contains("FSI"), ["Salt 1","Salt1 SMILES"]].drop_duplicates().head(8).to_string())
