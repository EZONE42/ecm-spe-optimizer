# -*- coding: utf-8 -*-
"""조성 → 숫자 특징 변환 (RDKit). 고분자는 반복단위, 염은 음이온 기준"""
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors, DataStructs
from rdkit.Chem import AllChem
RDLogger.DisableLog("rdApp.*")
def _mol(smi, repeat=False):
    if not isinstance(smi, str): return None
    s = smi.replace("[Cu]", "[*]").replace("[Au]", "[*]") if repeat else smi
    m = Chem.MolFromSmiles(s)
    return m
def poly_desc(smi):
    m = _mol(smi, True)
    if m is None: return [np.nan]*7
    heavy = m.GetNumHeavyAtoms() - sum(1 for a in m.GetAtoms() if a.GetAtomicNum()==0)
    nO = sum(1 for a in m.GetAtoms() if a.GetSymbol()=="O")
    patt_ether = Chem.MolFromSmarts("[OX2;!$(O-C=O)]")
    n_ether = len(m.GetSubstructMatches(patt_ether))
    n_carb = len(m.GetSubstructMatches(Chem.MolFromSmarts("OC(=O)O")))
    n_N = sum(1 for a in m.GetAtoms() if a.GetSymbol()=="N")
    return [heavy, nO/max(heavy,1), n_ether/max(heavy,1), n_carb, n_N/max(heavy,1),
            rdMolDescriptors.CalcTPSA(m)/max(heavy,1), rdMolDescriptors.CalcNumRotatableBonds(m)]
POLY_COLS = ["p_heavy","p_Ofrac","p_etherfrac","p_carbonate","p_Nfrac","p_tpsa_per","p_rot"]
def salt_desc(smi):
    if not isinstance(smi, str): return [np.nan]*5
    parts = [p for p in smi.split(".") if "Li" not in p]
    m = Chem.MolFromSmiles(parts[0]) if parts else None
    if m is None: return [np.nan]*5
    return [Descriptors.MolWt(m), sum(1 for a in m.GetAtoms() if a.GetSymbol()=="F"),
            m.GetNumHeavyAtoms(), int(m.HasSubstructMatch(Chem.MolFromSmarts("S(=O)(=O)[N-]S(=O)(=O)"))),
            rdMolDescriptors.CalcTPSA(m)]
SALT_COLS = ["a_mw","a_nF","a_heavy","a_imide","a_tpsa"]
def build(d):
    up = {s: poly_desc(s) for s in d["S1 SMILES"].dropna().unique()}
    us = {s: salt_desc(s) for s in d["Salt1 SMILES"].dropna().unique()}
    P = pd.DataFrame([up.get(s,[np.nan]*7) for s in d["S1 SMILES"]], columns=POLY_COLS, index=d.index)
    S = pd.DataFrame([us.get(s,[np.nan]*5) for s in d["Salt1 SMILES"]], columns=SALT_COLS, index=d.index)
    X = pd.concat([P,S], axis=1)
    X["log_m"] = np.log10(d["m"].clip(lower=1e-4))
    X["log_Mw"] = np.log10(d["Mw"].clip(lower=50))
    X["invT"] = 1000.0/d["T_K"]
    X["T"] = d["T_K"]
    return X
def similarity(smi_a, smi_b):
    a, b = _mol(smi_a, True), _mol(smi_b, True)
    fa = AllChem.GetMorganFingerprintAsBitVect(a, 2, 1024); fb = AllChem.GetMorganFingerprintAsBitVect(b, 2, 1024)
    return DataStructs.TanimotoSimilarity(fa, fb)
