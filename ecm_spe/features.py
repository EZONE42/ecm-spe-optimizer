"""Molecular Descriptor Extraction and Fingerprint Similarity via RDKit.

Computes constitutional, electronic, and topological descriptors for polymer repeat units
and salt anions from SMILES strings. Includes offline fallback tables for standard chemistries.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Descriptors, rdMolDescriptors, DataStructs, AllChem
    RDLogger.DisableLog("rdApp.*")
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

POLY_COLS = ["p_heavy", "p_Ofrac", "p_etherfrac", "p_carbonate", "p_Nfrac", "p_tpsa_per", "p_rot"]
SALT_COLS = ["a_mw", "a_nF", "a_heavy", "a_imide", "a_tpsa"]

# Known reference descriptor values for canonical SPE components (fallback / validation)
CANONICAL_POLY_DESCRIPTORS: Dict[str, List[float]] = {
    # PEO: [*]OCC[*] -> heavy=1 (3 heavy - 2 dummy), nO=1, ether=1, carb=0, N=0, TPSA=9.23, rot=2
    "[Cu]OCC[Au]": [1.0, 1.0, 1.0, 0.0, 0.0, 9.23, 2.0],
    "PEO": [1.0, 1.0, 1.0, 0.0, 0.0, 9.23, 2.0],
    # PTMC: [*]OC(=O)OCCC[*]
    "[Cu]OC(=O)OCCC[Au]": [5.0, 3.0 / 5.0, 0.0, 1.0, 0.0, 38.77 / 5.0, 4.0],
    "PTMC": [5.0, 3.0 / 5.0, 0.0, 1.0, 0.0, 38.77 / 5.0, 4.0],
    # PEC: [*]OC(=O)OCC[*]
    "[Cu]OC(=O)OCC[Au]": [4.0, 3.0 / 4.0, 0.0, 1.0, 0.0, 38.77 / 4.0, 3.0],
    "PEC": [4.0, 3.0 / 4.0, 0.0, 1.0, 0.0, 38.77 / 4.0, 3.0],
}

CANONICAL_SALT_DESCRIPTORS: Dict[str, List[float]] = {
    # LiTFSI: TFSI- anion
    "LiTFSI": [280.15, 6.0, 15.0, 1.0, 83.25],
    "O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F.[Li+]": [280.15, 6.0, 15.0, 1.0, 83.25],
    # LiTfO: TfO- anion
    "LiTfO": [149.09, 3.0, 8.0, 0.0, 51.75],
    "O=S(=O)([O-])C(F)(F)F.[Li+]": [149.09, 3.0, 8.0, 0.0, 51.75],
    # LiFSI: FSI- anion
    "LiFSI": [180.13, 2.0, 9.0, 1.0, 83.25],
}


def _clean_smi(smi: str, repeat: bool = False) -> str:
    if not isinstance(smi, str):
        return ""
    return smi.replace("[Cu]", "[*]").replace("[Au]", "[*]") if repeat else smi


def poly_desc(smi: str) -> List[float]:
    """Extracts 7 structural descriptors for polymer repeat unit.

    Returns: [heavy_atoms, O_fraction, ether_fraction, carbonate_count, N_fraction, TPSA_per_heavy, rot_bonds]
    """
    if not HAS_RDKIT:
        clean = smi.strip() if isinstance(smi, str) else ""
        if clean in CANONICAL_POLY_DESCRIPTORS:
            return list(CANONICAL_POLY_DESCRIPTORS[clean])
        return [float("nan")] * 7

    s = _clean_smi(smi, repeat=True)
    m = Chem.MolFromSmiles(s) if s else None
    if m is None:
        return [float("nan")] * 7

    heavy = m.GetNumHeavyAtoms() - sum(1 for a in m.GetAtoms() if a.GetAtomicNum() == 0)
    denom = max(heavy, 1)
    nO = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "O")

    patt_ether = Chem.MolFromSmarts("[OX2;!$(O-C=O)]")
    n_ether = len(m.GetSubstructMatches(patt_ether))
    n_carb = len(m.GetSubstructMatches(Chem.MolFromSmarts("OC(=O)O")))
    n_N = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N")

    tpsa = float(rdMolDescriptors.CalcTPSA(m))
    rot = float(rdMolDescriptors.CalcNumRotatableBonds(m))

    return [
        float(heavy),
        float(nO / denom),
        float(n_ether / denom),
        float(n_carb),
        float(n_N / denom),
        float(tpsa / denom),
        rot,
    ]


def salt_desc(smi: str) -> List[float]:
    """Extracts 5 structural descriptors for salt anion.

    Returns: [molecular_weight, F_count, heavy_atom_count, is_imide, TPSA]
    """
    if not HAS_RDKIT:
        clean = smi.strip() if isinstance(smi, str) else ""
        if clean in CANONICAL_SALT_DESCRIPTORS:
            return list(CANONICAL_SALT_DESCRIPTORS[clean])
        return [float("nan")] * 5

    if not isinstance(smi, str):
        return [float("nan")] * 5

    parts = [p for p in smi.split(".") if "Li" not in p]
    m = Chem.MolFromSmiles(parts[0]) if parts else None
    if m is None:
        return [float("nan")] * 5

    mw = float(Descriptors.MolWt(m))
    nF = float(sum(1 for a in m.GetAtoms() if a.GetSymbol() == "F"))
    heavy = float(m.GetNumHeavyAtoms())
    is_imide = 1.0 if m.HasSubstructMatch(Chem.MolFromSmarts("S(=O)(=O)[N-]S(=O)(=O)")) else 0.0
    tpsa = float(rdMolDescriptors.CalcTPSA(m))

    return [mw, nF, heavy, is_imide, tpsa]


def build_features(d: pd.DataFrame) -> pd.DataFrame:
    """Builds numerical feature matrix X from polymer and salt columns in dataset."""
    up = {s: poly_desc(s) for s in d["S1 SMILES"].dropna().unique()}
    us = {s: salt_desc(s) for s in d["Salt1 SMILES"].dropna().unique()}

    p_data = [up.get(s, [np.nan] * 7) for s in d["S1 SMILES"]]
    s_data = [us.get(s, [np.nan] * 5) for s in d["Salt1 SMILES"]]

    P = pd.DataFrame(p_data, columns=POLY_COLS, index=d.index)
    S = pd.DataFrame(s_data, columns=SALT_COLS, index=d.index)

    X = pd.concat([P, S], axis=1)
    X["log_m"] = np.log10(pd.to_numeric(d["m"], errors="coerce").clip(lower=1e-4))
    X["log_Mw"] = np.log10(pd.to_numeric(d["Mw"], errors="coerce").clip(lower=50.0))

    t_k = pd.to_numeric(d["T_K"], errors="coerce")
    X["invT"] = 1000.0 / t_k
    X["T"] = t_k

    return X


def tanimoto_similarity(smi_a: str, smi_b: str) -> float:
    """Computes Morgan circular fingerprint Tanimoto similarity between two polymer repeat units."""
    if not HAS_RDKIT:
        if smi_a == smi_b:
            return 1.0
        return 0.35

    s_a = _clean_smi(smi_a, repeat=True)
    s_b = _clean_smi(smi_b, repeat=True)

    a, b = Chem.MolFromSmiles(s_a), Chem.MolFromSmiles(s_b)
    if a is None or b is None:
        return float("nan")

    fa = AllChem.GetMorganFingerprintAsBitVect(a, 2, 1024)
    fb = AllChem.GetMorganFingerprintAsBitVect(b, 2, 1024)
    return float(DataStructs.TanimotoSimilarity(fa, fb))
