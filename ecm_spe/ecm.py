"""4-Resistance Equivalent Circuit Modeling (ECM) and Cell Overpotential Scorer.

Evaluates solid polymer electrolyte performance across 4 physical resistance components:
1. R_sep: Bulk membrane ionic resistance (L_sep / sigma)
2. R_cath: Cathode porous transmission line model (de Levie TLM with Bruggeman tortuosity)
3. R_conc: Steady-state concentration polarization resistance (R_sep * (1/t+ - 1))
4. R_int: Electrode-electrolyte interfacial charge-transfer ASR
Total cell overpotential: Delta_V = i * (R_sep + R_cath + R_conc + R_int)
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union

import numpy as np


@dataclass
class CellParameters:
    """Standard lithium polymer cell geometry and operational parameters."""
    L_sep: float = 20e-4        # Separator thickness in cm (20 um)
    L_cath: float = 60e-4       # Cathode thickness in cm (60 um)
    eps: float = 0.30           # Cathode composite electrolyte volume fraction (porosity)
    R_int: float = 40.0         # Total interfacial ASR in Ohm*cm^2 (assumed / measured)
    i: float = 1.5e-3           # Current density in A/cm^2 (0.5C at 3 mAh/cm^2)
    bruggeman_exp: float = 1.5  # MacMullin / Bruggeman tortuosity exponent
    t_plus: float = 0.147       # Default Li+ transference number (BPE-LiTFSI 15:1 55C)


DEFAULT_CELL = CellParameters()
EO_MW_KG = 44.05e-3  # kg/mol (Ethylene oxide repeat unit)


def de_levie_tlm(
    sigma: float,
    cell: Optional[CellParameters] = None,
) -> float:
    """Calculates effective cathode porous ionic resistance via de Levie transmission line model.

    R_ion = L_cath / (sigma * eps^1.5)
    nu = sqrt(R_ion / R_ct)
    R_cath = sqrt(R_ion * R_ct) * coth(nu)

    Args:
        sigma: Bulk ionic conductivity in S/cm.
        cell: Cell design parameters (defaults to DEFAULT_CELL).

    Returns:
        Effective cathode transmission line resistance in Ohm*cm^2.
    """
    if cell is None:
        cell = DEFAULT_CELL
    if sigma <= 0.0 or math.isnan(sigma):
        return float("inf")

    sig_eff = sigma * (cell.eps ** cell.bruggeman_exp)
    Ri = cell.L_cath / sig_eff
    Rct = cell.R_int

    if Rct <= 0.0:
        return Ri

    nu = math.sqrt(Ri / Rct)
    if nu > 50.0:  # coth(nu) -> 1.0 for large arguments to avoid overflow
        return math.sqrt(Ri * Rct)
    if nu < 1e-6:  # Taylor expansion: coth(x) ~ 1/x + x/3
        return Rct + Ri / 3.0

    return math.sqrt(Ri * Rct) / math.tanh(nu)


def concentration_polarization_resistance(
    R_sep: float,
    t_plus: float,
) -> float:
    """Computes steady-state concentration polarization resistance from Li+ transference number.

    Under DC polarization, immobile anions accumulate, amplifying bulk membrane resistance by:
    R_conc = R_sep * (1 / t+ - 1)

    Args:
        R_sep: Membrane bulk resistance in Ohm*cm^2.
        t_plus: Lithium cation transference number (0 < t+ <= 1.0).

    Returns:
        Concentration polarization resistance in Ohm*cm^2.
    """
    if t_plus <= 0.0 or math.isnan(t_plus):
        return float("inf")
    if t_plus >= 1.0:
        return 0.0
    return R_sep * (1.0 / t_plus - 1.0)


def cell_resistance_breakdown(
    sigma: float,
    t_plus: Optional[float] = None,
    cell: Optional[CellParameters] = None,
    include_conc: bool = True,
) -> Dict[str, float]:
    """Evaluates the complete cell ASR budget and voltage loss at specified conductivity.

    Args:
        sigma: Bulk ionic conductivity in S/cm.
        t_plus: Lithium transference number (defaults to cell.t_plus if include_conc=True).
        cell: Cell design parameters.
        include_conc: Whether to include concentration polarization in the total ASR.

    Returns:
        Dictionary containing individual resistances, total ASR, voltage loss, and percentage shares.
    """
    if cell is None:
        cell = DEFAULT_CELL

    R_sep = cell.L_sep / sigma if sigma > 0 else float("inf")
    R_cath = de_levie_tlm(sigma, cell)
    R_int = cell.R_int

    tp = t_plus if t_plus is not None else cell.t_plus
    R_conc = concentration_polarization_resistance(R_sep, tp) if include_conc else 0.0

    ASR_total = R_sep + R_cath + R_conc + R_int
    dV = cell.i * ASR_total

    shares = {
        "share_sep": R_sep / ASR_total if ASR_total > 0 else 0.0,
        "share_cath": R_cath / ASR_total if ASR_total > 0 else 0.0,
        "share_conc": R_conc / ASR_total if ASR_total > 0 else 0.0,
        "share_int": R_int / ASR_total if ASR_total > 0 else 0.0,
    }

    return {
        "sigma": float(sigma),
        "t_plus": float(tp),
        "R_sep": float(R_sep),
        "R_cath_TLM": float(R_cath),
        "R_conc": float(R_conc),
        "R_int": float(R_int),
        "ASR_total": float(ASR_total),
        "dV_0p5C": float(dV),
        **shares,
    }


def optimize_composition(
    conductivity_func,
    grid: Optional[np.ndarray] = None,
    temperatures: Tuple[float, ...] = (25.0, 60.0),
    cell: Optional[CellParameters] = None,
    include_conc: bool = False,
) -> Dict[str, Union[dict, list]]:
    """Sweeps salt ratio (EO:Li) grid across temperatures to find minimum dV composition.

    Args:
        conductivity_func: Callable (EO_Li: float, temp_c: float) -> sigma (S/cm).
        grid: Array of EO:Li ratios (default 4.0 to 40.0 by 1.0).
        temperatures: Temperatures in Celsius to evaluate.
        cell: Cell parameters.
        include_conc: Whether to include concentration polarization in the ranking score.

    Returns:
        Structured dictionary with best composition per temperature and full sweep rows.
    """
    if grid is None:
        grid = np.round(np.arange(4.0, 41.0, 1.0), 1)
    if cell is None:
        cell = DEFAULT_CELL

    all_rows = []
    summary = {}

    for T in temperatures:
        temp_rows = []
        for r in grid:
            sig = float(conductivity_func(r, T))
            res = cell_resistance_breakdown(sig, cell=cell, include_conc=include_conc)
            row = {
                "T_C": float(T),
                "EO_Li": float(r),
                "molality": float(1.0 / (r * EO_MW_KG)),
                "sigma_pred": sig,
                "R_sep": res["R_sep"],
                "R_cath_TLM": res["R_cath_TLM"],
                "R_conc": res["R_conc"],
                "R_int": res["R_int"],
                "ASR_total": res["ASR_total"],
                "dV_0p5C": res["dV_0p5C"],
            }
            temp_rows.append(row)
            all_rows.append(row)

        # Find best row by lowest dV (or lowest ASR_total)
        best = min(temp_rows, key=lambda x: x["dV_0p5C"])
        t_key = f"{int(T)}C" if T.is_integer() else f"{T}C"
        summary[f"{t_key}_best_EO_Li"] = best["EO_Li"]
        summary[f"{t_key}_best_sigma"] = best["sigma_pred"]
        summary[f"{t_key}_best_dV"] = round(best["dV_0p5C"], 4)
        summary[f"{t_key}_sep_share_at_best"] = round(best["R_sep"] / best["ASR_total"], 4)

    return {"summary": summary, "rows": all_rows}
