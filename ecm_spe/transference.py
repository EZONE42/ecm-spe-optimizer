"""Electrochemical Impedance Spectroscopy (EIS) and Transference Number Analysis.

Implements:
1. Kasa geometric circle fitting for interfacial and bulk EIS arcs.
2. Bruce-Vincent steady-state chronoamperometry + EIS lithium transference number (t+).
3. Symmetric cell interfacial area-specific resistance (ASR_int).
"""
from __future__ import annotations

import math
from typing import Dict, Optional, Tuple

import numpy as np


def circle_fit(
    x: np.ndarray,
    y: np.ndarray,
) -> Tuple[float, float, float, float, float]:
    """Fits a circle to 2D impedance points (Re(Z), -Im(Z)) using the Kasa algebraic method.

    Minimizes algebraic distance: sum ( (x - xc)^2 + (y - yc)^2 - r^2 )^2
    Equation: x^2 + y^2 = c0 * x + c1 * y + c2
    where xc = c0 / 2, yc = c1 / 2, r = sqrt(c2 + xc^2 + yc^2)

    Args:
        x: Real impedance array Re(Z).
        y: Imaginary impedance array -Im(Z).

    Returns:
        Tuple (x_left_intercept, x_right_intercept, xc, yc, radius).
        Intercepts are the points where the fitted circle crosses the real axis (y=0).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3:
        raise ValueError("Circle fitting requires at least 3 points.")

    A = np.c_[x, y, np.ones_like(x)]
    b = x ** 2 + y ** 2
    c = np.linalg.lstsq(A, b, rcond=None)[0]

    xc = float(c[0] / 2.0)
    yc = float(c[1] / 2.0)
    r_sq = float(c[2] + xc ** 2 + yc ** 2)

    if r_sq < 0.0:
        return float("nan"), float("nan"), xc, yc, float("nan")

    r = math.sqrt(r_sq)
    h_sq = r ** 2 - yc ** 2
    h = math.sqrt(max(h_sq, 0.0))

    return xc - h, xc + h, xc, yc, r


def bruce_vincent_tplus(
    Rb0: float,
    Ri0: float,
    Ri_ss: float,
    I0: float,
    Iss: float,
    dV: float = 0.05,
) -> float:
    """Calculates cation transference number via the Bruce-Vincent-Evans steady-state method.

    t+ = Iss * (dV - I0 * Ri0) / (I0 * (dV - Iss * Ri_ss))

    Accounting for initial (pre-polarization) and steady-state (post-polarization)
    interfacial resistance growth during chronoamperometry.

    Args:
        Rb0: Initial bulk electrolyte resistance (Ohm).
        Ri0: Initial interfacial charge-transfer resistance (Ohm).
        Ri_ss: Steady-state interfacial charge-transfer resistance (Ohm).
        I0: Initial polarization current immediately after step (A).
        Iss: Steady-state polarization current (A).
        dV: Applied potentiostatic voltage step (V, typically 10 to 50 mV = 0.01~0.05 V).

    Returns:
        Lithium cation transference number t+ (dimensionless, typically 0.05 ~ 0.50 for SPE).
    """
    if I0 <= 0.0 or Iss <= 0.0 or dV <= 0.0:
        return float("nan")

    num = Iss * (dV - I0 * Ri0)
    den = I0 * (dV - Iss * Ri_ss)

    if den == 0.0:
        return float("nan")

    return float(num / den)


def interfacial_asr(
    Ri: float,
    area_cm2: float,
    is_symmetric: bool = True,
) -> float:
    """Calculates single-interface area-specific resistance (ASR) from cell impedance.

    For symmetric cell (Li | SPE | Li), total Ri reflects two identical electrode interfaces:
    ASR_int = Ri * area_cm2 / 2

    Args:
        Ri: Interfacial charge-transfer resistance from EIS arc (Ohm).
        area_cm2: Active electrode geometric surface area (cm^2).
        is_symmetric: True if two identical interfaces in series.

    Returns:
        Interfacial area-specific resistance in Ohm*cm^2.
    """
    divisor = 2.0 if is_symmetric else 1.0
    return float(Ri * area_cm2 / divisor)


def extract_eis_intercepts(
    re: np.ndarray,
    mim: np.ndarray,
) -> Dict[str, float]:
    """Extracts bulk resistance (Rb) and interfacial resistance (Ri) from a single Nyquist spectrum.

    Rb is identified at the high-frequency / bulk minimum phase angle valley.
    Rtot is estimated via Kasa circle fitting of the interfacial semicircle.
    """
    re = np.asarray(re, dtype=float)
    mim = np.asarray(mim, dtype=float)

    if len(re) < 5:
        raise ValueError("Nyquist spectrum too short for intercept extraction.")

    # Find valley (minimum phase angle atan(-Im / Re) or first local minimum of -Im)
    valleys = [i for i in range(1, len(mim) - 1) if mim[i] <= mim[i - 1] and mim[i] <= mim[i + 1]]
    i_valley = valleys[0] if valleys else int(np.argmin(mim[len(mim) // 4:]))

    Rb = float(re[i_valley])

    # Interfacial arc apex and circle fit
    apex = i_valley + int(np.argmax(mim[i_valley:])) if i_valley < len(mim) - 1 else i_valley
    sel = slice(i_valley + 1, min(len(re), apex + 15))

    if sel.stop - sel.start >= 3:
        lo, hi, xc, yc, r = circle_fit(re[sel], mim[sel])
        Rtot = hi if not math.isnan(hi) else float(re[-1])
    else:
        Rtot = float(re[-1])

    Ri = max(0.0, Rtot - Rb)

    return {
        "Rb": Rb,
        "Rtot": Rtot,
        "Ri": Ri,
    }
