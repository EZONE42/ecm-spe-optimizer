#!/usr/bin/env python3
"""Quickstart Tutorial: 4-Resistance Equivalent Circuit Modeling & Composition Optimization.

Demonstrates:
1. Translating bulk ionic conductivity (sigma) and transference number (t+) to cell ASR.
2. Revealing the 44% concentration polarization bottleneck.
3. Finding optimal EO:Li ratio at 25°C and 60°C.
"""
import sys
from pathlib import Path

# Add project root to path
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ecm_spe.ecm import (
    CellParameters,
    cell_resistance_breakdown,
    de_levie_tlm,
    optimize_composition,
)
from ecm_spe.ml import peo_litfsi_empirical_conductivity


def main():
    print("=" * 72)
    print("  ECM-SPE: 4-Resistance Equivalent Circuit Modeling Tutorial")
    print("=" * 72)

    # 1. Define standard cell parameters
    cell = CellParameters(
        L_sep=20e-4,    # 20 um separator
        L_cath=60e-4,   # 60 um composite cathode
        eps=0.30,       # 30% electrolyte volume fraction in cathode
        R_int=40.0,     # 40 Ohm*cm^2 Li interface ASR
        i=1.5e-3,       # 1.5 mA/cm^2 (0.5C at 3 mAh/cm^2)
        t_plus=0.147,   # Bruce-Vincent measured Li+ transference number
    )

    print("\n[Step 1] Baseline Cell Design Parameters:")
    print(f"  Separator Thickness (L_sep)  : {cell.L_sep * 1e4:.0f} um")
    print(f"  Cathode Thickness   (L_cath) : {cell.L_cath * 1e4:.0f} um")
    print(f"  Cathode Porosity    (eps)    : {cell.eps * 100:.0f} %")
    print(f"  Interface Resistance(R_int)  : {cell.R_int:.1f} Ohm*cm^2")
    print(f"  Transference Number (t+)     : {cell.t_plus:.3f}")

    # 2. Evaluate at optimal 60°C composition (EO:Li = 18:1)
    sigma_60C = peo_litfsi_empirical_conductivity(18.0, 60.0)
    res_60C = cell_resistance_breakdown(sigma_60C, cell=cell, include_conc=True)

    print("\n[Step 2] 4-Resistance Breakdown @ 60°C (EO:Li = 18:1, sigma = {:.2e} S/cm):".format(sigma_60C))
    print(f"  1. Membrane Bulk Resistance    (R_sep) : {res_60C['R_sep']:6.2f} Ohm*cm^2 ({res_60C['share_sep']*100:4.1f}%)")
    print(f"  2. Cathode Porous TLM          (R_cath): {res_60C['R_cath_TLM']:6.2f} Ohm*cm^2 ({res_60C['share_cath']*100:4.1f}%)")
    print(f"  3. Concentration Polarization  (R_conc): {res_60C['R_conc']:6.2f} Ohm*cm^2 ({res_60C['share_conc']*100:4.1f}%)  <-- 1순위 병목!")
    print(f"  4. Interface Charge Transfer   (R_int) : {res_60C['R_int']:6.2f} Ohm*cm^2 ({res_60C['share_int']*100:4.1f}%)")
    print("  " + "-" * 56)
    print(f"  Total Area-Specific Resistance (ASR)   : {res_60C['ASR_total']:6.2f} Ohm*cm^2")
    print(f"  Cell Overpotential Loss (Delta V)      : {res_60C['dV_0p5C']:6.3f} V")

    print("\n  >> Core Physical Insight:")
    print("     Bulk membrane resistance accounts for only {:.1f}% of total resistance,".format(res_60C['share_sep'] * 100))
    print("     while concentration polarization (due to low t+ = {:.3f}) accounts for {:.1f}%.".format(cell.t_plus, res_60C['share_conc'] * 100))
    print("     Therefore, doubling t+ reduces cell overpotential far more than doubling bulk conductivity!")

    # 3. Grid optimization sweep
    print("\n[Step 3] Salt Ratio Grid Sweep (EO:Li = 4:1 ~ 40:1):")
    sweep = optimize_composition(
        peo_litfsi_empirical_conductivity,
        temperatures=(25.0, 60.0),
        cell=cell,
        include_conc=True,
    )
    for k, v in sweep["summary"].items():
        print(f"  {k}: {v}")

    print("\n" + "=" * 72)
    print("  Tutorial completed successfully.")
    print("=" * 72)


if __name__ == "__main__":
    main()
