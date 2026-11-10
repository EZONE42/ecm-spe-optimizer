"""Command-line interface for ECM-SPE.

Usage:
  python -m ecm_spe <command> [options]
  ecm-spe <command> [options]

Commands:
  optimize       Runs 4-resistance ECM composition optimization sweep
  transference   Runs Bruce-Vincent t+ and Kasa circle fitting demonstration
  report         Prints the 4-resistance cell breakdown and bottleneck analysis
  all            Runs full analysis pipeline
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from .ecm import (
    DEFAULT_CELL,
    CellParameters,
    cell_resistance_breakdown,
    optimize_composition,
)
from .ml import ConductivityModel, peo_litfsi_empirical_conductivity
from .transference import bruce_vincent_tplus, circle_fit, interfacial_asr

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "results"


def cmd_optimize(args) -> dict:
    """Executes EO:Li ratio optimization sweep across 25C and 60C."""
    cell = CellParameters(
        L_sep=args.L_sep * 1e-4,
        L_cath=args.L_cath * 1e-4,
        eps=args.eps,
        R_int=args.R_int,
        i=args.current_density * 1e-3,
        t_plus=args.t_plus,
    )

    model = ConductivityModel()
    cond_func = lambda r, T: model.predict_peo_litfsi(r, T)

    grid = np.round(np.arange(args.grid_min, args.grid_max + 0.1, args.grid_step), 1)
    res = optimize_composition(
        cond_func,
        grid=grid,
        temperatures=(25.0, 60.0),
        cell=cell,
        include_conc=args.include_conc,
    )

    summary = res["summary"]
    print("\n=== ECM-SPE Composition Optimization Results ===")
    for k, v in summary.items():
        print(f"  {k}: {v}")

    out_csv = getattr(args, "out", None)
    if out_csv is None and getattr(args, "save", False):
        out_csv = RESULTS_DIR / "optimize_demo.csv"

    if out_csv is not None:
        try:
            Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(res["rows"]).to_csv(out_csv, index=False)
            print(f"\nSaved sweep details -> {out_csv}")
        except (PermissionError, OSError) as e:
            print(f"\n[Note] Could not write CSV ({e}); proceeding with memory results.")
    return res


def cmd_transference(args) -> dict:
    """Demonstrates Bruce-Vincent transference number calculation and Kasa circle fitting."""
    # Synthetic / reference values matching KIT Zenodo 7428858 BPE-LiTFSI (15:1, 55C)
    Rb0 = 120.5   # Ohm
    Ri0 = 85.2    # Ohm
    Ri_ss = 112.4 # Ohm
    I0 = 0.05 / (Rb0 + Ri0) # 0.05 V / 205.7 Ohm = 243.1 uA
    Iss = 38.6e-6 # 38.6 uA
    dV = 0.05    # 50 mV

    t_plus = bruce_vincent_tplus(Rb0, Ri0, Ri_ss, I0, Iss, dV=dV)
    asr = interfacial_asr(Ri0, area_cm2=0.5026, is_symmetric=True)

    out = {
        "Rb0_Ohm": Rb0,
        "Ri0_Ohm": Ri0,
        "Ri_ss_Ohm": Ri_ss,
        "I0_uA": round(I0 * 1e6, 2),
        "Iss_uA": round(Iss * 1e6, 2),
        "t_plus": round(t_plus, 4),
        "ASR_int_Ohm_cm2": round(asr, 2),
    }

    print("\n=== Bruce-Vincent Transference Number & Interface ASR ===")
    for k, v in out.items():
        print(f"  {k}: {v}")
    return out


def cmd_report(args=None) -> dict:
    """Prints the 4-resistance cell breakdown demonstrating the 44% concentration polarization bottleneck."""
    # Condition A: KIT BPE-LiTFSI experimental benchmark (55°C, EO:Li 15:1, L_sep 30 um, t+ 0.147)
    cell_bpe = CellParameters(L_sep=30e-4, L_cath=60e-4, eps=0.30, R_int=116.8, t_plus=0.147, i=1.5e-3)
    res_bpe = cell_resistance_breakdown(6.1e-5, cell=cell_bpe, include_conc=True)

    # Condition B: High-temperature PEO-LiTFSI design target (60°C, EO:Li 18:1, L_sep 20 um)
    sigma_60C = peo_litfsi_empirical_conductivity(18.0, 60.0)
    res_60C = cell_resistance_breakdown(sigma_60C, t_plus=0.147, include_conc=True)

    print("\n" + "=" * 76)
    print("  Condition 1: KIT BPE-LiTFSI Experimental Benchmark @ 55°C (15:1, 0.5C)")
    print("=" * 76)
    print(f"  Membrane Bulk Resistance   (R_sep) : {res_bpe['R_sep']:6.2f} Ω·cm² ({res_bpe['share_sep']*100:4.1f}%)")
    print(f"  Cathode Porous TLM         (R_cath): {res_bpe['R_cath_TLM']:6.2f} Ω·cm² ({res_bpe['share_cath']*100:4.1f}%)")
    print(f"  Concentration Polarization (R_conc): {res_bpe['R_conc']:6.2f} Ω·cm² ({res_bpe['share_conc']*100:4.1f}%)  <-- 1순위 병목 (43~44%)")
    print(f"  Interfacial Charge Transfer(R_int) : {res_bpe['R_int']:6.2f} Ω·cm² ({res_bpe['share_int']*100:4.1f}%)")
    print("  " + "-" * 72)
    print(f"  Total Area-Specific Resistance     : {res_bpe['ASR_total']:6.2f} Ω·cm² (100.0%)")
    print(f"  Total Overpotential Loss (Delta V) : {res_bpe['dV_0p5C']:6.3f} V")

    print("\n" + "=" * 76)
    print("  Condition 2: PEO-LiTFSI Design Optimum @ 60°C (18:1, 0.5C)")
    print("=" * 76)
    print(f"  Membrane Bulk Resistance   (R_sep) : {res_60C['R_sep']:6.2f} Ω·cm² ({res_60C['share_sep']*100:4.1f}%)")
    print(f"  Cathode Porous TLM         (R_cath): {res_60C['R_cath_TLM']:6.2f} Ω·cm² ({res_60C['share_cath']*100:4.1f}%)")
    print(f"  Concentration Polarization (R_conc): {res_60C['R_conc']:6.2f} Ω·cm² ({res_60C['share_conc']*100:4.1f}%)")
    print(f"  Interfacial Charge Transfer(R_int) : {res_60C['R_int']:6.2f} Ω·cm² ({res_60C['share_int']*100:4.1f}%)")
    print("  " + "-" * 72)
    print(f"  Total Area-Specific Resistance     : {res_60C['ASR_total']:6.2f} Ω·cm² (100.0%)")
    print(f"  Total Overpotential Loss (Delta V) : {res_60C['dV_0p5C']:6.3f} V")
    print("=" * 76)

    return {"bpe_55C": res_bpe, "peo_60C": res_60C}


def cmd_all(args) -> None:
    """Runs all analyses and updates results."""
    cmd_optimize(args)
    cmd_transference(args)
    cmd_report(args)


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        prog="ecm-spe",
        description="ECM-SPE: Solid Polymer Electrolyte Composition Optimization",
    )
    sub = ap.add_subparsers(dest="command")

    # optimize
    p_opt = sub.add_parser("optimize", help="Run 4-resistance composition sweep")
    p_opt.add_argument("--grid-min", type=float, default=4.0)
    p_opt.add_argument("--grid-max", type=float, default=40.0)
    p_opt.add_argument("--grid-step", type=float, default=1.0)
    p_opt.add_argument("--L-sep", type=float, default=20.0, help="Separator thickness in um")
    p_opt.add_argument("--L-cath", type=float, default=60.0, help="Cathode thickness in um")
    p_opt.add_argument("--eps", type=float, default=0.30, help="Cathode electrolyte volume fraction")
    p_opt.add_argument("--R-int", type=float, default=40.0, help="Interfacial ASR in Ohm*cm^2")
    p_opt.add_argument("--current-density", type=float, default=1.5, help="Current density in mA/cm^2")
    p_opt.add_argument("--t-plus", type=float, default=0.147, help="Transference number")
    p_opt.add_argument("--include-conc", action="store_true", help="Include concentration polarization in ranking")
    p_opt.add_argument("--out", type=str, default=None, help="Output CSV path")
    p_opt.add_argument("--save", action="store_true", help="Save results to results/optimize_demo.csv")

    # transference
    sub.add_parser("transference", help="Bruce-Vincent transference number analysis")

    # report
    sub.add_parser("report", help="4-resistance bottleneck report")

    # all
    p_all = sub.add_parser("all", help="Run full pipeline")
    p_all.add_argument("--grid-min", type=float, default=4.0)
    p_all.add_argument("--grid-max", type=float, default=40.0)
    p_all.add_argument("--grid-step", type=float, default=1.0)
    p_all.add_argument("--L-sep", type=float, default=20.0)
    p_all.add_argument("--L-cath", type=float, default=60.0)
    p_all.add_argument("--eps", type=float, default=0.30)
    p_all.add_argument("--R-int", type=float, default=40.0)
    p_all.add_argument("--current-density", type=float, default=1.5)
    p_all.add_argument("--t-plus", type=float, default=0.147)
    p_all.add_argument("--include-conc", action="store_true")

    args = ap.parse_args(argv)
    if args.command == "optimize":
        cmd_optimize(args)
    elif args.command == "transference":
        cmd_transference(args)
    elif args.command == "report":
        cmd_report(args)
    elif args.command == "all":
        cmd_all(args)
    else:
        ap.print_help()


if __name__ == "__main__":
    main()
