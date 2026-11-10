"""Tests for 4-resistance equivalent circuit model and de Levie transmission line model."""
import math
import pytest
import numpy as np

from ecm_spe.ecm import (
    CellParameters,
    cell_resistance_breakdown,
    concentration_polarization_resistance,
    de_levie_tlm,
    optimize_composition,
)
from ecm_spe.ml import peo_litfsi_empirical_conductivity


def test_concentration_polarization_limits():
    R_sep = 20.0
    # Ideal single-ion conductor (t+ = 1.0) -> zero concentration overpotential
    assert concentration_polarization_resistance(R_sep, t_plus=1.0) == 0.0

    # Symmetric binary electrolyte (t+ = 0.5) -> R_conc equals R_sep
    assert concentration_polarization_resistance(R_sep, t_plus=0.5) == pytest.approx(20.0)

    # Typical SPE (t+ = 0.1) -> 9x amplification of membrane resistance
    assert concentration_polarization_resistance(R_sep, t_plus=0.1) == pytest.approx(180.0)

    # Invalid t+
    assert math.isinf(concentration_polarization_resistance(R_sep, t_plus=0.0))


def test_de_levie_tlm_properties(standard_cell):
    sig1 = 1e-4
    sig2 = 5e-4

    r1 = de_levie_tlm(sig1, standard_cell)
    r2 = de_levie_tlm(sig2, standard_cell)

    # Higher conductivity must strictly reduce cathode ionic resistance
    assert r1 > r2 > 0.0

    # Non-physical conductivity handling
    assert math.isinf(de_levie_tlm(0.0, standard_cell))
    assert math.isinf(de_levie_tlm(-1e-4, standard_cell))


def test_de_levie_tlm_infinite_charge_transfer_limit():
    """In the limit of very high Rct (blocking limit), R_cath approaches Rct + R_ion / 3."""
    cell = CellParameters(L_cath=60e-4, eps=0.30, R_int=1e8, bruggeman_exp=1.5)
    sig = 1e-3
    sig_eff = sig * (0.30 ** 1.5)
    R_ion = 60e-4 / sig_eff

    r_tlm = de_levie_tlm(sig, cell)
    # R_cath ~ R_int + R_ion / 3
    expected = cell.R_int + R_ion / 3.0
    assert r_tlm == pytest.approx(expected, rel=1e-3)


def test_cell_resistance_budget_and_shares(standard_cell):
    sig = 5e-4
    b = cell_resistance_breakdown(sig, cell=standard_cell, include_conc=True)

    # Total resistance must equal sum of 4 parts
    expected_tot = b["R_sep"] + b["R_cath_TLM"] + b["R_conc"] + b["R_int"]
    assert b["ASR_total"] == pytest.approx(expected_tot)

    # Shares must sum to 1.0 (100%)
    tot_share = b["share_sep"] + b["share_cath"] + b["share_conc"] + b["share_int"]
    assert tot_share == pytest.approx(1.0)

    # Ohm's law: dV = i * ASR_total
    assert b["dV_0p5C"] == pytest.approx(standard_cell.i * b["ASR_total"])


def test_concentration_polarization_dominance():
    """Verifies that at BPE measured t+ = 0.147 and 55°C, concentration polarization is the primary bottleneck (40~45%)."""
    cell_55C = CellParameters(L_sep=30e-4, L_cath=60e-4, eps=0.30, R_int=116.8, t_plus=0.147, i=1.5e-3)
    sig_55C = 6.1e-5  # BPE-LiTFSI 15:1 measured conductivity at 55°C (S/cm)
    b = cell_resistance_breakdown(sig_55C, cell=cell_55C, include_conc=True)

    # R_conc / R_sep = (1/0.147 - 1) ~ 5.8
    assert b["R_conc"] / b["R_sep"] == pytest.approx(1.0 / 0.147 - 1.0, rel=1e-3)

    # Concentration polarization share should be around 40% ~ 46% (matches 44% in slides/README)
    assert 0.38 <= b["share_conc"] <= 0.48
    assert b["share_conc"] > b["share_sep"]
    assert b["share_conc"] > b["share_cath"]


def test_optimize_composition_finds_optimum(standard_cell):
    sweep = optimize_composition(
        peo_litfsi_empirical_conductivity,
        temperatures=(25.0, 60.0),
        cell=standard_cell,
        include_conc=False,
    )
    summary = sweep["summary"]
    # Optimal EO:Li ratio should be in 15~21 at 60°C
    best_60C = summary["60C_best_EO_Li"]
    assert 14.0 <= best_60C <= 22.0
    assert summary["60C_best_dV"] > 0.0
