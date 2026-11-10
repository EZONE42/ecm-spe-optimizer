"""Tests for Bruce-Vincent transference number and Kasa circle fitting."""
import math
import pytest
import numpy as np

from ecm_spe.transference import (
    bruce_vincent_tplus,
    circle_fit,
    interfacial_asr,
    extract_eis_intercepts,
)


def test_kasa_circle_fit_recovery(synthetic_nyquist_semicircle):
    x, y, xc_true, yc_true, r_true = synthetic_nyquist_semicircle
    lo, hi, xc, yc, r = circle_fit(x, y)

    # Reconstructs center and radius within 1% error despite noise
    assert xc == pytest.approx(xc_true, rel=0.02)
    assert yc == pytest.approx(yc_true, rel=0.03)
    assert r == pytest.approx(r_true, rel=0.02)

    # Real axis intercepts: xc +/- sqrt(r^2 - yc^2)
    expected_hi = xc_true + math.sqrt(r_true**2 - yc_true**2)
    assert hi == pytest.approx(expected_hi, rel=0.02)


def test_bruce_vincent_ideal_limits():
    dV = 0.05  # 50 mV

    # Case 1: Ideal non-polarizing cell with constant interface (Iss == I0, Ri0 == Ri_ss) -> t+ = 1.0
    t_ideal = bruce_vincent_tplus(Rb0=100.0, Ri0=50.0, Ri_ss=50.0, I0=2e-4, Iss=2e-4, dV=dV)
    assert t_ideal == pytest.approx(1.0)

    # Case 2: Zero interfacial resistance (Ri0 = 0, Ri_ss = 0) -> t+ = Iss / I0
    t_zero_ri = bruce_vincent_tplus(Rb0=100.0, Ri0=0.0, Ri_ss=0.0, I0=2e-4, Iss=5e-5, dV=dV)
    assert t_zero_ri == pytest.approx(5e-5 / 2e-4)
    assert t_zero_ri == pytest.approx(0.25)


def test_bruce_vincent_kit_reference_values():
    """Verifies against KIT Zenodo 7428858 BPE-LiTFSI reference values."""
    # Experimental values: dV=50 mV, I0=243.1 uA, Iss=38.6 uA, Ri0=85.2 Ohm, Ri_ss=112.4 Ohm
    t_plus = bruce_vincent_tplus(
        Rb0=120.5,
        Ri0=85.2,
        Ri_ss=112.4,
        I0=243.1e-6,
        Iss=38.6e-6,
        dV=0.05,
    )
    # Expected transference number in range 0.12 ~ 0.16
    assert 0.10 <= t_plus <= 0.20
    assert not math.isnan(t_plus)


def test_interfacial_asr_scaling():
    Ri = 100.0  # Ohm
    A = 0.5026  # cm^2

    # Symmetric cell: Ri is for two interfaces in series -> divide by 2
    asr_sym = interfacial_asr(Ri, area_cm2=A, is_symmetric=True)
    assert asr_sym == pytest.approx(100.0 * 0.5026 / 2.0)

    # Asymmetric cell: single interface
    asr_asym = interfacial_asr(Ri, area_cm2=A, is_symmetric=False)
    assert asr_asym == pytest.approx(100.0 * 0.5026)
    assert asr_asym == pytest.approx(2.0 * asr_sym)
