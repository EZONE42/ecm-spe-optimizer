"""Pytest configuration and shared test fixtures for ECM-SPE."""
import pytest
import numpy as np

from ecm_spe.ecm import CellParameters, DEFAULT_CELL


@pytest.fixture
def standard_cell() -> CellParameters:
    return CellParameters(
        L_sep=20e-4,
        L_cath=60e-4,
        eps=0.30,
        R_int=40.0,
        i=1.5e-3,
        bruggeman_exp=1.5,
        t_plus=0.147,
    )


@pytest.fixture
def synthetic_nyquist_semicircle():
    """Generates a perfect noisy semicircle for testing Kasa circle fitting."""
    xc_true = 150.0  # Center real
    yc_true = -30.0  # Depressed circle center
    r_true = 80.0    # Radius
    theta = np.linspace(0.1, np.pi - 0.1, 40)
    x = xc_true + r_true * np.cos(theta)
    y = yc_true + r_true * np.sin(theta)
    # Add minimal reproducible jitter
    rng = np.random.default_rng(42)
    noise = rng.normal(0, 0.2, len(x))
    return x + noise, y + noise, xc_true, yc_true, r_true
