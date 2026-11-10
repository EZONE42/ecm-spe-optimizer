"""Tests for RDKit molecular descriptor extraction and similarity."""
import math
import pytest

from ecm_spe.features import poly_desc, salt_desc, tanimoto_similarity


def test_poly_descriptors_peo():
    # PEO: [*]OCC[*] -> heavy=1 (3 heavy - 2 dummy atoms), nO=1, ether=1, carbonate=0
    d = poly_desc("[Cu]OCC[Au]")
    assert len(d) == 7
    heavy, Ofrac, etherfrac, carb, Nfrac, tpsa_per, rot = d

    assert heavy == pytest.approx(1.0)
    assert Ofrac == pytest.approx(1.0)
    assert etherfrac == pytest.approx(1.0)
    assert carb == pytest.approx(0.0)
    assert Nfrac == pytest.approx(0.0)
    assert tpsa_per == pytest.approx(9.23)
    assert rot == pytest.approx(2.0)


def test_poly_descriptors_carbonate():
    # PTMC: [*]OC(=O)OCCC[*] -> carbonate count must be 1
    d_ptmc = poly_desc("[Cu]OC(=O)OCCC[Au]")
    assert d_ptmc[3] == pytest.approx(1.0)  # carbonate

    # PEC: [*]OC(=O)OCC[*] -> carbonate count must be 1
    d_pec = poly_desc("[Cu]OC(=O)OCC[Au]")
    assert d_pec[3] == pytest.approx(1.0)  # carbonate


def test_salt_descriptors_litfsi():
    # LiTFSI: TFSI- anion (CF3SO2NSO2CF3)-
    # MW ~ 280.15 g/mol, 6 Fluorines, 15 heavy atoms, is_imide = 1
    d = salt_desc("O=S(=O)([N-]S(=O)(=O)C(F)(F)F)C(F)(F)F.[Li+]")
    assert len(d) == 5
    mw, nF, heavy, is_imide, tpsa = d

    assert mw == pytest.approx(280.15, abs=0.5)
    assert nF == pytest.approx(6.0)
    assert heavy == pytest.approx(15.0)
    assert is_imide == pytest.approx(1.0)


def test_salt_descriptors_litfo():
    # LiTfO: TfO- anion (CF3SO3)-
    # MW ~ 149.09 g/mol, 3 Fluorines, 8 heavy atoms, is_imide = 0
    d = salt_desc("O=S(=O)([O-])C(F)(F)F.[Li+]")
    mw, nF, heavy, is_imide, tpsa = d

    assert mw == pytest.approx(149.09, abs=0.5)
    assert nF == pytest.approx(3.0)
    assert heavy == pytest.approx(8.0)
    assert is_imide == pytest.approx(0.0)


def test_polymer_tanimoto_similarity():
    smi_peo = "[Cu]OCC[Au]"
    # Self similarity must be 1.0
    assert tanimoto_similarity(smi_peo, smi_peo) == pytest.approx(1.0)

    # Similarity to PEC / PTMC should be in range 0.2 ~ 0.5
    sim_pec = tanimoto_similarity(smi_peo, "[Cu]OC(=O)OCC[Au]")
    sim_ptmc = tanimoto_similarity(smi_peo, "[Cu]OC(=O)OCCC[Au]")
    assert 0.15 <= sim_pec <= 0.60
    assert 0.15 <= sim_ptmc <= 0.60
