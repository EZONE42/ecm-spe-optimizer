"""Regression tests verifying ground truth quantitative results against results/*.csv."""
from pathlib import Path
import pytest
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "results"


def test_regression_li_tplus_ground_truth():
    p = RESULTS_DIR / "li_tplus.csv"
    assert p.exists(), "results/li_tplus.csv missing"
    df = pd.read_csv(p)
    d = dict(zip(df["quantity"], df["value"]))

    # Transference numbers
    assert float(d["t_plus_peakI0"]) == pytest.approx(0.147, abs=0.005)
    assert float(d["t_plus_calcI0"]) == pytest.approx(0.111, abs=0.005)
    assert float(d["t_plus_firstI0"]) == pytest.approx(0.172, abs=0.005)

    # Interfacial ASR per side: Ri0 * 0.5026 / 2
    ri0 = float(d["R_i0"])
    asr = float(d["ASR_int_per_side"])
    assert asr == pytest.approx(ri0 * 0.5026 / 2.0, rel=1e-2)


def test_regression_transfer_results_extrapolation_bounds():
    p = RESULTS_DIR / "transfer_results.csv"
    assert p.exists(), "results/transfer_results.csv missing"
    df = pd.read_csv(p)
    assert len(df) == 6

    # Verify salt replacement cases (A:C) have modest fold error (~3.7 - 4.6x)
    salt_cases = df[df["case"].str.contains("A:C")]
    for _, r in salt_cases.iterrows():
        assert 3.0 <= float(r["fold_err"]) <= 5.5

    # Verify novel polymer replacement cases (PEC, PTMC) exhibit severe extrapolation failure (> 20x)
    poly_fail_cases = df[df["case"].str.contains("PEC|PTMC")]
    for _, r in poly_fail_cases.iterrows():
        assert float(r["fold_err"]) > 20.0
        assert float(r["MAE_log"]) > 1.20


def test_regression_optimize_demo_sweep():
    import io
    p = RESULTS_DIR / "optimize_demo.csv"
    assert p.exists(), "results/optimize_demo.csv missing"
    with open(p, "r", encoding="utf-8", errors="ignore") as f:
        content = f.read()
    df = pd.read_csv(io.StringIO(content))
    assert "EO_Li" in df.columns
    assert "T_C" in df.columns
    assert "dV_0p5C" in df.columns

    # 60°C best composition should be around EO:Li = 15~21
    df_60 = df[df["T_C"] == 60]
    best_row = df_60.loc[df_60["dV_0p5C"].idxmin()]
    assert 14.0 <= float(best_row["EO_Li"]) <= 22.0
