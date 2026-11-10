"""CLI smoke tests for ecm_spe commands."""
import sys
import pytest

from ecm_spe.cli import main, cmd_optimize, cmd_transference, cmd_report


def test_cli_help(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    captured = capsys.readouterr()
    assert "ECM-SPE" in captured.out


def test_cli_transference(capsys):
    class Args:
        pass

    out = cmd_transference(Args())
    assert "t_plus" in out
    assert 0.10 <= out["t_plus"] <= 0.20
    captured = capsys.readouterr()
    assert "Bruce-Vincent Transference Number" in captured.out


def test_cli_report(capsys):
    class Args:
        pass

    res = cmd_report(Args())
    assert "bpe_55C" in res
    assert "peo_60C" in res
    assert 0.38 <= res["bpe_55C"]["share_conc"] <= 0.46  # ~40-44% concentration polarization bottleneck
    captured = capsys.readouterr()
    assert "Condition 1: KIT BPE-LiTFSI" in captured.out


def test_cli_optimize_short_grid(capsys):
    class Args:
        grid_min = 15.0
        grid_max = 25.0
        grid_step = 5.0
        L_sep = 20.0
        L_cath = 60.0
        eps = 0.30
        R_int = 40.0
        current_density = 1.5
        t_plus = 0.147
        include_conc = True

    res = cmd_optimize(Args())
    assert "summary" in res
    captured = capsys.readouterr()
    assert "ECM-SPE Composition Optimization Results" in captured.out
