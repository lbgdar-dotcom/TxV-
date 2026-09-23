"""Negative controls for the in-silico validation run.

`scripts/validate_in_silico.py` reports "no findings" on the shipped design.
That is only worth anything if the run can report findings, so each threshold
is moved past the panel here and the run must fail on it.
"""

import importlib.util
import sys
from pathlib import Path

import pytest

pytest.importorskip("RNA", reason="ViennaRNA not installed")

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "validate_in_silico.py"
pytestmark = pytest.mark.skipif(
    not (ROOT / "orders" / "pjet12_ivt_units" / "gblocks.fasta").exists(),
    reason="pJET package not generated",
)


@pytest.fixture(scope="module")
def validator():
    spec = importlib.util.spec_from_file_location("validate_in_silico", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["validate_in_silico"] = module
    spec.loader.exec_module(module)
    return module


def test_the_shipped_design_passes(validator, capsys):
    assert validator.main() == 0
    assert "no findings" in capsys.readouterr().out


@pytest.mark.parametrize("threshold,value,expected", [
    ("MIN_START_ACCESSIBILITY", 0.95, "accessibility"),
    ("MIN_CAI", 0.99, "CAI"),
    ("MAX_CPG_OE", 0.01, "CpG"),
    ("MIN_CAP_MFE", -1.0, "cap-proximal"),
    ("MAX_MFE_PER_NT_SPREAD", 0.01, "MFE/nt varies"),
])
def test_each_threshold_can_fail(validator, monkeypatch, capsys,
                                 threshold, value, expected):
    monkeypatch.setattr(validator, threshold, value)
    assert validator.main() == 1
    out = capsys.readouterr().out
    assert "FINDING(S)" in out
    assert expected in out


def test_a_failing_signal_peptide_call_is_surfaced(validator, monkeypatch, capsys):
    """Does the script act on the verdict, or just print it?

    The peptide itself cannot be swapped out here: the design library validates
    pinned DNA against its protein and refuses an inconsistent module, which is
    correct of it. So the wiring is what gets tested -- given a failing call,
    the run must report and exit non-zero. Whether the call is itself correct
    is covered by tests/test_insilico.py.
    """
    from txv.insilico import SignalPeptideCall

    broken = SignalPeptideCall(
        length=24, n_region_charge=-3.0, h_region_start=0, h_region_length=0,
        h_region_max=-2.0, minus_one="W", minus_three="W",
        problems=["no hydrophobic h-region (max window KD -2.00)"])
    monkeypatch.setattr(validator, "signal_peptide_call", lambda _: broken)
    assert validator.main() == 1
    out = capsys.readouterr().out
    assert "h-region" in out
    assert "FINDING(S)" in out


def test_a_tm_module_with_no_membrane_span_is_surfaced(validator, monkeypatch,
                                                       capsys):
    monkeypatch.setattr(validator, "tm_segments", lambda _: [])
    assert validator.main() == 1
    assert "cannot span a bilayer" in capsys.readouterr().out


def test_a_construct_whose_lamp1_tail_is_not_terminal_is_caught(
        validator, monkeypatch, capsys):
    """GYQTI binds AP-3 only as the last five residues; internal is silent."""
    from txv import pvax1_ag
    real = pvax1_ag.build_panel_construct

    def shifted(name, *a, **k):
        construct = real(name, *a, **k)
        if name == "P6":
            object.__setattr__(construct, "orf_protein",
                               construct.orf_protein + "AAA")
        return construct

    monkeypatch.setattr(validator, "build_panel_construct", shifted)
    assert validator.main() == 1
    assert "GYQTI is not the C-terminus" in capsys.readouterr().out
