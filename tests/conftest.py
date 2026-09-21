import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from txv.epitopes import Antigen  # noqa: E402


@pytest.fixture
def antigens():
    return [
        Antigen("KRAS_G12D", "LVVVGADGVGKSALTIQLIQNHFVD", mutation_offset=6),
        Antigen("TP53_R175H", "SSCMGGMNRRPILTIITLEDSSGNL", mutation_offset=12),
        Antigen("NYESO1", "SLLMWITQC", kind="tumor_associated"),
        Antigen("MAGEA3", "KVAELVHFL", kind="tumor_associated"),
        Antigen("PADRE", "AKFVAAWTLKAAA", kind="helper"),
    ]


@pytest.fixture
def construct(antigens):
    from txv.constructs import build_construct

    return build_construct("TEST-001", antigens)
