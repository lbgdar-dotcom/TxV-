"""Self-contained IVT units for blunt cloning into pJET1.2."""

import pytest

from txv.ivt_unit import (
    FWD_HANDLE, FWD_PRIMER, LEADER, REV_HANDLE, T7_CORE, UTR3,
    make_unit, reverse_primer, screen_handles, simulate_pcr,
    verify_orientation_independence, verify_unit,
)
from txv.pvax1_ag import PANEL, build_panel_construct
from txv.seqops import find_all, revcomp, translate

NAMES = [p.name for p in PANEL]


@pytest.fixture(scope="module")
def units():
    return {n: make_unit(n, build_panel_construct(n).orf) for n in NAMES}


@pytest.mark.parametrize("name", NAMES)
def test_unit_verifies(units, name):
    assert verify_unit(units[name], build_panel_construct(name).orf_protein) == []


@pytest.mark.parametrize("name", NAMES)
def test_orientation_does_not_matter(units, name):
    """The whole point: blunt ligation cannot control orientation."""
    assert verify_orientation_independence(units[name]) == []


@pytest.mark.parametrize("name", NAMES)
def test_unit_is_a_complete_transcription_unit(units, name):
    s = units[name].sequence
    for element in (T7_CORE, LEADER, UTR3):
        assert s.count(element) == 1
    assert s.index(T7_CORE) < s.index(LEADER) < s.index(UTR3)


@pytest.mark.parametrize("name", NAMES)
def test_transcript_starts_AGG_and_carries_the_whole_message(units, name):
    unit = units[name]
    assert unit.transcript.startswith("AGG")
    assert unit.orf in unit.transcript
    assert unit.transcript.endswith(UTR3)
    # Nothing from the primer handles reaches the RNA.
    assert FWD_HANDLE not in unit.transcript
    assert REV_HANDLE not in unit.transcript


@pytest.mark.parametrize("name", NAMES)
def test_translation_from_the_kozak_gives_the_designed_protein(units, name):
    unit = units[name]
    start = unit.sequence.index("GCCACCATG") + 6
    protein = translate(unit.sequence[start:], stop_at_stop=True).rstrip("*")
    assert protein == build_panel_construct(name).orf_protein


def test_no_upstream_aug_before_the_start(units):
    for unit in units.values():
        leader = unit.sequence[unit.sequence.index(T7_CORE) + len(T7_CORE):
                               unit.sequence.index("GCCACCATG") + 6]
        assert "ATG" not in leader


# --- the primers -----------------------------------------------------------

def test_handles_are_sane(units):
    body = "".join(u.orf for u in units.values()) + LEADER + UTR3
    assert screen_handles(context=body) == []


def test_reverse_primer_is_the_complement_of_its_handle():
    """The classic slip: using the handle itself as the primer gives no product."""
    assert reverse_primer(0) == revcomp(REV_HANDLE)
    assert reverse_primer(0) != REV_HANDLE


def test_reverse_primer_carries_the_tail_at_its_five_prime_end():
    primer = reverse_primer(120)
    assert len(primer) == 140
    assert primer.startswith("T" * 120)
    assert primer.endswith(revcomp(REV_HANDLE))


def test_pcr_simulation_finds_one_product_on_a_circular_template(units):
    unit = units["P3"]
    vector = "ACGTTGCA" * 40
    plasmid = vector[:80] + unit.sequence + vector[80:]
    products = simulate_pcr(plasmid, FWD_PRIMER, reverse_primer(0))
    assert len(products) == 1
    assert products[0] == unit.sequence


def test_pcr_simulation_finds_nothing_without_the_handles():
    assert simulate_pcr("ACGT" * 200, FWD_PRIMER, reverse_primer(0)) == []


def test_handles_are_absent_from_every_orf(units):
    for name, unit in units.items():
        assert find_all(unit.orf, FWD_HANDLE, both_strands=True) == []
        assert find_all(unit.orf, REV_HANDLE, both_strands=True) == []


# --- fragment properties ---------------------------------------------------

def test_fragments_are_orderable(units):
    for unit in units.values():
        assert 125 <= len(unit) <= 3000


def test_no_polya_is_encoded(units):
    """It comes from the reverse primer; a 120mer A-tract is unsynthesisable."""
    for unit in units.values():
        assert "A" * 20 not in unit.sequence


def test_verifier_catches_a_broken_kozak():
    unit = make_unit("bad", "ATGGCCTGA")
    unit.sequence = unit.sequence.replace("GCCACCATG", "GCCAAAATG")
    assert verify_unit(unit)


def test_verifier_catches_a_lost_t7_promoter():
    unit = make_unit("bad", "ATGGCCTGA")
    unit.sequence = unit.sequence.replace(T7_CORE, "A" * len(T7_CORE))
    assert any("T7" in p for p in verify_unit(unit))
