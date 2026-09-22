import csv
import io

import pytest

from txv.order import (
    PVAX1_AG_LEFT_OVERLAP,
    PVAX1_AG_RIGHT_OVERLAP,
    FragmentMode,
    PolyAMode,
    assess_synthesis,
    order_fasta,
    order_table,
    synthesis_fragment,
)
from txv.pvax1_ag import build_panel_construct


@pytest.fixture
def panel_construct():
    return build_panel_construct("P3")


def test_orf_mode_orders_only_the_coding_sequence(panel_construct):
    fragment = synthesis_fragment(panel_construct, FragmentMode.ORF)
    assert fragment.insert == panel_construct.orf
    assert fragment.sequence.startswith(PVAX1_AG_LEFT_OVERLAP)
    assert fragment.sequence.endswith(PVAX1_AG_RIGHT_OVERLAP)
    assert fragment.length == len(panel_construct.orf) + 40


def test_orf_fragment_begins_at_atg_without_repeating_the_kozak(panel_construct):
    fragment = synthesis_fragment(panel_construct, FragmentMode.ORF)
    assert fragment.insert.startswith("ATG")
    # The Kozak appears once, at the end of the left overlap.
    assert fragment.sequence.count("GCCACCATG") == 1


def test_cassette_mode_includes_the_utrs(panel_construct):
    fragment = synthesis_fragment(
        panel_construct, FragmentMode.CASSETTE, PolyAMode.PCR_ADDED
    )
    utr5 = panel_construct.features_of_kind("utr5")[0]
    utr3 = panel_construct.features_of_kind("utr3")[0]
    assert fragment.insert == panel_construct.template[utr5.start : utr3.end]
    assert panel_construct.orf in fragment.sequence


def test_encoded_polya_is_included_only_in_cassette_mode(panel_construct):
    encoded = synthesis_fragment(
        panel_construct, FragmentMode.CASSETTE, PolyAMode.ENCODED
    )
    pcr = synthesis_fragment(
        panel_construct, FragmentMode.CASSETTE, PolyAMode.PCR_ADDED
    )
    assert "A" * 100 in encoded.sequence
    assert "A" * 100 not in pcr.sequence
    assert encoded.length == pcr.length + 120


def test_encoded_polya_is_flagged_as_a_synthesis_risk(panel_construct):
    encoded = synthesis_fragment(
        panel_construct, FragmentMode.CASSETTE, PolyAMode.ENCODED
    )
    assert not encoded.risk.ok
    assert any("A-tract" in flag for flag in encoded.risk.flags)
    # Dropping it clears the flag, which is the whole point of the mode.
    assert synthesis_fragment(
        panel_construct, FragmentMode.CASSETTE, PolyAMode.PCR_ADDED
    ).risk.ok


def test_cassette_mode_warns_that_it_moves_the_primer_rails(panel_construct):
    fragment = synthesis_fragment(panel_construct, FragmentMode.CASSETTE)
    assert any("primer" in note for note in fragment.notes)


def test_orf_mode_says_the_backbone_utrs_are_what_is_used(panel_construct):
    fragment = synthesis_fragment(panel_construct, FragmentMode.ORF)
    assert any("backbone" in note for note in fragment.notes)


def test_polya_encoded_is_a_noop_in_orf_mode(panel_construct):
    encoded = synthesis_fragment(panel_construct, FragmentMode.ORF, PolyAMode.ENCODED)
    pcr = synthesis_fragment(panel_construct, FragmentMode.ORF, PolyAMode.PCR_ADDED)
    assert encoded.sequence == pcr.sequence
    assert any("no effect" in note for note in encoded.notes)


def test_custom_overlaps_are_honoured(panel_construct):
    fragment = synthesis_fragment(
        panel_construct, FragmentMode.ORF,
        left_overlap="AAAAAA", right_overlap="TTTTTT",
    )
    assert fragment.sequence.startswith("AAAAAA")
    assert fragment.sequence.endswith("TTTTTT")
    assert fragment.insert == panel_construct.orf


def _nonrepetitive(n: int, seed: int = 7) -> str:
    """A pseudo-random sequence: balanced GC, no long runs, no long repeats."""
    import random

    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_assess_synthesis_passes_a_clean_fragment():
    assert assess_synthesis(_nonrepetitive(600)).ok


def test_assess_synthesis_flags_a_long_homopolymer():
    flags = assess_synthesis(_nonrepetitive(200) + "A" * 30).flags
    assert any("A-tract" in f for f in flags)


def test_assess_synthesis_flags_extreme_gc():
    assert any("GC" in f for f in assess_synthesis("AT" * 200).flags)
    assert any("GC" in f for f in assess_synthesis("GC" * 200).flags)


def test_assess_synthesis_flags_tandem_repeats():
    """A tandem repeat is exactly what breaks synthesis and assembly."""
    assert any("repeat" in f for f in assess_synthesis("ATGC" * 50).flags)


def test_assess_synthesis_flags_oversize_fragments():
    assert any("exceeds" in f for f in assess_synthesis(_nonrepetitive(3600)).flags)


def test_order_table_is_valid_csv(panel_construct):
    fragments = [synthesis_fragment(panel_construct, m) for m in FragmentMode]
    rows = list(csv.DictReader(io.StringIO(order_table(fragments))))
    assert len(rows) == 2
    assert rows[0]["name"] == "P3"
    assert int(rows[0]["length_bp"]) == len(fragments[0].sequence)
    assert rows[0]["sequence"] == fragments[0].sequence


def test_order_fasta_round_trips(panel_construct):
    fragment = synthesis_fragment(panel_construct, FragmentMode.ORF)
    text = order_fasta([fragment])
    body = "".join(text.splitlines()[1:])
    assert body == fragment.sequence
    assert text.startswith(">P3|orf|")


def test_every_panel_member_orders_cleanly_in_the_default_mode():
    for name in [f"P{i}" for i in range(9)]:
        fragment = synthesis_fragment(build_panel_construct(name))
        assert fragment.risk.ok, (name, fragment.risk.flags)
