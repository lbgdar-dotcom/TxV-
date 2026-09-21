import pytest

from txv.epitopes import (
    AnchorMotifScorer,
    Antigen,
    Cassette,
    choose_linker,
    dedupe_antigens,
    junction_risk,
    optimize_order,
    scan_junctions,
)


def test_antigen_validation():
    with pytest.raises(ValueError):
        Antigen("bad", "ACDEFX")
    with pytest.raises(ValueError):
        Antigen("empty", "")
    with pytest.raises(ValueError):
        Antigen("off", "ACDEF", mutation_offset=99)
    assert Antigen("ok", " acdef ").sequence == "ACDEF"


def test_cassette_layout_is_contiguous_and_matches_the_sequence(antigens):
    cassette = Cassette(antigens)
    sequence = cassette.sequence()
    spans = cassette.layout()
    assert spans[0][1] == 0
    assert spans[-1][2] == len(sequence)
    for previous, current in zip(spans, spans[1:]):
        assert previous[2] == current[1], "layout has a gap or overlap"
    for name, start, end, _ in spans:
        if not name.startswith("linker_"):
            match = next(a for a in antigens if a.name == name)
            assert sequence[start:end] == match.sequence


def test_per_junction_linker_override():
    cassette = Cassette(
        [Antigen("a", "AAAA"), Antigen("b", "CCCC"), Antigen("c", "DDDD")],
        linker="GG",
        linker_overrides={1: "PPPP"},
    )
    assert cassette.sequence() == "AAAAGGCCCCPPPPDDDD"
    assert cassette.linker_at(0) == "GG"
    assert cassette.linker_at(1) == "PPPP"


def test_scan_ignores_peptides_contained_in_one_bead():
    # NYESO1's own epitope SLLMWITQC must never be reported as junctional.
    cassette = Cassette(
        [Antigen("x", "AAAAAAAAAAAA"), Antigen("nyeso", "SLLMWITQC")], linker="GG"
    )
    peptides = {h.peptide for h in scan_junctions(cassette, threshold=0.0)}
    assert "SLLMWITQC" not in peptides
    # And every reported peptide really does cross the boundary.
    junction = cassette.junctions()[0]
    for hit in scan_junctions(cassette, threshold=0.0):
        assert hit.peptide not in junction.left
        assert hit.peptide not in junction.right


def test_scan_ignores_peptides_wholly_inside_the_linker():
    cassette = Cassette(
        [Antigen("a", "AAAA"), Antigen("b", "CCCC")], linker="LLLLLLLLLLLLLLLL"
    )
    for hit in scan_junctions(cassette, threshold=0.0):
        assert "A" in hit.peptide or "C" in hit.peptide


def test_anchor_scorer_is_bounded_and_length_gated():
    scorer = AnchorMotifScorer()
    assert scorer("SLLMWITQV") > 0.0
    assert 0.0 <= scorer("SLLMWITQV") <= 1.0
    assert scorer("SL") == 0.0          # too short
    assert scorer("A" * 30) == 0.0      # outside class I lengths
    with pytest.raises(ValueError):
        AnchorMotifScorer(supertypes=["A99"])


def test_long_flexible_linker_beats_short_linkers(antigens):
    cassette = Cassette(antigens)
    chosen, table = choose_linker(cassette, ["GGSGGGGSGG", "AAY", "RAKR"])
    assert chosen == "GGSGGGGSGG"
    assert table["GGSGGGGSGG"] < table["AAY"]


def test_reordering_never_increases_risk(antigens):
    cassette = Cassette(antigens, linker="AAY")
    reordered, before, after = optimize_order(cassette)
    assert after <= before
    assert sorted(a.name for a in reordered.antigens) == sorted(
        a.name for a in antigens
    )
    assert junction_risk(reordered) == pytest.approx(after)


def test_reordering_respects_pinned_ends(antigens):
    cassette = Cassette(antigens, linker="AAY")
    reordered, _, _ = optimize_order(cassette, fixed_first=True, fixed_last=True)
    assert reordered.antigens[0].name == antigens[0].name
    assert reordered.antigens[-1].name == antigens[-1].name


def test_reordering_is_a_noop_below_three_beads():
    cassette = Cassette([Antigen("a", "AAAA"), Antigen("b", "CCCC")])
    reordered, before, after = optimize_order(cassette)
    assert before == after
    assert [a.name for a in reordered.antigens] == ["a", "b"]


def test_dedupe_drops_duplicates_and_contained_sequences():
    antigens = [
        Antigen("long", "AAAASLLMWITQCDDDD"),
        Antigen("short", "SLLMWITQC"),
        Antigen("exact", "AAAASLLMWITQCDDDD"),
        Antigen("other", "KVAELVHFL"),
    ]
    kept, dropped = dedupe_antigens(antigens)
    names = [a.name for a in kept]
    assert "long" in names and "other" in names
    assert "short" not in names and "exact" not in names
    assert len(dropped) == 2
    # Input order of survivors is preserved.
    assert names == ["long", "other"]


def test_dedupe_is_stable_when_nothing_is_redundant(antigens):
    kept, dropped = dedupe_antigens(antigens)
    assert dropped == []
    assert [a.name for a in kept] == [a.name for a in antigens]
