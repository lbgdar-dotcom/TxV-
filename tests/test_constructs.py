import pytest

from txv.constructs import ConstructSpec, build_construct
from txv.epitopes import Antigen, Cassette
from txv.parts import default_registry
from txv.seqops import translate


def test_feature_coordinates_tile_the_top_level_template(construct):
    top = [f for f in construct.features
           if f.kind in ("promoter", "utr5", "kozak", "orf", "utr3", "polya",
                         "linearization")]
    top.sort(key=lambda f: f.start)
    assert top[0].start == 0
    assert top[-1].end == len(construct.template)
    for previous, current in zip(top, top[1:]):
        assert previous.end == current.start, "top-level features do not tile"


def test_orf_subfeatures_lie_inside_the_orf(construct):
    orf = construct.feature("ORF")
    for feature in construct.features:
        if feature.kind in ("signal_peptide", "linker", "trafficking", "stop",
                            "neoepitope", "tumor_associated", "helper", "start"):
            assert orf.start <= feature.start < feature.end <= orf.end


def test_orf_translation_matches_the_declared_protein(construct):
    assert translate(construct.orf, stop_at_stop=True).rstrip("*") == construct.orf_protein
    assert construct.orf.startswith("ATG")
    assert len(construct.orf) % 3 == 0


def test_every_antigen_is_present_in_the_protein(construct, antigens):
    for antigen in antigens:
        assert antigen.sequence in construct.orf_protein


def test_transcript_starts_at_the_t7_plus_one_g(construct):
    promoter = construct.feature("T7_promoter")
    assert construct.transcript_start == promoter.end - 1
    assert construct.transcript.startswith("G")
    # The linearisation site sits outside the transcript.
    assert construct.transcript_end <= construct.feature("BspQI_site").start


def test_transcript_ends_with_the_polya_tail(construct):
    assert construct.transcript.endswith("A" * 50)


def test_tandem_stop_is_added_and_can_be_switched_off(antigens):
    with_tandem = build_construct("T1", antigens)
    spec = ConstructSpec(name="T2", tandem_stop=False)
    without = build_construct("T2", antigens, spec=spec)
    assert len(with_tandem.feature("stop")) == 6
    assert len(without.feature("stop")) == 3


def test_signal_peptide_and_trafficking_can_be_omitted(antigens):
    spec = ConstructSpec(name="bare", signal_peptide=None, trafficking=None)
    construct = build_construct("bare", antigens, spec=spec)
    assert construct.orf_protein.startswith("M")
    assert not construct.features_of_kind("signal_peptide")
    assert not construct.features_of_kind("trafficking")


def test_initiator_methionine_is_added_when_nothing_supplies_one(antigens):
    spec = ConstructSpec(name="noM", signal_peptide=None, trafficking=None)
    construct = build_construct("noM", [Antigen("x", "ACDEFGH")], spec=spec)
    assert construct.orf_protein.startswith("M")
    assert any("initiator methionine" in n for n in construct.design_notes)


def test_default_construct_contains_no_placeholders(construct):
    """The shipped defaults are real sequences, so a default build is clean."""
    assert construct.placeholders_used == []
    assert not any("Placeholder parts in use" in n for n in construct.design_notes)


def test_placeholders_are_reported_when_explicitly_requested(antigens):
    spec = ConstructSpec(name="ph", utr5="UTR5_placeholder")
    construct = build_construct("ph", antigens, spec=spec)
    assert "UTR5_placeholder" in construct.placeholders_used
    assert any("Placeholder parts in use" in n for n in construct.design_notes)


def test_default_utrs_are_the_validated_ones(construct):
    assert construct.feature("UTR5_hAg")
    assert construct.feature("UTR3_AES_mtRNR1")
    # 5' UTR + Kozak must reconstitute the validated 54-nt leader exactly.
    leader = construct.template[
        construct.feature("UTR5_hAg").start : construct.feature("ORF").start
    ]
    assert leader == "GAGAATAAACTAGTATTCTTCTGGTCCCCACAGACTCAGAGAGAACCCGCCACC"


def test_default_polya_is_120_nt(construct):
    assert len(construct.feature("polyA_120")) == 120


def test_empty_cassette_is_rejected():
    from txv.constructs import ConstructBuilder

    with pytest.raises(ValueError):
        ConstructBuilder().build(ConstructSpec(name="x"), Cassette([]))


def test_bad_stop_codon_is_rejected(antigens):
    with pytest.raises(ValueError):
        build_construct("x", antigens, spec=ConstructSpec(name="x", stop_codon="ATG"))


def test_unknown_part_name_is_rejected(antigens):
    with pytest.raises(KeyError):
        build_construct("x", antigens, spec=ConstructSpec(name="x", utr5="nope"))


def test_summary_is_self_consistent(construct):
    summary = construct.summary()
    assert summary["orf_nt"] == len(construct.orf)
    assert summary["transcript_nt"] == construct.transcript_length
    assert summary["protein_aa"] == len(construct.orf_protein)
    assert 0 < summary["gc"] < 1


def test_registry_is_honoured(antigens):
    registry = default_registry()
    construct = build_construct("x", antigens, registry=registry)
    assert construct.feature("SP_tPA")
