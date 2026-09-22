"""Tests against the real pVax1_AG_eGFP vector map."""

from pathlib import Path

import pytest

from txv.genbank_io import parse_genbank, read_genbank, write_genbank
from txv.plasmid import build_plasmid, verify_plasmid
from txv.pvax1_ag import PANEL, build_panel_construct
from txv.seqops import find_all, translate

BACKBONE = Path(__file__).resolve().parents[1] / "backbone" / "pvax1_ag_egfp.gb"
pytestmark = pytest.mark.skipif(
    not BACKBONE.exists(), reason="backbone vector map not present"
)


@pytest.fixture(scope="module")
def backbone():
    return read_genbank(BACKBONE)


# --- the parser and the vector itself -------------------------------------

def test_backbone_parses_to_the_expected_vector(backbone):
    assert backbone.name == "pVax1_AG_eGFP"
    assert len(backbone.sequence) == 3855
    assert backbone.is_circular


def test_genbank_round_trips(backbone):
    again = parse_genbank(write_genbank(backbone))
    assert again.sequence == backbone.sequence
    assert len(again.features) == len(backbone.features)
    assert {f.label for f in again.features} == {f.label for f in backbone.features}


def test_egfp_cds_includes_its_own_stop(backbone):
    egfp = backbone.find("eGFP")
    assert len(egfp) == 720
    protein = translate(egfp.slice(backbone.sequence), stop_at_stop=False)
    assert protein.startswith("MVSKGEEL")
    assert protein.endswith("*")
    assert "*" not in protein[:-1]


def test_homology_rails_bracket_the_egfp_cds_exactly(backbone):
    """The rails must abut the CDS, or the insert lands out of frame."""
    left, right = "GAAGAAATATAAGAGCCACC", "GCTGCCTTCTGCGGGGCTTG"
    lh, rh = find_all(backbone.sequence, left), find_all(backbone.sequence, right)
    assert len(lh) == 1 and len(rh) == 1
    egfp = backbone.find("eGFP")
    assert lh[0] + len(left) == egfp.start
    assert rh[0] == egfp.end


def test_transcript_starts_AGG_for_cleancap(backbone):
    core = "TAATACGACTCACTATA"
    hits = find_all(backbone.sequence, core)
    assert len(hits) == 1
    plus_one = hits[0] + len(core)
    assert backbone.sequence[plus_one : plus_one + 3] == "AGG"


def test_backbone_has_exactly_two_flanking_bsai_sites(backbone):
    sites = find_all(backbone.sequence, "GGTCTC", both_strands=True)
    assert len(sites) == 2
    t7 = find_all(backbone.sequence, "TAATACGACTCACTATA")[0]
    polya = backbone.find("poly(A) 120")
    assert sites[0] < t7, "upstream site must precede the T7 promoter"
    assert sites[1] >= polya.end, "downstream site must follow the poly(A)"


def test_five_prime_utr_is_kozak_terminated_and_has_no_uaug(backbone):
    utr5 = backbone.find("5' UTR").slice(backbone.sequence)
    assert len(utr5) == 44
    assert utr5.endswith("GCCACC")
    assert "ATG" not in utr5


def test_polya_tract_is_120_pure_adenines(backbone):
    assert backbone.find("poly(A) 120").slice(backbone.sequence) == "A" * 120


# --- assembly -------------------------------------------------------------

@pytest.mark.parametrize("name", [p.name for p in PANEL])
def test_every_panel_member_assembles_and_verifies(backbone, name):
    construct = build_panel_construct(name)
    build = build_plasmid(backbone, construct, f"pVax1_AG_{name}")
    failures = [c for c in verify_plasmid(build, construct, backbone) if not c.ok]
    assert not failures, [(c.name, c.detail) for c in failures]


def test_insert_replaces_exactly_the_egfp_cds(backbone):
    construct = build_panel_construct("P3")
    build = build_plasmid(backbone, construct, "pVax1_AG_P3")
    egfp = backbone.find("eGFP")
    assert build.removed_span == (egfp.start, egfp.end)
    assert build.removed_sequence == egfp.slice(backbone.sequence)
    assert build.insert == construct.orf


def test_plasmid_size_arithmetic(backbone):
    construct = build_panel_construct("P7")
    build = build_plasmid(backbone, construct, "pVax1_AG_P7")
    assert len(build) == 3855 - 720 + len(construct.orf)


def test_egfp_specific_annotations_are_dropped_not_left_lying(backbone):
    """A primer annotation that no longer anneals is a map that lies."""
    build = build_plasmid(backbone, build_panel_construct("P3"), "x")
    assert "eGFP" in build.dropped_features
    assert "eGFP_qPCR_F" in build.dropped_features
    assert "ivtGFP_HiFi_F" in build.dropped_features
    labels = {f.label for f in build.record.features}
    assert "eGFP" not in labels


def test_downstream_features_are_shifted_and_still_correct(backbone):
    construct = build_panel_construct("P3")
    build = build_plasmid(backbone, construct, "x")
    for label in ("NeoR/KanR", "ori", "poly(A) 120", "bGH poly(A) signal"):
        original = backbone.find(label)
        moved = build.record.find(label)
        assert moved.slice(build.record.sequence) == original.slice(backbone.sequence)


def test_upstream_features_keep_their_coordinates(backbone):
    build = build_plasmid(backbone, build_panel_construct("P3"), "x")
    for label in ("CMV enhancer", "CMV promoter", "T7 Promoter", "5' UTR"):
        assert build.record.find(label).start == backbone.find(label).start


def test_module_annotation_is_carried_into_the_plasmid(backbone):
    build = build_plasmid(backbone, build_panel_construct("P7"), "pVax1_AG_P7")
    labels = {f.label for f in build.record.features}
    assert "pVax1_AG_P7_ORF" in labels
    for module in ("HA", "CL1", "P2A", "LAMP1_SP", "FLAG", "LAMP1_TMT"):
        assert module in labels or any(l.startswith(f"{module}_") for l in labels)


def test_assembled_plasmid_round_trips_through_genbank(backbone):
    construct = build_panel_construct("P3")
    build = build_plasmid(backbone, construct, "pVax1_AG_P3")
    again = parse_genbank(write_genbank(build.record))
    assert again.sequence == build.record.sequence
    orf = again.sequence[build.insert_start : build.insert_end]
    assert translate(orf, stop_at_stop=True).rstrip("*") == construct.orf_protein


def test_assembly_refuses_a_backbone_without_unique_rails(backbone):
    doubled = parse_genbank(write_genbank(backbone))
    doubled.sequence = backbone.sequence + backbone.sequence
    with pytest.raises(ValueError, match="occurs 2 times"):
        build_plasmid(doubled, build_panel_construct("P0"), "x")
