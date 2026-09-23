"""Golden Gate design: the assembly is checked by simulating it, not asserted."""

import pytest

from txv.goldengate import (
    BSAI,
    KOZAK_OVERHANG,
    KOZAK_REMAINDER,
    UTR3_OVERHANG,
    VECTOR_LEFT_FLANK,
    VECTOR_RIGHT_FLANK,
    assemble,
    digest_bsai,
    make_cassette_part,
    make_part,
    overhang_fidelity,
    verify_part,
)
from txv.pvax1_ag import PANEL, build_panel_construct
from txv.seqops import find_all, revcomp, translate


@pytest.fixture(scope="module")
def parts():
    return {
        p.name: make_cassette_part(p.name, build_panel_construct(p.name).orf)
        for p in PANEL
    }


# --- fragment layout -------------------------------------------------------

@pytest.mark.parametrize("name", [p.name for p in PANEL])
def test_every_part_verifies_by_simulated_assembly(parts, name):
    assert verify_part(parts[name]) == []


def test_fragment_carries_exactly_one_bsai_site_per_strand(parts):
    for part in parts.values():
        assert part.sequence.count(BSAI) == 1
        assert part.sequence.count(revcomp(BSAI)) == 1


def test_bsai_sites_face_inward_and_are_lost_on_assembly(parts):
    """The product must be dead to the enzyme, or the reaction undoes itself."""
    for part in parts.values():
        digested = digest_bsai(part.sequence)[0]
        product = assemble(VECTOR_LEFT_FLANK, digested, VECTOR_RIGHT_FLANK)
        assert BSAI not in product
        assert revcomp(BSAI) not in product


def test_digest_releases_exactly_the_payload(parts):
    for name, part in parts.items():
        digested = digest_bsai(part.sequence)[0]
        assert digested.core == part.payload
        assert digested.left_overhang == KOZAK_OVERHANG
        assert digested.right_overhang == UTR3_OVERHANG


def test_payload_is_the_orf_plus_the_kozak_remainder(parts):
    for name, part in parts.items():
        orf = build_panel_construct(name).orf
        assert part.payload == KOZAK_REMAINDER + orf


# --- the product -----------------------------------------------------------

@pytest.mark.parametrize("name", [p.name for p in PANEL])
def test_product_reconstitutes_the_kozak_junction_exactly(parts, name):
    digested = digest_bsai(parts[name].sequence)[0]
    product = assemble(VECTOR_LEFT_FLANK, digested, VECTOR_RIGHT_FLANK)
    assert product.startswith("GAAGAAATATAAGAGCCACCATG")
    assert find_all(product, "GCCACCATG") == [14]


@pytest.mark.parametrize("name", [p.name for p in PANEL])
def test_product_orf_translates_to_the_designed_protein(parts, name):
    construct = build_panel_construct(name)
    digested = digest_bsai(parts[name].sequence)[0]
    product = assemble(VECTOR_LEFT_FLANK, digested, VECTOR_RIGHT_FLANK)
    start = product.index("ATG", product.index("GCCACC"))
    protein = translate(product[start:], stop_at_stop=True).rstrip("*")
    assert protein == construct.orf_protein


@pytest.mark.parametrize("name", [p.name for p in PANEL])
def test_stop_codon_abuts_the_three_prime_utr(parts, name):
    digested = digest_bsai(parts[name].sequence)[0]
    product = assemble(VECTOR_LEFT_FLANK, digested, VECTOR_RIGHT_FLANK)
    assert product.endswith(UTR3_OVERHANG + VECTOR_RIGHT_FLANK)
    stop_and_utr = product[-len(UTR3_OVERHANG) - len(VECTOR_RIGHT_FLANK) - 3:]
    assert stop_and_utr[:3] in ("TGA", "TAA", "TAG")


# --- overhang fidelity -----------------------------------------------------

def test_chosen_overhangs_are_clean():
    assert overhang_fidelity(KOZAK_OVERHANG, UTR3_OVERHANG) == []


def test_the_obvious_overhang_choice_is_correctly_rejected():
    """CACC is the tempting choice and it lets the insert self-circularise."""
    problems = overhang_fidelity("CACC", UTR3_OVERHANG)
    assert problems
    assert any("circularise" in p for p in problems)


def test_palindromic_overhangs_are_rejected():
    assert any("palindromic" in p for p in overhang_fidelity("GATC", "GCTG"))
    assert any("palindromic" in p for p in overhang_fidelity("GCCA", "AATT"))


def test_fidelity_check_sees_a_near_complement_pair():
    # GCCA and TGGT: revcomp(GCCA) = TGGC, one mismatch from TGGT.
    assert overhang_fidelity("GCCA", "TGGT")


# --- what the fragments must not contain -----------------------------------

@pytest.mark.parametrize("name", [p.name for p in PANEL])
def test_no_internal_bsai_site_in_any_orf(name):
    """An internal site would cut the cassette apart during assembly."""
    orf = build_panel_construct(name).orf
    assert find_all(orf, BSAI, both_strands=True) == []


def test_fragments_are_shorter_than_the_homology_arm_version(parts):
    """22 nt of Type IIS overhead against 40 nt of homology arms."""
    for name, part in parts.items():
        gibson = len(build_panel_construct(name).orf) + 40
        assert len(part) < gibson


def test_fragments_stay_inside_gene_fragment_limits(parts):
    for part in parts.values():
        assert 125 <= len(part) <= 3000


# --- the simulator itself --------------------------------------------------

def test_digest_rejects_a_sequence_with_no_bsai_pair():
    with pytest.raises(ValueError):
        digest_bsai("ACGT" * 30)


def test_make_part_round_trips_an_arbitrary_payload():
    part = make_part("x", "ATGAAACCCTGA", left_overhang="GCCA",
                     right_overhang="GCTG")
    assert digest_bsai(part.sequence)[0].core == "ATGAAACCCTGA"
