"""In-silico validation primitives, and proof that each one can fail.

Every check here is paired with a negative control. A validator that cannot
report a problem is worse than no validator, because it produces a clean run
that reads like evidence.
"""

import pytest

from txv.codon import HUMAN_CODON_USAGE, gc_floor, uridine_floor
from txv.insilico import (
    KD, TM_MIN_LENGTH, cai, cpg_stats, hydropathy, signal_peptide_call,
    start_accessibility, tm_segments,
)
from txv.pvax1_ag import MODULES
from txv.seqops import AA_TO_CODONS, CODON_TO_AA

RNA = pytest.importorskip("RNA", reason="ViennaRNA not installed")


# --- codon adaptation index ------------------------------------------------

def _best_codon(aa: str) -> str:
    return max((c for c in AA_TO_CODONS[aa] if CODON_TO_AA[c] != "*"),
               key=lambda c: HUMAN_CODON_USAGE[c])


def _worst_codon(aa: str) -> str:
    return min((c for c in AA_TO_CODONS[aa] if CODON_TO_AA[c] != "*"),
               key=lambda c: HUMAN_CODON_USAGE[c])


def test_cai_is_one_for_the_most_frequent_codon_everywhere():
    protein = "LRSAVPTGKEQD"
    assert cai("".join(_best_codon(a) for a in protein)) == pytest.approx(1.0)


def test_cai_falls_when_the_rarest_codons_are_used():
    protein = "LRSAVPTGKEQD"
    best = cai("".join(_best_codon(a) for a in protein))
    worst = cai("".join(_worst_codon(a) for a in protein))
    assert worst < 0.5 < best


def test_cai_ignores_the_stop_codon():
    assert cai("ATGTGA") == cai("ATG")


# --- CpG -------------------------------------------------------------------

def test_cpg_counts_only_the_cg_dinucleotide():
    count, _ = cpg_stats("CGCGCG")
    assert count == 3
    assert cpg_stats("GCGCGC")[0] == 2


def test_cpg_observed_over_expected_normalises_out_composition():
    """A GC-rich sequence with CpG avoided must not look CpG-rich."""
    depleted = "GCGCCAGGCCAGGCCAGGCCAGGCCAGGCCAGG"
    dense = "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"
    assert cpg_stats(depleted)[1] < cpg_stats(dense)[1]


# --- hydropathy ------------------------------------------------------------

def test_hydropathy_separates_hydrophobic_from_charged():
    assert hydropathy("L" * 19)[0] == pytest.approx(KD["L"])
    assert hydropathy("D" * 19)[0] == pytest.approx(KD["D"])
    assert hydropathy("L" * 19)[0] > hydropathy("D" * 19)[0]


def test_hydropathy_is_empty_below_the_window():
    assert hydropathy("LLL", window=19) == []


def test_tm_segments_finds_a_membrane_span():
    segments = tm_segments("D" * 10 + "L" * 25 + "D" * 10)
    assert segments
    lo, hi, _ = segments[0]
    assert hi - lo + 1 >= TM_MIN_LENGTH


def test_tm_segments_finds_nothing_in_a_soluble_protein():
    assert tm_segments("DEKRDEKRDEKRDEKRDEKRDEKRDEKRDEKR") == []


def test_tm_segments_rejects_a_span_that_is_too_short_to_cross():
    """Nine leucines are hydrophobic but cannot span a bilayer."""
    assert tm_segments("D" * 20 + "L" * 9 + "D" * 20) == []


@pytest.mark.parametrize("module", ["CTLA4_TMT", "LAMP1_TMT"])
def test_the_designed_tm_modules_contain_a_real_span(module):
    segments = tm_segments(MODULES[module])
    assert any(hi - lo + 1 >= TM_MIN_LENGTH for lo, hi, _ in segments)


# --- signal peptides -------------------------------------------------------

@pytest.mark.parametrize("module", ["CTLA4_SP", "LAMP1_SP"])
def test_the_designed_signal_peptides_have_the_right_architecture(module):
    call = signal_peptide_call(MODULES[module])
    assert call.ok, call.problems


def test_signal_peptide_call_rejects_a_peptide_with_no_h_region():
    call = signal_peptide_call("MKKDEDEDEDEDEDEDEDEDEKKA")
    assert not call.ok
    assert any("h-region" in p for p in call.problems)


def test_signal_peptide_call_rejects_a_bulky_minus_one_residue():
    """Signal peptidase I needs a small residue at -1."""
    peptide = MODULES["LAMP1_SP"][:-1] + "W"
    call = signal_peptide_call(peptide)
    assert not call.ok
    assert any("-1 residue" in p for p in call.problems)


def test_signal_peptide_call_rejects_a_negatively_charged_n_region():
    call = signal_peptide_call("MDDDD" + "LLLLLLLLLLLL" + "AHGASA")
    assert not call.ok
    assert any("n-region" in p for p in call.problems)


def test_signal_peptide_call_survives_a_peptide_too_short_to_score():
    call = signal_peptide_call("MA")
    assert not call.ok


# --- start-codon accessibility --------------------------------------------

def test_a_hairpin_over_the_start_codon_lowers_accessibility():
    """The measure must actually respond to structure on the AUG."""
    stem = "GGGGCCCCGGGGCCCC"
    open_seq = "A" * 40 + "ATG" + "A" * 40
    closed = "A" * 24 + stem + "ATG" + "GGGGCCCCGGGGCCCC"[::-1] + "A" * 24
    assert start_accessibility(open_seq, 40) > start_accessibility(closed, 40)


def test_a_polya_context_is_almost_fully_accessible():
    assert start_accessibility("A" * 40 + "ATG" + "A" * 40, 40) > 0.9


# --- the composition floors ------------------------------------------------

def test_gc_floor_is_bounded_by_amino_acid_composition():
    """Ala, Pro and Gly are GCN/CCN/GGN: no codon has fewer than two G+C."""
    for aa in "APG":
        total, codons = gc_floor(aa * 10)
        assert total / (codons * 3) >= 2 / 3


def test_arginine_escapes_the_gc_floor_through_AGA():
    """Arg looks GC-locked as CGN but AGA carries only one G+C."""
    total, codons = gc_floor("R" * 10)
    assert total / (codons * 3) == pytest.approx(1 / 3)


def test_gc_floor_is_low_for_residues_with_at_rich_codons():
    total, codons = gc_floor("KKKKNNNN")
    assert total / (codons * 3) == 0.0


def test_gc_floor_never_exceeds_what_the_optimiser_achieves():
    from txv.pvax1_ag import canonical_encodings
    from txv.seqops import gc_fraction
    canon = canonical_encodings()
    for module, dna in canon.items():
        if module in ("L", "E5"):
            continue
        total, codons = gc_floor(MODULES[module], min_codon_usage=0.10)
        assert gc_fraction(dna) >= total / (codons * 3) - 1e-9, module


def test_uridine_floor_and_gc_floor_disagree_as_they_should():
    """They pull in opposite directions; that tension is the design problem."""
    protein = "FFFFYYYY"          # no U-free codon, but AT-rich codons exist
    u_total, codons = uridine_floor(protein)
    g_total, _ = gc_floor(protein)
    assert u_total > 0
    assert g_total / (codons * 3) < 0.5
