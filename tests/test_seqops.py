import pytest

from txv import seqops as s


def test_genetic_code_is_complete_and_correct():
    assert len(s.CODON_TO_AA) == 64
    assert s.CODON_TO_AA["ATG"] == "M"
    assert s.CODON_TO_AA["TGG"] == "W"
    assert {s.CODON_TO_AA[c] for c in s.STOP_CODONS} == {"*"}
    # Every sense codon maps back through AA_TO_CODONS.
    for codon, aa in s.CODON_TO_AA.items():
        assert codon in s.AA_TO_CODONS[aa]


def test_translate_and_revcomp():
    assert s.translate("ATGGCCTGTTAA") == "MAC*"
    assert s.translate("ATGGCCTGTTAAGGG") == "MAC*"
    assert s.translate("ATGGCCTGTTAAGGG", stop_at_stop=False) == "MAC*G"
    assert s.revcomp("ATGC") == "GCAT"
    assert s.revcomp(s.revcomp("ATTGCA")) == "ATTGCA"


def test_clean_normalises_and_rejects_junk():
    assert s.clean(" au gc\n") == "ATGC"
    with pytest.raises(ValueError):
        s.clean("ATGX")


def test_gc_and_base_fractions():
    assert s.gc_fraction("GGCC") == 1.0
    assert s.gc_fraction("ATAT") == 0.0
    assert s.gc_fraction("") == 0.0
    assert s.base_fraction("ATTT", "T") == 0.75


def test_find_all_handles_iupac_and_both_strands():
    assert s.find_all("AAGAATTCAA", "GAATTC") == [2]
    # GGTRAGT matches GGTAAGT and GGTGAGT.
    assert s.find_all("CCGGTAAGTCC", "GGTRAGT") == [2]
    assert s.find_all("CCGGTGAGTCC", "GGTRAGT") == [2]
    # GAATTC is its own reverse complement, so both-strand search is stable.
    assert s.find_all("AAGAATTCAA", "GAATTC", both_strands=True) == [2]
    # A non-palindromic motif is found only on the strand it is on.
    assert s.find_all("AAGGATCCAA", "GGATCC") == [2]
    assert s.find_all("TTTGCTCTTCTTT", "GCTCTTC", both_strands=True) == [3]
    assert s.find_all("TTTGAAGAGCTTT", "GCTCTTC", both_strands=True) == [3]


def test_homopolymer_helpers():
    assert s.longest_homopolymer("ATTTTTGC") == ("T", 5, 1)
    assert s.longest_homopolymer("") == ("", 0, 0)
    runs = s.homopolymer_runs("AAAACGTTTTTTT", 5)
    assert runs == [("T", 7, 6)]


def test_has_internal_stop():
    assert not s.has_internal_stop("ATGGCCTAA")
    assert s.has_internal_stop("ATGTAAGCCTAA")


def test_gc_windows_shorter_than_window():
    assert s.gc_windows("GGCC", window=100) == [(0, 1.0)]
    assert s.gc_windows("", window=100) == []
