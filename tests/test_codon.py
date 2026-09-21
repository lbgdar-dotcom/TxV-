import pytest

from txv.codon import (
    DEFAULT_FORBIDDEN,
    CodonOptimizer,
    OptimizerConfig,
    codon_adaptation_index,
)
from txv.seqops import AA_TO_CODONS, translate

PROTEIN = (
    "MDAMKRGLCCVLLLCGAVFVSPSGGSGGGGSGGKLFEKIRSQLYNPEAAYSLLMWITQCFLPVF"
    "HGDTWKENRQLPSAVIYQYMDDLYVGSDLEIGQHRTKIEELRQHLLRWGLTTPDKKHQKEPPFL"
)


def test_optimized_dna_encodes_the_input_protein():
    result = CodonOptimizer().optimize(PROTEIN)
    assert translate(result.dna, stop_at_stop=False).rstrip("*") == PROTEIN
    assert len(result.dna) == (len(PROTEIN) + 1) * 3
    assert result.dna[-3:] == "TGA"


def test_stop_codon_handling():
    opt = CodonOptimizer()
    assert opt.optimize("MAC", add_stop="TAA").dna.endswith("TAA")
    # A trailing '*' in the protein means "terminate", not "encode a residue".
    assert opt.optimize("MAC*").dna.endswith("TGA")
    # A translational fusion gets no stop at all.
    fusion = opt.optimize("MAC", add_stop=None)
    assert len(fusion.dna) == 9


def test_forbidden_motifs_are_avoided():
    result = CodonOptimizer().optimize(PROTEIN)
    assert result.residual_forbidden == []
    assert result.ok


def test_find_forbidden_detects_both_strands():
    opt = CodonOptimizer()
    assert ("GAATTC", 3) in opt.find_forbidden("AAAGAATTCAAA")
    # GCTCTTC on the reverse strand appears as GAAGAGC on the forward one.
    hits = opt.find_forbidden("AAAGAAGAGCAAA")
    assert any(m.endswith("(rc)") for m, _ in hits)


def test_uridine_weight_reduces_uridine():
    heavy = CodonOptimizer(OptimizerConfig(w_uridine=2.0)).optimize(PROTEIN)
    none = CodonOptimizer(OptimizerConfig(w_uridine=0.0)).optimize(PROTEIN)
    assert heavy.uridine_fraction < none.uridine_fraction
    # Uridine depletion is bought with codon adaptation; that is the trade.
    assert heavy.cai <= none.cai + 1e-9


def test_gc_target_moves_gc():
    low = CodonOptimizer(OptimizerConfig(target_gc=0.40, w_gc=400.0)).optimize(PROTEIN)
    high = CodonOptimizer(OptimizerConfig(target_gc=0.70, w_gc=400.0)).optimize(PROTEIN)
    assert low.gc < high.gc


def test_cai_beats_a_worst_case_encoding():
    result = CodonOptimizer().optimize(PROTEIN)
    worst = "".join(
        min(AA_TO_CODONS[aa], key=lambda c: CodonOptimizer().config.usage.get(c, 0.0))
        for aa in PROTEIN
    )
    assert result.cai > codon_adaptation_index(worst)
    assert 0.0 < result.cai <= 1.0


def test_rare_codons_are_excluded_but_never_orphan_a_residue():
    opt = CodonOptimizer(OptimizerConfig(min_codon_usage=0.99))
    for aa in AA_TO_CODONS:
        if aa == "*":
            continue
        assert opt.choices_for(aa), f"{aa} left with no codon"
    # Leucine's rarest codons are dropped at a normal threshold.
    assert "TTA" not in CodonOptimizer().choices_for("L")


def test_rejects_bad_input():
    opt = CodonOptimizer()
    with pytest.raises(ValueError):
        opt.optimize("")
    with pytest.raises(ValueError):
        opt.optimize("MAZ")
    with pytest.raises(ValueError):
        OptimizerConfig(beam_width=0)
    with pytest.raises(ValueError):
        OptimizerConfig(target_gc=1.5)


def test_repair_preserves_protein_when_a_motif_is_forced():
    # 'EF' can only be encoded GA(A|G) TT(T|C); GAATTC (EcoRI) is reachable but
    # avoidable. Whatever the optimiser does, the protein must survive.
    protein = "MEFEFEFEFEFEFK"
    result = CodonOptimizer().optimize(protein)
    assert translate(result.dna, stop_at_stop=False).rstrip("*") == protein


def test_forbidden_list_is_screened_as_given():
    custom = CodonOptimizer(OptimizerConfig(forbidden=("GGGGGG",)))
    assert custom.find_forbidden("AAGGGGGGAA")
    assert not custom.find_forbidden("AAGAATTCAA")
    assert "AATAAA" in DEFAULT_FORBIDDEN
