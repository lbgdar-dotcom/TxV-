"""Two codon tables, two panels, and the ways they must not mix."""

import pytest

from txv.codon import CODON_USAGE, HUMAN_CODON_USAGE, MOUSE_CODON_USAGE
from txv.pvax1_ag import (
    DEFAULT_SPECIES, MODULES, PANEL, build_panel_construct,
    canonical_encodings, linker_variants, panel_optimizer_config,
)
from txv.seqops import AA_TO_CODONS, translate


# --- the table itself ------------------------------------------------------

def test_the_mouse_table_covers_every_codon():
    assert len(MOUSE_CODON_USAGE) == 64
    assert set(MOUSE_CODON_USAGE) == set(HUMAN_CODON_USAGE)


@pytest.mark.parametrize("aa", sorted(AA_TO_CODONS))
def test_each_amino_acids_fractions_sum_to_one(aa):
    """A usage table is per-amino-acid fractions; it must normalise."""
    total = sum(MOUSE_CODON_USAGE[c] for c in AA_TO_CODONS[aa])
    assert total == pytest.approx(1.0, abs=0.021)


def test_the_two_tables_are_close_but_not_identical():
    """Two mammals. If this ever fails, one table was replaced wholesale."""
    diffs = [abs(MOUSE_CODON_USAGE[c] - HUMAN_CODON_USAGE[c])
             for c in MOUSE_CODON_USAGE]
    assert max(diffs) < 0.10
    assert sum(diffs) / len(diffs) < 0.03
    assert any(d > 0 for d in diffs)


# --- the panels ------------------------------------------------------------

@pytest.mark.parametrize("name", [p.name for p in PANEL])
def test_both_species_give_the_same_protein(name):
    """Only synonymous codons may move. Anything else is a design change."""
    human = build_panel_construct(name)
    mouse = build_panel_construct(name, species="mouse")
    assert human.orf_protein == mouse.orf_protein
    assert len(human.orf) == len(mouse.orf)
    assert translate(mouse.orf, stop_at_stop=True).rstrip("*") == mouse.orf_protein


@pytest.mark.parametrize("name", [p.name for p in PANEL])
def test_the_two_species_differ_somewhere(name):
    """If they were identical the re-optimisation did nothing."""
    assert build_panel_construct(name).orf != \
           build_panel_construct(name, species="mouse").orf


def test_each_species_keeps_one_encoding_per_module():
    """The canonical-encoding rule has to hold *within* a species.

    This is the bug a second table invites: antigen A encoded from the human
    table in one construct and the mouse table in the next would reintroduce
    exactly the confound the whole panel is built to avoid.
    """
    for species in ("human", "mouse"):
        canon = canonical_encodings(species=species)
        for name in [p.name for p in PANEL]:
            orf = build_panel_construct(name, species=species).orf
            for module in PANEL_MODULES(name):
                if module in ("L", "E5"):
                    continue
                assert canon[module] in orf, f"{module} in {name} ({species})"


def PANEL_MODULES(name):
    return next(p for p in PANEL if p.name == name).modules


def test_the_two_panels_do_not_share_encodings_for_changed_modules():
    """A module that differs between tables must differ in every construct."""
    human = canonical_encodings()
    mouse = canonical_encodings(species="mouse")
    changed = [m for m in human if human[m] != mouse[m]]
    assert changed, "no module changed; the mouse table is not being applied"
    for module in changed:
        for name in [p.name for p in PANEL]:
            orf = build_panel_construct(name, species="mouse").orf
            if module in PANEL_MODULES(name):
                assert human[module] not in orf, \
                    f"{name} carries the human encoding of {module}"


# --- configuration ---------------------------------------------------------

def test_the_config_carries_the_right_table():
    assert panel_optimizer_config("mouse").usage == MOUSE_CODON_USAGE
    assert panel_optimizer_config("human").usage == HUMAN_CODON_USAGE
    assert panel_optimizer_config().usage == CODON_USAGE[DEFAULT_SPECIES]


def test_an_unknown_species_is_refused():
    with pytest.raises(ValueError):
        panel_optimizer_config("tardigrade")


def test_the_default_is_still_human():
    """The shipped human package must stay reproducible byte for byte."""
    assert DEFAULT_SPECIES == "human"
    assert build_panel_construct("P3").orf == \
           build_panel_construct("P3", species="human").orf


def test_linkers_are_derived_per_species():
    assert len(linker_variants("mouse" and None, species="mouse")) >= 2
    for variants in (linker_variants(), linker_variants(species="mouse")):
        for dna in variants:
            assert translate(dna, stop_at_stop=False) == MODULES["L"]
