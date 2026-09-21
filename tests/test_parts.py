import json

import pytest

from txv.parts import (
    Part,
    PartKind,
    PartRegistry,
    Provenance,
    default_registry,
    load_parts_file,
    segmented_polya,
)


def test_part_requires_exactly_one_sequence_form():
    with pytest.raises(ValueError):
        Part("x", PartKind.LINKER, Provenance.CANONICAL)
    with pytest.raises(ValueError):
        Part("x", PartKind.LINKER, Provenance.CANONICAL, dna="ATG", protein="M")


def test_part_normalises_and_validates():
    dna = Part("x", PartKind.UTR5, Provenance.CANONICAL, dna=" a tg c\n")
    assert dna.dna == "ATGC"
    protein = Part("y", PartKind.LINKER, Provenance.CANONICAL, protein=" ggs ")
    assert protein.protein == "GGS"
    with pytest.raises(ValueError):
        Part("z", PartKind.LINKER, Provenance.CANONICAL, protein="GGZ")
    with pytest.raises(ValueError):
        Part("z", PartKind.UTR5, Provenance.CANONICAL, dna="ATGX")


def test_rna_input_is_accepted_as_dna():
    assert Part("x", PartKind.UTR5, Provenance.CANONICAL, dna="AUGC").dna == "ATGC"


def test_t7_promoter_ends_in_the_plus_one_g():
    part = default_registry().get("T7_promoter")
    assert part.dna == "TAATACGACTCACTATAG"
    assert part.provenance is Provenance.CANONICAL


def test_placeholders_are_flagged():
    registry = default_registry()
    names = {p.name for p in registry.placeholders()}
    assert {"UTR5_placeholder", "UTR3_placeholder", "MITD_placeholder"} <= names
    assert registry.get("UTR5_placeholder").is_placeholder
    assert not registry.get("T7_promoter").is_placeholder


def test_placeholder_utrs_carry_no_upstream_aug():
    """A stand-in must not itself break the construct it stands in for."""
    assert "ATG" not in default_registry().get("UTR5_placeholder").dna


def test_registry_lookup_and_membership():
    registry = default_registry()
    assert "T7_promoter" in registry
    assert len(registry) == len(registry.names())
    with pytest.raises(KeyError):
        registry.get("nope")


def test_registry_rejects_silent_overwrite_when_asked():
    registry = PartRegistry([])
    part = Part("x", PartKind.LINKER, Provenance.CANONICAL, protein="GGS")
    registry.add(part)
    with pytest.raises(KeyError):
        registry.add(part, overwrite=False)
    registry.add(part)  # overwrite=True is the default


def test_of_kind_filters():
    registry = default_registry()
    linkers = registry.of_kind(PartKind.LINKER)
    assert linkers and all(p.kind is PartKind.LINKER for p in linkers)


def test_load_parts_file_overrides_builtins(tmp_path):
    path = tmp_path / "parts.json"
    path.write_text(json.dumps({"parts": [{
        "name": "UTR5_placeholder", "kind": "utr5", "provenance": "literature",
        "dna": "GGGAAATAAGAGAGAAAAGAAGAGTAAGAAG", "source": "internal",
    }]}))
    registry = load_parts_file(path)
    part = registry.get("UTR5_placeholder")
    assert not part.is_placeholder
    assert part.source == "internal"
    assert "T7_promoter" in registry, "built-ins survive the merge"


def test_segmented_polya_structure():
    part = segmented_polya(first=30, linker="GCATATGACT", second=70)
    assert part.dna == "A" * 30 + "GCATATGACT" + "A" * 70
    assert part.is_placeholder, "the spacer is programme-specific"
    assert part.kind is PartKind.POLYA
