import re

from txv.export import feature_table, to_fasta, to_genbank
import pytest


def test_genbank_has_the_required_records(construct):
    text = to_genbank(construct)
    assert text.startswith("LOCUS")
    assert text.rstrip().endswith("//")
    for token in ("DEFINITION", "FEATURES", "ORIGIN", "source"):
        assert token in text
    assert f"{len(construct.template)} bp" in text


def test_genbank_coordinates_are_one_based_inclusive(construct):
    text = to_genbank(construct)
    promoter = construct.feature("T7_promoter")
    assert f"{promoter.start + 1}..{promoter.end}" in text
    assert "1..%d" % len(construct.template) in text


def test_genbank_origin_block_reproduces_the_sequence(construct):
    text = to_genbank(construct)
    origin = text.split("ORIGIN\n", 1)[1].rsplit("//", 1)[0]
    bases = "".join(re.findall(r"[acgtn]+", origin)).upper()
    assert bases == construct.template


def test_genbank_carries_the_translation(construct):
    text = to_genbank(construct)
    flat = "".join(
        line.strip() for line in text.splitlines() if line.startswith(" " * 21)
    )
    assert construct.orf_protein[:30] in flat.replace('"', "")


def test_fasta_targets(construct):
    for target in ("template", "transcript", "mrna", "orf", "protein"):
        text = to_fasta(construct, target)
        assert text.startswith(f">{construct.name}|{target}")
        body = "".join(text.splitlines()[1:])
        assert body
        if target == "mrna":
            assert "T" not in body and "U" in body
        elif target != "protein":
            assert set(body) <= set("ACGTN")


def test_fasta_line_wrapping(construct):
    lines = to_fasta(construct, "template", width=60).splitlines()[1:]
    assert all(len(line) <= 60 for line in lines)
    assert "".join(lines) == construct.template


def test_unknown_fasta_target_rejected(construct):
    with pytest.raises(ValueError):
        to_fasta(construct, "nonsense")


def test_feature_table_lists_every_feature(construct):
    table = feature_table(construct)
    assert len(table.splitlines()) == len(construct.features) + 2
    assert "T7_promoter" in table
