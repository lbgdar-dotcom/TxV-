"""Serialise constructs to GenBank and FASTA.

GenBank is the practical interchange format here: Benchling, SnapGene and
every sequence-analysis tool ingest it, and unlike FASTA it carries the feature
annotation, which is the part of a construct design that is actually worth
reviewing.
"""

from __future__ import annotations

import textwrap
from datetime import date

from .constructs import Construct
from .seqops import to_rna

#: construct feature kind -> GenBank feature key
_FEATURE_KEY = {
    "promoter": "promoter",
    "utr5": "5'UTR",
    "utr3": "3'UTR",
    "kozak": "regulatory",
    "orf": "CDS",
    "signal_peptide": "sig_peptide",
    "trafficking": "misc_feature",
    "linker": "misc_feature",
    "neoepitope": "misc_feature",
    "tumor_associated": "misc_feature",
    "full_length": "misc_feature",
    "helper": "misc_feature",
    "start": "misc_feature",
    "stop": "misc_feature",
    "polya": "polyA_site",
    "linearization": "misc_feature",
}


def _wrap_qualifier(key: str, value: str) -> list[str]:
    text = str(value).replace('"', "'")
    body = f'/{key}="{text}"'
    return textwrap.wrap(body, width=58, subsequent_indent=" " * 21,
                         initial_indent=" " * 21, break_long_words=True) or [
        " " * 21 + body
    ]


def to_genbank(construct: Construct, molecule: str = "DNA",
               topology: str = "linear", division: str = "SYN") -> str:
    """Render a construct as an annotated GenBank record.

    Coordinates are converted from the 0-based half-open internal
    representation to GenBank's 1-based inclusive convention.
    """
    seq = construct.template
    stamp = date.today().strftime("%d-%b-%Y").upper()
    lines = [
        f"LOCUS       {construct.name[:16]:<16} {len(seq)} bp    {molecule}"
        f"     {topology}   {division} {stamp}",
        f"DEFINITION  IVT mRNA template: {construct.name}, "
        f"{len(construct.cassette.antigens)}-antigen cassette.",
        f"ACCESSION   {construct.name}",
        "VERSION     .",
        "KEYWORDS    IVT mRNA; cancer vaccine; multi-antigen; LNP.",
        "SOURCE      synthetic construct",
        "  ORGANISM  synthetic construct",
        "FEATURES             Location/Qualifiers",
    ]

    source = f"     {'source':<16}1..{len(seq)}"
    lines.append(source)
    lines += _wrap_qualifier("organism", "synthetic construct")
    lines += _wrap_qualifier("mol_type", "other DNA")

    for feature in construct.features:
        key = _FEATURE_KEY.get(feature.kind, "misc_feature")
        lines.append(f"     {key:<16}{feature.start + 1}..{feature.end}")
        lines += _wrap_qualifier("label", feature.name)
        if feature.kind == "orf":
            lines += _wrap_qualifier("translation", construct.orf_protein)
            lines += _wrap_qualifier("codon_start", "1")
        if feature.note:
            lines += _wrap_qualifier("note", feature.note)

    for note in construct.design_notes:
        lines.append(f"COMMENT     {note}")
    lines.append(
        f"COMMENT     Transcript (run-off from T7 +1): "
        f"{construct.transcript_start + 1}..{construct.transcript_end}, "
        f"{construct.transcript_length} nt."
    )

    lines.append("ORIGIN")
    for i in range(0, len(seq), 60):
        chunk = seq[i : i + 60].lower()
        groups = " ".join(chunk[j : j + 10] for j in range(0, len(chunk), 10))
        lines.append(f"{i + 1:>9} {groups}")
    lines.append("//")
    return "\n".join(lines) + "\n"


def to_fasta(construct: Construct, what: str = "template", width: int = 60) -> str:
    """FASTA for ``template``, ``transcript``, ``mrna``, ``orf`` or ``protein``."""
    bodies = {
        "template": (construct.template, "DNA template for IVT"),
        "transcript": (construct.transcript, "run-off transcript, DNA alphabet"),
        "mrna": (to_rna(construct.transcript), "run-off transcript, RNA alphabet"),
        "orf": (construct.orf, "open reading frame"),
        "protein": (construct.orf_protein, "encoded polypeptide"),
    }
    try:
        seq, description = bodies[what]
    except KeyError:
        raise ValueError(
            f"unknown target {what!r}; choose from {sorted(bodies)}"
        ) from None
    header = f">{construct.name}|{what} {description} len={len(seq)}"
    wrapped = "\n".join(seq[i : i + width] for i in range(0, len(seq), width))
    return f"{header}\n{wrapped}\n"


def feature_table(construct: Construct) -> str:
    """Plain-text feature table for a design review or a lab notebook entry."""
    rows = [f"{'feature':<24}{'kind':<18}{'start':>7}{'end':>8}{'len':>7}"]
    rows.append("-" * 64)
    for feature in construct.features:
        rows.append(
            f"{feature.name[:23]:<24}{feature.kind[:17]:<18}"
            f"{feature.start + 1:>7}{feature.end:>8}{len(feature):>7}"
        )
    return "\n".join(rows)


__all__ = ["to_genbank", "to_fasta", "feature_table"]
