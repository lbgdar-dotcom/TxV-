"""Turn a construct into something you can paste into a vendor order form.

Two decisions have to be made before a construct becomes an order, and both
are easy to get wrong silently.

**1. What are you replacing in the backbone?**

``FragmentMode.ORF`` orders only the coding sequence and drops it between the
backbone's existing UTR rails. Minimal change, minimal risk, and the UTRs stay
whatever the backbone already has.

``FragmentMode.CASSETTE`` orders the whole transcription unit -- 5' UTR, ORF,
3' UTR -- and replaces the backbone's UTRs with the ones in the construct. This
is what you want if the point of the build is to *upgrade* the UTRs, but it
moves the primer-annealing rails, so the IVT primers must be re-designed
against the new sequence. It is not a drop-in.

**2. Is the poly(A) tail encoded or added by PCR?**

A 120-nt A-tract is genuinely hard to synthesise: most vendors will fail,
refuse, or silently deliver a contracted tract, and it contracts again during
propagation in *E. coli*. If your IVT forward/reverse primers already add the
tail by PCR, the encoded tract is doing no work for the transcript and
``PolyAMode.PCR_ADDED`` drops it from the order -- which is usually the single
biggest improvement to synthesis success on a construct like this.
``PolyAMode.ENCODED`` keeps it, for when the plasmid itself must carry it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from .constructs import Construct
from .seqops import clean, find_all, gc_fraction, gc_windows, homopolymer_runs, revcomp

#: pVax1_AG assembly overlaps, for dropping an ORF between the existing rails.
PVAX1_AG_LEFT_OVERLAP = "GAAGAAATATAAGAGCCACC"
PVAX1_AG_RIGHT_OVERLAP = "GCTGCCTTCTGCGGGGCTTG"


class FragmentMode(str, Enum):
    ORF = "orf"
    CASSETTE = "cassette"


class PolyAMode(str, Enum):
    ENCODED = "encoded"
    PCR_ADDED = "pcr_added"


def longest_hairpin(sequence: str, min_loop: int = 3) -> int:
    """Longest perfect inverted repeat (hairpin stem) in a sequence.

    Synthesis vendors fail on strong hairpins because the oligos that build the
    fragment self-anneal instead of extending. A perfect stem is a crude but
    directionally correct proxy; anything at or above ~12 bp is worth redesign.
    """
    best = 0
    n = len(sequence)
    for i in range(n - 8):
        for j in range(i + min_loop + 9, n):
            length = 0
            while (i + length < j - length - min_loop
                   and sequence[i + length] == revcomp(sequence[j - length])):
                length += 1
            best = max(best, length)
    return best


@dataclass
class SynthesisRisk:
    """Why a vendor might reject or fail a fragment."""

    length: int
    gc: float
    longest_homopolymer: tuple[str, int, int] | None
    long_repeats: int
    flags: list[str] = field(default_factory=list)
    min_window_gc: float = 0.5
    max_window_gc: float = 0.5
    hairpin_bp: int = 0

    @property
    def ok(self) -> bool:
        return not self.flags


@dataclass
class OrderFragment:
    name: str
    sequence: str
    mode: FragmentMode
    polya: PolyAMode
    left_overlap: str
    right_overlap: str
    risk: SynthesisRisk
    notes: list[str] = field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.sequence)

    @property
    def insert(self) -> str:
        """The fragment without the assembly overlaps."""
        start = len(self.left_overlap)
        end = len(self.sequence) - len(self.right_overlap) or None
        return self.sequence[start:end]


def assess_synthesis(sequence: str, homopolymer_limit: int = 9,
                     repeat_k: int = 40, gblock: bool = True) -> SynthesisRisk:
    """Flag the properties that make a synthetic fragment fail.

    With ``gblock=True`` the thresholds are the ones a double-stranded gene
    fragment is actually judged against: overall length, overall and *windowed*
    GC, homopolymer tracts, internal repeats and hairpins. Windowed GC matters
    independently of the overall figure -- a fragment can average 60% and still
    be rejected for a single 100-bp window at 80%.
    """
    sequence = clean(sequence)
    gc = gc_fraction(sequence)
    runs = homopolymer_runs(sequence, homopolymer_limit)
    longest = max(runs, key=lambda r: r[1]) if runs else None

    seen: set[str] = set()
    repeats = 0
    for i in range(len(sequence) - repeat_k + 1):
        kmer = sequence[i : i + repeat_k]
        if kmer in seen:
            repeats += 1
        seen.add(kmer)

    windows = [v for _, v in gc_windows(sequence, 100, 10)] or [gc]
    hairpin = longest_hairpin(sequence) if gblock else 0

    flags: list[str] = []
    if gblock:
        if len(sequence) < 125:
            flags.append(
                f"{len(sequence)} nt is below the 125-nt gene-fragment minimum; "
                "order as an oligo pair or pad the fragment"
            )
        if min(windows) < 0.25 or max(windows) > 0.75:
            flags.append(
                f"a 100-nt window reaches {max(windows):.0%} / {min(windows):.0%} GC, "
                "outside the 25-75% window limit"
            )
        if hairpin >= 12:
            flags.append(f"{hairpin}-bp perfect hairpin stem: synthesis oligos "
                         "self-anneal rather than extend")
    if longest and longest[1] >= 20:
        flags.append(
            f"{longest[1]}-nt {longest[0]}-tract at {longest[2]}: most vendors "
            "cannot synthesise this reliably. Consider PolyAMode.PCR_ADDED, or "
            "a segmented tail"
        )
    elif longest:
        flags.append(f"{longest[1]}-nt {longest[0]}-tract at {longest[2]}")
    if not 0.25 <= gc <= 0.75:
        flags.append(f"GC {gc:.1%} is outside the 25-75% most vendors accept")
    if repeats:
        flags.append(f"{repeats} exact repeat(s) of {repeat_k}+ nt")
    if len(sequence) > 3000:
        flags.append(
            f"{len(sequence)} nt exceeds the usual single-fragment limit; "
            "expect it to be split and assembled"
        )
    return SynthesisRisk(
        len(sequence), round(gc, 4), longest, repeats, flags,
        min_window_gc=round(min(windows), 4), max_window_gc=round(max(windows), 4),
        hairpin_bp=hairpin,
    )


def synthesis_fragment(
    construct: Construct,
    mode: FragmentMode | str = FragmentMode.ORF,
    polya: PolyAMode | str = PolyAMode.PCR_ADDED,
    left_overlap: str = PVAX1_AG_LEFT_OVERLAP,
    right_overlap: str = PVAX1_AG_RIGHT_OVERLAP,
) -> OrderFragment:
    """Build the orderable fragment for one construct."""
    mode = FragmentMode(mode)
    polya = PolyAMode(polya)
    notes: list[str] = []

    if mode is FragmentMode.ORF:
        body = construct.orf
        notes.append(
            "ORF only: drops between the backbone's existing UTR rails. The "
            "construct's own UTRs are NOT part of this order, so the backbone's "
            "UTRs are what the transcript will carry."
        )
        if left_overlap.endswith("GCCACC"):
            notes.append(
                "Left overlap ends in the Kozak (GCCACC), so the fragment "
                "starts at ATG and must not repeat it."
            )
        if polya is PolyAMode.ENCODED:
            notes.append(
                "PolyAMode.ENCODED has no effect in ORF mode -- the tail lives "
                "in the backbone, not in this fragment."
            )
    else:
        utr5 = construct.features_of_kind("utr5")[0]
        end_feature = (
            construct.features_of_kind("polya")[0]
            if polya is PolyAMode.ENCODED
            else construct.features_of_kind("utr3")[0]
        )
        body = construct.template[utr5.start : end_feature.end]
        notes.append(
            "Full cassette: 5' UTR + ORF + 3' UTR"
            + (" + poly(A)" if polya is PolyAMode.ENCODED else "")
            + ". This REPLACES the backbone's UTRs, which moves the IVT "
              "primer-annealing rails -- re-design the IVT primers against the "
              "new sequence before ordering."
        )
        if polya is PolyAMode.PCR_ADDED:
            notes.append(
                "Poly(A) omitted from the order: the IVT primer adds the tail, "
                "so the encoded tract does no work for the transcript and only "
                "makes synthesis and propagation harder."
            )

    sequence = left_overlap + body + right_overlap
    return OrderFragment(
        name=construct.name,
        sequence=sequence,
        mode=mode,
        polya=polya,
        left_overlap=left_overlap,
        right_overlap=right_overlap,
        risk=assess_synthesis(sequence),
        notes=notes,
    )


def order_table(fragments: list[OrderFragment]) -> str:
    """CSV ready for a vendor upload or a purchasing spreadsheet."""
    import csv
    import io

    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(
        ["name", "length_bp", "gc", "mode", "polya", "synthesis_flags", "sequence"]
    )
    for fragment in fragments:
        writer.writerow([
            fragment.name, fragment.length, f"{fragment.risk.gc:.4f}",
            fragment.mode.value, fragment.polya.value,
            "; ".join(fragment.risk.flags) or "none",
            fragment.sequence,
        ])
    return buffer.getvalue()


def order_fasta(fragments: list[OrderFragment], width: int = 60) -> str:
    out = []
    for fragment in fragments:
        out.append(
            f">{fragment.name}|{fragment.mode.value}|{fragment.length}bp "
            f"gc={fragment.risk.gc:.3f}"
        )
        out += [
            fragment.sequence[i : i + width]
            for i in range(0, fragment.length, width)
        ]
    return "\n".join(out) + "\n"


__all__ = [
    "FragmentMode", "PolyAMode", "OrderFragment", "SynthesisRisk",
    "AssembledPlasmid", "PVAX1_AG_LEFT_OVERLAP", "PVAX1_AG_RIGHT_OVERLAP",
    "assess_synthesis", "longest_hairpin", "synthesis_fragment",
    "order_table", "order_fasta", "fragment_genbank", "assemble_into_backbone",
]


# ---------------------------------------------------------------------------
# Annotated records and in-silico assembly
# ---------------------------------------------------------------------------

def fragment_genbank(fragment: OrderFragment, construct: Construct) -> str:
    """Annotated GenBank of the fragment *as ordered*.

    This is the record to import into Benchling for an order, and it
    deliberately does **not** show the construct's own UTRs: in
    :attr:`FragmentMode.ORF` those come from the backbone, and a record showing
    them would misrepresent what you are buying.
    """
    from datetime import date
    import textwrap

    sequence = fragment.sequence
    offset = len(fragment.left_overlap) - construct.feature("ORF").start

    features: list[tuple[str, int, int, str, str]] = [
        ("left_homology_arm", 0, len(fragment.left_overlap), "misc_feature",
         "20-bp homology to the linearised backbone; ends in the Kozak."),
        ("right_homology_arm", len(sequence) - len(fragment.right_overlap),
         len(sequence), "misc_feature",
         "20-bp homology to the linearised backbone (3' UTR start)."),
    ]
    keep = {
        "orf": "CDS", "signal_peptide": "sig_peptide", "trafficking": "misc_feature",
        "linker": "misc_feature", "neoepitope": "misc_feature", "tag": "misc_feature",
        "degron": "misc_feature", "skip_peptide": "misc_feature",
        "start": "misc_feature", "stop": "terminator",
        "tumor_associated": "misc_feature", "helper": "misc_feature",
        "full_length": "misc_feature",
    }
    for feature in construct.features:
        if feature.kind not in keep:
            continue
        start, end = feature.start + offset, feature.end + offset
        if start < 0 or end > len(sequence):
            continue
        features.append((feature.name, start, end, keep[feature.kind], feature.note))

    stamp = date.today().strftime("%d-%b-%Y").upper()
    lines = [
        f"LOCUS       {fragment.name[:16]:<16} {len(sequence)} bp    DNA"
        f"     linear   SYN {stamp}",
        f"DEFINITION  {fragment.name} synthesis fragment ({fragment.mode.value} mode), "
        f"{len(sequence)} bp, for homology assembly into the pVax1_AG backbone.",
        f"ACCESSION   {fragment.name}",
        "VERSION     .",
        "KEYWORDS    gene fragment; homology assembly; IVT mRNA.",
        "SOURCE      synthetic construct",
        "  ORGANISM  synthetic construct",
        "FEATURES             Location/Qualifiers",
        f"     {'source':<16}1..{len(sequence)}",
        '                     /organism="synthetic construct"',
        '                     /mol_type="other DNA"',
    ]
    for name, start, end, key, note in sorted(features, key=lambda f: (f[1], -(f[2] - f[1]))):
        lines.append(f"     {key:<16}{start + 1}..{end}")
        lines.append(f'                     /label="{name}"')
        if key == "CDS":
            lines.append('                     /codon_start="1"')
            lines += textwrap.wrap(
                f'/translation="{construct.orf_protein}"', width=58,
                initial_indent=" " * 21, subsequent_indent=" " * 21,
            )
        if note:
            lines += textwrap.wrap(
                f'/note="{note[:220].replace(chr(34), chr(39))}"', width=58,
                initial_indent=" " * 21, subsequent_indent=" " * 21,
            )
    for note in fragment.notes:
        lines += textwrap.wrap(note, width=68, initial_indent="COMMENT     ",
                               subsequent_indent="            ")
    lines.append("ORIGIN")
    for i in range(0, len(sequence), 60):
        chunk = sequence[i : i + 60].lower()
        groups = " ".join(chunk[j : j + 10] for j in range(0, len(chunk), 10))
        lines.append(f"{i + 1:>9} {groups}")
    lines.append("//")
    return "\n".join(lines) + "\n"


@dataclass
class AssembledPlasmid:
    name: str
    sequence: str
    insert_start: int
    insert_end: int
    backbone_length: int

    def __len__(self) -> int:
        return len(self.sequence)


def assemble_into_backbone(
    backbone: str,
    insert: str,
    name: str,
    left_rail: str = "GAAGAAATATAAGAGCCACC",
    right_rail: str = PVAX1_AG_RIGHT_OVERLAP,
) -> AssembledPlasmid:
    """Simulate the assembly, returning the full circular plasmid sequence.

    ``backbone`` is the *whole circular* parent plasmid (e.g. pVax1_AG_eGFP).
    The region between the end of ``left_rail`` and the start of ``right_rail``
    -- that is, the existing ORF -- is replaced by ``insert``.

    Both rails must occur exactly once, or the assembly is ambiguous and this
    raises rather than guessing.
    """
    backbone = clean(backbone)
    insert = clean(insert)

    left_hits = find_all(backbone, left_rail)
    right_hits = find_all(backbone, right_rail)
    for label, hits, rail in (("left", left_hits, left_rail),
                              ("right", right_hits, right_rail)):
        if len(hits) != 1:
            raise ValueError(
                f"{label} rail {rail} occurs {len(hits)} times in the backbone "
                "(expected exactly 1); assembly would be ambiguous"
            )

    start = left_hits[0] + len(left_rail)
    end = right_hits[0]
    if end < start:
        raise ValueError(
            "the right rail precedes the left rail in this linear representation; "
            "rotate the backbone so the insertion site does not span the origin"
        )
    return AssembledPlasmid(
        name=name,
        sequence=backbone[:start] + insert + backbone[end:],
        insert_start=start,
        insert_end=start + len(insert),
        backbone_length=len(backbone),
    )
