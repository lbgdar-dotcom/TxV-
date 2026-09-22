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
from .seqops import clean, gc_fraction, homopolymer_runs

#: pVax1_AG assembly overlaps, for dropping an ORF between the existing rails.
PVAX1_AG_LEFT_OVERLAP = "GAAGAAATATAAGAGCCACC"
PVAX1_AG_RIGHT_OVERLAP = "GCTGCCTTCTGCGGGGCTTG"


class FragmentMode(str, Enum):
    ORF = "orf"
    CASSETTE = "cassette"


class PolyAMode(str, Enum):
    ENCODED = "encoded"
    PCR_ADDED = "pcr_added"


@dataclass
class SynthesisRisk:
    """Why a vendor might reject or fail a fragment."""

    length: int
    gc: float
    longest_homopolymer: tuple[str, int, int] | None
    long_repeats: int
    flags: list[str] = field(default_factory=list)

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
                     repeat_k: int = 40) -> SynthesisRisk:
    """Flag the properties that make a synthetic fragment fail."""
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

    flags: list[str] = []
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
    return SynthesisRisk(len(sequence), round(gc, 4), longest, repeats, flags)


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
    "PVAX1_AG_LEFT_OVERLAP", "PVAX1_AG_RIGHT_OVERLAP",
    "assess_synthesis", "synthesis_fragment", "order_table", "order_fasta",
]
