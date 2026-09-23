"""Golden Gate assembly of the antigen cassette into an IVT destination vector.

Why this replaces the homology-arm route
----------------------------------------
The pVax1_AG backbone is a dual-purpose vector: it drives the same ORF from a
CMV promoter in mammalian cells *and* from T7 for IVT. For a template that only
ever has to grow in *E. coli* and be transcribed in a tube, the eukaryotic half
is dead weight -- 809 bp of CMV enhancer, CMV promoter and bGH poly(A) signal,
21% of the plasmid, none of it transcribed by T7.

Dropping it does more than shrink the plasmid. The cryptic-splice and internal
``AATAAA`` screens exist because a nuclear transcript would be spliced and
prematurely polyadenylated. With no CMV promoter there is no nuclear transcript,
so those become sequence hygiene rather than functional constraints.

The mechanics
-------------
BsaI cuts ``GGTCTC(1/5)``: one nucleotide past its recognition site on the top
strand and five on the bottom, leaving a **four-nucleotide 5' overhang** whose
sequence you choose. Both the vector and the insert carry their BsaI sites
*facing inward*, so the sites leave with the dropout and the stuffer and the
product contains no BsaI site at all. That is what makes the reaction one-pot:
cut and ligate run together, and every correct ligation is dead to the enzyme.

The two overhangs here are chosen to be invisible in the product:

``GCCA``
    the first four bases of the Kozak. The fragment therefore carries ``CC``
    before its ``ATG``; the product still reads ``GCCACC|ATG`` exactly.
``GCTG``
    the first four bases of the 3' UTR, so the stop codon runs straight into it.

The obvious choice for the left overhang is ``CACC``, the four bases
immediately before the ATG, because then the fragment's payload *is* the ORF.
It was rejected: ``CACC`` is a single mismatch from the reverse complement of
``GCTG``, so the insert's two ends can ligate to each other and the fragment
self-circularises. ``GCCA`` differs from it at all four positions. Two extra
bases in the fragment is a cheap price for removing a background pathway.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from .seqops import clean, revcomp, translate

#: BsaI. Golden Gate here uses BsaI because IVT linearisation is done by PCR,
#: so no Type IIS site has to survive in the finished plasmid.
BSAI = "GGTCTC"
BSAI_SPACER = "A"          # one base between recognition site and overhang

#: Overhangs, chosen so both junctions are scarless (see module docstring).
KOZAK_OVERHANG = "GCCA"
UTR3_OVERHANG = "GCTG"

#: What the destination vector presents at each arm, i.e. the sequence flanking
#: the overhangs. Left arm is the 5' UTR up to the Kozak; right arm is the 3'
#: UTR from just after its first four bases.
VECTOR_LEFT_FLANK = "GAAGAAATATAAGA"        # 5' UTR, ends just before GCCA
VECTOR_RIGHT_FLANK = "CCTTCTGCGGGGCTTG"     # 3' UTR, continues after GCTG

#: Bases the fragment must carry between the left overhang and its ATG, so that
#: overhang + this + ATG reconstitutes GCCACCATG.
KOZAK_REMAINDER = "CC"


@dataclass
class GoldenGatePart:
    """One orderable fragment, with its BsaI sites facing inward."""

    name: str
    sequence: str
    left_overhang: str
    right_overhang: str
    payload: str
    notes: list[str] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def overhead(self) -> int:
        """Bases added to the payload purely to make the assembly work."""
        return len(self.sequence) - len(self.payload)


def make_cassette_part(name: str, orf: str) -> "GoldenGatePart":
    """The orderable fragment for one antigen cassette.

    The payload is ``CC`` + the ORF, because the left overhang sits two bases
    earlier than the ATG (see the module docstring).
    """
    return make_part(name, KOZAK_REMAINDER + clean(orf))


def make_part(
    name: str,
    payload: str,
    left_overhang: str = KOZAK_OVERHANG,
    right_overhang: str = UTR3_OVERHANG,
    spacer: str = BSAI_SPACER,
) -> GoldenGatePart:
    """Wrap a payload in inward-facing BsaI sites.

    ``payload`` is what survives into the product *after* the overhangs, i.e.
    for an antigen cassette it starts at ``ATG`` and ends at the stop codon.
    """
    payload = clean(payload)
    sequence = (
        BSAI + spacer + left_overhang
        + payload
        + right_overhang + spacer + revcomp(BSAI)
    )
    notes = [
        f"BsaI sites face inward and are removed by the reaction; the product "
        f"contains no BsaI site.",
        f"Left overhang {left_overhang} plus the fragment's leading "
        f"{KOZAK_REMAINDER} reconstitutes GCCACC|ATG with no scar.",
        f"Right overhang {right_overhang} is the start of the 3' UTR, so the "
        f"stop codon abuts it directly.",
    ]
    return GoldenGatePart(name, sequence, left_overhang, right_overhang,
                          payload, notes)


# ---------------------------------------------------------------------------
# Simulation — the only honest way to check a Golden Gate design
# ---------------------------------------------------------------------------

@dataclass
class Digested:
    """A fragment's sticky ends after a Type IIS cut."""

    left_overhang: str
    core: str
    right_overhang: str

    @property
    def insert(self) -> str:
        """What this fragment contributes to the product, overhangs included."""
        return self.left_overhang + self.core + self.right_overhang


def digest_bsai(sequence: str, circular: bool = False) -> list[Digested]:
    """Cut a linear or circular sequence at every BsaI site.

    Returns the piece(s) lying between inward-facing sites. This models the
    geometry that matters -- where the overhangs fall -- rather than every
    product of the reaction.
    """
    sequence = clean(sequence)
    search = sequence + sequence[:10] if circular else sequence

    top = [m.start() for m in re.finditer(BSAI, search)]
    bottom = [m.start() for m in re.finditer(revcomp(BSAI), search)]
    if not top or not bottom:
        raise ValueError("sequence lacks an inward-facing BsaI pair")

    # Top-strand site cuts 1 nt downstream of its recognition sequence; the
    # 4-nt overhang begins there. Bottom-strand site cuts the mirror image.
    left_site = min(top)
    right_site = max(bottom)
    left_overhang_start = left_site + len(BSAI) + 1
    right_overhang_end = right_site - 1

    left_overhang = search[left_overhang_start:left_overhang_start + 4]
    right_overhang = search[right_overhang_end - 4:right_overhang_end]
    core = search[left_overhang_start + 4:right_overhang_end - 4]
    return [Digested(left_overhang, core, right_overhang)]


def assemble(vector_left: str, part: Digested, vector_right: str) -> str:
    """Ligate a digested part between two vector arms.

    The overhangs are contributed once, not twice -- which is exactly the
    property that makes the junction scarless, and exactly the one that is easy
    to get wrong on paper.
    """
    return vector_left + part.insert + vector_right


def verify_part(part: GoldenGatePart,
                vector_left: str = VECTOR_LEFT_FLANK,
                vector_right: str = VECTOR_RIGHT_FLANK) -> list[str]:
    """Re-derive the assembly product and report anything wrong with it.

    An empty list means the part assembles as intended.
    """
    problems: list[str] = []

    if len(re.findall(BSAI, part.sequence)) != 1:
        problems.append("expected exactly one top-strand BsaI site")
    if len(re.findall(revcomp(BSAI), part.sequence)) != 1:
        problems.append("expected exactly one bottom-strand BsaI site")

    digested = digest_bsai(part.sequence)[0]
    if digested.left_overhang != part.left_overhang:
        problems.append(
            f"digest yields left overhang {digested.left_overhang}, "
            f"not {part.left_overhang}"
        )
    if digested.right_overhang != part.right_overhang:
        problems.append(
            f"digest yields right overhang {digested.right_overhang}, "
            f"not {part.right_overhang}"
        )
    if digested.core != part.payload:
        problems.append("digest does not release the intended payload")

    product = assemble(vector_left, digested, vector_right)
    if "GCCACCATG" not in product:
        problems.append("product does not reconstitute the Kozak/ATG junction")
    if BSAI in product or revcomp(BSAI) in product:
        problems.append("product still contains a BsaI site; it would be re-cut")

    # Read the ORF the way a ribosome would -- from the Kozak ATG to the first
    # in-frame stop -- rather than by looking for a delimiter that also occurs
    # inside the coding sequence.
    orf_start = product.index("ATG", product.index("GCCACC"))
    protein = translate(product[orf_start:], stop_at_stop=True)
    if not protein.endswith("*"):
        problems.append("no in-frame stop codon downstream of the start")
    else:
        orf_end = orf_start + len(protein) * 3
        expected_end = (len(vector_left) + len(part.left_overhang)
                        + len(part.payload))
        if orf_end != expected_end:
            problems.append(
                f"the first in-frame stop is at {orf_end}, but the payload ends "
                f"at {expected_end}: the ORF does not span the whole insert"
            )
        if product[orf_end:orf_end + 4] != part.right_overhang:
            problems.append("the stop codon does not abut the right overhang")
    return problems


def overhang_fidelity(left: str, right: str) -> list[str]:
    """Flag overhang pairs that can mis-ligate.

    A two-part Golden Gate reaction puts four overhangs in the pot: the
    insert's two, and their reverse complements on the vector arms. Two
    overhangs ligate when one is the reverse complement of the other, so a
    *near* reverse complement is a background pathway. The two that matter here
    are the insert circularising on itself and the vector re-closing empty.
    """
    def mismatches(a: str, b: str) -> int:
        return sum(x != y for x, y in zip(a, b))

    problems = []
    for name, oh in (("left", left), ("right", right)):
        if oh == revcomp(oh):
            problems.append(f"{name} overhang {oh} is palindromic and self-ligates")

    insert_self = mismatches(right, revcomp(left))
    if insert_self <= 1:
        problems.append(
            f"{left} and {right} are {insert_self} mismatch(es) from being "
            "ligation partners: the insert can circularise on itself"
        )
    vector_self = mismatches(revcomp(right), left)
    if vector_self <= 1:
        problems.append(
            f"the vector arms are {vector_self} mismatch(es) from closing empty"
        )
    return problems


__all__ = [
    "BSAI", "KOZAK_OVERHANG", "UTR3_OVERHANG",
    "VECTOR_LEFT_FLANK", "VECTOR_RIGHT_FLANK",
    "GoldenGatePart", "Digested",
    "make_part", "make_cassette_part", "KOZAK_REMAINDER", "digest_bsai", "assemble", "verify_part", "overhang_fidelity",
]
